"""Neutral-atom X-ray structure factors from canonical crystal structures."""

from __future__ import annotations

import cmath
import math
from typing import Any

from cristma.crystallography import SpaceGroupSetting
from cristma.reference_data import ElementCatalog
from cristma.structure import CrystalStructure, DisplacementParameters
from cristma.symmetry import expand_orbit

from .diagnostics import (
    DiffractionInvariantError,
    STRUCTURE_FACTOR_CELL_MISMATCH,
    STRUCTURE_FACTOR_EXTINCT_NONZERO,
    STRUCTURE_FACTOR_SETTING_MISMATCH,
    STRUCTURE_FACTOR_SYMMETRY_MISMATCH,
    STRUCTURE_FACTOR_UNSUPPORTED_ANISOTROPIC_ADP,
    STRUCTURE_FACTOR_UNSUPPORTED_SPECIES,
)
from .models import ReflectionSet
from .reciprocal import cell_fingerprint
from .structure_factor_models import (
    StructureFactor,
    StructureFactorProvenance,
    StructureFactorSet,
    XRayScatteringContext,
)


CALCULATOR_METHOD = "independent_site_symmetry_expansion"
CALCULATOR_VERSION = "1"


def _operation_key(operation: Any) -> tuple[object, object]:
    normalized = operation.normalized()
    return normalized.rotation, normalized.translation


def _operation_set(operations: tuple[Any, ...]) -> tuple[tuple[object, object], ...]:
    return tuple(sorted(_operation_key(operation) for operation in operations))


def _raise(code: str, message: str, **evidence: object) -> None:
    raise DiffractionInvariantError(code, message, evidence)


def _isotropic_factor(
    displacement: DisplacementParameters | None,
    s: float,
    *,
    site_id: str,
) -> float:
    if displacement is None:
        return 1.0
    if displacement.kind == "U_aniso" or displacement.tensor is not None:
        _raise(
            STRUCTURE_FACTOR_UNSUPPORTED_ANISOTROPIC_ADP,
            "anisotropic displacement parameters are not supported in v1",
            site_id=site_id,
            displacement_kind=displacement.kind,
        )
    if displacement.kind not in {"U_iso", "B_iso"}:
        raise ValueError(f"unsupported displacement kind: {displacement.kind!r}")
    if displacement.isotropic is None or displacement.isotropic.value is None:
        raise ValueError("isotropic displacement value is required")
    value = float(displacement.isotropic.value)
    if not math.isfinite(value) or value < 0:
        raise ValueError("isotropic displacement value must be finite and non-negative")
    exponent = -8.0 * math.pi**2 * value * s**2 if displacement.kind == "U_iso" else -value * s**2
    return math.exp(exponent)


class StructureFactorCalculator:
    """Calculate neutral-atom X-ray amplitudes without powder corrections."""

    def __init__(
        self,
        *,
        absolute_tolerance: float = 1e-12,
        relative_tolerance: float = 1e-12,
    ) -> None:
        self.absolute_tolerance = self._tolerance(
            absolute_tolerance, "absolute_tolerance"
        )
        self.relative_tolerance = self._tolerance(
            relative_tolerance, "relative_tolerance"
        )

    @staticmethod
    def _tolerance(value: float, name: str) -> float:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise TypeError(f"{name} must be a real number")
        normalized = float(value)
        if not math.isfinite(normalized) or normalized < 0:
            raise ValueError(f"{name} must be finite and non-negative")
        return normalized

    def get_config(self) -> dict[str, float]:
        return {
            "absolute_tolerance": self.absolute_tolerance,
            "relative_tolerance": self.relative_tolerance,
        }

    def clone(self, **overrides: float) -> "StructureFactorCalculator":
        unknown = set(overrides) - set(self.get_config())
        if unknown:
            raise TypeError(f"unknown configuration field: {sorted(unknown)[0]}")
        return type(self)(**(self.get_config() | overrides))

    def calculate(
        self,
        structure: CrystalStructure,
        space_group: SpaceGroupSetting,
        reflections: ReflectionSet,
        context: XRayScatteringContext,
    ) -> StructureFactorSet:
        if not isinstance(structure, CrystalStructure):
            raise TypeError("structure must be CrystalStructure")
        if not isinstance(space_group, SpaceGroupSetting):
            raise TypeError("space_group must be SpaceGroupSetting")
        if not isinstance(reflections, ReflectionSet):
            raise TypeError("reflections must be ReflectionSet")
        if not isinstance(context, XRayScatteringContext):
            raise TypeError("context must be XRayScatteringContext")

        if reflections.space_group_setting_id != space_group.setting_id:
            _raise(
                STRUCTURE_FACTOR_SETTING_MISMATCH,
                "ReflectionSet belongs to another space-group setting",
                reflection_setting_id=reflections.space_group_setting_id,
                supplied_setting_id=space_group.setting_id,
            )
        fingerprint = cell_fingerprint(structure.cell)
        if reflections.provenance.cell_fingerprint != fingerprint:
            _raise(
                STRUCTURE_FACTOR_CELL_MISMATCH,
                "ReflectionSet belongs to another unit cell",
                reflection_cell_fingerprint=reflections.provenance.cell_fingerprint,
                structure_cell_fingerprint=fingerprint,
            )
        if structure.space_group is None or _operation_set(
            structure.space_group.operations
        ) != _operation_set(space_group.symmetry_operations):
            _raise(
                STRUCTURE_FACTOR_SYMMETRY_MISMATCH,
                "structure symmetry operations do not match the supplied setting",
                supplied_setting_id=space_group.setting_id,
                structure_operation_count=(
                    0 if structure.space_group is None else len(structure.space_group.operations)
                ),
                setting_operation_count=len(space_group.symmetry_operations),
            )

        for site in structure.sites:
            if site.displacement is not None and (
                site.displacement.kind == "U_aniso" or site.displacement.tensor is not None
            ):
                _raise(
                    STRUCTURE_FACTOR_UNSUPPORTED_ANISOTROPIC_ADP,
                    "anisotropic displacement parameters are not supported in v1",
                    site_id=site.id,
                    site_label=site.label,
                    displacement_kind=site.displacement.kind,
                )

        atoms = tuple(
            atom
            for site in structure.sites
            for atom in expand_orbit(
                site,
                space_group.symmetry_operations,
                cell=structure.cell,
                structure_id=structure.id,
            )
        )
        elements = ElementCatalog.default()
        results: list[StructureFactor] = []
        for reflection in reflections.reflections:
            h, k, l = reflection.representative_hkl.as_tuple()
            s = 1.0 / (2.0 * reflection.d_spacing)
            contributions: list[complex] = []
            for atom in atoms:
                amplitudes: list[float] = []
                for component in atom.components:
                    element = component.element
                    if element is None:
                        _raise(
                            STRUCTURE_FACTOR_UNSUPPORTED_SPECIES,
                            "species has no supported parent element",
                            site_id=atom.source_site_id,
                            species=component.species.label,
                        )
                    atomic_number = elements.by_symbol(element).atomic_number
                    try:
                        f0 = context.form_factors.evaluate(atomic_number, s)
                    except ValueError as exc:
                        _raise(
                            STRUCTURE_FACTOR_UNSUPPORTED_SPECIES,
                            "species lies outside the neutral-atom form-factor table",
                            site_id=atom.source_site_id,
                            species=component.species.label,
                            atomic_number=atomic_number,
                            reason=str(exc),
                        )
                    amplitudes.append(float(component.occupancy.value) * f0)
                site_amplitude = math.fsum(amplitudes)
                thermal = _isotropic_factor(
                    atom.displacement,
                    s,
                    site_id=atom.source_site_id,
                )
                x, y, z = atom.fractional
                phase = cmath.exp(2j * math.pi * (h * x + k * y + l * z))
                contributions.append(site_amplitude * thermal * phase)

            raw = complex(
                math.fsum(value.real for value in contributions),
                math.fsum(value.imag for value in contributions),
            )
            contribution_scale = math.fsum(abs(value) for value in contributions)
            extinction_tolerance = max(
                self.absolute_tolerance,
                self.relative_tolerance * contribution_scale,
            )
            normalized_to_zero = reflection.extinction.absent
            if normalized_to_zero and abs(raw) > extinction_tolerance:
                _raise(
                    STRUCTURE_FACTOR_EXTINCT_NONZERO,
                    "systematically absent reflection has a nonzero structure factor",
                    reflection_id=reflection.reflection_id,
                    hkl=reflection.representative_hkl.as_tuple(),
                    space_group_setting_id=space_group.setting_id,
                    raw_f_complex=(raw.real, raw.imag),
                    raw_amplitude=abs(raw),
                    contribution_scale=contribution_scale,
                    tolerance_used=extinction_tolerance,
                    absolute_tolerance=self.absolute_tolerance,
                    relative_tolerance=self.relative_tolerance,
                    extinction_evidence=reflection.extinction.evidence,
                )
            published = 0j if normalized_to_zero else raw
            provenance = StructureFactorProvenance(
                method=CALCULATOR_METHOD,
                version=CALCULATOR_VERSION,
                space_group_setting_id=space_group.setting_id,
                cell_fingerprint=fingerprint,
                table_id=context.table_id,
                table_version=context.table_version,
                raw_f_complex=raw,
                contribution_scale=contribution_scale,
                extinction_tolerance=extinction_tolerance,
                normalized_to_zero=normalized_to_zero,
            )
            results.append(
                StructureFactor(
                    reflection_id=reflection.reflection_id,
                    representative_hkl=reflection.representative_hkl,
                    f_complex=published,
                    provenance=provenance,
                )
            )
        return StructureFactorSet(tuple(results), reflections, context)


__all__ = ["StructureFactorCalculator"]
