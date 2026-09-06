"""Intrinsic powder lines from validated structure factors and radiation."""

from __future__ import annotations

import hashlib
import math

from .diagnostics import (
    DiffractionInvariantError,
    POWDER_FRIEDEL_D_SPACING_MISMATCH,
    POWDER_MISSING_FRIEDEL_MATE,
    POWDER_NONRECIPROCAL_FRIEDEL_LINK,
    POWDER_SCATTERING_PROBE_MISMATCH,
)
from .powder_models import (
    PowderLine,
    PowderLineProvenance,
    PowderLineSet,
    PowderReflectionFamily,
    RadiationProbe,
    RadiationSpectrum,
)
from .structure_factor_models import StructureFactorSet


CALCULATOR_METHOD = "friedel_grouped_bragg_lines"
CALCULATOR_VERSION = "1"
FRIEDEL_D_RELATIVE_TOLERANCE = 1e-10


def _stable_id(prefix: str, values: tuple[str, ...]) -> str:
    digest = hashlib.sha256()
    for value in values:
        encoded = value.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, byteorder="big"))
        digest.update(encoded)
    return f"{prefix}:{digest.hexdigest()}"


def _raise(code: str, message: str, **evidence: object) -> None:
    raise DiffractionInvariantError(code, message, evidence)


class PowderLineCalculator:
    """Convert structure factors into uncorrected radiation-dependent lines."""

    def calculate(
        self,
        structure_factors: StructureFactorSet,
        spectrum: RadiationSpectrum,
    ) -> PowderLineSet:
        if not isinstance(structure_factors, StructureFactorSet):
            raise TypeError("structure_factors must be StructureFactorSet")
        if not isinstance(spectrum, RadiationSpectrum):
            raise TypeError("spectrum must be RadiationSpectrum")
        if spectrum.probe is not RadiationProbe.XRAY:
            _raise(
                POWDER_SCATTERING_PROBE_MISMATCH,
                "X-ray structure factors require an X-ray wavelength spectrum",
                structure_factor_context=type(structure_factors.context).__name__,
                radiation_probe=spectrum.probe.value,
                radiation_source_id=spectrum.source_id,
            )

        reflections = structure_factors.reflection_set.reflections
        by_reflection_id = {
            reflection.reflection_id: reflection for reflection in reflections
        }
        factors_by_id = {
            factor.reflection_id: factor
            for factor in structure_factors.structure_factors
        }
        normalized_weights = spectrum.normalized_weights
        consumed: set[str] = set()
        families: list[PowderReflectionFamily] = []
        skipped_components = 0

        for reflection in reflections:
            if reflection.reflection_id in consumed:
                continue
            mate_id = reflection.friedel_mate_id
            if mate_id is None:
                members = (reflection,)
            else:
                mate = by_reflection_id.get(mate_id)
                if mate is None or mate_id not in factors_by_id:
                    _raise(
                        POWDER_MISSING_FRIEDEL_MATE,
                        "Friedel mate is missing from the source StructureFactorSet",
                        reflection_id=reflection.reflection_id,
                        friedel_mate_id=mate_id,
                    )
                if mate.friedel_mate_id != reflection.reflection_id:
                    _raise(
                        POWDER_NONRECIPROCAL_FRIEDEL_LINK,
                        "Friedel links are not reciprocal",
                        reflection_id=reflection.reflection_id,
                        friedel_mate_id=mate_id,
                        mate_backlink=mate.friedel_mate_id,
                    )
                members = tuple(
                    sorted((reflection, mate), key=lambda item: item.reflection_id)
                )
                if not math.isclose(
                    members[0].d_spacing,
                    members[1].d_spacing,
                    rel_tol=FRIEDEL_D_RELATIVE_TOLERANCE,
                    abs_tol=0.0,
                ):
                    _raise(
                        POWDER_FRIEDEL_D_SPACING_MISMATCH,
                        "Friedel mates have inconsistent d spacings",
                        reflection_ids=tuple(item.reflection_id for item in members),
                        d_spacings=tuple(item.d_spacing for item in members),
                        relative_tolerance=FRIEDEL_D_RELATIVE_TOLERANCE,
                    )

            member_ids = tuple(item.reflection_id for item in members)
            consumed.update(member_ids)
            if any(item.extinction.absent for item in members):
                continue

            family_id = _stable_id("powder-family", member_ids)
            d_spacing = members[0].d_spacing
            family_strength = math.fsum(
                member.multiplicity_crystallographic
                * factors_by_id[member.reflection_id].f_squared
                for member in members
            )
            family_multiplicity = sum(
                item.multiplicity_crystallographic for item in members
            )
            lines: list[PowderLine] = []
            for component, normalized_weight in zip(
                spectrum.components,
                normalized_weights,
                strict=True,
            ):
                bragg_argument = component.wavelength_angstrom / (2.0 * d_spacing)
                if bragg_argument > 1.0:
                    skipped_components += 1
                    continue
                two_theta_deg = math.degrees(2.0 * math.asin(bragg_argument))
                lines.append(
                    PowderLine(
                        line_id=_stable_id(
                            "powder-line",
                            (family_id, component.component_id),
                        ),
                        family_id=family_id,
                        radiation_component_id=component.component_id,
                        wavelength_angstrom=component.wavelength_angstrom,
                        normalized_radiation_weight=normalized_weight,
                        two_theta_deg=two_theta_deg,
                        intrinsic_line_intensity=(
                            normalized_weight * family_strength
                        ),
                    )
                )
            if not lines:
                continue
            families.append(
                PowderReflectionFamily(
                    family_id=family_id,
                    reflection_ids=member_ids,
                    representative_hkls=tuple(
                        item.representative_hkl for item in members
                    ),
                    d_spacing=d_spacing,
                    multiplicity_crystallographic=family_multiplicity,
                    family_strength=family_strength,
                    lines=tuple(lines),
                )
            )

        ordered_families = tuple(
            sorted(
                families,
                key=lambda item: (item.family_sort_angle, item.family_id),
            )
        )
        return PowderLineSet(
            families=ordered_families,
            structure_factors=structure_factors,
            spectrum=spectrum,
            provenance=PowderLineProvenance(
                method=CALCULATOR_METHOD,
                version=CALCULATOR_VERSION,
                reflections_considered=len(reflections),
                families_emitted=len(ordered_families),
                radiation_components_skipped=skipped_components,
            ),
        )


__all__ = ["PowderLineCalculator"]
