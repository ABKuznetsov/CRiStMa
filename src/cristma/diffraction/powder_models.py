"""Immutable radiation inputs and intrinsic powder-line results."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
from importlib.resources import files
import json
import math
from numbers import Real

from cristma.diagnostics import Diagnostic

from .models import MillerIndex, ReflectionSetStatus
from .structure_factor_models import StructureFactorSet


def _nonempty(value: str, name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must not be empty")


def _positive_finite(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number")
    normalized = float(value)
    if not math.isfinite(normalized) or normalized <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return normalized


def _nonnegative_finite(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number")
    normalized = float(value)
    if not math.isfinite(normalized) or normalized < 0.0:
        raise ValueError(f"{name} must be finite and non-negative")
    return normalized


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")


@dataclass(frozen=True, slots=True)
class RadiationComponent:
    """One monochromatic component of an X-ray spectrum."""

    component_id: str
    label: str
    wavelength_angstrom: float
    relative_weight: float

    def __post_init__(self) -> None:
        _nonempty(self.component_id, "radiation component ID")
        _nonempty(self.label, "radiation component label")
        object.__setattr__(
            self,
            "wavelength_angstrom",
            _positive_finite(self.wavelength_angstrom, "wavelength"),
        )
        object.__setattr__(
            self,
            "relative_weight",
            _positive_finite(self.relative_weight, "relative weight"),
        )


@dataclass(frozen=True, slots=True)
class RadiationSpectrumProvenance:
    """Origin and conversion metadata for a radiation spectrum."""

    dataset_id: str
    dataset_version: str
    source: str
    energy_source: str | None = None
    radiative_rate_source: str | None = None
    energy_to_wavelength_formula: str | None = None
    hc_value: float | None = None
    hc_units: str | None = None
    xraylib_version: str | None = None
    xraylib_commit: str | None = None
    resource_checksum: str | None = None
    generator: str | None = None

    def __post_init__(self) -> None:
        _nonempty(self.dataset_id, "radiation dataset ID")
        _nonempty(self.dataset_version, "radiation dataset version")
        _nonempty(self.source, "radiation source")
        preset_values = (
            self.energy_source,
            self.radiative_rate_source,
            self.energy_to_wavelength_formula,
            self.hc_value,
            self.hc_units,
            self.xraylib_version,
            self.xraylib_commit,
            self.resource_checksum,
            self.generator,
        )
        present = tuple(value is not None for value in preset_values)
        if any(present) and not all(present):
            raise ValueError("preset provenance must be complete")
        if not any(present):
            return
        for value, name in (
            (self.energy_source, "energy source"),
            (self.radiative_rate_source, "radiative-rate source"),
            (self.energy_to_wavelength_formula, "energy conversion formula"),
            (self.hc_units, "hc units"),
            (self.xraylib_version, "xraylib version"),
            (self.generator, "radiation resource generator"),
        ):
            _nonempty(value, name)  # type: ignore[arg-type]
        _positive_finite(self.hc_value, "hc value")  # type: ignore[arg-type]
        assert self.xraylib_commit is not None
        if len(self.xraylib_commit) != 40 or any(
            character not in "0123456789abcdef" for character in self.xraylib_commit
        ):
            raise ValueError("xraylib commit must be a lowercase Git SHA-1")
        assert self.resource_checksum is not None
        if len(self.resource_checksum) != 64 or any(
            character not in "0123456789abcdef"
            for character in self.resource_checksum
        ):
            raise ValueError("resource checksum must be a lowercase SHA-256")

    @classmethod
    def user_supplied(cls, dataset_id: str) -> "RadiationSpectrumProvenance":
        """Describe an explicit caller-provided spectrum."""

        return cls(
            dataset_id=dataset_id,
            dataset_version="1",
            source="user_supplied",
        )


@dataclass(frozen=True, slots=True)
class RadiationSpectrum:
    """An ordered X-ray spectrum with scale-independent relative weights."""

    components: tuple[RadiationComponent, ...]
    provenance: RadiationSpectrumProvenance

    def __post_init__(self) -> None:
        if not isinstance(self.components, tuple):
            raise TypeError("radiation components must be a tuple")
        if not self.components:
            raise ValueError("radiation spectrum must not be empty")
        if not all(isinstance(item, RadiationComponent) for item in self.components):
            raise TypeError("radiation spectrum must contain RadiationComponent values")
        component_ids = tuple(item.component_id for item in self.components)
        if len(set(component_ids)) != len(component_ids):
            raise ValueError("radiation component IDs must be unique")
        if not isinstance(self.provenance, RadiationSpectrumProvenance):
            raise TypeError("radiation provenance must be RadiationSpectrumProvenance")

    @property
    def normalized_weights(self) -> tuple[float, ...]:
        """Return component weights normalized without changing stored inputs."""

        total = math.fsum(item.relative_weight for item in self.components)
        return tuple(item.relative_weight / total for item in self.components)

    @classmethod
    @lru_cache(maxsize=1)
    def copper_k_alpha(cls) -> "RadiationSpectrum":
        """Load the packaged Cu K-alpha1/K-alpha2 reference spectrum."""

        resource = files("cristma.reference_data").joinpath(
            "resources", "xray", "xray_radiation.json"
        )
        payload = json.loads(resource.read_text(encoding="ascii"))
        metadata = payload["metadata"]
        components = payload["components"]
        checksum = hashlib.sha256(_canonical_bytes(components)).hexdigest()
        if checksum != metadata["resource_checksum"]:
            raise ValueError("X-ray radiation resource checksum mismatch")
        provenance = RadiationSpectrumProvenance(
            dataset_id=metadata["dataset_id"],
            dataset_version=metadata["version"],
            source=metadata["source"],
            energy_source=metadata["energy_source"],
            radiative_rate_source=metadata["radiative_rate_source"],
            energy_to_wavelength_formula=metadata["energy_to_wavelength_formula"],
            hc_value=metadata["hc_value"],
            hc_units=metadata["hc_units"],
            xraylib_version=metadata["xraylib_version"],
            xraylib_commit=metadata["xraylib_commit"],
            resource_checksum=metadata["resource_checksum"],
            generator=metadata["generator"],
        )
        return cls(
            tuple(
                RadiationComponent(
                    component_id=item["component_id"],
                    label=item["label"],
                    wavelength_angstrom=item["wavelength_angstrom"],
                    relative_weight=item["relative_weight"],
                )
                for item in components
            ),
            provenance,
        )


@dataclass(frozen=True, slots=True)
class PowderLine:
    """One radiation-dependent angular line with intrinsic intensity."""

    line_id: str
    family_id: str
    radiation_component_id: str
    wavelength_angstrom: float
    normalized_radiation_weight: float
    two_theta_deg: float
    intrinsic_line_intensity: float

    def __post_init__(self) -> None:
        _nonempty(self.line_id, "powder line ID")
        _nonempty(self.family_id, "powder family ID")
        _nonempty(self.radiation_component_id, "radiation component ID")
        object.__setattr__(
            self,
            "wavelength_angstrom",
            _positive_finite(self.wavelength_angstrom, "wavelength"),
        )
        weight = _positive_finite(
            self.normalized_radiation_weight,
            "normalized radiation weight",
        )
        if weight > 1.0:
            raise ValueError("normalized radiation weight must not exceed one")
        object.__setattr__(self, "normalized_radiation_weight", weight)
        if isinstance(self.two_theta_deg, bool) or not isinstance(
            self.two_theta_deg, Real
        ):
            raise TypeError("two-theta must be a real number")
        angle = float(self.two_theta_deg)
        if not math.isfinite(angle) or not 0.0 <= angle <= 180.0:
            raise ValueError("two-theta must be finite and within [0, 180] degrees")
        object.__setattr__(self, "two_theta_deg", angle)
        object.__setattr__(
            self,
            "intrinsic_line_intensity",
            _nonnegative_finite(
                self.intrinsic_line_intensity,
                "intrinsic line intensity",
            ),
        )


@dataclass(frozen=True, slots=True)
class PowderReflectionFamily:
    """One Friedel-grouped powder family and its radiation lines."""

    family_id: str
    reflection_ids: tuple[str, ...]
    representative_hkls: tuple[MillerIndex, ...]
    d_spacing: float
    multiplicity_crystallographic: int
    family_strength: float
    lines: tuple[PowderLine, ...]

    def __post_init__(self) -> None:
        _nonempty(self.family_id, "powder family ID")
        if not self.reflection_ids:
            raise ValueError("powder family reflection IDs must not be empty")
        if self.reflection_ids != tuple(sorted(set(self.reflection_ids))):
            raise ValueError("powder family reflection IDs must be unique and sorted")
        if len(self.representative_hkls) != len(self.reflection_ids) or not all(
            isinstance(item, MillerIndex) for item in self.representative_hkls
        ):
            raise ValueError("representative Miller indices must match reflection IDs")
        object.__setattr__(
            self,
            "d_spacing",
            _positive_finite(self.d_spacing, "powder family d spacing"),
        )
        if (
            isinstance(self.multiplicity_crystallographic, bool)
            or not isinstance(self.multiplicity_crystallographic, int)
            or self.multiplicity_crystallographic <= 0
        ):
            raise ValueError("powder family multiplicity must be a positive integer")
        object.__setattr__(
            self,
            "family_strength",
            _nonnegative_finite(self.family_strength, "powder family strength"),
        )
        if not self.lines or not all(isinstance(item, PowderLine) for item in self.lines):
            raise ValueError("powder family must contain PowderLine values")
        if any(item.family_id != self.family_id for item in self.lines):
            raise ValueError("powder line belongs to another family")
        line_ids = tuple(item.line_id for item in self.lines)
        component_ids = tuple(item.radiation_component_id for item in self.lines)
        if len(set(line_ids)) != len(line_ids) or len(set(component_ids)) != len(
            component_ids
        ):
            raise ValueError("powder family line and component IDs must be unique")

    @property
    def family_sort_angle(self) -> float:
        return min(item.two_theta_deg for item in self.lines)


@dataclass(frozen=True, slots=True)
class PowderLineProvenance:
    """Deterministic counters for a powder-line calculation."""

    method: str
    version: str
    reflections_considered: int
    families_emitted: int
    radiation_components_skipped: int

    def __post_init__(self) -> None:
        _nonempty(self.method, "powder calculation method")
        _nonempty(self.version, "powder calculation version")
        for value, name in (
            (self.reflections_considered, "reflections_considered"),
            (self.families_emitted, "families_emitted"),
            (self.radiation_components_skipped, "radiation_components_skipped"),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a non-negative integer")


@dataclass(frozen=True, slots=True)
class PowderLineSet:
    """Intrinsic powder families derived from one structure-factor set."""

    families: tuple[PowderReflectionFamily, ...]
    structure_factors: StructureFactorSet
    spectrum: RadiationSpectrum
    provenance: PowderLineProvenance
    diagnostics: tuple[Diagnostic, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.families, tuple) or not all(
            isinstance(item, PowderReflectionFamily) for item in self.families
        ):
            raise TypeError("powder families must contain PowderReflectionFamily values")
        if not isinstance(self.structure_factors, StructureFactorSet):
            raise TypeError("structure_factors must be StructureFactorSet")
        if not isinstance(self.spectrum, RadiationSpectrum):
            raise TypeError("spectrum must be RadiationSpectrum")
        if not isinstance(self.provenance, PowderLineProvenance):
            raise TypeError("provenance must be PowderLineProvenance")
        if not all(isinstance(item, Diagnostic) for item in self.diagnostics):
            raise TypeError("powder diagnostics must contain Diagnostic values")
        family_ids = tuple(item.family_id for item in self.families)
        if len(set(family_ids)) != len(family_ids):
            raise ValueError("powder family IDs must be unique")
        ordering = tuple(
            (item.family_sort_angle, item.family_id) for item in self.families
        )
        if ordering != tuple(sorted(ordering)):
            raise ValueError("powder families must be ordered by family_sort_angle")
        if self.provenance.families_emitted != len(self.families):
            raise ValueError("families_emitted must match the powder family count")
        if self.provenance.reflections_considered != len(
            self.structure_factors.structure_factors
        ):
            raise ValueError("reflections_considered must match StructureFactorSet")
        declared_order = {
            item.component_id: index for index, item in enumerate(self.spectrum.components)
        }
        line_ids: list[str] = []
        for family in self.families:
            try:
                indices = tuple(
                    declared_order[item.radiation_component_id] for item in family.lines
                )
            except KeyError as exc:
                raise ValueError("powder line uses an unknown radiation component") from exc
            if indices != tuple(sorted(indices)):
                raise ValueError("powder lines must follow spectrum component order")
            line_ids.extend(item.line_id for item in family.lines)
        if len(set(line_ids)) != len(line_ids):
            raise ValueError("powder line IDs must be unique within a set")

    @property
    def status(self) -> ReflectionSetStatus:
        return self.structure_factors.status

    @property
    def lines(self) -> tuple[PowderLine, ...]:
        return tuple(line for family in self.families for line in family.lines)

    @property
    def lines_by_angle(self) -> tuple[PowderLine, ...]:
        return tuple(sorted(self.lines, key=lambda item: (item.two_theta_deg, item.line_id)))


__all__ = [
    "PowderLine",
    "PowderLineProvenance",
    "PowderLineSet",
    "PowderReflectionFamily",
    "RadiationComponent",
    "RadiationSpectrum",
    "RadiationSpectrumProvenance",
]
