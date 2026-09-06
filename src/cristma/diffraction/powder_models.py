"""Immutable radiation inputs and intrinsic powder-line results."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
from importlib.resources import files
import json
import math
from numbers import Real


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


__all__ = [
    "RadiationComponent",
    "RadiationSpectrum",
    "RadiationSpectrumProvenance",
]
