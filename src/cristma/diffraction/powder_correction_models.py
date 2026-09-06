"""Immutable physical correction inputs and corrected powder lines."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from numbers import Real

from .models import ReflectionSetStatus
from .powder_models import PowderLineSet, RadiationProbe


def _finite_nonnegative(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number")
    normalized = float(value)
    if not math.isfinite(normalized) or normalized < 0.0:
        raise ValueError(f"{name} must be finite and non-negative")
    return normalized


@dataclass(frozen=True, slots=True)
class BraggBrentanoGeometry:
    """Symmetric reflection geometry with explicit X-ray polarization."""

    perpendicular_polarization_fraction: float = 0.5
    geometry_id: str = field(default="bragg_brentano_symmetric_reflection", init=False)

    def __post_init__(self) -> None:
        fraction = _finite_nonnegative(
            self.perpendicular_polarization_fraction,
            "perpendicular polarization fraction",
        )
        if fraction > 1.0:
            raise ValueError("perpendicular polarization fraction must not exceed one")
        object.__setattr__(self, "perpendicular_polarization_fraction", fraction)


@dataclass(frozen=True, slots=True)
class CorrectedPowderLine:
    """One powder line after an explicit geometrical intensity correction."""

    powder_line_id: str
    two_theta_deg: float
    intrinsic_line_intensity: float
    lorentz_factor: float
    polarization_factor: float
    corrected_line_intensity: float

    def __post_init__(self) -> None:
        if not isinstance(self.powder_line_id, str) or not self.powder_line_id.strip():
            raise ValueError("powder line ID must not be empty")
        normalized = (
            _finite_nonnegative(self.two_theta_deg, "two-theta"),
            _finite_nonnegative(
                self.intrinsic_line_intensity,
                "intrinsic line intensity",
            ),
            _finite_nonnegative(self.lorentz_factor, "Lorentz factor"),
            _finite_nonnegative(self.polarization_factor, "polarization factor"),
            _finite_nonnegative(
                self.corrected_line_intensity,
                "corrected line intensity",
            ),
        )
        for field_name, value in zip(
            (
                "two_theta_deg",
                "intrinsic_line_intensity",
                "lorentz_factor",
                "polarization_factor",
                "corrected_line_intensity",
            ),
            normalized,
            strict=True,
        ):
            object.__setattr__(self, field_name, value)

    @property
    def lorentz_polarization_factor(self) -> float:
        return self.lorentz_factor * self.polarization_factor


@dataclass(frozen=True, slots=True)
class PowderCorrectionProvenance:
    """Exact convention used for geometrical powder corrections."""

    method: str
    version: str
    geometry_id: str
    radiation_probe: RadiationProbe
    radiation_source_id: str
    lorentz_formula: str
    polarization_formula: str
    lines_corrected: int

    def __post_init__(self) -> None:
        for value, name in (
            (self.method, "method"),
            (self.version, "version"),
            (self.geometry_id, "geometry ID"),
            (self.radiation_source_id, "radiation source ID"),
            (self.lorentz_formula, "Lorentz formula"),
            (self.polarization_formula, "polarization formula"),
        ):
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must not be empty")
        if not isinstance(self.radiation_probe, RadiationProbe):
            raise TypeError("radiation probe must be RadiationProbe")
        if (
            isinstance(self.lines_corrected, bool)
            or not isinstance(self.lines_corrected, int)
            or self.lines_corrected < 0
        ):
            raise ValueError("lines_corrected must be a non-negative integer")


@dataclass(frozen=True, slots=True)
class CorrectedPowderLineSet:
    """Corrected one-to-one projection of an intrinsic PowderLineSet."""

    lines: tuple[CorrectedPowderLine, ...]
    powder_lines: PowderLineSet
    geometry: BraggBrentanoGeometry
    provenance: PowderCorrectionProvenance

    def __post_init__(self) -> None:
        if not isinstance(self.lines, tuple) or not all(
            isinstance(item, CorrectedPowderLine) for item in self.lines
        ):
            raise TypeError("corrected lines must contain CorrectedPowderLine values")
        if not isinstance(self.powder_lines, PowderLineSet):
            raise TypeError("powder_lines must be PowderLineSet")
        if not isinstance(self.geometry, BraggBrentanoGeometry):
            raise TypeError("geometry must be BraggBrentanoGeometry")
        if not isinstance(self.provenance, PowderCorrectionProvenance):
            raise TypeError("provenance must be PowderCorrectionProvenance")
        expected = tuple(item.line_id for item in self.powder_lines.lines)
        observed = tuple(item.powder_line_id for item in self.lines)
        if observed != expected:
            raise ValueError("corrected lines must exactly match PowderLineSet order")
        for corrected, source in zip(
            self.lines,
            self.powder_lines.lines,
            strict=True,
        ):
            if (
                corrected.two_theta_deg != source.two_theta_deg
                or corrected.intrinsic_line_intensity
                != source.intrinsic_line_intensity
            ):
                raise ValueError(
                    "corrected line must retain its source angle and intrinsic intensity"
                )
        if self.provenance.lines_corrected != len(self.lines):
            raise ValueError("lines_corrected must match the corrected line count")
        if self.provenance.geometry_id != self.geometry.geometry_id:
            raise ValueError("correction provenance belongs to another geometry")
        if self.provenance.radiation_probe is not self.powder_lines.spectrum.probe:
            raise ValueError("correction provenance belongs to another radiation probe")
        if self.provenance.radiation_source_id != self.powder_lines.spectrum.source_id:
            raise ValueError("correction provenance belongs to another radiation source")

    @property
    def status(self) -> ReflectionSetStatus:
        return self.powder_lines.status

    @property
    def lines_by_angle(self) -> tuple[CorrectedPowderLine, ...]:
        return tuple(
            sorted(self.lines, key=lambda item: (item.two_theta_deg, item.powder_line_id))
        )


__all__ = [
    "BraggBrentanoGeometry",
    "CorrectedPowderLine",
    "CorrectedPowderLineSet",
    "PowderCorrectionProvenance",
]
