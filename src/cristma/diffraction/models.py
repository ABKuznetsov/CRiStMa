"""Immutable scientific values returned by diffraction calculations."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
import math

from cristma.diagnostics import Diagnostic

from .diagnostics import SEARCH_LIMIT_REACHED


def _nonempty(value: str, name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must not be empty")


def _positive_finite(value: float, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a real number")
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")


def _setting_id(value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError("space-group setting ID must be an integer")
    if not 1 <= value <= 530:
        raise ValueError("space-group setting ID must be between 1 and 530")


@dataclass(frozen=True, slots=True, order=True)
class MillerIndex:
    h: int
    k: int
    l: int

    def __post_init__(self) -> None:
        for value in (self.h, self.k, self.l):
            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError("Miller-index components must be integers")

    def negated(self) -> MillerIndex:
        return MillerIndex(-self.h, -self.k, -self.l)

    @property
    def is_zero(self) -> bool:
        return self.h == self.k == self.l == 0

    def as_tuple(self) -> tuple[int, int, int]:
        return (self.h, self.k, self.l)


class ExtinctionCauseKind(str, Enum):
    CENTERING = "centering"
    SCREW_AXIS = "screw_axis"
    GLIDE_PLANE = "glide_plane"
    COMBINED = "combined"


@dataclass(frozen=True, slots=True)
class PhaseBucketEvidence:
    transformed_hkl: MillerIndex
    operation_ids: tuple[str, ...]
    translation_parts: tuple[tuple[Fraction, Fraction, Fraction], ...]
    exact_phases: tuple[Fraction, ...]
    relative_phases: tuple[Fraction, ...]
    cancels: bool

    def __post_init__(self) -> None:
        count = len(self.operation_ids)
        if count == 0 or len(set(self.operation_ids)) != count:
            raise ValueError("phase bucket operation IDs must be non-empty and unique")
        if not all(isinstance(value, str) and value.strip() for value in self.operation_ids):
            raise ValueError("phase bucket operation IDs must not be empty")
        if not (
            len(self.translation_parts)
            == len(self.exact_phases)
            == len(self.relative_phases)
            == count
        ):
            raise ValueError("phase bucket evidence sequences must have equal lengths")
        for translation in self.translation_parts:
            if len(translation) != 3 or not all(isinstance(value, Fraction) for value in translation):
                raise TypeError("translation parts must contain three Fractions")
        for phase in self.exact_phases + self.relative_phases:
            if not isinstance(phase, Fraction):
                raise TypeError("exact phases must be Fractions")
            if not Fraction(0) <= phase < Fraction(1):
                raise ValueError("exact phases must be normalized into [0, 1)")
        if type(self.cancels) is not bool:
            raise TypeError("bucket cancellation verdict must be bool")
        expected_relative = tuple(
            (phase - self.exact_phases[0]) % 1 for phase in self.exact_phases
        )
        if self.relative_phases != expected_relative:
            raise ValueError("relative phases must be derived from the reference phase")
        if self.cancels != any(phase != 0 for phase in self.relative_phases):
            raise ValueError("bucket verdict must match its exact relative phases")


@dataclass(frozen=True, slots=True)
class ExtinctionCause:
    kind: ExtinctionCauseKind
    operation_ids: tuple[str, ...]
    evidence: tuple[PhaseBucketEvidence, ...]
    condition: str

    def __post_init__(self) -> None:
        if not isinstance(self.kind, ExtinctionCauseKind):
            raise TypeError("extinction cause kind must be ExtinctionCauseKind")
        if not self.operation_ids or len(set(self.operation_ids)) != len(self.operation_ids):
            raise ValueError("extinction cause operation IDs must be non-empty and unique")
        if not all(isinstance(value, str) and value.strip() for value in self.operation_ids):
            raise ValueError("extinction cause operation IDs must not be empty")
        if not self.evidence:
            raise ValueError("extinction cause must retain exact evidence")
        _nonempty(self.condition, "extinction condition")


@dataclass(frozen=True, slots=True)
class ExtinctionResult:
    absent: bool
    causes: tuple[ExtinctionCause, ...]
    evidence: tuple[PhaseBucketEvidence, ...]

    def __post_init__(self) -> None:
        if type(self.absent) is not bool:
            raise TypeError("absent must be bool")
        if self.absent and not self.causes:
            raise ValueError("absent reflection must have an extinction cause")
        if not self.absent and self.causes:
            raise ValueError("allowed reflection must not have extinction causes")
        if not self.evidence:
            raise ValueError("extinction result must retain phase-bucket evidence")


class ReflectionSetStatus(str, Enum):
    COMPLETE = "complete"
    INCOMPLETE = "incomplete"


@dataclass(frozen=True, slots=True)
class ReflectionProvenance:
    space_group_setting_id: int
    method: str

    def __post_init__(self) -> None:
        _setting_id(self.space_group_setting_id)
        _nonempty(self.method, "reflection provenance method")


@dataclass(frozen=True, slots=True)
class ReflectionGenerationProvenance:
    method: str
    version: str
    space_group_setting_id: int
    cell_fingerprint: str
    d_min: float
    reciprocal_convention: str
    boundary_tolerance: float
    metric_compatibility_tolerance: float
    max_candidates: int | None
    integer_points_tested: int
    reflections_within_d_min: int
    orbits_created: int
    status: ReflectionSetStatus

    def __post_init__(self) -> None:
        _nonempty(self.method, "generation method")
        _nonempty(self.version, "generation version")
        _setting_id(self.space_group_setting_id)
        if len(self.cell_fingerprint) != 64 or any(c not in "0123456789abcdef" for c in self.cell_fingerprint):
            raise ValueError("cell fingerprint must be a lowercase SHA-256 digest")
        _positive_finite(self.d_min, "d_min")
        _nonempty(self.reciprocal_convention, "reciprocal convention")
        _positive_finite(self.boundary_tolerance, "boundary tolerance")
        _positive_finite(self.metric_compatibility_tolerance, "metric compatibility tolerance")
        if self.max_candidates is not None and (
            isinstance(self.max_candidates, bool)
            or not isinstance(self.max_candidates, int)
            or self.max_candidates <= 0
        ):
            raise ValueError("max_candidates must be a positive integer or None")
        for value, name in (
            (self.integer_points_tested, "integer_points_tested"),
            (self.reflections_within_d_min, "reflections_within_d_min"),
            (self.orbits_created, "orbits_created"),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a non-negative integer")
        if self.reflections_within_d_min > self.integer_points_tested:
            raise ValueError("accepted reflection candidates cannot exceed tested candidates")
        if not isinstance(self.status, ReflectionSetStatus):
            raise TypeError("generation status must be ReflectionSetStatus")


@dataclass(frozen=True, slots=True)
class Reflection:
    reflection_id: str
    representative_hkl: MillerIndex
    equivalent_hkls: tuple[MillerIndex, ...]
    d_spacing: float
    reciprocal_norm: float
    multiplicity_crystallographic: int
    friedel_mate_id: str | None
    extinction: ExtinctionResult
    provenance: ReflectionProvenance

    def __post_init__(self) -> None:
        _nonempty(self.reflection_id, "reflection ID")
        if not self.equivalent_hkls:
            raise ValueError("reflection orbit must not be empty")
        if self.equivalent_hkls != tuple(sorted(set(self.equivalent_hkls))):
            raise ValueError("reflection orbit must be unique and sorted")
        if any(index.is_zero for index in self.equivalent_hkls):
            raise ValueError("reflection orbit must not contain the zero index")
        if self.representative_hkl != self.equivalent_hkls[0]:
            raise ValueError("representative must be the first orbit index")
        if self.multiplicity_crystallographic != len(self.equivalent_hkls):
            raise ValueError("crystallographic multiplicity must equal orbit length")
        _positive_finite(self.d_spacing, "d spacing")
        _positive_finite(self.reciprocal_norm, "reciprocal norm")
        if not math.isclose(self.d_spacing * self.reciprocal_norm, 1.0, rel_tol=1e-10):
            raise ValueError("d spacing and reciprocal norm must be reciprocal")
        if self.friedel_mate_id is not None:
            _nonempty(self.friedel_mate_id, "Friedel mate ID")
            if self.friedel_mate_id == self.reflection_id:
                raise ValueError("self-Friedel reflection uses None, not its own ID")
        self_friedel = self.representative_hkl.negated() in self.equivalent_hkls
        if (self.friedel_mate_id is None) != self_friedel:
            raise ValueError("Friedel mate ID contradicts the reciprocal orbit")
        if not isinstance(self.extinction, ExtinctionResult):
            raise TypeError("extinction must be ExtinctionResult")
        if self.provenance.space_group_setting_id < 1:
            raise ValueError("reflection provenance must identify a setting")


@dataclass(frozen=True, slots=True)
class ReflectionSet:
    reflections: tuple[Reflection, ...]
    space_group_setting_id: int
    d_min: float
    status: ReflectionSetStatus
    diagnostics: tuple[Diagnostic, ...]
    provenance: ReflectionGenerationProvenance

    def __post_init__(self) -> None:
        _setting_id(self.space_group_setting_id)
        _positive_finite(self.d_min, "d_min")
        if not isinstance(self.status, ReflectionSetStatus):
            raise TypeError("reflection-set status must be ReflectionSetStatus")
        representatives = tuple(item.representative_hkl for item in self.reflections)
        if representatives != tuple(sorted(representatives)):
            raise ValueError("reflections must be sorted by representative index")
        ids = tuple(item.reflection_id for item in self.reflections)
        if len(set(ids)) != len(ids):
            raise ValueError("reflection IDs must be unique within a set")
        by_id = {item.reflection_id: item for item in self.reflections}
        for item in self.reflections:
            if item.provenance.space_group_setting_id != self.space_group_setting_id:
                raise ValueError("reflection belongs to another setting")
            if item.friedel_mate_id is not None:
                mate = by_id.get(item.friedel_mate_id)
                if mate is None or mate.friedel_mate_id != item.reflection_id:
                    raise ValueError("Friedel links must exist and be reciprocal")
        has_limit = any(item.code == SEARCH_LIMIT_REACHED for item in self.diagnostics)
        if self.status is ReflectionSetStatus.INCOMPLETE and not has_limit:
            raise ValueError("incomplete result requires a search-limit diagnostic")
        if self.status is ReflectionSetStatus.COMPLETE and has_limit:
            raise ValueError("complete result must not contain a search-limit diagnostic")
        if self.provenance.status is not self.status:
            raise ValueError("provenance status must match reflection-set status")
        if self.provenance.space_group_setting_id != self.space_group_setting_id:
            raise ValueError("provenance setting must match reflection-set setting")

    @property
    def allowed(self) -> tuple[Reflection, ...]:
        return tuple(item for item in self.reflections if not item.extinction.absent)

    @property
    def systematically_absent(self) -> tuple[Reflection, ...]:
        return tuple(item for item in self.reflections if item.extinction.absent)


__all__ = [
    "ExtinctionCause",
    "ExtinctionCauseKind",
    "ExtinctionResult",
    "MillerIndex",
    "PhaseBucketEvidence",
    "Reflection",
    "ReflectionGenerationProvenance",
    "ReflectionProvenance",
    "ReflectionSet",
    "ReflectionSetStatus",
]
