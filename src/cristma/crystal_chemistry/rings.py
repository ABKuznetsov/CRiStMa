"""Immutable orbit-first results of structural-ring analysis."""
from __future__ import annotations

from dataclasses import dataclass, replace
from enum import StrEnum
from cristma.crystallography import PeriodicSymmetryRelation
from cristma.diagnostics import Diagnostic


class RingAnalysisStatus(StrEnum):
    COMPLETE = "complete"
    INCOMPLETE = "incomplete"


class StructuralRingScope(StrEnum):
    LOCAL = "local"
    FRAMEWORK = "framework"


@dataclass(frozen=True, slots=True)
class RingSearchPolicy:
    maximum_ring_size: int = 12
    maximum_states_per_connection: int = 50_000
    maximum_paths_per_connection: int = 4_096

    def __post_init__(self) -> None:
        values = (self.maximum_ring_size, self.maximum_states_per_connection, self.maximum_paths_per_connection)
        if any(isinstance(x, bool) or not isinstance(x, int) or x <= 0 for x in values):
            raise ValueError("ring-search limits must be positive integers")
        if self.maximum_ring_size < 3:
            raise ValueError("maximum ring size must be at least three")

    def get_config(self) -> dict[str, int]:
        return {
            "maximum_ring_size": self.maximum_ring_size,
            "maximum_states_per_connection": self.maximum_states_per_connection,
            "maximum_paths_per_connection": self.maximum_paths_per_connection,
        }

    def clone(self, **changes: object) -> "RingSearchPolicy":
        return replace(self, **changes)


@dataclass(frozen=True, slots=True, order=True)
class PeriodicUnitOrbitRef:
    unit_orbit_id: str
    periodic_relation: PeriodicSymmetryRelation

    def __post_init__(self) -> None:
        if not self.unit_orbit_id:
            raise ValueError("periodic unit-orbit reference ID must not be empty")


@dataclass(frozen=True, slots=True)
class StructuralRingOrbit:
    ring_orbit_id: str
    parent_block_id: str
    representation_id: str
    unit_orbit_refs: tuple[PeriodicUnitOrbitRef, ...]
    connection_orbit_ids: tuple[str, ...]
    connector_site_refs: tuple[str, ...]
    multiplicity_in_reference_cell: int
    scope: StructuralRingScope = StructuralRingScope.LOCAL
    provenance: tuple[tuple[str, object], ...] = ()

    def __post_init__(self) -> None:
        if not self.ring_orbit_id or not self.parent_block_id or not self.representation_id:
            raise ValueError("structural ring-orbit identities must not be empty")
        if len(self.unit_orbit_refs) < 3 or len(set(self.unit_orbit_refs)) < 3:
            raise ValueError("structural ring requires three distinct periodic unit-orbit images")
        if len(self.connection_orbit_ids) != len(self.unit_orbit_refs):
            raise ValueError("structural ring requires one connection orbit per unit")
        if tuple(sorted(set(self.connector_site_refs))) != self.connector_site_refs:
            raise ValueError("ring connector site references must be unique and sorted")
        if self.multiplicity_in_reference_cell <= 0:
            raise ValueError("ring multiplicity must be positive")

    @property
    def size(self) -> int:
        return len(self.unit_orbit_refs)

@dataclass(frozen=True, slots=True)
class RingAnalysisResult:
    ring_orbits: tuple[StructuralRingOrbit, ...]
    status: RingAnalysisStatus
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()

    def __post_init__(self) -> None:
        ids = tuple(x.ring_orbit_id for x in self.ring_orbits)
        if tuple(sorted(ids)) != ids or len(set(ids)) != len(ids):
            raise ValueError("ring-orbit IDs must be unique and sorted")
        if self.status is RingAnalysisStatus.INCOMPLETE and not self.diagnostics:
            raise ValueError("incomplete ring analysis requires diagnostics")

__all__ = ["PeriodicUnitOrbitRef", "RingAnalysisResult", "RingAnalysisStatus",
           "RingSearchPolicy", "StructuralRingOrbit", "StructuralRingScope"]
