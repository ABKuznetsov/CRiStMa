"""Immutable results of structural-ring analysis."""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import StrEnum

from cristma.chemistry import Composition
from cristma.diagnostics import Diagnostic
from cristma.structure import PeriodicAtomRef


Translation = tuple[int, int, int]
_ZERO_TRANSLATION: Translation = (0, 0, 0)


def _validate_translation(value: Translation, *, name: str) -> None:
    if len(value) != 3 or any(
        isinstance(item, bool) or not isinstance(item, int) for item in value
    ):
        raise ValueError(f"{name} must contain three integers")


@dataclass(frozen=True, slots=True, order=True)
class PeriodicUnitRef:
    """Reference to one lattice-translated image of a structural unit."""

    unit_id: str
    cell_translation: Translation

    def __post_init__(self) -> None:
        if not self.unit_id:
            raise ValueError("periodic unit reference ID must not be empty")
        _validate_translation(self.cell_translation, name="cell translation")


class RingAnalysisStatus(StrEnum):
    """Whether every eligible connection was searched without truncation."""

    COMPLETE = "complete"
    INCOMPLETE = "incomplete"


class StructuralRingScope(StrEnum):
    """Topological scale of a ring inside its parent structural block."""

    LOCAL = "local"
    FRAMEWORK = "framework"


@dataclass(frozen=True, slots=True)
class RingSearchPolicy:
    """Explicit deterministic safety limits for periodic ring traversal."""

    maximum_ring_size: int = 12
    maximum_states_per_connection: int = 50_000
    maximum_paths_per_connection: int = 4_096

    def __post_init__(self) -> None:
        if any(
            isinstance(value, bool) or not isinstance(value, int) or value <= 0
            for value in (
                self.maximum_ring_size,
                self.maximum_states_per_connection,
                self.maximum_paths_per_connection,
            )
        ):
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


@dataclass(frozen=True, slots=True)
class StructuralRing:
    """One finite locally-shortest cycle in a structural block."""

    ring_id: str
    parent_block_id: str
    representation_id: str
    unit_refs: tuple[PeriodicUnitRef, ...]
    connection_ids: tuple[str, ...]
    connector_atom_refs: tuple[PeriodicAtomRef, ...]
    composition: Composition
    translation_sum: Translation
    scope: StructuralRingScope = StructuralRingScope.LOCAL
    provenance: tuple[tuple[str, object], ...] = ()
    parent_block_orbit_id: str = ""

    def __post_init__(self) -> None:
        if not self.ring_id or not self.parent_block_id or not self.representation_id:
            raise ValueError("structural ring identities must not be empty")
        if len(self.unit_refs) < 3 or len(set(self.unit_refs)) < 3:
            raise ValueError("structural ring requires three distinct periodic units")
        if len(self.connection_ids) != len(self.unit_refs):
            raise ValueError("structural ring must have one connection per unit")
        if len(set(self.connector_atom_refs)) != len(self.connector_atom_refs):
            raise ValueError("structural ring connector atoms must be unique")
        _validate_translation(self.translation_sum, name="ring translation sum")
        if self.translation_sum != _ZERO_TRANSLATION:
            raise ValueError("structural ring must have zero translation sum")
        if self.parent_block_orbit_id and not self.parent_block_orbit_id.strip():
            raise ValueError("parent block orbit ID must not be blank")

    @property
    def size(self) -> int:
        return len(self.unit_refs)


@dataclass(frozen=True, slots=True)
class StructuralRingOrbit:
    """Crystallographically equivalent ring instances in one context."""

    orbit_id: str
    parent_block_id: str
    representation_id: str
    representative_ring_id: str
    ring_ids: tuple[str, ...]
    multiplicity: int
    composition: Composition
    size: int
    scope: StructuralRingScope = StructuralRingScope.LOCAL
    parent_block_orbit_id: str = ""

    def __post_init__(self) -> None:
        if not self.orbit_id or not self.parent_block_id or not self.representation_id:
            raise ValueError("ring-orbit identities must not be empty")
        if not self.representative_ring_id or not self.ring_ids:
            raise ValueError("ring orbit requires a representative and members")
        if len(set(self.ring_ids)) != len(self.ring_ids):
            raise ValueError("ring-orbit member IDs must be unique")
        if self.representative_ring_id not in self.ring_ids:
            raise ValueError("ring-orbit representative must be a member")
        if self.multiplicity != len(self.ring_ids):
            raise ValueError("ring-orbit multiplicity must equal member count")
        if self.size < 3:
            raise ValueError("ring-orbit size must be at least three")
        if self.parent_block_orbit_id and not self.parent_block_orbit_id.strip():
            raise ValueError("parent block orbit ID must not be blank")


@dataclass(frozen=True, slots=True)
class RingAnalysisResult:
    """All ring instances and symmetry orbits for one representation."""

    rings: tuple[StructuralRing, ...]
    orbits: tuple[StructuralRingOrbit, ...]
    status: RingAnalysisStatus
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()

    def __post_init__(self) -> None:
        ring_by_id = {ring.ring_id: ring for ring in self.rings}
        if len(ring_by_id) != len(self.rings):
            raise ValueError("ring IDs must be unique")
        if self.status is RingAnalysisStatus.INCOMPLETE and not self.diagnostics:
            raise ValueError("incomplete ring analysis requires a diagnostic")
        assigned: set[str] = set()
        for orbit in self.orbits:
            if any(ring_id not in ring_by_id for ring_id in orbit.ring_ids):
                raise ValueError("ring orbit references an unknown ring")
            members = tuple(ring_by_id[ring_id] for ring_id in orbit.ring_ids)
            representative = ring_by_id[orbit.representative_ring_id]
            if (
                representative.parent_block_id != orbit.parent_block_id
                or representative.parent_block_orbit_id
                != orbit.parent_block_orbit_id
            ):
                raise ValueError("ring-orbit parent identity must follow its representative")
            if any(
                ring.representation_id != orbit.representation_id
                or ring.size != orbit.size
                or ring.composition != orbit.composition
                or ring.scope is not orbit.scope
                or (
                    orbit.parent_block_orbit_id
                    and ring.parent_block_orbit_id != orbit.parent_block_orbit_id
                )
                for ring in members
            ):
                raise ValueError("ring-orbit members must share context and identity")
            if assigned.intersection(orbit.ring_ids):
                raise ValueError("ring instance occurs in more than one orbit")
            assigned.update(orbit.ring_ids)
        if self.orbits and assigned != set(ring_by_id):
            raise ValueError("every structural ring must belong to exactly one orbit")


__all__ = [
    "PeriodicUnitRef",
    "RingAnalysisResult",
    "RingAnalysisStatus",
    "RingSearchPolicy",
    "StructuralRing",
    "StructuralRingOrbit",
    "StructuralRingScope",
]
