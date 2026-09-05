"""Space-group action on periodic structural rings."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math

import numpy as np

from cristma.diagnostics import Diagnostic, Severity
from cristma.structure import AtomicView, ExpandedAtom, PeriodicAtomRef
from cristma.symmetry import AffineOperation
from cristma.symmetry.orbit import DEFAULT_FRACTIONAL_TOLERANCE

from .representation import StructuralRepresentation
from .ring_finder import _canonical_cycle_key, _connection_joins, _first_unit_image
from .rings import PeriodicUnitRef, StructuralRing, StructuralRingOrbit
from .structural_graph import StructuralConnection


Translation = tuple[int, int, int]


def _add(left: Translation, right: Translation) -> Translation:
    return (left[0] + right[0], left[1] + right[1], left[2] + right[2])


def _subtract(left: Translation, right: Translation) -> Translation:
    return (left[0] - right[0], left[1] - right[1], left[2] - right[2])


def _periodically_equal(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
    tolerance: float,
) -> bool:
    return all(
        abs((a - b + 0.5) % 1.0 - 0.5) <= tolerance + 1e-12
        for a, b in zip(left, right, strict=True)
    )


def map_periodic_atom_ref(
    operation: AffineOperation,
    atom_ref: PeriodicAtomRef,
    view: AtomicView[ExpandedAtom],
    tolerance: float = DEFAULT_FRACTIONAL_TOLERANCE,
) -> PeriodicAtomRef:
    """Map one full periodic atom coordinate through a space-group operation."""

    atom_by_id = {atom.id: atom for atom in view.atoms}
    if atom_ref.atom_id not in atom_by_id:
        raise ValueError("periodic atom reference is absent from the atomic view")
    source = atom_by_id[atom_ref.atom_id]
    full = np.asarray(source.fractional, dtype=float) + np.asarray(
        atom_ref.cell_translation, dtype=float
    )
    rotation = np.asarray(operation.rotation, dtype=float)
    offset = np.asarray(operation.translation, dtype=float)
    raw = rotation @ full + offset
    wrapped = tuple(float(value % 1.0) for value in raw)
    candidates = tuple(
        atom for atom in view.atoms
        if atom.source_site_id == source.source_site_id
        and _periodically_equal(atom.fractional, wrapped, tolerance)
    )
    if len(candidates) != 1:
        raise ValueError("space-group operation has no unique expanded-atom image")
    target = candidates[0]
    image = tuple(
        int(round(float(raw[index]) - target.fractional[index]))
        for index in range(3)
    )
    if any(
        not math.isclose(
            float(raw[index]), target.fractional[index] + image[index], abs_tol=tolerance
        )
        for index in range(3)
    ):
        raise ValueError("mapped atom image is not an integer lattice translation")
    return PeriodicAtomRef(target.id, image)  # type: ignore[arg-type]


def _absolute_unit_atoms(unit, unit_ref: PeriodicUnitRef) -> frozenset[PeriodicAtomRef]:
    return frozenset(
        PeriodicAtomRef(atom_ref.atom_id, _add(atom_ref.cell_translation, unit_ref.cell_translation))
        for atom_ref in unit.atom_refs
    )


def map_periodic_unit_ref(
    operation: AffineOperation,
    unit_ref: PeriodicUnitRef,
    view: AtomicView[ExpandedAtom],
    representation: StructuralRepresentation,
) -> PeriodicUnitRef:
    """Map a unit by transformed atom membership, never by its textual ID."""

    unit_by_id = {unit.unit_id: unit for unit in representation.units}
    if unit_ref.unit_id not in unit_by_id:
        raise ValueError("periodic unit reference is absent from the representation")
    source = unit_by_id[unit_ref.unit_id]
    transformed = frozenset(
        map_periodic_atom_ref(operation, atom_ref, view)
        for atom_ref in _absolute_unit_atoms(source, unit_ref)
    )
    matches: set[PeriodicUnitRef] = set()
    for candidate in representation.units:
        if (
            candidate.kind is not source.kind
            or candidate.interaction_layers != source.interaction_layers
            or candidate.contact_classifications != source.contact_classifications
            or len(candidate.atom_refs) != len(transformed)
        ):
            continue
        for transformed_ref in transformed:
            for candidate_ref in candidate.atom_refs:
                if candidate_ref.atom_id != transformed_ref.atom_id:
                    continue
                shift = _subtract(
                    transformed_ref.cell_translation,
                    candidate_ref.cell_translation,
                )
                candidate_members = frozenset(
                    PeriodicAtomRef(
                        atom_ref.atom_id,
                        _add(atom_ref.cell_translation, shift),
                    )
                    for atom_ref in candidate.atom_refs
                )
                if candidate_members == transformed:
                    matches.add(PeriodicUnitRef(candidate.unit_id, shift))
    if len(matches) != 1:
        raise ValueError("space-group operation has no unique structural-unit image")
    return next(iter(matches))


def _absolute_shared_atoms(
    connection: StructuralConnection,
    first: PeriodicUnitRef,
    second: PeriodicUnitRef,
) -> frozenset[PeriodicAtomRef]:
    first_image = _first_unit_image(connection, first, second)
    return frozenset(
        PeriodicAtomRef(atom_ref.atom_id, _add(atom_ref.cell_translation, first_image))
        for atom_ref in connection.shared_atom_refs
    )


def _mapped_connection_id(
    operation: AffineOperation,
    source_connection: StructuralConnection,
    source_first: PeriodicUnitRef,
    source_second: PeriodicUnitRef,
    target_first: PeriodicUnitRef,
    target_second: PeriodicUnitRef,
    view: AtomicView[ExpandedAtom],
    representation: StructuralRepresentation,
) -> str:
    transformed_shared = frozenset(
        map_periodic_atom_ref(operation, atom_ref, view)
        for atom_ref in _absolute_shared_atoms(
            source_connection, source_first, source_second
        )
    )
    matches = tuple(
        candidate for candidate in representation.connections
        if candidate.connection_kind is source_connection.connection_kind
        and candidate.interaction_layers == source_connection.interaction_layers
        and candidate.contact_classifications == source_connection.contact_classifications
        and _connection_joins(candidate, target_first, target_second)
        and _absolute_shared_atoms(candidate, target_first, target_second)
        == transformed_shared
    )
    if len(matches) != 1:
        raise ValueError("space-group operation has no unique structural-connection image")
    return matches[0].connection_id


def _mapped_ring_key(
    operation: AffineOperation,
    ring: StructuralRing,
    view: AtomicView[ExpandedAtom],
    representation: StructuralRepresentation,
) -> tuple[tuple[str, Translation, str], ...]:
    target_units = tuple(
        map_periodic_unit_ref(operation, unit_ref, view, representation)
        for unit_ref in ring.unit_refs
    )
    connection_by_id = {
        connection.connection_id: connection for connection in representation.connections
    }
    target_connection_ids = tuple(
        _mapped_connection_id(
            operation,
            connection_by_id[connection_id],
            ring.unit_refs[index],
            ring.unit_refs[(index + 1) % ring.size],
            target_units[index],
            target_units[(index + 1) % ring.size],
            view,
            representation,
        )
        for index, connection_id in enumerate(ring.connection_ids)
    )
    return _canonical_cycle_key(target_units, target_connection_ids)


@dataclass(slots=True)
class _DisjointSet:
    parents: dict[str, str]

    def find(self, item: str) -> str:
        parent = self.parents[item]
        if parent != item:
            self.parents[item] = self.find(parent)
        return self.parents[item]

    def union(self, first: str, second: str) -> None:
        first_root = self.find(first)
        second_root = self.find(second)
        if first_root != second_root:
            self.parents[max(first_root, second_root)] = min(first_root, second_root)


def build_ring_orbits(
    structure,
    view: AtomicView[ExpandedAtom],
    representation: StructuralRepresentation,
    rings: tuple[StructuralRing, ...],
) -> tuple[tuple[StructuralRingOrbit, ...], tuple[Diagnostic, ...]]:
    """Group instances only when an actual space-group operation maps them."""

    if structure.space_group is None:
        raise ValueError("ring symmetry grouping requires a space group")
    key_to_ring = {
        _canonical_cycle_key(ring.unit_refs, ring.connection_ids): ring
        for ring in rings
    }
    disjoint = _DisjointSet({ring.ring_id: ring.ring_id for ring in rings})
    diagnostics: list[Diagnostic] = []
    for ring in rings:
        for operation in structure.space_group.operations:
            try:
                mapped_key = _mapped_ring_key(
                    operation, ring, view, representation
                )
            except ValueError as exc:
                diagnostics.append(Diagnostic(
                    Severity.WARNING,
                    "crystal_chemistry.rings.symmetry_mapping_incomplete",
                    f"Could not map {ring.ring_id} by {operation.id or operation.source}: {exc}",
                ))
                continue
            target = key_to_ring.get(mapped_key)
            if (
                target is not None
                and target.parent_block_id == ring.parent_block_id
                and target.representation_id == ring.representation_id
            ):
                disjoint.union(ring.ring_id, target.ring_id)

    groups: dict[str, list[str]] = {}
    for ring in rings:
        groups.setdefault(disjoint.find(ring.ring_id), []).append(ring.ring_id)
    ring_by_id = {ring.ring_id: ring for ring in rings}
    orbits: list[StructuralRingOrbit] = []
    for member_ids in groups.values():
        ordered_ids = tuple(sorted(member_ids))
        representative = ring_by_id[ordered_ids[0]]
        digest = hashlib.sha256("|".join(ordered_ids).encode("utf-8")).hexdigest()[:24]
        orbits.append(StructuralRingOrbit(
            f"ring-orbit:{digest}",
            representative.parent_block_id,
            representative.representation_id,
            representative.ring_id,
            ordered_ids,
            len(ordered_ids),
            representative.composition,
            representative.size,
            representative.scope,
        ))
    return tuple(sorted(orbits, key=lambda item: item.orbit_id)), tuple(diagnostics)


__all__ = [
    "build_ring_orbits",
    "map_periodic_atom_ref",
    "map_periodic_unit_ref",
]
