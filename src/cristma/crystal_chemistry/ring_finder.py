"""Structural-ring discovery over existing crystal-chemistry blocks."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import hashlib

from cristma.chemistry import Composition
from cristma.diagnostics import Diagnostic, Severity
from cristma.structure import AtomicView, CrystalStructure, ExpandedAtom, PeriodicAtomRef

from ._ring_search import _LiftedPath, find_shortest_return_paths
from .representation import StructuralRepresentation
from .rings import (
    PeriodicUnitRef,
    RingAnalysisResult,
    RingAnalysisStatus,
    RingSearchPolicy,
    StructuralRing,
    StructuralRingScope,
)
from .structural_blocks import StructuralBlock, StructuralBlockResult
from .structural_graph import StructuralConnection, StructuralConnectionKind


Translation = tuple[int, int, int]
_ZERO_TRANSLATION: Translation = (0, 0, 0)


def _add(left: Translation, right: Translation) -> Translation:
    return (left[0] + right[0], left[1] + right[1], left[2] + right[2])


def _subtract(left: Translation, right: Translation) -> Translation:
    return (left[0] - right[0], left[1] - right[1], left[2] - right[2])


def _eligible(connection: StructuralConnection) -> bool:
    return connection.connection_kind in {
        StructuralConnectionKind.SHARED_VERTEX,
        StructuralConnectionKind.SHARED_EDGE,
        StructuralConnectionKind.SHARED_FACE,
    }


def _normalized_tokens(
    unit_refs: tuple[PeriodicUnitRef, ...],
    connection_ids: tuple[str, ...],
) -> tuple[tuple[str, Translation, str], ...]:
    origin = unit_refs[0].cell_translation
    return tuple(
        (unit_ref.unit_id, _subtract(unit_ref.cell_translation, origin), connection_id)
        for unit_ref, connection_id in zip(unit_refs, connection_ids, strict=True)
    )


def _canonical_cycle_key(
    unit_refs: tuple[PeriodicUnitRef, ...],
    connection_ids: tuple[str, ...],
) -> tuple[tuple[str, Translation, str], ...]:
    """Canonicalize rotation, reversal, and a common periodic translation."""

    count = len(unit_refs)
    candidates: list[tuple[tuple[str, Translation, str], ...]] = []
    for offset in range(count):
        vertices = unit_refs[offset:] + unit_refs[:offset]
        edges = connection_ids[offset:] + connection_ids[:offset]
        candidates.append(_normalized_tokens(vertices, edges))

    reversed_vertices = (unit_refs[0], *reversed(unit_refs[1:]))
    reversed_edges = (connection_ids[-1], *reversed(connection_ids[:-1]))
    for offset in range(count):
        vertices = reversed_vertices[offset:] + reversed_vertices[:offset]
        edges = reversed_edges[offset:] + reversed_edges[:offset]
        candidates.append(_normalized_tokens(vertices, edges))
    return min(candidates)


def _connection_joins(
    connection: StructuralConnection,
    first: PeriodicUnitRef,
    second: PeriodicUnitRef,
) -> bool:
    return (
        connection.first_unit_id == first.unit_id
        and connection.second_unit_id == second.unit_id
        and _add(first.cell_translation, connection.lattice_translation)
        == second.cell_translation
    ) or (
        connection.second_unit_id == first.unit_id
        and connection.first_unit_id == second.unit_id
        and _subtract(first.cell_translation, connection.lattice_translation)
        == second.cell_translation
    )


def _is_chordless(
    unit_refs: tuple[PeriodicUnitRef, ...],
    connections: tuple[StructuralConnection, ...],
) -> bool:
    count = len(unit_refs)
    for first_index in range(count):
        for second_index in range(first_index + 1, count):
            if second_index == first_index + 1 or (
                first_index == 0 and second_index == count - 1
            ):
                continue
            if any(
                _connection_joins(connection, unit_refs[first_index], unit_refs[second_index])
                for connection in connections
            ):
                return False
    return True


def _first_unit_image(
    connection: StructuralConnection,
    first: PeriodicUnitRef,
    second: PeriodicUnitRef,
) -> Translation:
    if (
        connection.first_unit_id == first.unit_id
        and connection.second_unit_id == second.unit_id
        and _add(first.cell_translation, connection.lattice_translation)
        == second.cell_translation
    ):
        return first.cell_translation
    if (
        connection.second_unit_id == first.unit_id
        and connection.first_unit_id == second.unit_id
        and _subtract(first.cell_translation, connection.lattice_translation)
        == second.cell_translation
    ):
        return second.cell_translation
    raise ValueError("ring connection does not join its adjacent periodic units")


def _translated_atom_ref(reference: PeriodicAtomRef, shift: Translation) -> PeriodicAtomRef:
    return PeriodicAtomRef(reference.atom_id, _add(reference.cell_translation, shift))


def _materialize_ring(
    key: tuple[tuple[str, Translation, str], ...],
    block: StructuralBlock,
    representation: StructuralRepresentation,
    atomic_view: AtomicView[ExpandedAtom],
) -> StructuralRing:
    unit_by_id = {unit.unit_id: unit for unit in representation.units}
    connection_by_id = {
        connection.connection_id: connection for connection in representation.connections
    }
    unit_refs = tuple(PeriodicUnitRef(unit_id, translation) for unit_id, translation, _ in key)
    connection_ids = tuple(connection_id for _, _, connection_id in key)

    atom_refs = {
        _translated_atom_ref(atom_ref, unit_ref.cell_translation)
        for unit_ref in unit_refs
        for atom_ref in unit_by_id[unit_ref.unit_id].atom_refs
    }
    atom_by_id = {atom.id: atom for atom in atomic_view.atoms}
    amounts: dict[str, float] = {}
    for atom_ref in atom_refs:
        atom = atom_by_id[atom_ref.atom_id]
        for component in atom.components:
            amount = float(component.occupancy.value)
            if amount > 0:
                symbol = component.species.require_element()
                amounts[symbol] = amounts.get(symbol, 0.0) + amount

    connector_refs: set[PeriodicAtomRef] = set()
    for index, connection_id in enumerate(connection_ids):
        first = unit_refs[index]
        second = unit_refs[(index + 1) % len(unit_refs)]
        connection = connection_by_id[connection_id]
        first_image = _first_unit_image(connection, first, second)
        connector_refs.update(
            _translated_atom_ref(atom_ref, first_image)
            for atom_ref in connection.shared_atom_refs
        )

    digest = hashlib.sha256(repr(key).encode("utf-8")).hexdigest()[:24]
    return StructuralRing(
        ring_id=f"ring:{digest}",
        parent_block_id=block.block_id,
        representation_id=representation.representation_id,
        unit_refs=unit_refs,
        connection_ids=connection_ids,
        connector_atom_refs=tuple(sorted(
            connector_refs, key=lambda item: (item.atom_id, item.cell_translation)
        )),
        composition=Composition.from_mapping(amounts),
        translation_sum=_ZERO_TRANSLATION,
        provenance=(("method", "cristma.ring_finder:1"),),
    )


def _limit_diagnostic(
    connection: StructuralConnection,
    limit_name: str | None,
) -> Diagnostic:
    return Diagnostic(
        Severity.WARNING,
        "crystal_chemistry.rings.search_limit_reached",
        f"Ring search reached {limit_name or 'an unknown limit'} at {connection.connection_id}.",
        recovery="Increase the explicit RingSearchPolicy limit and repeat the analysis.",
    )


def _classify_ring_scopes(
    rings: tuple[StructuralRing, ...],
) -> tuple[StructuralRing, ...]:
    """Separate local rings from larger circuits using their own edge evidence."""

    shortest_by_connection: dict[str, int] = {}
    for ring in rings:
        for connection_id in ring.connection_ids:
            shortest_by_connection[connection_id] = min(
                ring.size,
                shortest_by_connection.get(connection_id, ring.size),
            )
    return tuple(
        replace(
            ring,
            scope=(
                StructuralRingScope.FRAMEWORK
                if any(
                    shortest_by_connection[connection_id] < ring.size
                    for connection_id in ring.connection_ids
                )
                else StructuralRingScope.LOCAL
            ),
        )
        for ring in rings
    )


def _validate_inputs(
    structure: CrystalStructure,
    atomic_view: AtomicView[ExpandedAtom],
    representation: StructuralRepresentation,
    blocks: StructuralBlockResult,
) -> None:
    if blocks.representation_id != representation.representation_id:
        raise ValueError("ring blocks belong to another structural representation")
    if atomic_view.cell is None or not any(atomic_view.periodic):
        raise ValueError("structural-ring analysis requires a periodic atomic view")
    if structure.cell != atomic_view.cell:
        raise ValueError("atomic view belongs to another unit cell")
    unit_ids = {unit.unit_id for unit in representation.units}
    connection_ids = {item.connection_id for item in representation.connections}
    atom_ids = set(atomic_view.ids)
    if any(
        not set(block.unit_ids) <= unit_ids
        or not set(block.connection_ids) <= connection_ids
        for block in blocks.blocks
    ):
        raise ValueError("structural block references unknown representation objects")
    if any(
        atom_ref.atom_id not in atom_ids
        for unit in representation.units
        for atom_ref in unit.atom_refs
    ):
        raise ValueError("structural unit references an atom absent from the atomic view")


@dataclass(frozen=True, slots=True)
class RingFinder:
    """Find locally shortest finite rings without repeating chemistry inference."""

    policy: RingSearchPolicy = field(default_factory=RingSearchPolicy)

    def get_config(self) -> dict[str, int]:
        return self.policy.get_config()

    def clone(self, **changes: object) -> "RingFinder":
        return replace(self, **changes)

    def find_instances(
        self,
        structure: CrystalStructure,
        atomic_view: AtomicView[ExpandedAtom],
        representation: StructuralRepresentation,
        blocks: StructuralBlockResult,
    ) -> RingAnalysisResult:
        _validate_inputs(structure, atomic_view, representation, blocks)
        all_connections = {
            item.connection_id: item for item in representation.connections
        }
        eligible_connections = tuple(
            item for item in representation.connections if _eligible(item)
        )
        eligible_ids = {item.connection_id for item in eligible_connections}
        eligible_representation = replace(
            representation,
            connections=eligible_connections,
        )
        found: dict[tuple[tuple[str, Translation, str], ...], StructuralRing] = {}
        diagnostics: list[Diagnostic] = []
        complete = True

        for block in blocks.blocks:
            eligible_block = replace(
                block,
                connection_ids=tuple(
                    item for item in block.connection_ids if item in eligible_ids
                ),
            )
            block_connections = tuple(
                all_connections[item]
                for item in eligible_block.connection_ids
            )
            for connection_id in eligible_block.connection_ids:
                removed = all_connections[connection_id]
                search = find_shortest_return_paths(
                    eligible_representation,
                    eligible_block,
                    removed,
                    self.policy,
                )
                if not search.complete:
                    complete = False
                    diagnostics.append(_limit_diagnostic(removed, search.limit_name))
                for path in search.paths:
                    unit_refs = path.states
                    connection_ids = tuple(step.connection_id for step in path.steps) + (
                        removed.connection_id,
                    )
                    if len(unit_refs) < 3 or len(set(unit_refs)) != len(unit_refs):
                        continue
                    if not _is_chordless(unit_refs, block_connections):
                        continue
                    key = _canonical_cycle_key(unit_refs, connection_ids)
                    if key not in found:
                        found[key] = _materialize_ring(
                            key,
                            block,
                            eligible_representation,
                            atomic_view,
                        )

        status = (
            RingAnalysisStatus.COMPLETE
            if complete
            else RingAnalysisStatus.INCOMPLETE
        )
        rings = _classify_ring_scopes(
            tuple(found[key] for key in sorted(found))
        )
        return RingAnalysisResult(
            rings,
            (),
            status,
            tuple(diagnostics),
            (("method", "cristma.ring_finder:1"),),
        )

    def find(
        self,
        structure: CrystalStructure,
        atomic_view: AtomicView[ExpandedAtom],
        representation: StructuralRepresentation,
        blocks: StructuralBlockResult,
    ) -> RingAnalysisResult:
        """Find instances and group them by the structure's actual symmetry."""

        from ._ring_symmetry import build_ring_orbits

        instances = self.find_instances(
            structure, atomic_view, representation, blocks
        )
        block_orbit_by_id = {
            block.block_id: orbit.block_orbit_id
            for orbit in blocks.block_orbits
            for block in orbit.blocks
        }
        if block_orbit_by_id:
            instances = replace(
                instances,
                rings=tuple(
                    replace(
                        ring,
                        parent_block_orbit_id=block_orbit_by_id[ring.parent_block_id],
                    )
                    for ring in instances.rings
                ),
            )
        orbits, symmetry_diagnostics = build_ring_orbits(
            structure,
            atomic_view,
            representation,
            instances.rings,
            blocks,
        )
        diagnostics = instances.diagnostics + symmetry_diagnostics
        status = (
            RingAnalysisStatus.INCOMPLETE
            if diagnostics
            else instances.status
        )
        return replace(
            instances,
            orbits=orbits,
            status=status,
            diagnostics=diagnostics,
        )


__all__ = ["RingFinder"]
