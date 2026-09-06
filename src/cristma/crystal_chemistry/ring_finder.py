"""Bounded structural-ring discovery directly on the symmetry quotient graph."""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from cristma.chemistry import InteractionLayer
from cristma.crystallography import PeriodicSymmetryRelation
from cristma.crystallography.symmetry_context import _digest
from cristma.diagnostics import Diagnostic, Severity
from ._ring_search import find_shortest_return_paths
from .representation import StructuralRepresentation
from .rings import (
    PeriodicUnitOrbitRef,
    RingAnalysisResult,
    RingAnalysisStatus,
    RingSearchPolicy,
    StructuralRingOrbit,
    StructuralRingScope,
)
from .structural_blocks import StructuralBlock, StructuralBlockResult
from .structural_graph import StructuralConnectionKind, StructuralConnectionOrbit


def _eligible(edge: StructuralConnectionOrbit) -> bool:
    return edge.relation_type in {
        StructuralConnectionKind.SHARED_VERTEX,
        StructuralConnectionKind.SHARED_EDGE,
        StructuralConnectionKind.SHARED_FACE,
    } and bool(set(edge.interaction_layers) & {
        InteractionLayer.STRUCTURAL,
        InteractionLayer.INTRA_SUBSYSTEM,
        InteractionLayer.INTRAMOLECULAR,
    })


def _has_cycle_candidate(edges) -> bool:
    parents: dict[str, str] = {}

    def find(item):
        parents.setdefault(item, item)
        if parents[item] != item:
            parents[item] = find(parents[item])
        return parents[item]

    for edge in edges:
        first, second = find(edge.first_unit_orbit_id), find(edge.second_unit_orbit_id)
        if first == second:
            return True
        parents[max(first, second)] = min(first, second)
    return False


def _joins(edge, first, second, context) -> bool:
    return (
        edge.first_unit_orbit_id == first.unit_orbit_id
        and edge.second_unit_orbit_id == second.unit_orbit_id
        and first.periodic_relation.compose(edge.periodic_relation, context) == second.periodic_relation
    ) or (
        edge.second_unit_orbit_id == first.unit_orbit_id
        and edge.first_unit_orbit_id == second.unit_orbit_id
        and first.periodic_relation.compose(edge.periodic_relation.inverse(context), context) == second.periodic_relation
    )


def _is_chordless(refs, edges, context) -> bool:
    count = len(refs)
    for first in range(count):
        for second in range(first + 1, count):
            if second == first + 1 or (first == 0 and second == count - 1):
                continue
            if any(_joins(edge, refs[first], refs[second], context) for edge in edges):
                return False
    return True


def _action_relation(operation_key, context):
    return PeriodicSymmetryRelation(operation_key, (0, 0, 0))


def _lattice_normalized_tokens(refs, edge_ids, context):
    shift = tuple(-x for x in refs[0].periodic_relation.lattice_translation)
    translation = PeriodicSymmetryRelation(context.identity_operation_key, shift)
    normalized = tuple(
        PeriodicUnitOrbitRef(ref.unit_orbit_id, translation.compose(ref.periodic_relation, context))
        for ref in refs
    )
    return tuple(
        (ref.unit_orbit_id, ref.periodic_relation.operation_key,
         ref.periodic_relation.lattice_translation, edge_id)
        for ref, edge_id in zip(normalized, edge_ids, strict=True)
    )


def _orientation_keys(refs, edge_ids, context):
    count = len(refs)
    keys = []
    for direction in (1, -1):
        for start in range(count):
            indices = tuple((start + direction * step) % count for step in range(count))
            ordered_refs = tuple(refs[index] for index in indices)
            ordered_edges = tuple(
                edge_ids[index if direction == 1 else (index - 1) % count]
                for index in indices
            )
            keys.append(_lattice_normalized_tokens(ordered_refs, ordered_edges, context))
    return tuple(keys)


def _canonical_cycle_key(refs, edge_ids, context, *, quotient_symmetry: bool):
    actions = context.operation_keys if quotient_symmetry else (context.identity_operation_key,)
    candidates = []
    for key in actions:
        action = _action_relation(key, context)
        transformed = tuple(
            PeriodicUnitOrbitRef(ref.unit_orbit_id, action.compose(ref.periodic_relation, context))
            for ref in refs
        )
        candidates.extend(_orientation_keys(transformed, edge_ids, context))
    return min(candidates)


def _refs_from_key(key):
    return tuple(
        PeriodicUnitOrbitRef(unit_id, PeriodicSymmetryRelation(operation_key, translation))
        for unit_id, operation_key, translation, _ in key
    )


def _multiplicity(refs, edge_ids, representation):
    context = representation._symmetry_context
    images = set()
    for key in context.operation_keys:
        action = _action_relation(key, context)
        transformed = tuple(
            PeriodicUnitOrbitRef(ref.unit_orbit_id, action.compose(ref.periodic_relation, context))
            for ref in refs
        )
        images.add(min(_orientation_keys(transformed, edge_ids, context)))
    return len(images)


def _scope(rings):
    shortest: dict[str, int] = {}
    for ring in rings:
        for edge_id in ring.connection_orbit_ids:
            shortest[edge_id] = min(shortest.get(edge_id, ring.size), ring.size)
    return tuple(
        replace(ring, scope=(
            StructuralRingScope.FRAMEWORK
            if any(shortest[x] < ring.size for x in ring.connection_orbit_ids)
            else StructuralRingScope.LOCAL
        ))
        for ring in rings
    )


def _limit(edge, name):
    return Diagnostic(
        Severity.WARNING,
        "crystal_chemistry.rings.search_limit_reached",
        f"Ring search reached {name or 'an unknown limit'} at {edge.connection_orbit_id}.",
        recovery="Increase the explicit RingSearchPolicy limit and repeat the analysis.",
    )


@dataclass(frozen=True, slots=True)
class RingFinder:
    policy: RingSearchPolicy = field(default_factory=RingSearchPolicy)

    def get_config(self) -> dict[str, int]:
        return self.policy.get_config()

    def clone(self, **changes: object) -> "RingFinder":
        return replace(self, **changes)

    def find(self, representation: StructuralRepresentation,
             blocks: StructuralBlockResult) -> RingAnalysisResult:
        if blocks.representation_id != representation.representation_id:
            raise ValueError("ring blocks belong to another representation")
        all_edges = {x.connection_orbit_id: x for x in representation.connection_orbits}
        eligible_ids = {key for key, edge in all_edges.items() if _eligible(edge)}
        found: dict[tuple, tuple[StructuralBlock, tuple]] = {}
        diagnostics: list[Diagnostic] = []
        context = representation._symmetry_context
        for block in blocks.blocks:
            edge_ids = tuple(x for x in block.connection_orbit_ids if x in eligible_ids)
            edges = tuple(all_edges[x] for x in edge_ids)
            if not _has_cycle_candidate(edges):
                continue
            eligible_block = replace(block, connection_orbit_ids=edge_ids)
            eligible_representation = replace(representation, connection_orbits=edges)
            for edge_id in edge_ids:
                removed = all_edges[edge_id]
                search = find_shortest_return_paths(
                    eligible_representation, eligible_block, removed, self.policy
                )
                if not search.complete:
                    diagnostics.append(_limit(removed, search.limit_name))
                for path in search.paths:
                    refs = path.states
                    cycle_edges = tuple(x.connection_orbit_id for x in path.steps) + (edge_id,)
                    if len(refs) < 3 or len(set(refs)) != len(refs) or not _is_chordless(refs, edges, context):
                        continue
                    key = _canonical_cycle_key(refs, cycle_edges, context, quotient_symmetry=True)
                    found.setdefault(key, (block, cycle_edges))
        rings = []
        for key, (block, _) in sorted(found.items()):
            refs = _refs_from_key(key)
            edge_ids = tuple(item[3] for item in key)
            connector_refs = tuple(sorted({
                site
                for edge_id in edge_ids
                for site in all_edges[edge_id].connector_site_refs
            }))
            ring_id = "structural-ring-orbit:" + _digest({
                "representation_id": representation.representation_id,
                "parent_block_id": block.block_id,
                "cycle": key,
            })
            rings.append(StructuralRingOrbit(
                ring_id, block.block_id, representation.representation_id,
                refs, edge_ids, connector_refs,
                _multiplicity(refs, edge_ids, representation),
                provenance=(("method", "cristma.ring_finder:3"),),
            ))
        ordered = tuple(sorted(_scope(tuple(rings)), key=lambda x: x.ring_orbit_id))
        return RingAnalysisResult(
            ordered,
            RingAnalysisStatus.INCOMPLETE if diagnostics else RingAnalysisStatus.COMPLETE,
            tuple(dict.fromkeys(diagnostics)),
            (("method", "cristma.ring_finder:3"), ("policy", self.get_config())),
        )


__all__ = ["RingFinder"]
