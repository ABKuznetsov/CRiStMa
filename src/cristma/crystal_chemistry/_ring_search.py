"""Bounded shortest-return traversal over an exact affine quotient graph."""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from cristma.crystallography import identity_relation
from .representation import StructuralRepresentation
from .rings import PeriodicUnitOrbitRef, RingSearchPolicy
from .structural_blocks import StructuralBlock
from .structural_graph import StructuralConnectionOrbit


@dataclass(frozen=True, slots=True)
class _LiftedStep:
    connection_orbit_id: str
    source: PeriodicUnitOrbitRef
    target: PeriodicUnitOrbitRef


@dataclass(frozen=True, slots=True)
class _LiftedPath:
    states: tuple[PeriodicUnitOrbitRef, ...]
    steps: tuple[_LiftedStep, ...]


@dataclass(frozen=True, slots=True)
class _ReturnPathResult:
    paths: tuple[_LiftedPath, ...]
    complete: bool
    limit_name: str | None = None


def _steps_from(state: PeriodicUnitOrbitRef, connections, context) -> tuple[_LiftedStep, ...]:
    steps: list[_LiftedStep] = []
    for edge in connections:
        if edge.first_unit_orbit_id == state.unit_orbit_id:
            target = PeriodicUnitOrbitRef(
                edge.second_unit_orbit_id,
                state.periodic_relation.compose(edge.periodic_relation, context),
            )
            steps.append(_LiftedStep(edge.connection_orbit_id, state, target))
        if edge.second_unit_orbit_id == state.unit_orbit_id:
            target = PeriodicUnitOrbitRef(
                edge.first_unit_orbit_id,
                state.periodic_relation.compose(edge.periodic_relation.inverse(context), context),
            )
            steps.append(_LiftedStep(edge.connection_orbit_id, state, target))
    return tuple(sorted(steps, key=lambda x: (x.connection_orbit_id, x.target)))


def _backtrack(predecessors, source, target, maximum_paths):
    paths: list[_LiftedPath] = []
    truncated = False

    def visit(state, reversed_steps):
        nonlocal truncated
        if len(paths) >= maximum_paths:
            truncated = True
            return
        if state == source:
            steps = tuple(reversed(reversed_steps))
            paths.append(_LiftedPath((source, *(x.target for x in steps)), steps))
            return
        for step in predecessors.get(state, ()):
            visit(step.source, (*reversed_steps, step))
            if truncated:
                return

    visit(target, ())
    paths.sort(key=lambda x: x.states)
    return tuple(paths), truncated


def find_shortest_return_paths(
    representation: StructuralRepresentation,
    block: StructuralBlock,
    removed: StructuralConnectionOrbit,
    policy: RingSearchPolicy,
) -> _ReturnPathResult:
    if block.representation_id != representation.representation_id:
        raise ValueError("structural block belongs to another representation")
    if removed.connection_orbit_id not in block.connection_orbit_ids:
        raise ValueError("removed connection does not belong to structural block")
    by_id = {x.connection_orbit_id: x for x in representation.connection_orbits}
    if not set(block.connection_orbit_ids) <= set(by_id):
        raise ValueError("structural block references an unknown connection orbit")
    connections = tuple(by_id[x] for x in block.connection_orbit_ids)
    context = representation._symmetry_context
    source = PeriodicUnitOrbitRef(removed.first_unit_orbit_id, identity_relation(context))
    target = PeriodicUnitOrbitRef(removed.second_unit_orbit_id, removed.periodic_relation)
    distance = {source: 0}
    predecessors: dict[PeriodicUnitOrbitRef, list[_LiftedStep]] = {}
    queue = deque((source,))
    target_depth = None
    while queue:
        current = queue.popleft()
        depth = distance[current]
        if target_depth is not None and depth >= target_depth:
            continue
        if depth >= policy.maximum_ring_size - 1:
            continue
        for step in _steps_from(current, connections, context):
            if step.connection_orbit_id == removed.connection_orbit_id and (
                (step.source == source and step.target == target)
                or (step.source == target and step.target == source)
            ):
                continue
            candidate_depth = depth + 1
            known = distance.get(step.target)
            if known is None:
                if len(distance) >= policy.maximum_states_per_connection:
                    return _ReturnPathResult((), False, "maximum_states_per_connection")
                distance[step.target] = candidate_depth
                predecessors[step.target] = [step]
                queue.append(step.target)
                if step.target == target:
                    target_depth = candidate_depth
            elif known == candidate_depth:
                predecessors.setdefault(step.target, []).append(step)
    if target_depth is None:
        return _ReturnPathResult((), True)
    paths, truncated = _backtrack(predecessors, source, target, policy.maximum_paths_per_connection)
    return _ReturnPathResult(paths, not truncated, "maximum_paths_per_connection" if truncated else None)


__all__ = ["find_shortest_return_paths"]
