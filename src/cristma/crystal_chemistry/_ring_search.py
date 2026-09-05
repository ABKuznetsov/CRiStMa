"""Bounded shortest-return traversal over a lifted periodic unit graph."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass

from .rings import PeriodicUnitRef, RingSearchPolicy
from .structural_blocks import StructuralBlock
from .structural_graph import StructuralConnection
from .representation import StructuralRepresentation


Translation = tuple[int, int, int]


def _add(left: Translation, right: Translation) -> Translation:
    return (left[0] + right[0], left[1] + right[1], left[2] + right[2])


def _negate(value: Translation) -> Translation:
    return (-value[0], -value[1], -value[2])


@dataclass(frozen=True, slots=True)
class _LiftedStep:
    connection_id: str
    source: PeriodicUnitRef
    target: PeriodicUnitRef
    first_unit_image: Translation


@dataclass(frozen=True, slots=True)
class _LiftedPath:
    states: tuple[PeriodicUnitRef, ...]
    steps: tuple[_LiftedStep, ...]


@dataclass(frozen=True, slots=True)
class _ReturnPathResult:
    paths: tuple[_LiftedPath, ...]
    complete: bool
    limit_name: str | None = None


def _steps_from(
    state: PeriodicUnitRef,
    connections: tuple[StructuralConnection, ...],
) -> tuple[_LiftedStep, ...]:
    steps: list[_LiftedStep] = []
    for connection in connections:
        if connection.first_unit_id == state.unit_id:
            target = PeriodicUnitRef(
                connection.second_unit_id,
                _add(state.cell_translation, connection.lattice_translation),
            )
            steps.append(_LiftedStep(
                connection.connection_id,
                state,
                target,
                state.cell_translation,
            ))
        if connection.second_unit_id == state.unit_id:
            target = PeriodicUnitRef(
                connection.first_unit_id,
                _add(state.cell_translation, _negate(connection.lattice_translation)),
            )
            steps.append(_LiftedStep(
                connection.connection_id,
                state,
                target,
                target.cell_translation,
            ))
    return tuple(sorted(
        steps,
        key=lambda step: (
            step.connection_id,
            step.target.unit_id,
            step.target.cell_translation,
        ),
    ))


def _is_removed_instance(
    step: _LiftedStep,
    removed_id: str,
    source: PeriodicUnitRef,
    target: PeriodicUnitRef,
) -> bool:
    return (
        step.connection_id == removed_id
        and (
            (step.source == source and step.target == target)
            or (step.source == target and step.target == source)
        )
    )


def _backtrack_paths(
    predecessors: dict[PeriodicUnitRef, list[_LiftedStep]],
    source: PeriodicUnitRef,
    target: PeriodicUnitRef,
    maximum_paths: int,
) -> tuple[tuple[_LiftedPath, ...], bool]:
    paths: list[_LiftedPath] = []
    truncated = False

    def visit(state: PeriodicUnitRef, reversed_steps: tuple[_LiftedStep, ...]) -> None:
        nonlocal truncated
        if len(paths) >= maximum_paths:
            truncated = True
            return
        if state == source:
            steps = tuple(reversed(reversed_steps))
            paths.append(_LiftedPath(
                (source, *(step.target for step in steps)),
                steps,
            ))
            return
        for step in predecessors.get(state, ()):
            visit(step.source, (*reversed_steps, step))
            if truncated:
                return

    visit(target, ())
    paths.sort(key=lambda path: tuple(
        (state.unit_id, state.cell_translation) for state in path.states
    ))
    return tuple(paths), truncated


def find_shortest_return_paths(
    representation: StructuralRepresentation,
    block: StructuralBlock,
    removed: StructuralConnection,
    policy: RingSearchPolicy,
) -> _ReturnPathResult:
    """Find every shortest lifted path replacing one exact edge instance."""

    if block.representation_id != representation.representation_id:
        raise ValueError("structural block belongs to another representation")
    if removed.connection_id not in block.connection_ids:
        raise ValueError("removed connection does not belong to structural block")
    connection_by_id = {
        connection.connection_id: connection for connection in representation.connections
    }
    if not set(block.connection_ids) <= set(connection_by_id):
        raise ValueError("structural block references an unknown connection")
    connections = tuple(connection_by_id[item] for item in block.connection_ids)

    source = PeriodicUnitRef(removed.first_unit_id, (0, 0, 0))
    target = PeriodicUnitRef(removed.second_unit_id, removed.lattice_translation)
    distance: dict[PeriodicUnitRef, int] = {source: 0}
    predecessors: dict[PeriodicUnitRef, list[_LiftedStep]] = {}
    queue = deque((source,))
    target_depth: int | None = None
    maximum_return_edges = policy.maximum_ring_size - 1

    while queue:
        current = queue.popleft()
        depth = distance[current]
        if target_depth is not None and depth >= target_depth:
            continue
        if depth >= maximum_return_edges:
            continue
        for step in _steps_from(current, connections):
            if _is_removed_instance(
                step, removed.connection_id, source, target
            ):
                continue
            candidate_depth = depth + 1
            known_depth = distance.get(step.target)
            if known_depth is None:
                if len(distance) >= policy.maximum_states_per_connection:
                    return _ReturnPathResult(
                        (), False, "maximum_states_per_connection"
                    )
                distance[step.target] = candidate_depth
                predecessors[step.target] = [step]
                queue.append(step.target)
                if step.target == target:
                    target_depth = candidate_depth
            elif known_depth == candidate_depth:
                predecessors.setdefault(step.target, []).append(step)

    if target_depth is None:
        return _ReturnPathResult((), True)
    paths, truncated = _backtrack_paths(
        predecessors,
        source,
        target,
        policy.maximum_paths_per_connection,
    )
    if truncated:
        return _ReturnPathResult(
            paths, False, "maximum_paths_per_connection"
        )
    return _ReturnPathResult(paths, True)


__all__ = ["find_shortest_return_paths"]
