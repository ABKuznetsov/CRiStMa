"""Exact periodic connectivity of integer-labelled structural unit graphs."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass

from cristma.diagnostics import Diagnostic

from .representation import StructuralRepresentation
from .structural_graph import StructuralConnection


Translation = tuple[int, int, int]
_ZERO_TRANSLATION: Translation = (0, 0, 0)


@dataclass(frozen=True, slots=True)
class PeriodicComponent:
    """One connected component and its exact translation subgroup rank."""

    component_id: str
    unit_ids: tuple[str, ...]
    connection_ids: tuple[str, ...]
    image_offsets: tuple[tuple[str, Translation], ...]
    closure_translations: tuple[Translation, ...]
    periodic_rank: int
    periodic_generators: tuple[Translation, ...]

    def __post_init__(self) -> None:
        if not self.component_id or not self.unit_ids:
            raise ValueError("periodic component requires an ID and at least one unit")
        if self.periodic_rank not in range(4):
            raise ValueError("periodic rank must lie between zero and three")
        if self.periodic_rank != len(self.periodic_generators):
            raise ValueError("periodic generator count must equal periodic rank")


@dataclass(frozen=True, slots=True)
class PeriodicConnectivityResult:
    """Exact connectivity analysis for one selected representation."""

    representation_id: str
    components: tuple[PeriodicComponent, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()

    def __post_init__(self) -> None:
        if not self.representation_id:
            raise ValueError("connectivity result requires a representation ID")


def _negated(value: Translation) -> Translation:
    return (-value[0], -value[1], -value[2])


def _added(first: Translation, second: Translation) -> Translation:
    return (
        first[0] + second[0],
        first[1] + second[1],
        first[2] + second[2],
    )


def _subtracted(first: Translation, second: Translation) -> Translation:
    return _added(first, _negated(second))


def _canonical_vector(value: Translation) -> Translation:
    for component in value:
        if component < 0:
            return _negated(value)
        if component > 0:
            return value
    return value


def _extended_gcd(first: int, second: int) -> tuple[int, int, int]:
    """Return positive gcd and Bezout coefficients for two integers."""

    old_r, remainder = abs(first), abs(second)
    old_s, coefficient_s = 1, 0
    old_t, coefficient_t = 0, 1
    while remainder:
        quotient = old_r // remainder
        old_r, remainder = remainder, old_r - quotient * remainder
        old_s, coefficient_s = coefficient_s, old_s - quotient * coefficient_s
        old_t, coefficient_t = coefficient_t, old_t - quotient * coefficient_t
    return (
        old_r,
        old_s if first >= 0 else -old_s,
        old_t if second >= 0 else -old_t,
    )


def _integer_lattice_basis(
    closures: tuple[Translation, ...],
) -> tuple[Translation, ...]:
    """Return a row-Hermite basis of the exact generated Z-subgroup."""

    rows = [list(vector) for vector in closures if vector != _ZERO_TRANSLATION]
    pivot_row = 0
    pivots: list[int] = []
    for column in range(3):
        nonzero = next(
            (index for index in range(pivot_row, len(rows)) if rows[index][column]),
            None,
        )
        if nonzero is None:
            continue
        rows[pivot_row], rows[nonzero] = rows[nonzero], rows[pivot_row]
        for index in range(pivot_row + 1, len(rows)):
            if not rows[index][column]:
                continue
            first = rows[pivot_row][:]
            second = rows[index][:]
            gcd, first_coefficient, second_coefficient = _extended_gcd(
                first[column], second[column]
            )
            rows[pivot_row] = [
                first_coefficient * left + second_coefficient * right
                for left, right in zip(first, second, strict=True)
            ]
            rows[index] = [
                -(second[column] // gcd) * left
                + (first[column] // gcd) * right
                for left, right in zip(first, second, strict=True)
            ]
        if rows[pivot_row][column] < 0:
            rows[pivot_row] = [-value for value in rows[pivot_row]]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows) or pivot_row == 3:
            break

    basis = rows[:pivot_row]
    for index in range(len(basis) - 1, -1, -1):
        column = pivots[index]
        pivot = basis[index][column]
        for earlier in range(index):
            quotient = basis[earlier][column] // pivot
            basis[earlier] = [
                value - quotient * pivot_value
                for value, pivot_value in zip(
                    basis[earlier], basis[index], strict=True
                )
            ]
    return tuple(tuple(row) for row in basis)  # type: ignore[return-value]


def _component_connections(
    unit_ids: set[str],
    connections: tuple[StructuralConnection, ...],
) -> tuple[StructuralConnection, ...]:
    return tuple(
        item for item in connections
        if item.first_unit_id in unit_ids and item.second_unit_id in unit_ids
    )


@dataclass(frozen=True, slots=True)
class PeriodicConnectivityAnalyzer:
    """Calculate exact gain-graph closures without using coordinates."""

    def get_config(self) -> dict[str, object]:
        return {}

    def clone(self, **changes: object) -> "PeriodicConnectivityAnalyzer":
        if changes:
            names = ", ".join(sorted(changes))
            raise TypeError(f"unknown PeriodicConnectivityAnalyzer configuration: {names}")
        return self

    def analyze(
        self,
        representation: StructuralRepresentation,
    ) -> PeriodicConnectivityResult:
        unit_ids = {unit.unit_id for unit in representation.units}
        connection_ids = {
            connection.connection_id for connection in representation.connections
        }
        if len(connection_ids) != len(representation.connections):
            raise ValueError("structural connection IDs must be unique")
        if any(
            connection.first_unit_id not in unit_ids
            or connection.second_unit_id not in unit_ids
            for connection in representation.connections
        ):
            raise ValueError("structural connection references an unknown unit")

        adjacency: dict[str, list[tuple[str, Translation, str]]] = {
            unit_id: [] for unit_id in unit_ids
        }
        for connection in representation.connections:
            adjacency[connection.first_unit_id].append((
                connection.second_unit_id,
                connection.lattice_translation,
                connection.connection_id,
            ))
            adjacency[connection.second_unit_id].append((
                connection.first_unit_id,
                _negated(connection.lattice_translation),
                connection.connection_id,
            ))
        for neighbors in adjacency.values():
            neighbors.sort(key=lambda item: (item[2], item[0], item[1]))

        components: list[PeriodicComponent] = []
        remaining = set(unit_ids)
        while remaining:
            root = min(remaining)
            offsets: dict[str, Translation] = {root: _ZERO_TRANSLATION}
            queue = deque((root,))
            while queue:
                current = queue.popleft()
                for neighbor, gain, _ in adjacency[current]:
                    if neighbor in offsets:
                        continue
                    offsets[neighbor] = _added(offsets[current], gain)
                    queue.append(neighbor)

            component_unit_ids = tuple(sorted(offsets))
            component_unit_set = set(component_unit_ids)
            remaining -= component_unit_set
            component_connections = _component_connections(
                component_unit_set,
                representation.connections,
            )
            closures = tuple(sorted({
                _canonical_vector(_subtracted(
                    _added(
                        offsets[connection.first_unit_id],
                        connection.lattice_translation,
                    ),
                    offsets[connection.second_unit_id],
                ))
                for connection in component_connections
                if _subtracted(
                    _added(
                        offsets[connection.first_unit_id],
                        connection.lattice_translation,
                    ),
                    offsets[connection.second_unit_id],
                ) != _ZERO_TRANSLATION
            }))
            generators = _integer_lattice_basis(closures)
            components.append(PeriodicComponent(
                component_id=f"component:{root}",
                unit_ids=component_unit_ids,
                connection_ids=tuple(sorted(
                    item.connection_id for item in component_connections
                )),
                image_offsets=tuple(sorted(offsets.items())),
                closure_translations=closures,
                periodic_rank=len(generators),
                periodic_generators=generators,
            ))

        return PeriodicConnectivityResult(
            representation_id=representation.representation_id,
            components=tuple(sorted(components, key=lambda item: item.component_id)),
            provenance=(("method", "cristma.periodic_connectivity_analyzer:1"),),
        )


__all__ = [
    "PeriodicComponent",
    "PeriodicConnectivityAnalyzer",
    "PeriodicConnectivityResult",
]
