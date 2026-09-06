"""Exact periodic connectivity of affine-relation-labelled quotient graphs."""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from cristma.crystallography import PeriodicSymmetryRelation, identity_relation
from cristma.diagnostics import Diagnostic
from .representation import StructuralRepresentation
from .structural_graph import StructuralConnectionOrbit

Translation = tuple[int, int, int]
ZERO: Translation = (0, 0, 0)


@dataclass(frozen=True, slots=True)
class PeriodicComponent:
    component_id: str
    unit_orbit_ids: tuple[str, ...]
    connection_orbit_ids: tuple[str, ...]
    image_relations: tuple[tuple[str, PeriodicSymmetryRelation], ...]
    closure_relations: tuple[PeriodicSymmetryRelation, ...]
    rank: int
    periodic_generators: tuple[Translation, ...]

    def __post_init__(self) -> None:
        if not self.component_id or not self.unit_orbit_ids:
            raise ValueError("periodic component requires identity and unit orbits")
        if self.rank not in range(4) or self.rank != len(self.periodic_generators):
            raise ValueError("periodic rank and exact generator count disagree")

    @property
    def unit_ids(self):
        return self.unit_orbit_ids

    @property
    def connection_ids(self):
        return self.connection_orbit_ids

    @property
    def periodic_rank(self):
        return self.rank


@dataclass(frozen=True, slots=True)
class PeriodicConnectivityResult:
    representation_id: str
    components: tuple[PeriodicComponent, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()

    def __post_init__(self) -> None:
        if not self.representation_id:
            raise ValueError("connectivity result requires a representation ID")


def _negative(vector: Translation) -> Translation:
    return tuple(-x for x in vector)  # type: ignore[return-value]


def _canonical_vector(vector: Translation) -> Translation:
    for value in vector:
        if value < 0:
            return _negative(vector)
        if value > 0:
            return vector
    return vector


def _extended_gcd(first: int, second: int) -> tuple[int, int, int]:
    old_r, remainder = abs(first), abs(second)
    old_s, coefficient_s, old_t, coefficient_t = 1, 0, 0, 1
    while remainder:
        quotient = old_r // remainder
        old_r, remainder = remainder, old_r - quotient * remainder
        old_s, coefficient_s = coefficient_s, old_s - quotient * coefficient_s
        old_t, coefficient_t = coefficient_t, old_t - quotient * coefficient_t
    return old_r, old_s if first >= 0 else -old_s, old_t if second >= 0 else -old_t


def integer_translation_lattice_basis(closures: tuple[Translation, ...]) -> tuple[Translation, ...]:
    """Canonical exact basis of the generated integer translation subgroup."""
    rows = [list(x) for x in closures if x != ZERO]
    pivot_row = 0
    pivots: list[int] = []
    for column in range(3):
        nonzero = next((i for i in range(pivot_row, len(rows)) if rows[i][column]), None)
        if nonzero is None:
            continue
        rows[pivot_row], rows[nonzero] = rows[nonzero], rows[pivot_row]
        for index in range(pivot_row + 1, len(rows)):
            if not rows[index][column]:
                continue
            first, second = rows[pivot_row][:], rows[index][:]
            gcd, a, b = _extended_gcd(first[column], second[column])
            rows[pivot_row] = [a*x + b*y for x, y in zip(first, second, strict=True)]
            rows[index] = [-(second[column] // gcd)*x + (first[column] // gcd)*y for x, y in zip(first, second, strict=True)]
        if rows[pivot_row][column] < 0:
            rows[pivot_row] = [-x for x in rows[pivot_row]]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows) or pivot_row == 3:
            break
    basis = rows[:pivot_row]
    for index in range(len(basis) - 1, -1, -1):
        column, pivot = pivots[index], basis[index][pivots[index]]
        for earlier in range(index):
            quotient = basis[earlier][column] // pivot
            basis[earlier] = [x - quotient*y for x, y in zip(basis[earlier], basis[index], strict=True)]
    return tuple(tuple(x) for x in basis)  # type: ignore[return-value]


def _translation_kernel(relations, context) -> tuple[Translation, ...]:
    """Schreier generators of the pure-translation kernel of an affine subgroup."""
    if not relations:
        return ()
    generators = tuple(sorted(set((*relations, *(x.inverse(context) for x in relations)))))
    identity = identity_relation(context)
    representatives = {identity.operation_key: identity}
    queue = deque((identity.operation_key,))
    translations: set[Translation] = set()
    while queue:
        key = queue.popleft()
        representative = representatives[key]
        for generator in generators:
            candidate = representative.compose(generator, context)
            existing = representatives.get(candidate.operation_key)
            if existing is None:
                representatives[candidate.operation_key] = candidate
                queue.append(candidate.operation_key)
                continue
            delta = candidate.compose(existing.inverse(context), context)
            if delta.operation_key != context.identity_operation_key:
                raise ValueError("affine subgroup reduction did not produce a translation")
            if delta.lattice_translation != ZERO:
                translations.add(_canonical_vector(delta.lattice_translation))
    return tuple(sorted(translations))


def _component_connections(unit_ids: set[str], connections: tuple[StructuralConnectionOrbit, ...]):
    return tuple(x for x in connections if x.first_unit_orbit_id in unit_ids and x.second_unit_orbit_id in unit_ids)


@dataclass(frozen=True, slots=True)
class PeriodicConnectivityAnalyzer:
    def get_config(self) -> dict[str, object]:
        return {}

    def clone(self, **changes: object) -> "PeriodicConnectivityAnalyzer":
        if changes:
            raise TypeError("unknown PeriodicConnectivityAnalyzer configuration: " + ", ".join(sorted(changes)))
        return self

    def analyze(self, representation: StructuralRepresentation) -> PeriodicConnectivityResult:
        context = representation._symmetry_context
        unit_ids = {x.unit_orbit_id for x in representation.unit_orbits}
        adjacency: dict[str, list[tuple[str, PeriodicSymmetryRelation, str]]] = {x: [] for x in unit_ids}
        for edge in representation.connection_orbits:
            adjacency[edge.first_unit_orbit_id].append((edge.second_unit_orbit_id, edge.periodic_relation, edge.connection_orbit_id))
            adjacency[edge.second_unit_orbit_id].append((edge.first_unit_orbit_id, edge.periodic_relation.inverse(context), edge.connection_orbit_id))
        for values in adjacency.values():
            values.sort(key=lambda x: (x[2], x[0], x[1]))
        components: list[PeriodicComponent] = []
        remaining = set(unit_ids)
        while remaining:
            root = min(remaining)
            images = {root: identity_relation(context)}
            queue = deque((root,))
            closures: set[PeriodicSymmetryRelation] = set()
            while queue:
                current = queue.popleft()
                for neighbor, gain, _ in adjacency[current]:
                    proposed = images[current].compose(gain, context)
                    known = images.get(neighbor)
                    if known is None:
                        images[neighbor] = proposed
                        queue.append(neighbor)
                    else:
                        closure = proposed.compose(known.inverse(context), context)
                        if closure != identity_relation(context):
                            closures.add(closure)
            component_units = tuple(sorted(images))
            remaining -= set(component_units)
            edges = _component_connections(set(component_units), representation.connection_orbits)
            translations = _translation_kernel(tuple(sorted(closures)), context)
            basis = integer_translation_lattice_basis(translations)
            components.append(PeriodicComponent(
                "component:" + root, component_units, tuple(sorted(x.connection_orbit_id for x in edges)),
                tuple(sorted(images.items())), tuple(sorted(closures)), len(basis), basis,
            ))
        return PeriodicConnectivityResult(representation.representation_id,
                                          tuple(sorted(components, key=lambda x: x.component_id)),
                                          provenance=(("method", "cristma.periodic_connectivity_analyzer:2"),))


__all__ = ["PeriodicComponent", "PeriodicConnectivityAnalyzer", "PeriodicConnectivityResult",
           "integer_translation_lattice_basis"]
