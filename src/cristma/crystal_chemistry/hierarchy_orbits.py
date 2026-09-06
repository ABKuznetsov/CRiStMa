"""Exact space-group orbits and scientific geometry for hierarchy objects."""

from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from types import SimpleNamespace

import numpy as np

from cristma.diagnostics import Diagnostic, Severity
from cristma.structure import AtomicView, CrystalStructure, ExpandedAtom, PeriodicAtomRef
from cristma.symmetry.orbit import DEFAULT_FRACTIONAL_TOLERANCE

from .periodic_connectivity import PeriodicConnectivityResult
from .polyhedra import CoordinationPolyhedron, _convex_hull_faces
from .representation import StructuralRepresentation
from .rings import PeriodicUnitRef
from .structural_blocks import StructuralBlock, StructuralBlockOrbit
from .structural_units import (
    StructuralUnit,
    StructuralUnitGeometry,
    StructuralUnitGeometryKind,
    StructuralUnitKind,
    StructuralUnitOrbit,
)


Translation = tuple[int, int, int]
_ZERO: Translation = (0, 0, 0)
UNIT_SYMMETRY_INCOMPLETE = "crystal_chemistry.unit_orbit.symmetry_incomplete"
BLOCK_SYMMETRY_INCOMPLETE = "crystal_chemistry.block_orbit.symmetry_incomplete"
UNIT_GEOMETRY_INCOMPLETE = "crystal_chemistry.unit_geometry.incomplete"


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


def _stable_id(prefix: str, values: tuple[str, ...]) -> str:
    payload = json.dumps(values, ensure_ascii=True, separators=(",", ":")).encode(
        "ascii"
    )
    return f"{prefix}:{hashlib.sha256(payload).hexdigest()}"


def build_unit_orbits(
    structure: CrystalStructure,
    view: AtomicView[ExpandedAtom],
    units: tuple[StructuralUnit, ...],
) -> tuple[
    tuple[StructuralUnit, ...],
    tuple[StructuralUnitOrbit, ...],
    tuple[Diagnostic, ...],
]:
    """Group units only when an actual space-group operation maps membership."""

    if not units:
        return (), (), ()
    if structure.space_group is None:
        raise ValueError("structural-unit orbit construction requires a space group")
    from ._ring_symmetry import map_periodic_unit_ref

    unit_ids = {item.unit_id for item in units}
    collection = SimpleNamespace(units=units)
    disjoint = _DisjointSet({item.unit_id: item.unit_id for item in units})
    incomplete: set[str] = set()
    diagnostics: list[Diagnostic] = []
    for unit in units:
        source = PeriodicUnitRef(unit.unit_id, _ZERO)
        for operation in structure.space_group.operations:
            try:
                mapped = map_periodic_unit_ref(operation, source, view, collection)
            except ValueError:
                incomplete.add(unit.unit_id)
                continue
            if mapped.unit_id not in unit_ids:
                incomplete.add(unit.unit_id)
                continue
            disjoint.union(unit.unit_id, mapped.unit_id)

    groups: dict[str, list[StructuralUnit]] = {}
    for unit in units:
        groups.setdefault(disjoint.find(unit.unit_id), []).append(unit)
    replacements: dict[str, StructuralUnit] = {}
    orbits: list[StructuralUnitOrbit] = []
    for rows in groups.values():
        ordered = tuple(sorted(rows, key=lambda item: item.unit_id))
        orbit_id = _stable_id("unit-orbit", tuple(item.unit_id for item in ordered))
        members = tuple(replace(item, unit_orbit_id=orbit_id) for item in ordered)
        orbit_diagnostics: tuple[Diagnostic, ...] = ()
        if any(item.unit_id in incomplete for item in members):
            orbit_diagnostics = (
                Diagnostic(
                    Severity.WARNING,
                    UNIT_SYMMETRY_INCOMPLETE,
                    "A space-group image of the structural unit is missing",
                ),
            )
            diagnostics.extend(orbit_diagnostics)
        for item in members:
            replacements[item.unit_id] = item
        orbits.append(
            StructuralUnitOrbit(
                unit_orbit_id=orbit_id,
                representative_unit_id=members[0].unit_id,
                units=members,
                diagnostics=orbit_diagnostics,
            )
        )
    return (
        tuple(replacements[item.unit_id] for item in units),
        tuple(sorted(orbits, key=lambda item: item.unit_orbit_id)),
        tuple(dict.fromkeys(diagnostics)),
    )


def _integer_rotation(operation) -> np.ndarray:
    rows: list[list[int]] = []
    for row in operation.rotation:
        values: list[int] = []
        for value in row:
            exact = Fraction(value)
            if exact.denominator != 1:
                raise ValueError("space-group rotation is not integral")
            values.append(exact.numerator)
        rows.append(values)
    return np.asarray(rows, dtype=int)


def _transform_translation(rotation: np.ndarray, value: Translation) -> Translation:
    transformed = rotation @ np.asarray(value, dtype=int)
    return tuple(int(item) for item in transformed)  # type: ignore[return-value]


def _lattice_contains(value: Translation, generators: tuple[Translation, ...]) -> bool:
    if not generators:
        return value == _ZERO
    matrix = np.asarray(generators, dtype=int).T
    target = np.asarray(value, dtype=int)
    rank = len(generators)
    for axes in combinations(range(3), rank):
        square = [[Fraction(int(matrix[axis, column])) for column in range(rank)] for axis in axes]
        determinant = round(float(np.linalg.det(np.asarray(square, dtype=float))))
        if determinant == 0:
            continue
        augmented = [row + [Fraction(int(target[axis]))] for row, axis in zip(square, axes, strict=True)]
        for column in range(rank):
            pivot = next((row for row in range(column, rank) if augmented[row][column]), None)
            if pivot is None:
                break
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
            divisor = augmented[column][column]
            augmented[column] = [item / divisor for item in augmented[column]]
            for row in range(rank):
                if row == column:
                    continue
                factor = augmented[row][column]
                augmented[row] = [
                    left - factor * right
                    for left, right in zip(augmented[row], augmented[column], strict=True)
                ]
        coefficients = tuple(augmented[index][-1] for index in range(rank))
        if any(item.denominator != 1 for item in coefficients):
            return False
        reconstructed = matrix @ np.asarray([int(item) for item in coefficients], dtype=int)
        return bool(np.array_equal(reconstructed, target))
    return False


def _same_lattice(
    first: tuple[Translation, ...],
    second: tuple[Translation, ...],
) -> bool:
    return len(first) == len(second) and all(
        _lattice_contains(item, second) for item in first
    ) and all(_lattice_contains(item, first) for item in second)


def _block_unit_refs(
    block: StructuralBlock,
    connectivity: PeriodicConnectivityResult,
) -> tuple[PeriodicUnitRef, ...]:
    component_id = dict(block.provenance).get("component_id")
    component = next(
        (item for item in connectivity.components if item.component_id == component_id),
        None,
    )
    if component is None:
        raise ValueError("structural block lacks its periodic component provenance")
    offsets = dict(component.image_offsets)
    return tuple(
        PeriodicUnitRef(unit_id, offsets[unit_id]) for unit_id in component.unit_ids
    )


def _mapped_block_id(
    operation,
    source: StructuralBlock,
    blocks: tuple[StructuralBlock, ...],
    refs_by_id: dict[str, tuple[PeriodicUnitRef, ...]],
    view: AtomicView[ExpandedAtom],
    representation: StructuralRepresentation,
) -> str | None:
    from ._ring_symmetry import map_periodic_unit_ref

    mapped_refs = tuple(
        map_periodic_unit_ref(operation, item, view, representation)
        for item in refs_by_id[source.block_id]
    )
    rotation = _integer_rotation(operation)
    mapped_generators = tuple(
        _transform_translation(rotation, item) for item in source.periodic_generators
    )
    mapped_by_unit = {item.unit_id: item.cell_translation for item in mapped_refs}
    matches: list[str] = []
    for candidate in blocks:
        if (
            candidate.periodic_rank != source.periodic_rank
            or candidate.classification is not source.classification
            or not _same_lattice(mapped_generators, candidate.periodic_generators)
        ):
            continue
        candidate_refs = refs_by_id[candidate.block_id]
        candidate_by_unit = {
            item.unit_id: item.cell_translation for item in candidate_refs
        }
        if set(candidate_by_unit) != set(mapped_by_unit):
            continue
        anchor = min(mapped_by_unit)
        shift = tuple(
            mapped_by_unit[anchor][index] - candidate_by_unit[anchor][index]
            for index in range(3)
        )
        if all(
            _lattice_contains(
                tuple(
                    mapped_by_unit[unit_id][index]
                    - candidate_by_unit[unit_id][index]
                    - shift[index]
                    for index in range(3)
                ),
                candidate.periodic_generators,
            )
            for unit_id in mapped_by_unit
        ):
            matches.append(candidate.block_id)
    return matches[0] if len(matches) == 1 else None


def build_block_orbits(
    structure: CrystalStructure,
    view: AtomicView[ExpandedAtom],
    representation: StructuralRepresentation,
    connectivity: PeriodicConnectivityResult,
    blocks: tuple[StructuralBlock, ...],
) -> tuple[
    tuple[StructuralBlock, ...],
    tuple[StructuralBlockOrbit, ...],
    tuple[Diagnostic, ...],
]:
    """Group periodic components by exact symmetry action modulo translations."""

    if not blocks:
        return (), (), ()
    if structure.space_group is None:
        raise ValueError("structural-block orbit construction requires a space group")
    refs_by_id = {
        item.block_id: _block_unit_refs(item, connectivity) for item in blocks
    }
    disjoint = _DisjointSet({item.block_id: item.block_id for item in blocks})
    incomplete: set[str] = set()
    for block in blocks:
        for operation in structure.space_group.operations:
            try:
                target_id = _mapped_block_id(
                    operation,
                    block,
                    blocks,
                    refs_by_id,
                    view,
                    representation,
                )
            except ValueError:
                target_id = None
            if target_id is None:
                incomplete.add(block.block_id)
            else:
                disjoint.union(block.block_id, target_id)

    groups: dict[str, list[StructuralBlock]] = {}
    for block in blocks:
        groups.setdefault(disjoint.find(block.block_id), []).append(block)
    replacements: dict[str, StructuralBlock] = {}
    orbits: list[StructuralBlockOrbit] = []
    diagnostics: list[Diagnostic] = []
    for rows in groups.values():
        ordered = tuple(sorted(rows, key=lambda item: item.block_id))
        orbit_id = _stable_id("block-orbit", tuple(item.block_id for item in ordered))
        members = tuple(replace(item, block_orbit_id=orbit_id) for item in ordered)
        orbit_diagnostics: tuple[Diagnostic, ...] = ()
        if any(item.block_id in incomplete for item in members):
            orbit_diagnostics = (
                Diagnostic(
                    Severity.WARNING,
                    BLOCK_SYMMETRY_INCOMPLETE,
                    "A space-group image of the structural block is missing",
                ),
            )
            diagnostics.extend(orbit_diagnostics)
        for item in members:
            replacements[item.block_id] = item
        orbits.append(
            StructuralBlockOrbit(
                block_orbit_id=orbit_id,
                representative_block_id=members[0].block_id,
                blocks=members,
                diagnostics=orbit_diagnostics,
            )
        )
    return (
        tuple(replacements[item.block_id] for item in blocks),
        tuple(sorted(orbits, key=lambda item: item.block_orbit_id)),
        tuple(dict.fromkeys(diagnostics)),
    )


def _coordinates(
    refs: tuple[PeriodicAtomRef, ...],
    view: AtomicView[ExpandedAtom],
) -> np.ndarray:
    if view.cell_matrix is None:
        raise ValueError("structural-unit geometry requires a unit cell")
    atom_by_id = {item.id: item for item in view.atoms}
    rows = []
    for ref in refs:
        atom = atom_by_id.get(ref.atom_id)
        if atom is None:
            raise ValueError("structural-unit atom is absent from the atomic view")
        fractional = np.asarray(atom.fractional, dtype=float) + np.asarray(
            ref.cell_translation, dtype=float
        )
        rows.append(fractional @ view.cell_matrix)
    return np.asarray(rows, dtype=float)


def _canonical_polygon(points: np.ndarray) -> tuple[int, ...]:
    center = points.mean(axis=0)
    _, _, vectors = np.linalg.svd(points - center)
    first, second = vectors[0], vectors[1]
    angles = np.arctan2((points - center) @ second, (points - center) @ first)
    row = tuple(int(item) for item in np.argsort(angles))
    candidates = []
    for candidate in (row, tuple(reversed(row))):
        candidates.extend(
            candidate[index:] + candidate[:index] for index in range(len(candidate))
        )
    return min(candidates)


def build_unit_geometry(
    unit: StructuralUnit,
    view: AtomicView[ExpandedAtom],
    polyhedra: dict[str, CoordinationPolyhedron],
) -> StructuralUnitGeometry | None:
    """Derive only affine/convex geometry fixed by the unit's atom positions."""

    center_ref: PeriodicAtomRef | None = None
    if unit.source_polyhedron_id is not None:
        polyhedron = polyhedra.get(unit.source_polyhedron_id)
        if polyhedron is None:
            return None
        refs = tuple(item.atom_ref for item in polyhedron.vertices)
        if not refs:
            return None
        return StructuralUnitGeometry(
            kind=StructuralUnitGeometryKind.POLYHEDRON,
            affine_dimension=3,
            vertex_atom_refs=refs,
            faces=polyhedron.faces,
            center_atom_ref=polyhedron.center_atom_ref,
            diagnostics=polyhedron.diagnostics,
        ) if polyhedron.faces else None
    if unit.source_coordination_id is not None and len(unit.atom_refs) > 1:
        center_ref = unit.atom_refs[0]
        refs = tuple(
            sorted(
                unit.atom_refs[1:],
                key=lambda item: (item.atom_id, item.cell_translation),
            )
        )
    else:
        refs = tuple(
            sorted(
                unit.atom_refs,
                key=lambda item: (item.atom_id, item.cell_translation),
            )
        )
    if not refs:
        return None
    points = _coordinates(refs, view)
    dimension = 0 if len(points) == 1 else int(
        np.linalg.matrix_rank(points - points.mean(axis=0), tol=1e-8)
    )
    if dimension == 0:
        return StructuralUnitGeometry(
            StructuralUnitGeometryKind.POINT, 0, refs, center_atom_ref=center_ref
        )
    if dimension == 1:
        return StructuralUnitGeometry(
            StructuralUnitGeometryKind.LINEAR, 1, refs, center_atom_ref=center_ref
        )
    if dimension == 2 and len(refs) >= 3:
        return StructuralUnitGeometry(
            StructuralUnitGeometryKind.PLANAR_POLYGON,
            2,
            refs,
            (_canonical_polygon(points),),
            center_ref,
        )
    if dimension == 3 and len(refs) >= 4:
        faces = _convex_hull_faces(points, 1e-9)
        if faces:
            return StructuralUnitGeometry(
                StructuralUnitGeometryKind.POLYHEDRON,
                3,
                refs,
                faces,
                center_ref,
            )
    return None


__all__ = [
    "BLOCK_SYMMETRY_INCOMPLETE",
    "UNIT_GEOMETRY_INCOMPLETE",
    "UNIT_SYMMETRY_INCOMPLETE",
    "build_block_orbits",
    "build_unit_geometry",
    "build_unit_orbits",
]
