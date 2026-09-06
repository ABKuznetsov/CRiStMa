"""Space-group orbits of coordination polyhedra."""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json

from cristma.diagnostics import Diagnostic, Severity
from cristma.structure import AtomicView, ExpandedAtom
from cristma.symmetry import AffineOperation

from .contact_orbits import _match_atom, _transform
from .contacts import ResolutionStatus
from .polyhedra import CoordinationPolyhedron, CoordinationPolyhedronOrbit


SYMMETRY_INCOMPLETE = "crystal_chemistry.polyhedron_orbit.symmetry_incomplete"
GEOMETRY_INCONSISTENT = "crystal_chemistry.polyhedron_orbit.geometry_inconsistent"

Translation = tuple[int, int, int]
VertexKey = tuple[str, Translation, str]
PolyhedronKey = tuple[str, tuple[VertexKey, ...]]


@dataclass(frozen=True, slots=True)
class PolyhedronOrbitBuildResult:
    polyhedra: tuple[CoordinationPolyhedron, ...]
    polyhedron_orbits: tuple[CoordinationPolyhedronOrbit, ...]
    diagnostics: tuple[Diagnostic, ...]
    complete: bool


def _polyhedron_key(polyhedron: CoordinationPolyhedron) -> PolyhedronKey:
    return (
        polyhedron.center_atom_id,
        tuple(sorted(
            (
                vertex.atom_ref.atom_id,
                vertex.atom_ref.cell_translation,
                vertex.contact_orbit_id,
            )
            for vertex in polyhedron.vertices
        )),
    )


def _transformed_key(
    polyhedron: CoordinationPolyhedron,
    operation: AffineOperation,
    atoms: dict[str, ExpandedAtom],
    atoms_by_site: dict[str, tuple[ExpandedAtom, ...]],
    tolerance: float,
) -> PolyhedronKey | None:
    center = atoms[polyhedron.center_atom_ref.atom_id]
    center_fractional, center_lattice = _transform(
        operation,
        center.fractional,
        polyhedron.center_atom_ref.cell_translation,
        tolerance,
    )
    center_image = _match_atom(
        atoms_by_site,
        center.source_site_id,
        center_fractional,
        tolerance,
    )
    if center_image is None:
        return None

    vertices: list[VertexKey] = []
    for vertex in polyhedron.vertices:
        ligand = atoms[vertex.atom_ref.atom_id]
        ligand_fractional, ligand_lattice = _transform(
            operation,
            ligand.fractional,
            vertex.atom_ref.cell_translation,
            tolerance,
        )
        ligand_image = _match_atom(
            atoms_by_site,
            ligand.source_site_id,
            ligand_fractional,
            tolerance,
        )
        if ligand_image is None:
            return None
        relative_translation = tuple(
            ligand_lattice[index] - center_lattice[index]
            for index in range(3)
        )
        vertices.append((
            ligand_image.id,
            relative_translation,
            vertex.contact_orbit_id,
        ))
    return center_image.id, tuple(sorted(vertices))


def _orbit_id(keys: tuple[PolyhedronKey, ...]) -> str:
    payload = json.dumps(keys, ensure_ascii=True, separators=(",", ":")).encode("ascii")
    return f"polyhedron-orbit:{hashlib.sha256(payload).hexdigest()}"


def _aggregate_status(polyhedra: tuple[CoordinationPolyhedron, ...]) -> ResolutionStatus:
    statuses = {item.status for item in polyhedra}
    if ResolutionStatus.INCOMPLETE in statuses:
        return ResolutionStatus.INCOMPLETE
    if ResolutionStatus.AMBIGUOUS in statuses:
        return ResolutionStatus.AMBIGUOUS
    return ResolutionStatus.RESOLVED


def build_polyhedron_orbits(
    view: AtomicView[ExpandedAtom],
    polyhedra: tuple[CoordinationPolyhedron, ...],
    operations: tuple[AffineOperation, ...],
    tolerance: float,
) -> PolyhedronOrbitBuildResult:
    """Group polyhedra by the exact space-group action on oriented vertices."""

    if not polyhedra:
        return PolyhedronOrbitBuildResult((), (), (), True)
    if view.fractional is None or not all(view.periodic):
        raise ValueError("polyhedron orbit construction requires a 3D periodic atomic view")
    if not operations:
        raise ValueError("polyhedron orbit construction requires symmetry operations")
    if tolerance <= 0:
        raise ValueError("polyhedron orbit tolerance must be positive")

    atoms = {atom.id: atom for atom in view.atoms}
    atoms_by_site_rows: dict[str, list[ExpandedAtom]] = {}
    for atom in view.atoms:
        atoms_by_site_rows.setdefault(atom.source_site_id, []).append(atom)
    atoms_by_site = {
        site_id: tuple(sorted(rows, key=lambda item: item.id))
        for site_id, rows in atoms_by_site_rows.items()
    }
    by_key: dict[PolyhedronKey, CoordinationPolyhedron] = {}
    for polyhedron in polyhedra:
        key = _polyhedron_key(polyhedron)
        if key in by_key:
            raise ValueError("polyhedra must identify unique oriented coordination shells")
        if key[0] not in atoms or any(vertex[0] not in atoms for vertex in key[1]):
            raise ValueError("polyhedron references an atom outside the atomic view")
        by_key[key] = polyhedron

    parent = {key: key for key in by_key}

    def find(key: PolyhedronKey) -> PolyhedronKey:
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def union(left: PolyhedronKey, right: PolyhedronKey) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    diagnostics: list[Diagnostic] = []
    incomplete_keys: set[PolyhedronKey] = set()
    for key, polyhedron in by_key.items():
        for operation in operations:
            transformed = _transformed_key(
                polyhedron, operation, atoms, atoms_by_site, tolerance
            )
            if transformed is None or transformed not in by_key:
                incomplete_keys.add(key)
                continue
            union(key, transformed)

    groups: dict[PolyhedronKey, list[tuple[PolyhedronKey, CoordinationPolyhedron]]] = {}
    for key, polyhedron in by_key.items():
        groups.setdefault(find(key), []).append((key, polyhedron))

    replacements: dict[str, CoordinationPolyhedron] = {}
    orbits: list[CoordinationPolyhedronOrbit] = []
    for rows in groups.values():
        ordered = tuple(sorted(rows, key=lambda item: (item[0], item[1].polyhedron_id)))
        keys = tuple(item[0] for item in ordered)
        orbit_id = _orbit_id(keys)
        members = tuple(replace(item[1], polyhedron_orbit_id=orbit_id) for item in ordered)
        orbit_diagnostics: list[Diagnostic] = []
        if any(key in incomplete_keys for key in keys):
            orbit_diagnostics.append(Diagnostic(
                Severity.WARNING,
                SYMMETRY_INCOMPLETE,
                "A space-group image of the oriented polyhedron is missing",
            ))
        signatures = {
            (item.coordination_number, item.face_signature) for item in members
        }
        if len(signatures) != 1:
            orbit_diagnostics.append(Diagnostic(
                Severity.ERROR,
                GEOMETRY_INCONSISTENT,
                "Symmetry-equivalent polyhedra have different face incidence graphs",
            ))
        status = (
            ResolutionStatus.INCOMPLETE
            if orbit_diagnostics
            else _aggregate_status(members)
        )
        updated_members = tuple(
            replace(
                item,
                status=(ResolutionStatus.INCOMPLETE if orbit_diagnostics else item.status),
                diagnostics=tuple(dict.fromkeys((*item.diagnostics, *orbit_diagnostics))),
            )
            for item in members
        )
        for item in updated_members:
            replacements[item.polyhedron_id] = item
        diagnostics.extend(orbit_diagnostics)
        orbits.append(CoordinationPolyhedronOrbit(
            polyhedron_orbit_id=orbit_id,
            representative_polyhedron_id=updated_members[0].polyhedron_id,
            polyhedra=updated_members,
            status=status,
            diagnostics=tuple(orbit_diagnostics),
        ))

    output = tuple(replacements[item.polyhedron_id] for item in polyhedra)
    return PolyhedronOrbitBuildResult(
        output,
        tuple(sorted(orbits, key=lambda item: item.polyhedron_orbit_id)),
        tuple(dict.fromkeys(diagnostics)),
        not incomplete_keys and not diagnostics,
    )


__all__ = [
    "GEOMETRY_INCONSISTENT",
    "PolyhedronOrbitBuildResult",
    "SYMMETRY_INCOMPLETE",
    "build_polyhedron_orbits",
]
