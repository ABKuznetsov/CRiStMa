"""Space-group orbits of coordination polyhedra."""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
import math

import numpy as np

from cristma.crystallography import periodic_endpoint_instance
from cristma.crystallography.symmetry_context import _digest
from cristma.diagnostics import Diagnostic, Severity
from cristma.structure import AtomicView, ExpandedAtom, PeriodicAtomRef
from cristma.symmetry import AffineOperation

from .contact_orbits import _match_atom, _transform
from .contact_analysis import ContactAnalysisResult
from .contacts import ResolutionStatus
from .polyhedra import (
    DESCRIPTOR_UNAVAILABLE,
    GEOMETRY_DEGENERATE,
    MIXED_POSITION,
    PARTIAL_OCCUPANCY,
    CoordinationPolyhedron,
    CoordinationPolyhedronOrbit,
    PolyhedronVertex,
    _bond_length_metrics,
    _convex_hull_faces,
    _edge_angle_dispersion_deg,
    _volume_and_centroid,
    canonical_face_signature,
    unique_hull_edges,
)


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

    @property
    def diagnostic_codes(self) -> tuple[str, ...]:
        return tuple(item.code for item in self.diagnostics)


def _polyhedron_key(polyhedron: CoordinationPolyhedron) -> PolyhedronKey:
    return (
        polyhedron.center_atom_id,
        tuple(sorted(
            (
                vertex.atom_ref.atom_id,
                vertex.atom_ref.cell_translation,
                vertex.resolved_contact_orbit_id,
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
            vertex.resolved_contact_orbit_id,
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


def _source_image(contact_result: ContactAnalysisResult, site_id: str):
    site_mapping = contact_result._asymmetric_unit_mapping.by_site_id[site_id]
    matches = tuple(
        image
        for image in site_mapping.reference_cell_images
        if image.equivalent_relations == site_mapping.stabilizer_relations
    )
    if len(matches) != 1:
        raise ValueError("independent centre does not identify one source image")
    return matches[0]


def _image_by_id(contact_result: ContactAnalysisResult, site_id: str, image_id: str):
    return next(
        image
        for image in contact_result._asymmetric_unit_mapping.by_site_id[
            site_id
        ].reference_cell_images
        if image.image_id == image_id
    )


def _local_polyhedron_vertices(
    contact_result: ContactAnalysisResult,
    shell,
    alternative,
) -> tuple[PolyhedronVertex, ...]:
    incidence_by_id = {
        item.incidence_orbit_id: item
        for item in contact_result.contact_incidence_orbits
    }
    center_image = _source_image(contact_result, shell.center_independent_site_id)
    center_fractional = np.asarray(center_image.fractional_position, dtype=float)
    matrix = contact_result._structure.cell.matrix
    vertices: list[PolyhedronVertex] = []
    for incidence_id in alternative.primary_incidence_ids:
        incidence = incidence_by_id[incidence_id]
        for relation in incidence.equivalent_oriented_relations:
            _, image_id, cell_translation = periodic_endpoint_instance(
                incidence.ligand_independent_site_id,
                relation,
                contact_result._asymmetric_unit_mapping,
            )
            image = _image_by_id(
                contact_result,
                incidence.ligand_independent_site_id,
                image_id,
            )
            ligand_fractional = np.asarray(image.fractional_position, dtype=float)
            ligand_fractional += np.asarray(cell_translation, dtype=float)
            local = (ligand_fractional - center_fractional) @ matrix
            distance = float(np.linalg.norm(local))
            vertices.append(
                PolyhedronVertex(
                    PeriodicAtomRef(image_id, cell_translation),
                    incidence.incidence_orbit_id,
                    incidence.resolved_contact_orbit_id,
                    tuple(float(value) for value in local),
                    distance,
                    incidence.effective_neighbor_occupancy,
                )
            )
    ordered = tuple(
        sorted(
            vertices,
            key=lambda item: (
                item.atom_ref.atom_id,
                item.atom_ref.cell_translation,
                item.incidence_orbit_id,
            ),
        )
    )
    if len(ordered) != alternative.geometric_CN:
        raise ValueError("realized polyhedron vertex count disagrees with shell CN")
    if len({item.atom_ref for item in ordered}) != len(ordered):
        raise ValueError("shell incidences realize a duplicate periodic ligand")
    return ordered


def _orbit_ligand_composition(
    contact_result: ContactAnalysisResult,
    vertices: tuple[PolyhedronVertex, ...],
) -> tuple[tuple[str, float], ...]:
    incidence_by_id = {
        item.incidence_orbit_id: item
        for item in contact_result.contact_incidence_orbits
    }
    interpretation_by_id = {
        item.interpretation_id: item
        for orbit in contact_result.contact_orbits
        for item in orbit.interpretations
    }
    geometry_by_resolved = {
        resolved.resolved_contact_orbit_id: next(
            geometry
            for geometry in contact_result.pair_table.contact_orbits
            if geometry.geometry_orbit_id == resolved.geometry_orbit_id
        )
        for resolved in contact_result.contact_orbits
    }
    sites = {item.id: item for item in contact_result._structure.sites}
    totals: dict[str, list[float]] = {}
    for vertex in vertices:
        incidence = incidence_by_id[vertex.incidence_orbit_id]
        interpretation = interpretation_by_id[incidence.interpretation_id]
        geometry = geometry_by_resolved[incidence.resolved_contact_orbit_id]
        if (
            incidence.center_independent_site_id == geometry.first_independent_site_id
            and incidence.ligand_independent_site_id == geometry.second_independent_site_id
        ):
            participating = {
                item.second_species
                for item in interpretation.component_pair_interpretations
            }
        elif (
            incidence.center_independent_site_id == geometry.second_independent_site_id
            and incidence.ligand_independent_site_id == geometry.first_independent_site_id
        ):
            participating = {
                item.first_species
                for item in interpretation.component_pair_interpretations
            }
        else:
            participating = {
                species
                for item in interpretation.component_pair_interpretations
                for species in (item.first_species, item.second_species)
            }
        for component in sites[incidence.ligand_independent_site_id].components:
            if component.species not in participating:
                continue
            label = component.element or component.species.label
            totals.setdefault(label, []).append(float(component.occupancy.value))
    return tuple(
        (label, math.fsum(values)) for label, values in sorted(totals.items())
    )


def _representative_polyhedron(
    contact_result: ContactAnalysisResult,
    shell,
    tolerance: float,
) -> CoordinationPolyhedron:
    alternative = shell.selected
    if alternative is None:
        raise ValueError("resolved shell lacks a selected alternative")
    vertices = _local_polyhedron_vertices(contact_result, shell, alternative)
    local = np.asarray(
        tuple(item.local_cartesian for item in vertices), dtype=float
    ).reshape((-1, 3))
    distances = tuple(item.distance for item in vertices)
    mean, minimum, maximum, distortion = _bond_length_metrics(distances)
    diagnostics: list[Diagnostic] = []
    incidence_by_id = {
        item.incidence_orbit_id: item
        for item in contact_result.contact_incidence_orbits
    }
    sites = {item.id: item for item in contact_result._structure.sites}
    if any(item.occupancy < 1.0 - 1e-12 for item in vertices):
        diagnostics.append(Diagnostic(
            Severity.WARNING,
            PARTIAL_OCCUPANCY,
            "Polyhedron contains a partially occupied ligand vertex",
        ))
    if any(
        len(sites[incidence_by_id[item.incidence_orbit_id].ligand_independent_site_id].components) > 1
        for item in vertices
    ):
        diagnostics.append(Diagnostic(
            Severity.WARNING,
            MIXED_POSITION,
            "Polyhedron contains a mixed ligand position",
        ))

    faces: tuple[tuple[int, ...], ...] = ()
    face_signature = None
    angle_dispersion = None
    volume = None
    centroid = None
    center_offset = None
    status = ResolutionStatus.RESOLVED
    is_three_dimensional = (
        len(local) >= 4
        and np.linalg.matrix_rank(local[1:] - local[0], tol=tolerance) >= 3
    )
    if is_three_dimensional:
        try:
            faces = _convex_hull_faces(local, tolerance)
            face_signature = canonical_face_signature(len(vertices), faces)
            angle_dispersion = _edge_angle_dispersion_deg(
                local, unique_hull_edges(faces)
            )
            volume_value, centroid_value = _volume_and_centroid(local, faces)
            volume = volume_value
            centroid = tuple(float(value) for value in centroid_value)
            center_offset = float(np.linalg.norm(centroid_value))
        except ValueError as error:
            status = ResolutionStatus.INCOMPLETE
            diagnostics.append(Diagnostic(Severity.WARNING, GEOMETRY_DEGENERATE, str(error)))
    else:
        status = ResolutionStatus.INCOMPLETE
        diagnostics.append(Diagnostic(
            Severity.WARNING,
            GEOMETRY_DEGENERATE,
            "Coordination shell is not a closed three-dimensional polyhedron",
        ))
    if distortion is None or angle_dispersion is None:
        diagnostics.append(Diagnostic(
            Severity.WARNING,
            DESCRIPTOR_UNAVAILABLE,
            "One or more polyhedron descriptors are unavailable",
        ))

    polyhedron_orbit_id = "polyhedron-orbit:" + _digest(
        {"shell_orbit_id": shell.shell_orbit_id}
    )
    polyhedron_id = "polyhedron:" + _digest(
        {
            "shell_orbit_id": shell.shell_orbit_id,
            "alternative_id": alternative.alternative_id,
            "vertices": tuple(
                (
                    item.atom_ref.atom_id,
                    item.atom_ref.cell_translation,
                    item.incidence_orbit_id,
                )
                for item in vertices
            ),
        }
    )
    center_image = _source_image(contact_result, shell.center_independent_site_id)
    return CoordinationPolyhedron(
        polyhedron_id=polyhedron_id,
        polyhedron_orbit_id=polyhedron_orbit_id,
        source_site_id=shell.center_independent_site_id,
        center_atom_id=center_image.image_id,
        coordination_number=len(vertices),
        ligand_composition=_orbit_ligand_composition(contact_result, vertices),
        vertices=vertices,
        faces=faces,
        face_signature=face_signature,
        mean_bond_length=mean,
        min_bond_length=minimum,
        max_bond_length=maximum,
        bond_length_distortion=distortion,
        edge_angle_dispersion_deg=angle_dispersion,
        volume=volume,
        geometric_centroid=centroid,
        center_offset=center_offset,
        status=status,
        provenance=(
            *shell.provenance,
            ("method", "cristma.orbit_polyhedron_builder:1"),
            ("selected_alternative_id", alternative.alternative_id),
            ("tolerance", tolerance),
            ("bond_distortion", "Baur mean absolute relative deviation"),
            ("edge_angle_dispersion", "population standard deviation in degrees"),
        ),
        local_vertices=tuple(item.local_cartesian for item in vertices),
        shell_provenance=shell.provenance,
        diagnostics=tuple(dict.fromkeys(diagnostics)),
    )


@dataclass(frozen=True, slots=True)
class PolyhedronOrbitBuilder:
    tolerance: float = 1e-9

    def __post_init__(self) -> None:
        if not math.isfinite(self.tolerance) or self.tolerance <= 0:
            raise ValueError("polyhedron tolerance must be positive and finite")

    def build(self, contact_result: ContactAnalysisResult) -> PolyhedronOrbitBuildResult:
        if not isinstance(contact_result, ContactAnalysisResult):
            raise TypeError("contact_result must be ContactAnalysisResult")
        polyhedra: list[CoordinationPolyhedron] = []
        orbits: list[CoordinationPolyhedronOrbit] = []
        diagnostics: list[Diagnostic] = []
        skipped = False
        for shell in contact_result.coordination_shell_orbits:
            if shell.status is not ResolutionStatus.RESOLVED:
                skipped = True
                code = (
                    "crystal_chemistry.polyhedron.shell_ambiguous"
                    if shell.status is ResolutionStatus.AMBIGUOUS
                    else "crystal_chemistry.polyhedron.shell_incomplete"
                )
                diagnostic = Diagnostic(
                    Severity.WARNING,
                    code,
                    f"polyhedron is not selected from a {shell.status.value} shell",
                )
                diagnostics.append(diagnostic)
                continue
            polyhedron = _representative_polyhedron(
                contact_result, shell, self.tolerance
            )
            orbit = CoordinationPolyhedronOrbit(
                polyhedron.polyhedron_orbit_id,
                polyhedron.polyhedron_id,
                (polyhedron,),
                polyhedron.status,
                polyhedron.diagnostics,
            )
            polyhedra.append(polyhedron)
            orbits.append(orbit)
            diagnostics.extend(polyhedron.diagnostics)
        return PolyhedronOrbitBuildResult(
            tuple(sorted(polyhedra, key=lambda item: item.polyhedron_id)),
            tuple(sorted(orbits, key=lambda item: item.polyhedron_orbit_id)),
            tuple(dict.fromkeys(diagnostics)),
            not skipped and all(
                item.status is ResolutionStatus.RESOLVED for item in polyhedra
            ),
        )


__all__ = [
    "GEOMETRY_INCONSISTENT",
    "PolyhedronOrbitBuildResult",
    "PolyhedronOrbitBuilder",
    "SYMMETRY_INCOMPLETE",
    "build_polyhedron_orbits",
]
