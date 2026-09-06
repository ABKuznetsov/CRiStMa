"""Structural-unit orbits derived without materializing contact instances."""
from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
import numpy as np

from cristma.chemistry import InteractionLayer
from cristma.diagnostics import Diagnostic
from cristma.structure import PeriodicAtomRef
from .contact_analysis import ContactAnalysisResult
from .contacts import ResolutionStatus
from .polyhedron_orbits import PolyhedronOrbitBuildResult, _source_image
from .shell_orbits import ShellRole


class StructuralUnitKind(StrEnum):
    POLYHEDRON = "polyhedron"
    COORDINATION = "coordination"
    FINITE_GROUP = "finite_group"
    ATOM = "atom"


class StructuralUnitGeometryKind(StrEnum):
    POINT = "point"
    LINEAR = "linear"
    PLANAR_POLYGON = "planar_polygon"
    POLYHEDRON = "polyhedron"


@dataclass(frozen=True, slots=True)
class StructuralUnitGeometry:
    kind: StructuralUnitGeometryKind
    affine_dimension: int
    vertex_atom_refs: tuple[PeriodicAtomRef, ...]
    faces: tuple[tuple[int, ...], ...] = ()
    center_atom_ref: PeriodicAtomRef | None = None
    diagnostics: tuple[Diagnostic, ...] = ()

    def __post_init__(self) -> None:
        expected = {StructuralUnitGeometryKind.POINT: 0, StructuralUnitGeometryKind.LINEAR: 1,
                    StructuralUnitGeometryKind.PLANAR_POLYGON: 2, StructuralUnitGeometryKind.POLYHEDRON: 3}[self.kind]
        if self.affine_dimension != expected:
            raise ValueError("structural-unit geometry kind must match affine dimension")
        if not self.vertex_atom_refs or len(set(self.vertex_atom_refs)) != len(self.vertex_atom_refs):
            raise ValueError("structural-unit geometry requires unique vertices")
        for face in self.faces:
            if len(face) < 3 or len(set(face)) != len(face) or any(i < 0 or i >= len(self.vertex_atom_refs) for i in face):
                raise ValueError("invalid structural-unit face")
        if self.kind is StructuralUnitGeometryKind.PLANAR_POLYGON and len(self.faces) != 1:
            raise ValueError("planar structural-unit geometry requires one polygon face")
        if self.affine_dimension < 2 and self.faces:
            raise ValueError("zero- and one-dimensional geometry has no faces")


@dataclass(frozen=True, slots=True)
class StructuralUnitOrbit:
    """One scientific quotient-graph node, not expanded member state."""
    unit_orbit_id: str
    kind: StructuralUnitKind
    center_independent_site_id: str
    source_shell_orbit_id: str | None
    source_polyhedron_orbit_id: str | None
    constituent_site_refs: tuple[PeriodicAtomRef, ...]
    source_resolved_contact_orbit_ids: tuple[str, ...]
    multiplicity_in_reference_cell: int
    provenance: tuple[tuple[str, object], ...]
    interaction_layers: tuple[InteractionLayer, ...] = ()
    shell_roles: tuple[ShellRole, ...] = ()
    representative_geometry: StructuralUnitGeometry | None = None
    diagnostics: tuple[Diagnostic, ...] = ()

    def __post_init__(self) -> None:
        if not self.unit_orbit_id or not self.center_independent_site_id:
            raise ValueError("structural-unit orbit identities must not be empty")
        if not self.constituent_site_refs or len(set(self.constituent_site_refs)) != len(self.constituent_site_refs):
            raise ValueError("structural-unit orbit requires unique constituent references")
        if tuple(sorted(set(self.source_resolved_contact_orbit_ids))) != self.source_resolved_contact_orbit_ids:
            raise ValueError("source contact-orbit IDs must be unique and sorted")
        if self.multiplicity_in_reference_cell <= 0:
            raise ValueError("structural-unit multiplicity must be positive")
        if self.kind is StructuralUnitKind.POLYHEDRON and not self.source_polyhedron_orbit_id:
            raise ValueError("polyhedron unit requires a polyhedron-orbit source")
        if self.kind is StructuralUnitKind.COORDINATION and not self.source_shell_orbit_id:
            raise ValueError("coordination unit requires a shell-orbit source")
        if self.kind is StructuralUnitKind.ATOM and (self.source_shell_orbit_id or self.source_polyhedron_orbit_id):
            raise ValueError("atom unit cannot reference shell geometry")

    @property
    def unit_id(self) -> str:
        return self.unit_orbit_id

    @property
    def atom_refs(self) -> tuple[PeriodicAtomRef, ...]:
        return self.constituent_site_refs


StructuralUnit = StructuralUnitOrbit


@dataclass(frozen=True, slots=True)
class StructuralUnitBuildResult:
    unit_orbits: tuple[StructuralUnitOrbit, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()

    def __post_init__(self) -> None:
        ids = tuple(item.unit_orbit_id for item in self.unit_orbits)
        if tuple(sorted(ids)) != ids or len(set(ids)) != len(ids):
            raise ValueError("structural-unit orbit IDs must be unique and sorted")

    @property
    def units(self) -> tuple[StructuralUnitOrbit, ...]:
        return self.unit_orbits


def _ordered_planar_face(coordinates) -> tuple[int, ...]:
    """Return a deterministic non-self-intersecting boundary of coplanar points."""
    points = np.asarray(coordinates, dtype=float).reshape((-1, 3))
    centered = points - points.mean(axis=0)
    _, _, axes = np.linalg.svd(centered, full_matrices=False)
    projected = centered @ axes[:2].T
    angles = np.arctan2(projected[:, 1], projected[:, 0])
    cycle = tuple(int(index) for index in np.argsort(angles))
    candidates = []
    for ordered in (cycle, tuple(reversed(cycle))):
        for offset in range(len(ordered)):
            candidates.append(ordered[offset:] + ordered[:offset])
    return min(candidates)


def _geometry(polyhedron) -> StructuralUnitGeometry:
    vertices = polyhedron.vertices
    coordinates = np.asarray(polyhedron.local_vertices, dtype=float).reshape((-1, 3))
    dimension = 0 if len(coordinates) < 2 else int(np.linalg.matrix_rank(coordinates[1:] - coordinates[0], tol=1e-9))
    if dimension >= 3 and polyhedron.faces:
        kind, faces, dimension = StructuralUnitGeometryKind.POLYHEDRON, polyhedron.faces, 3
    elif dimension == 2 and len(vertices) >= 3:
        kind, faces = StructuralUnitGeometryKind.PLANAR_POLYGON, (_ordered_planar_face(coordinates),)
    elif dimension == 1:
        kind, faces = StructuralUnitGeometryKind.LINEAR, ()
    else:
        kind, faces, dimension = StructuralUnitGeometryKind.POINT, (), 0
    return StructuralUnitGeometry(kind, dimension, tuple(v.atom_ref for v in vertices), faces,
                                  polyhedron.center_atom_ref, polyhedron.diagnostics)


@dataclass(frozen=True, slots=True)
class StructuralUnitBuilder:
    """Project shell and atom orbits into quotient-graph nodes."""
    def get_config(self) -> dict[str, object]:
        return {}

    def clone(self, **changes: object) -> "StructuralUnitBuilder":
        if changes:
            raise TypeError("unknown StructuralUnitBuilder configuration: " + ", ".join(sorted(changes)))
        return self

    def build(self, contact_result: ContactAnalysisResult,
              polyhedra: PolyhedronOrbitBuildResult) -> StructuralUnitBuildResult:
        if not isinstance(contact_result, ContactAnalysisResult) or not isinstance(polyhedra, PolyhedronOrbitBuildResult):
            raise TypeError("orbit-first contact and polyhedron results are required")
        incidence_by_id = {x.incidence_orbit_id: x for x in contact_result.contact_incidence_orbits}
        interpretation_by_id = {x.interpretation_id: x for o in contact_result.contact_orbits for x in o.interpretations}
        polyhedron_by_shell = {o.representative.source_shell_orbit_id: o for o in polyhedra.polyhedron_orbits}
        represented_centers: set[str] = set()
        units: list[StructuralUnitOrbit] = []
        for shell in contact_result.coordination_shell_orbits:
            if shell.status is not ResolutionStatus.RESOLVED or shell.selected is None:
                continue
            poly_orbit = polyhedron_by_shell.get(shell.shell_orbit_id)
            if poly_orbit is None:
                raise ValueError("resolved shell has no polyhedron-orbit projection")
            poly = poly_orbit.representative
            incidences = tuple(incidence_by_id[x] for x in shell.selected.primary_incidence_ids)
            source_ids = tuple(sorted({x.resolved_contact_orbit_id for x in incidences}))
            layers = tuple(sorted({interpretation_by_id[x.interpretation_id].interaction_layer for x in incidences}, key=lambda x: x.value))
            geometry = _geometry(poly)
            kind = StructuralUnitKind.POLYHEDRON if geometry.kind is StructuralUnitGeometryKind.POLYHEDRON else StructuralUnitKind.COORDINATION
            mapping = contact_result._asymmetric_unit_mapping.by_site_id[shell.center_independent_site_id]
            units.append(StructuralUnitOrbit(
                "structural-unit-orbit:" + shell.shell_orbit_id, kind, shell.center_independent_site_id,
                shell.shell_orbit_id, poly_orbit.polyhedron_orbit_id if kind is StructuralUnitKind.POLYHEDRON else None,
                tuple(dict.fromkeys((poly.center_atom_ref, *(v.atom_ref for v in poly.vertices)))), source_ids,
                len(mapping.reference_cell_images), (("method", "cristma.structural_unit_builder:3"),),
                layers, (ShellRole.PRIMARY,), geometry, poly.diagnostics,
            ))
            represented_centers.add(shell.center_independent_site_id)
        for site in contact_result._structure.sites:
            if site.id in represented_centers:
                continue
            image = _source_image(contact_result, site.id)
            mapping = contact_result._asymmetric_unit_mapping.by_site_id[site.id]
            atom_ref = PeriodicAtomRef(image.image_id, (0, 0, 0))
            units.append(StructuralUnitOrbit(
                "structural-unit-orbit:atom:" + site.id, StructuralUnitKind.ATOM, site.id, None, None,
                (atom_ref,), (), len(mapping.reference_cell_images), (("method", "cristma.structural_unit_builder:3"),),
                representative_geometry=StructuralUnitGeometry(StructuralUnitGeometryKind.POINT, 0, (atom_ref,), center_atom_ref=atom_ref),
            ))
        return StructuralUnitBuildResult(tuple(sorted(units, key=lambda x: x.unit_orbit_id)),
                                         tuple(dict.fromkeys((*contact_result.diagnostics, *polyhedra.diagnostics))),
                                         (("method", "cristma.structural_unit_builder:3"),))


__all__ = ["StructuralUnit", "StructuralUnitBuildResult", "StructuralUnitBuilder",
           "StructuralUnitGeometry", "StructuralUnitGeometryKind", "StructuralUnitKind", "StructuralUnitOrbit"]
