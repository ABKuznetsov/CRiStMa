"""Native coordination-polyhedron construction from resolved shells."""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import combinations
import math

import numpy as np

from cristma.diagnostics import Diagnostic, Severity
from cristma.structure import AtomicView

from .contacts import CoordinationShell, ResolutionStatus, ResolvedContact


@dataclass(frozen=True, slots=True)
class CoordinationPolyhedron:
    polyhedron_id: str
    source_site_id: str
    center_atom_id: str
    shell_provenance: tuple[tuple[str, object], ...]
    vertex_contacts: tuple[ResolvedContact, ...]
    local_vertices: tuple[tuple[float, float, float], ...]
    faces: tuple[tuple[int, ...], ...]
    volume: float
    geometric_centroid: tuple[float, float, float]
    center_offset: float
    diagnostics: tuple[Diagnostic, ...] = ()


@dataclass(frozen=True, slots=True)
class PolyhedronBuildResult:
    status: ResolutionStatus
    polyhedron: CoordinationPolyhedron | None
    diagnostics: tuple[Diagnostic, ...]

    @classmethod
    def resolved(cls, polyhedron: CoordinationPolyhedron) -> "PolyhedronBuildResult":
        return cls(ResolutionStatus.RESOLVED, polyhedron, polyhedron.diagnostics)

    @classmethod
    def failure(
        cls,
        status: ResolutionStatus,
        code: str,
        message: str,
    ) -> "PolyhedronBuildResult":
        return cls(status, None, (Diagnostic(Severity.WARNING, code, message),))


def _local_vertices(
    shell: CoordinationShell,
    view: AtomicView,
) -> np.ndarray:
    atom_ids = set(view.ids)
    if shell.center_atom_id not in atom_ids:
        raise ValueError("coordination-shell centre is absent from the atomic view")
    rows = []
    for contact in shell.contacts:
        geometric = contact.geometric_contact
        if shell.center_atom_id == geometric.first_atom_id:
            ligand_id = geometric.second_atom_id
            vector = geometric.vector_cartesian
        elif shell.center_atom_id == geometric.second_atom_id:
            ligand_id = geometric.first_atom_id
            vector = tuple(-value for value in geometric.vector_cartesian)
        else:
            raise ValueError("coordination-shell contact does not contain its centre")
        if ligand_id not in atom_ids:
            raise ValueError("coordination-shell ligand is absent from the atomic view")
        rows.append(vector)
    return np.asarray(rows, dtype=float).reshape((-1, 3))


def _ordered_face(
    vertices: np.ndarray,
    indices: tuple[int, ...],
    outward_normal: np.ndarray,
) -> tuple[int, ...]:
    face_points = vertices[list(indices)]
    center = face_points.mean(axis=0)
    first = face_points[0] - center
    first /= np.linalg.norm(first)
    second = np.cross(outward_normal, first)
    angles = np.arctan2(
        (face_points - center) @ second,
        (face_points - center) @ first,
    )
    ordered = tuple(indices[index] for index in np.argsort(angles))
    points = vertices[list(ordered)]
    polygon_normal = np.zeros(3)
    for index in range(len(points)):
        polygon_normal += np.cross(points[index], points[(index + 1) % len(points)])
    if float(np.dot(polygon_normal, outward_normal)) < 0:
        ordered = tuple(reversed(ordered))
    return ordered


def _convex_hull_faces(vertices: np.ndarray, tolerance: float) -> tuple[tuple[int, ...], ...]:
    interior = vertices.mean(axis=0)
    planes: dict[frozenset[int], np.ndarray] = {}
    for first, second, third in combinations(range(len(vertices)), 3):
        normal = np.cross(vertices[second] - vertices[first], vertices[third] - vertices[first])
        norm = float(np.linalg.norm(normal))
        if norm <= tolerance:
            continue
        normal /= norm
        signed = (vertices - vertices[first]) @ normal
        if not (np.all(signed <= tolerance) or np.all(signed >= -tolerance)):
            continue
        if float(np.dot(interior - vertices[first], normal)) > 0:
            normal = -normal
            signed = -signed
        face = frozenset(int(index) for index in np.flatnonzero(np.abs(signed) <= tolerance))
        if len(face) >= 3:
            planes[face] = normal
    ordered = (
        _ordered_face(vertices, tuple(sorted(indices)), normal)
        for indices, normal in planes.items()
    )
    return tuple(sorted(ordered, key=lambda face: (len(face), face)))


def _volume_and_centroid(
    vertices: np.ndarray,
    faces: tuple[tuple[int, ...], ...],
) -> tuple[float, np.ndarray]:
    signed_volume = 0.0
    first_moment = np.zeros(3)
    for face in faces:
        anchor = vertices[face[0]]
        for index in range(1, len(face) - 1):
            second = vertices[face[index]]
            third = vertices[face[index + 1]]
            tetrahedron_volume = float(np.dot(anchor, np.cross(second, third))) / 6.0
            signed_volume += tetrahedron_volume
            first_moment += tetrahedron_volume * (anchor + second + third) / 4.0
    if math.isclose(signed_volume, 0.0, abs_tol=1e-15):
        raise ValueError("polyhedron hull has zero volume")
    return abs(signed_volume), first_moment / signed_volume


def polyhedron_face_signature(polyhedron: CoordinationPolyhedron) -> tuple[int, tuple[int, ...]]:
    """Return an identity-independent face-topology signature."""

    return len(polyhedron.local_vertices), tuple(sorted(len(face) for face in polyhedron.faces))


@dataclass(frozen=True, slots=True)
class PolyhedronBuilder:
    tolerance: float = 1e-9

    def __post_init__(self) -> None:
        if not math.isfinite(self.tolerance) or self.tolerance <= 0:
            raise ValueError("polyhedron tolerance must be positive and finite")

    def get_config(self) -> dict[str, float]:
        return {"tolerance": self.tolerance}

    def clone(self, **changes: float) -> "PolyhedronBuilder":
        return replace(self, **changes)

    def build(
        self,
        shell: CoordinationShell,
        view: AtomicView,
    ) -> PolyhedronBuildResult:
        if shell.status is not ResolutionStatus.RESOLVED:
            return PolyhedronBuildResult.failure(
                shell.status,
                "crystal_chemistry.polyhedron.shell_not_resolved",
                f"polyhedron cannot be built from a {shell.status.value} shell",
            )
        vertices = _local_vertices(shell, view)
        if len(vertices) < 4 or np.linalg.matrix_rank(
            vertices[1:] - vertices[0], tol=self.tolerance
        ) < 3:
            return PolyhedronBuildResult.failure(
                ResolutionStatus.NOT_APPLICABLE,
                "crystal_chemistry.polyhedron.not_three_dimensional",
                "coordination shell is not a three-dimensional polyhedron",
            )
        faces = _convex_hull_faces(vertices, self.tolerance)
        volume, centroid = _volume_and_centroid(vertices, faces)
        polyhedron = CoordinationPolyhedron(
            polyhedron_id=f"polyhedron:{shell.center_atom_id}",
            source_site_id=shell.source_site_id,
            center_atom_id=shell.center_atom_id,
            shell_provenance=shell.provenance,
            vertex_contacts=shell.contacts,
            local_vertices=tuple(tuple(float(value) for value in row) for row in vertices),
            faces=faces,
            volume=volume,
            geometric_centroid=tuple(float(value) for value in centroid),
            center_offset=float(np.linalg.norm(centroid)),
        )
        return PolyhedronBuildResult.resolved(polyhedron)

    def validate_orbit(
        self,
        polyhedra: tuple[CoordinationPolyhedron, ...],
    ) -> tuple[Diagnostic, ...]:
        signatures = {polyhedron_face_signature(item) for item in polyhedra}
        if len(signatures) <= 1:
            return ()
        return (Diagnostic(
            Severity.ERROR,
            "crystal_chemistry.polyhedron.symmetry_inconsistent",
            "symmetry-equivalent centres have different face topology",
        ),)


__all__ = [
    "CoordinationPolyhedron",
    "PolyhedronBuildResult",
    "PolyhedronBuilder",
    "polyhedron_face_signature",
]
