"""Native coordination-polyhedron construction from resolved shells."""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
import math

import numpy as np

from cristma.diagnostics import Diagnostic
from cristma.structure import PeriodicAtomRef
from .models import ResolutionStatus


FaceSignature = tuple[int, tuple[tuple[int, ...], ...]]
Edge = tuple[int, int]

PARTIAL_OCCUPANCY = "crystal_chemistry.polyhedron.partial_occupancy"
MIXED_POSITION = "crystal_chemistry.polyhedron.mixed_position"
GEOMETRY_DEGENERATE = "crystal_chemistry.polyhedron.geometry_degenerate"
DESCRIPTOR_UNAVAILABLE = "crystal_chemistry.polyhedron.descriptor_unavailable"


@dataclass(frozen=True, slots=True)
class PolyhedronVertex:
    """One oriented periodic ligand vertex in a coordination polyhedron."""

    atom_ref: PeriodicAtomRef
    incidence_orbit_id: str
    resolved_contact_orbit_id: str
    local_cartesian: tuple[float, float, float]
    distance: float
    occupancy: float

    def __post_init__(self) -> None:
        if not self.incidence_orbit_id or not self.resolved_contact_orbit_id:
            raise ValueError("polyhedron vertex incidence identities must not be empty")
        if len(self.local_cartesian) != 3 or not all(
            math.isfinite(value) for value in self.local_cartesian
        ):
            raise ValueError("polyhedron vertex requires three finite coordinates")
        if not math.isfinite(self.distance) or self.distance <= 0:
            raise ValueError("polyhedron vertex distance must be finite and positive")
        if not math.isfinite(self.occupancy) or not 0 <= self.occupancy <= 1:
            raise ValueError("polyhedron vertex occupancy must lie between zero and one")


@dataclass(frozen=True, slots=True)
class CoordinationPolyhedron:
    polyhedron_id: str
    source_site_id: str
    center_atom_id: str
    shell_provenance: tuple[tuple[str, object], ...]
    local_vertices: tuple[tuple[float, float, float], ...]
    faces: tuple[tuple[int, ...], ...]
    volume: float | None
    geometric_centroid: tuple[float, float, float] | None
    center_offset: float | None
    diagnostics: tuple[Diagnostic, ...] = ()
    polyhedron_orbit_id: str = field(default="", kw_only=True)
    source_shell_orbit_id: str | None = field(default=None, kw_only=True)
    center_atom_ref: PeriodicAtomRef = field(init=False)
    coordination_number: int | None = field(default=None, kw_only=True)
    ligand_composition: tuple[tuple[str, float], ...] = field(default=(), kw_only=True)
    vertices: tuple[PolyhedronVertex, ...] = field(default=(), kw_only=True)
    face_signature: FaceSignature | None = field(default=None, kw_only=True)
    mean_bond_length: float | None = field(default=None, kw_only=True)
    min_bond_length: float | None = field(default=None, kw_only=True)
    max_bond_length: float | None = field(default=None, kw_only=True)
    bond_length_distortion: float | None = field(default=None, kw_only=True)
    edge_angle_dispersion_deg: float | None = field(default=None, kw_only=True)
    status: ResolutionStatus = field(default=ResolutionStatus.RESOLVED, kw_only=True)
    provenance: tuple[tuple[str, object], ...] = field(default=(), kw_only=True)

    def __post_init__(self) -> None:
        if not self.polyhedron_id or not self.source_site_id or not self.center_atom_id:
            raise ValueError("polyhedron identities must not be empty")
        object.__setattr__(
            self,
            "center_atom_ref",
            PeriodicAtomRef(self.center_atom_id, (0, 0, 0)),
        )
        if self.coordination_number is None:
            object.__setattr__(self, "coordination_number", len(self.vertices))
        if self.coordination_number != len(self.vertices):
            raise ValueError("coordination number must equal the vertex count")
        atom_refs = tuple(item.atom_ref for item in self.vertices)
        if len(atom_refs) != len(set(atom_refs)):
            raise ValueError("polyhedron vertices must identify unique periodic atoms")
        if tuple(item.local_cartesian for item in self.vertices) != self.local_vertices:
            raise ValueError("polyhedron local vertices must agree with vertex records")
        labels = tuple(item[0] for item in self.ligand_composition)
        if labels != tuple(sorted(labels)) or len(labels) != len(set(labels)):
            raise ValueError("ligand composition must have unique sorted labels")
        if any(
            not label or not math.isfinite(amount) or amount <= 0
            for label, amount in self.ligand_composition
        ):
            raise ValueError("ligand composition amounts must be positive and finite")
        if self.faces and self.face_signature is None:
            object.__setattr__(
                self,
                "face_signature",
                canonical_face_signature(len(self.vertices), self.faces),
            )
        for face in self.faces:
            if len(face) < 3 or len(set(face)) != len(face):
                raise ValueError("polyhedron face must contain at least three unique vertices")
            if any(index < 0 or index >= self.coordination_number for index in face):
                raise ValueError("polyhedron face references an unknown vertex")


@dataclass(frozen=True, slots=True)
class CoordinationPolyhedronOrbit:
    """One space-group orbit of coordination polyhedra."""

    polyhedron_orbit_id: str
    representative_polyhedron_id: str
    polyhedra: tuple[CoordinationPolyhedron, ...]
    status: ResolutionStatus
    diagnostics: tuple[Diagnostic, ...] = ()

    def __post_init__(self) -> None:
        if not self.polyhedron_orbit_id or not self.representative_polyhedron_id:
            raise ValueError("polyhedron orbit identities must not be empty")
        if not self.polyhedra:
            raise ValueError("polyhedron orbit must not be empty")
        ids = tuple(item.polyhedron_id for item in self.polyhedra)
        if len(ids) != len(set(ids)):
            raise ValueError("polyhedron orbit members must have unique IDs")
        if ids[0] != self.representative_polyhedron_id:
            raise ValueError("representative polyhedron must be the first orbit member")
        if any(item.polyhedron_orbit_id != self.polyhedron_orbit_id for item in self.polyhedra):
            raise ValueError("polyhedron orbit member has a different orbit ID")

    @property
    def representative(self) -> CoordinationPolyhedron:
        return self.polyhedra[0]


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


def unique_hull_edges(faces: tuple[tuple[int, ...], ...]) -> tuple[Edge, ...]:
    """Return every unordered face-boundary edge exactly once."""

    return tuple(sorted({
        tuple(sorted((face[index], face[(index + 1) % len(face)])))
        for face in faces
        for index in range(len(face))
    }))


def _rotate_to_edge(face: tuple[int, ...], first: int, second: int) -> tuple[int, ...] | None:
    for candidate in (face, tuple(reversed(face))):
        for index, value in enumerate(candidate):
            if value == first and candidate[(index + 1) % len(candidate)] == second:
                return candidate[index:] + candidate[:index]
    return None


def _canonical_embedding_from_root(
    faces: tuple[tuple[int, ...], ...],
    root_index: int,
    root: tuple[int, ...],
) -> tuple[tuple[int, ...], ...] | None:
    edge_faces: dict[Edge, list[int]] = {}
    for face_index, face in enumerate(faces):
        for edge in unique_hull_edges((face,)):
            edge_faces.setdefault(edge, []).append(face_index)
    if any(len(indices) != 2 for indices in edge_faces.values()):
        return None

    oriented = {root_index: root}
    queue = [root_index]
    labels: dict[int, int] = {}
    while queue:
        face_index = queue.pop(0)
        face = oriented[face_index]
        for vertex in face:
            labels.setdefault(vertex, len(labels))
        for index, first in enumerate(face):
            second = face[(index + 1) % len(face)]
            adjacent = next(
                item for item in edge_faces[tuple(sorted((first, second)))]
                if item != face_index
            )
            candidate = _rotate_to_edge(faces[adjacent], second, first)
            if candidate is None:
                return None
            if adjacent not in oriented:
                oriented[adjacent] = candidate
                queue.append(adjacent)
    if len(oriented) != len(faces):
        return None
    encoded = []
    for face in oriented.values():
        row = tuple(labels[item] for item in face)
        rotations = tuple(row[index:] + row[:index] for index in range(len(row)))
        encoded.append(min(rotations))
    return tuple(sorted(encoded, key=lambda row: (len(row), row)))


def canonical_face_signature(
    vertex_count: int,
    faces: tuple[tuple[int, ...], ...],
) -> FaceSignature:
    """Canonically encode a closed vertex-edge-face incidence graph."""

    if vertex_count < 0 or any(
        len(face) < 3
        or len(set(face)) != len(face)
        or any(index < 0 or index >= vertex_count for index in face)
        for face in faces
    ):
        raise ValueError("invalid polyhedron face graph")
    if not faces:
        return vertex_count, ()
    candidates: list[tuple[tuple[int, ...], ...]] = []
    for face_index, face in enumerate(faces):
        for reverse in (False, True):
            row = tuple(reversed(face)) if reverse else face
            for offset in range(len(row)):
                root = row[offset:] + row[:offset]
                encoded = _canonical_embedding_from_root(faces, face_index, root)
                if encoded is not None:
                    candidates.append(encoded)
    if not candidates:
        raise ValueError("polyhedron face graph must be a connected closed manifold")
    return vertex_count, min(candidates)


def _bond_length_metrics(
    distances: tuple[float, ...],
) -> tuple[float | None, float | None, float | None, float | None]:
    if not distances:
        return None, None, None, None
    mean = math.fsum(distances) / len(distances)
    minimum = min(distances)
    maximum = max(distances)
    if not math.isfinite(mean) or mean <= 0:
        return mean, minimum, maximum, None
    distortion = math.fsum(abs(value - mean) / mean for value in distances) / len(distances)
    return mean, minimum, maximum, distortion


def _edge_angle_dispersion_deg(
    vertices: np.ndarray,
    edges: tuple[Edge, ...],
) -> float | None:
    if not edges:
        return None
    angles: list[float] = []
    for first, second in edges:
        first_norm = float(np.linalg.norm(vertices[first]))
        second_norm = float(np.linalg.norm(vertices[second]))
        if first_norm <= 0 or second_norm <= 0:
            return None
        cosine = float(np.dot(vertices[first], vertices[second])) / (first_norm * second_norm)
        angles.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
    mean = math.fsum(angles) / len(angles)
    return math.sqrt(math.fsum((value - mean) ** 2 for value in angles) / len(angles))


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


def polyhedron_face_signature(polyhedron: CoordinationPolyhedron) -> FaceSignature | None:
    """Return the stored canonical vertex-edge-face signature."""

    return polyhedron.face_signature


__all__ = [
    "CoordinationPolyhedron",
    "CoordinationPolyhedronOrbit",
    "FaceSignature",
    "PolyhedronVertex",
    "canonical_face_signature",
    "polyhedron_face_signature",
    "unique_hull_edges",
]
