"""Native coordination-polyhedron construction from resolved shells."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import hashlib
from itertools import combinations
import json
import math

import numpy as np

from cristma.diagnostics import Diagnostic, Severity
from cristma.structure import AtomicView, PeriodicAtomRef

from .contacts import CoordinationShell, ResolutionStatus, ResolvedContact


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


@dataclass(frozen=True, slots=True)
class PolyhedronBuildResult:
    status: ResolutionStatus
    polyhedron: CoordinationPolyhedron | None
    diagnostics: tuple[Diagnostic, ...]

    @classmethod
    def from_polyhedron(cls, polyhedron: CoordinationPolyhedron) -> "PolyhedronBuildResult":
        return cls(polyhedron.status, polyhedron, polyhedron.diagnostics)

    @classmethod
    def resolved(cls, polyhedron: CoordinationPolyhedron) -> "PolyhedronBuildResult":
        return cls.from_polyhedron(polyhedron)

    @classmethod
    def failure(
        cls,
        status: ResolutionStatus,
        code: str,
        message: str,
    ) -> "PolyhedronBuildResult":
        return cls(status, None, (Diagnostic(Severity.WARNING, code, message),))


def _oriented_contact(
    center_atom_id: str,
    contact: ResolvedContact,
) -> tuple[PeriodicAtomRef, tuple[float, float, float]]:
    geometric = contact.geometric_contact
    translation = geometric.cell_translation or (0, 0, 0)
    if center_atom_id == geometric.first_atom_id:
        return PeriodicAtomRef(geometric.second_atom_id, translation), geometric.vector_cartesian
    if center_atom_id == geometric.second_atom_id:
        return (
            PeriodicAtomRef(
                geometric.first_atom_id,
                tuple(-value for value in translation),
            ),
            tuple(-value for value in geometric.vector_cartesian),
        )
    raise ValueError("coordination-shell contact does not contain its centre")


def _polyhedron_vertices(
    shell: CoordinationShell,
    view: AtomicView,
) -> tuple[PolyhedronVertex, ...]:
    atoms = {atom.id: atom for atom in view.atoms}
    if shell.center_atom_id not in atoms:
        raise ValueError("coordination-shell centre is absent from the atomic view")
    vertices: list[PolyhedronVertex] = []
    for contact in shell.contacts:
        atom_ref, vector = _oriented_contact(shell.center_atom_id, contact)
        ligand = atoms.get(atom_ref.atom_id)
        if ligand is None:
            raise ValueError("coordination-shell ligand is absent from the atomic view")
        occupancy = math.fsum(
            float(component.occupancy.value) for component in ligand.components
        )
        vertices.append(PolyhedronVertex(
            atom_ref=atom_ref,
            incidence_orbit_id=(
                "legacy-incidence:" + contact.geometric_contact.contact_id
            ),
            resolved_contact_orbit_id=contact.contact_orbit_id,
            local_cartesian=tuple(float(value) for value in vector),
            distance=float(contact.geometric_contact.distance),
            occupancy=occupancy,
        ))
    return tuple(vertices)


def _local_vertices(vertices: tuple[PolyhedronVertex, ...]) -> np.ndarray:
    return np.asarray(
        tuple(item.local_cartesian for item in vertices), dtype=float
    ).reshape((-1, 3))


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


def _polyhedron_id(center_atom_id: str, vertices: tuple[PolyhedronVertex, ...]) -> str:
    payload = json.dumps(
        (
            center_atom_id,
            tuple(sorted(
                (
                    item.atom_ref.atom_id,
                    item.atom_ref.cell_translation,
                    item.incidence_orbit_id,
                )
                for item in vertices
            )),
        ),
        ensure_ascii=True,
        separators=(",", ":"),
    ).encode("ascii")
    return f"polyhedron:{hashlib.sha256(payload).hexdigest()}"


def _ligand_composition(
    vertices: tuple[PolyhedronVertex, ...],
    view: AtomicView,
) -> tuple[tuple[str, float], ...]:
    atoms = {atom.id: atom for atom in view.atoms}
    totals: dict[str, list[float]] = {}
    for vertex in vertices:
        for component in atoms[vertex.atom_ref.atom_id].components:
            if float(component.occupancy.value) <= 0:
                continue
            label = component.element or component.species.label
            totals.setdefault(label, []).append(float(component.occupancy.value))
    return tuple(
        (label, math.fsum(values)) for label, values in sorted(totals.items())
    )


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
                "crystal_chemistry.polyhedron.open_coordination_shell",
                f"polyhedron cannot be completed from a {shell.status.value} shell",
            )
        vertices = _polyhedron_vertices(shell, view)
        local = _local_vertices(vertices)
        distances = tuple(item.distance for item in vertices)
        mean, minimum, maximum, distortion = _bond_length_metrics(distances)
        diagnostics: list[Diagnostic] = []
        atoms = {atom.id: atom for atom in view.atoms}
        if any(item.occupancy < 1.0 - 1e-12 for item in vertices):
            diagnostics.append(Diagnostic(
                Severity.WARNING,
                PARTIAL_OCCUPANCY,
                "Polyhedron contains a partially occupied ligand vertex",
            ))
        if any(len(atoms[item.atom_ref.atom_id].components) > 1 for item in vertices):
            diagnostics.append(Diagnostic(
                Severity.WARNING,
                MIXED_POSITION,
                "Polyhedron contains a mixed ligand position",
            ))

        faces: tuple[tuple[int, ...], ...] = ()
        face_signature: FaceSignature | None = None
        angle_dispersion: float | None = None
        volume: float | None = None
        centroid: tuple[float, float, float] | None = None
        center_offset: float | None = None
        status = ResolutionStatus.RESOLVED
        is_three_dimensional = (
            len(local) >= 4
            and np.linalg.matrix_rank(local[1:] - local[0], tol=self.tolerance) >= 3
        )
        if is_three_dimensional:
            try:
                faces = _convex_hull_faces(local, self.tolerance)
                face_signature = canonical_face_signature(len(vertices), faces)
                angle_dispersion = _edge_angle_dispersion_deg(
                    local, unique_hull_edges(faces)
                )
                calculated_volume, calculated_centroid = _volume_and_centroid(local, faces)
                volume = calculated_volume
                centroid = tuple(float(value) for value in calculated_centroid)
                center_offset = float(np.linalg.norm(calculated_centroid))
            except ValueError as error:
                status = ResolutionStatus.INCOMPLETE
                diagnostics.append(Diagnostic(
                    Severity.WARNING,
                    GEOMETRY_DEGENERATE,
                    str(error),
                ))
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

        polyhedron = CoordinationPolyhedron(
            polyhedron_id=_polyhedron_id(shell.center_atom_id, vertices),
            polyhedron_orbit_id="",
            source_site_id=shell.source_site_id,
            center_atom_id=shell.center_atom_id,
            coordination_number=len(vertices),
            ligand_composition=_ligand_composition(vertices, view),
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
                ("method", "cristma.polyhedron_builder:2"),
                ("tolerance", self.tolerance),
                ("bond_distortion", "Baur mean absolute relative deviation"),
                ("edge_angle_dispersion", "population standard deviation in degrees"),
            ),
            local_vertices=tuple(item.local_cartesian for item in vertices),
            shell_provenance=shell.provenance,
            diagnostics=tuple(dict.fromkeys(diagnostics)),
        )
        return PolyhedronBuildResult.from_polyhedron(polyhedron)

    def validate_orbit(
        self,
        polyhedra: tuple[CoordinationPolyhedron, ...],
    ) -> tuple[Diagnostic, ...]:
        signatures = {item.face_signature for item in polyhedra}
        if len(signatures) <= 1:
            return ()
        return (Diagnostic(
            Severity.ERROR,
            "crystal_chemistry.polyhedron.symmetry_inconsistent",
            "symmetry-equivalent centres have different face incidence graphs",
        ),)


__all__ = [
    "CoordinationPolyhedron",
    "CoordinationPolyhedronOrbit",
    "FaceSignature",
    "PolyhedronBuildResult",
    "PolyhedronBuilder",
    "PolyhedronVertex",
    "canonical_face_signature",
    "polyhedron_face_signature",
    "unique_hull_edges",
]
