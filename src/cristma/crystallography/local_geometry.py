"""Canonical geometric contacts derived from finite neighbour graphs."""

from __future__ import annotations

from dataclasses import dataclass
import math

from cristma.geometry.neighbors import (
    Neighbor,
    NeighborGraph,
    PeriodicNeighbor,
    PeriodicNeighborGraph,
)
from cristma.structure.position import AtomicPosition
from cristma.structure.view import AtomicView


Translation = tuple[int, int, int]
CanonicalKey = tuple[str, str, Translation]


@dataclass(frozen=True, slots=True)
class GeometricContact:
    """One finite contact, independent of directed graph traversal."""

    contact_id: str
    first_atom_id: str
    second_atom_id: str
    cell_translation: Translation | None
    distance: float
    vector_cartesian: tuple[float, float, float]
    first_source_site_id: str | None
    second_source_site_id: str | None
    geometric_provenance: str

    def __post_init__(self) -> None:
        if not self.contact_id or not self.first_atom_id or not self.second_atom_id:
            raise ValueError("contact and endpoint IDs must not be empty")
        if not math.isfinite(self.distance) or self.distance <= 0:
            raise ValueError("contact distance must be finite and positive")
        if len(self.vector_cartesian) != 3 or not all(
            math.isfinite(value) for value in self.vector_cartesian
        ):
            raise ValueError("contact vector must contain three finite values")


def _negative(translation: Translation) -> Translation:
    return tuple(-value for value in translation)


def _contact_id(key: CanonicalKey, *, periodic: bool) -> str:
    first, second, translation = key
    suffix = ",".join(str(value) for value in translation) if periodic else "finite"
    return f"contact:{first}|{second}|{suffix}"


def geometric_contacts(
    view: AtomicView[AtomicPosition],
    graph: NeighborGraph | PeriodicNeighborGraph,
) -> tuple[GeometricContact, ...]:
    """Collapse reverse directed edges into deterministic physical contacts."""

    view_atoms = {atom.id: atom for atom in view.atoms}
    graph_ids = {atom.id for atom in graph.atoms}
    if set(view_atoms) != graph_ids or len(view_atoms) != len(graph.atoms):
        raise ValueError("geometric contact graph and atomic view must contain the same atoms")

    periodic = isinstance(graph, PeriodicNeighborGraph)
    contacts: dict[CanonicalKey, GeometricContact] = {}
    for edge in graph.edges:
        if isinstance(edge, PeriodicNeighbor):
            translation = edge.target.cell_translation
            target_atom_id = edge.target.atom_id
        elif isinstance(edge, Neighbor):
            translation = (0, 0, 0)
            target_atom_id = edge.target_atom_id
        else:
            raise TypeError(f"unsupported neighbour edge: {type(edge)!r}")

        forward: CanonicalKey = (edge.source_atom_id, target_atom_id, translation)
        reverse: CanonicalKey = (target_atom_id, edge.source_atom_id, _negative(translation))
        key = min(forward, reverse)
        forward_wins = key == forward
        first_atom_id, second_atom_id, canonical_translation = key
        vector = (
            edge.vector_cartesian
            if forward_wins
            else tuple(-value for value in edge.vector_cartesian)
        )
        first_atom = view_atoms[first_atom_id]
        second_atom = view_atoms[second_atom_id]
        contacts[key] = GeometricContact(
            contact_id=_contact_id(key, periodic=periodic),
            first_atom_id=first_atom_id,
            second_atom_id=second_atom_id,
            cell_translation=canonical_translation if periodic else None,
            distance=edge.distance,
            vector_cartesian=tuple(float(value) for value in vector),
            first_source_site_id=getattr(first_atom, "source_site_id", None),
            second_source_site_id=getattr(second_atom, "source_site_id", None),
            geometric_provenance=(
                "periodic_neighbor_graph" if periodic else "neighbor_graph"
            ),
        )
    return tuple(contacts[key] for key in sorted(contacts))


__all__ = ["GeometricContact", "geometric_contacts"]
