"""Immutable finite neighbor graph records."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from types import MappingProxyType
from typing import Generic, Mapping, Protocol, TypeVar

from cristma.diagnostics import Diagnostic
from cristma.structure.position import AtomicPosition
from cristma.structure.identity import PeriodicAtomRef


TAtom = TypeVar("TAtom", bound=AtomicPosition)
TNeighbor = TypeVar("TNeighbor", bound="NeighborLike", covariant=True)


class NeighborLike(Protocol):
    source_atom_id: str
    target_atom_id: str
    distance: float
    vector_cartesian: tuple[float, float, float]


class NeighborGraphLike(Protocol[TAtom, TNeighbor]):
    atoms: tuple[TAtom, ...]

    def neighbors(self, atom_id: str) -> tuple[TNeighbor, ...]: ...


@dataclass(frozen=True, slots=True)
class Neighbor:
    source_atom_id: str
    target_atom_id: str
    distance: float
    vector_cartesian: tuple[float, float, float]

    def __post_init__(self) -> None:
        if not self.source_atom_id or not self.target_atom_id:
            raise ValueError("neighbor endpoint IDs must not be empty")
        if not math.isfinite(self.distance) or self.distance <= 0:
            raise ValueError("neighbor distance must be finite and positive")
        if len(self.vector_cartesian) != 3 or not all(
            math.isfinite(value) for value in self.vector_cartesian
        ):
            raise ValueError("neighbor vector must contain three finite values")


@dataclass(frozen=True, slots=True)
class NeighborGraph(Generic[TAtom]):
    atoms: tuple[TAtom, ...]
    edges: tuple[Neighbor, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    _adjacency: Mapping[str, tuple[Neighbor, ...]] = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        atom_ids = tuple(atom.id for atom in self.atoms)
        if len(set(atom_ids)) != len(atom_ids):
            raise ValueError("neighbor graph atom IDs must be unique")
        known = set(atom_ids)
        rows: dict[str, list[Neighbor]] = {atom_id: [] for atom_id in atom_ids}
        for edge in self.edges:
            if edge.source_atom_id not in known or edge.target_atom_id not in known:
                raise ValueError("neighbor edge references an unknown atom")
            rows[edge.source_atom_id].append(edge)
        adjacency = {
            atom_id: tuple(sorted(row, key=lambda edge: (edge.distance, edge.target_atom_id)))
            for atom_id, row in rows.items()
        }
        object.__setattr__(self, "_adjacency", MappingProxyType(adjacency))

    def neighbors(self, atom_id: str) -> tuple[Neighbor, ...]:
        return self._adjacency[atom_id]


@dataclass(frozen=True, slots=True)
class PeriodicNeighbor:
    source_atom_id: str
    target: PeriodicAtomRef
    distance: float
    vector_cartesian: tuple[float, float, float]

    def __post_init__(self) -> None:
        Neighbor(
            self.source_atom_id,
            self.target.atom_id,
            self.distance,
            self.vector_cartesian,
        )

    @property
    def target_atom_id(self) -> str:
        return self.target.atom_id


@dataclass(frozen=True, slots=True)
class PeriodicNeighborGraph(Generic[TAtom]):
    atoms: tuple[TAtom, ...]
    edges: tuple[PeriodicNeighbor, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    _adjacency: Mapping[str, tuple[PeriodicNeighbor, ...]] = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        atom_ids = tuple(atom.id for atom in self.atoms)
        if len(set(atom_ids)) != len(atom_ids):
            raise ValueError("periodic neighbor graph atom IDs must be unique")
        known = set(atom_ids)
        rows: dict[str, list[PeriodicNeighbor]] = {atom_id: [] for atom_id in atom_ids}
        for edge in self.edges:
            if edge.source_atom_id not in known or edge.target.atom_id not in known:
                raise ValueError("periodic neighbor edge references an unknown atom")
            rows[edge.source_atom_id].append(edge)
        adjacency = {
            atom_id: tuple(
                sorted(
                    row,
                    key=lambda edge: (
                        edge.distance,
                        edge.target.atom_id,
                        edge.target.cell_translation,
                    ),
                )
            )
            for atom_id, row in rows.items()
        }
        object.__setattr__(self, "_adjacency", MappingProxyType(adjacency))

    def neighbors(self, atom_id: str) -> tuple[PeriodicNeighbor, ...]:
        return self._adjacency[atom_id]


__all__ = [
    "Neighbor",
    "NeighborGraph",
    "NeighborGraphLike",
    "NeighborLike",
    "PeriodicNeighbor",
    "PeriodicNeighborGraph",
]
