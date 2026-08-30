"""Graph-neutral coordination environments."""

from __future__ import annotations

from dataclasses import dataclass

from cristma.diagnostics import Diagnostic
from cristma.structure.occupation import SiteComponent
from cristma.structure.position import AtomicPosition
from cristma.structure.view import AtomicView

from .neighbors import Neighbor, NeighborGraphLike, PeriodicNeighbor


@dataclass(frozen=True, slots=True)
class CoordinationEnvironment:
    center_atom_id: str
    center_components: tuple[SiteComponent, ...]
    neighbors: tuple[Neighbor | PeriodicNeighbor, ...]

    @property
    def coordination_number(self) -> int:
        return len(self.neighbors)


@dataclass(frozen=True, slots=True)
class CoordinationResult:
    environments: tuple[CoordinationEnvironment, ...]
    diagnostics: tuple[Diagnostic, ...] = ()

    def by_atom(self, atom_id: str) -> CoordinationEnvironment:
        for environment in self.environments:
            if environment.center_atom_id == atom_id:
                return environment
        raise KeyError(atom_id)


@dataclass(frozen=True, slots=True)
class CoordinationAnalyzer:
    """Convert an already selected neighbor graph into coordination records."""

    def analyze(
        self,
        view: AtomicView[AtomicPosition],
        graph: NeighborGraphLike,
    ) -> CoordinationResult:
        view_ids = tuple(atom.id for atom in view.atoms)
        graph_ids = tuple(atom.id for atom in graph.atoms)
        if set(view_ids) != set(graph_ids) or len(view_ids) != len(graph_ids):
            raise ValueError("coordination graph and atomic view must contain the same atoms")
        environments = tuple(
            CoordinationEnvironment(
                center_atom_id=atom.id,
                center_components=atom.components,
                neighbors=graph.neighbors(atom.id),
            )
            for atom in view.atoms
        )
        return CoordinationResult(
            environments,
            tuple(getattr(graph, "diagnostics", ())),
        )


__all__ = [
    "CoordinationAnalyzer",
    "CoordinationEnvironment",
    "CoordinationResult",
]
