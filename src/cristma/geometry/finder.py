"""Configurable finite and periodic neighbor search."""

from __future__ import annotations

from dataclasses import dataclass, replace
import math

import numpy as np

from cristma.diagnostics import Diagnostic, Severity
from cristma.structure.position import AtomicPosition
from cristma.structure.view import AtomicView

from .neighbors import Neighbor, NeighborGraph


@dataclass(frozen=True, slots=True)
class NeighborFinder:
    cutoff: float
    tolerance: float = 1e-12

    def __post_init__(self) -> None:
        if not math.isfinite(self.cutoff) or self.cutoff <= 0:
            raise ValueError("cutoff must be positive and finite")
        if not math.isfinite(self.tolerance) or self.tolerance <= 0:
            raise ValueError("tolerance must be positive and finite")

    def get_config(self) -> dict[str, float]:
        return {"cutoff": self.cutoff, "tolerance": self.tolerance}

    def clone(self, **changes: float) -> NeighborFinder:
        return replace(self, **changes)

    def find(self, view: AtomicView[AtomicPosition]) -> NeighborGraph[AtomicPosition]:
        if any(view.periodic):
            raise NotImplementedError("periodic neighbor search is implemented separately")

        edges: list[Neighbor] = []
        diagnostics: list[Diagnostic] = []
        for left_index, left in enumerate(view.atoms):
            for right in view.atoms[left_index + 1 :]:
                vector_array = np.asarray(right.cartesian) - np.asarray(left.cartesian)
                distance = float(np.linalg.norm(vector_array))
                if distance <= self.tolerance:
                    diagnostics.append(
                        Diagnostic(
                            Severity.WARNING,
                            "geometry.coincident_positions",
                            f"positions {left.id!r} and {right.id!r} coincide",
                        )
                    )
                    continue
                if distance <= self.cutoff + self.tolerance:
                    vector = tuple(float(value) for value in vector_array)
                    reverse = tuple(-value for value in vector)
                    edges.extend(
                        (
                            Neighbor(left.id, right.id, distance, vector),
                            Neighbor(right.id, left.id, distance, reverse),
                        )
                    )
        return NeighborGraph(view.atoms, tuple(edges), tuple(diagnostics))


__all__ = ["NeighborFinder"]
