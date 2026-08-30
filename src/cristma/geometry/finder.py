"""Configurable finite and periodic neighbor search."""

from __future__ import annotations

from dataclasses import dataclass, replace
import itertools
import math
from typing import overload

import numpy as np

from cristma.diagnostics import Diagnostic, Severity
from cristma.structure.identity import ExpandedAtom, PeriodicAtomRef
from cristma.structure.molecular import MolecularAtom
from cristma.structure.position import AtomicPosition
from cristma.structure.view import AtomicView

from .neighbors import Neighbor, NeighborGraph, PeriodicNeighbor, PeriodicNeighborGraph


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

    @overload
    def find(self, view: AtomicView[MolecularAtom]) -> NeighborGraph[MolecularAtom]: ...

    @overload
    def find(self, view: AtomicView[ExpandedAtom]) -> PeriodicNeighborGraph[ExpandedAtom]: ...

    def find(
        self,
        view: AtomicView[AtomicPosition],
    ) -> NeighborGraph[AtomicPosition] | PeriodicNeighborGraph[AtomicPosition]:
        if any(view.periodic):
            return self._find_periodic(view)

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

    def _find_periodic(
        self,
        view: AtomicView[AtomicPosition],
    ) -> PeriodicNeighborGraph[AtomicPosition]:
        if view.cell is None or view.fractional is None:
            raise ValueError("periodic neighbor search requires cell and fractional coordinates")
        matrix = view.cell.matrix
        inverse = np.linalg.inv(matrix)
        component_bounds = tuple(
            self.cutoff * float(np.linalg.norm(inverse[:, axis]))
            for axis in range(3)
        )
        edges: dict[tuple[str, str, tuple[int, int, int]], PeriodicNeighbor] = {}
        diagnostics: list[Diagnostic] = []
        coincident: set[tuple[str, str, tuple[int, int, int]]] = set()

        for source_index, source in enumerate(view.atoms):
            source_fractional = view.fractional[source_index]
            for target_index, target in enumerate(view.atoms):
                delta = view.fractional[target_index] - source_fractional
                ranges: list[tuple[int, ...]] = []
                for axis, periodic in enumerate(view.periodic):
                    if not periodic:
                        ranges.append((0,))
                        continue
                    lower = math.ceil(
                        -float(delta[axis]) - component_bounds[axis] - self.tolerance
                    )
                    upper = math.floor(
                        -float(delta[axis]) + component_bounds[axis] + self.tolerance
                    )
                    ranges.append(tuple(range(lower, upper + 1)))

                for translation_values in itertools.product(*ranges):
                    translation = tuple(int(value) for value in translation_values)
                    if source.id == target.id and translation == (0, 0, 0):
                        continue
                    vector_fractional = delta + np.asarray(translation, dtype=float)
                    vector_array = vector_fractional @ matrix
                    distance = float(np.linalg.norm(vector_array))
                    key = (source.id, target.id, translation)
                    if distance <= self.tolerance:
                        if key not in coincident:
                            diagnostics.append(
                                Diagnostic(
                                    Severity.WARNING,
                                    "geometry.coincident_positions",
                                    f"periodic positions {source.id!r} and {target.id!r} coincide",
                                )
                            )
                            coincident.add(key)
                        continue
                    if distance <= self.cutoff + self.tolerance:
                        edges[key] = PeriodicNeighbor(
                            source_atom_id=source.id,
                            target=PeriodicAtomRef(target.id, translation),
                            distance=distance,
                            vector_cartesian=tuple(float(value) for value in vector_array),
                        )
        return PeriodicNeighborGraph(
            view.atoms,
            tuple(edges[key] for key in sorted(edges)),
            tuple(diagnostics),
        )


__all__ = ["NeighborFinder"]
