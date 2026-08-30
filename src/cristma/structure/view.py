"""Shared immutable numerical view over canonical structures."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cristma.chemistry.species import ChemicalSpecies

from .properties import AtomicPropertyTable


def _immutable_float_array(value: np.ndarray | object) -> np.ndarray:
    array = np.array(value, dtype=float, copy=True)
    array.flags.writeable = False
    return array


@dataclass(frozen=True, slots=True)
class AtomicView:
    """Application-neutral atom coordinates, identity, and properties."""

    ids: tuple[str, ...]
    species: tuple[ChemicalSpecies, ...]
    cartesian: np.ndarray
    fractional: np.ndarray | None
    cell: np.ndarray | None
    periodic: tuple[bool, bool, bool]
    properties: AtomicPropertyTable
    source_site_ids: tuple[str | None, ...]

    def __post_init__(self) -> None:
        count = len(self.ids)
        if len(self.species) != count or len(self.source_site_ids) != count:
            raise ValueError("atomic view identity arrays must have equal lengths")
        cartesian = _immutable_float_array(self.cartesian).reshape((-1, 3))
        if cartesian.shape != (count, 3):
            raise ValueError("Cartesian coordinates must have shape (atom_count, 3)")
        object.__setattr__(self, "cartesian", cartesian)
        if self.fractional is not None:
            fractional = _immutable_float_array(self.fractional).reshape((-1, 3))
            if fractional.shape != (count, 3):
                raise ValueError("fractional coordinates must have shape (atom_count, 3)")
            object.__setattr__(self, "fractional", fractional)
        if self.cell is not None:
            cell = _immutable_float_array(self.cell)
            if cell.shape != (3, 3):
                raise ValueError("cell matrix must have shape (3, 3)")
            object.__setattr__(self, "cell", cell)
        if self.properties.atom_count != count:
            raise ValueError("property table atom count does not match atomic view")


__all__ = ["AtomicView"]
