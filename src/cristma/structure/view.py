"""Shared immutable numerical view over canonical atomic positions."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Generic, TypeVar

import numpy as np

from cristma.core.cell import UnitCell

from .position import AtomicPosition
from .properties import AtomicPropertyTable


TAtom = TypeVar("TAtom", bound=AtomicPosition)


def _immutable_rows(value: object) -> np.ndarray:
    array = np.array(value, dtype=float, copy=True).reshape((-1, 3))
    array.flags.writeable = False
    return array


@dataclass(frozen=True, slots=True)
class AtomicView(Generic[TAtom]):
    """Finite atomic rows plus derived read-only numerical coordinate arrays."""

    atoms: tuple[TAtom, ...]
    cell: UnitCell | None
    periodic: tuple[bool, bool, bool]
    properties: AtomicPropertyTable
    cartesian: np.ndarray = field(init=False, repr=False, compare=False)
    fractional: np.ndarray | None = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        ids = tuple(atom.id for atom in self.atoms)
        if len(set(ids)) != len(ids):
            raise ValueError("atomic view atom IDs must be unique")
        if len(self.periodic) != 3:
            raise ValueError("periodicity must contain three axes")

        cartesian_rows = tuple(
            tuple(float(value) for value in atom.cartesian) for atom in self.atoms
        )
        if any(
            len(row) != 3 or not all(math.isfinite(value) for value in row)
            for row in cartesian_rows
        ):
            raise ValueError("Cartesian coordinates must contain three finite values per atom")
        object.__setattr__(self, "cartesian", _immutable_rows(cartesian_rows))

        fractional_rows = tuple(getattr(atom, "fractional", None) for atom in self.atoms)
        fractional: np.ndarray | None = None
        if self.atoms and all(row is not None for row in fractional_rows):
            numeric_rows = tuple(
                tuple(float(value) for value in row) for row in fractional_rows
            )
            if any(
                len(row) != 3 or not all(math.isfinite(value) for value in row)
                for row in numeric_rows
            ):
                raise ValueError("fractional coordinates must contain three finite values per atom")
            fractional = _immutable_rows(numeric_rows)
        elif not self.atoms and any(self.periodic) and self.cell is not None:
            fractional = _immutable_rows(())
        object.__setattr__(self, "fractional", fractional)

        if any(self.periodic) and (self.cell is None or fractional is None):
            raise ValueError("periodic atomic view requires a cell and fractional coordinates")
        if self.properties.atom_count != len(self.atoms):
            raise ValueError("property table atom count does not match atomic view")

    @property
    def ids(self) -> tuple[str, ...]:
        return tuple(atom.id for atom in self.atoms)

    @property
    def cell_matrix(self) -> np.ndarray | None:
        return None if self.cell is None else self.cell.matrix


__all__ = ["AtomicView"]
