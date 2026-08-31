"""Immutable source records shared by native VASP readers."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Literal

import numpy as np

from cristma.chemistry.species import ChemicalSpecies, as_species
from cristma.structure import SourceReference


CoordinateMode = Literal["direct", "cartesian"]


def _array(value: object, *, dtype: object = float) -> np.ndarray:
    result = np.array(value, dtype=dtype, copy=True)
    result.flags.writeable = False
    return result


@dataclass(frozen=True, slots=True)
class VaspScale:
    """One universal or three Cartesian-component VASP scale values."""

    values: tuple[float, ...]

    def __post_init__(self) -> None:
        values = tuple(float(value) for value in self.values)
        if len(values) not in {1, 3}:
            raise ValueError("VASP scale requires one or three numbers")
        if not all(math.isfinite(value) for value in values):
            raise ValueError("VASP scale values must be finite")
        if len(values) == 1 and values[0] == 0:
            raise ValueError("VASP scale must not be zero")
        if len(values) == 3 and not all(value > 0 for value in values):
            raise ValueError("three VASP scale values must be positive")
        object.__setattr__(self, "values", values)


@dataclass(frozen=True, slots=True)
class VaspAtomRow:
    """One reported coordinate row and optional relaxation flags."""

    coordinates: tuple[float, float, float]
    selective: tuple[bool, bool, bool] | None = None
    start_offset: int | None = None
    end_offset: int | None = None


@dataclass(frozen=True, slots=True)
class VaspHeader:
    """POSCAR-like lattice, populations, and coordinate convention."""

    title: str
    scale: VaspScale
    raw_lattice: np.ndarray
    species_labels: tuple[str, ...] | None
    counts: tuple[int, ...]
    coordinate_mode: CoordinateMode
    selective_dynamics: bool = False

    def __post_init__(self) -> None:
        lattice = _array(self.raw_lattice)
        if lattice.shape != (3, 3) or not np.isfinite(lattice).all():
            raise ValueError("raw VASP lattice must be a finite 3x3 array")
        if not self.counts or any(
            isinstance(value, bool) or not isinstance(value, int) or value < 0
            for value in self.counts
        ):
            raise ValueError("VASP counts must contain non-negative integers")
        if self.species_labels is not None and len(self.species_labels) != len(self.counts):
            raise ValueError("VASP species and count columns must have equal length")
        object.__setattr__(self, "raw_lattice", lattice)


@dataclass(frozen=True, slots=True)
class VaspFrameSpan:
    """Cheap index entry for one complete structure frame."""

    index: int
    start_offset: int
    end_offset: int
    reported_index: int | None = None

    def __post_init__(self) -> None:
        if self.index < 0 or self.start_offset < 0 or self.end_offset < self.start_offset:
            raise ValueError("invalid VASP frame span")


@dataclass(frozen=True, slots=True)
class VaspSnapshot:
    """Format-neutral numerical payload ready for canonical structure mapping."""

    name: str
    lattice: np.ndarray
    species: tuple[ChemicalSpecies | str, ...]
    fractional: np.ndarray
    frame_index: int
    source: SourceReference
    selective_dynamics: np.ndarray | None = None
    velocities: np.ndarray | None = None
    velocity_mode: CoordinateMode | None = None
    velocity_unit: str | None = None
    forces: np.ndarray | None = None
    force_unit: str | None = None

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("VASP snapshot name must not be empty")
        if isinstance(self.frame_index, bool) or self.frame_index < 0:
            raise ValueError("VASP frame index must be non-negative")
        lattice = _array(self.lattice)
        fractional = _array(self.fractional)
        if lattice.shape != (3, 3) or not np.isfinite(lattice).all():
            raise ValueError("snapshot lattice must be a finite 3x3 array")
        if abs(float(np.linalg.det(lattice))) <= 1e-15:
            raise ValueError("snapshot lattice must be non-singular")
        if fractional.ndim != 2 or fractional.shape[1:] != (3,) or not np.isfinite(fractional).all():
            raise ValueError("fractional coordinates must be a finite Nx3 array")
        species = tuple(as_species(value) for value in self.species)
        if len(species) != fractional.shape[0]:
            raise ValueError("species count must match fractional coordinate rows")
        object.__setattr__(self, "lattice", lattice)
        object.__setattr__(self, "fractional", fractional)
        object.__setattr__(self, "species", species)

        selective = self.selective_dynamics
        if selective is not None:
            selective = _array(selective, dtype=bool)
            if selective.shape != fractional.shape:
                raise ValueError("selective_dynamics must be an Nx3 array")
            object.__setattr__(self, "selective_dynamics", selective)

        velocities_present = self.velocities is not None
        convention_present = self.velocity_mode is not None and self.velocity_unit is not None
        if velocities_present != convention_present:
            raise ValueError("velocity mode and unit must accompany velocity values")
        if velocities_present:
            if self.velocity_mode not in {"direct", "cartesian"}:
                raise ValueError("unknown velocity mode")
            velocities = _array(self.velocities)
            if velocities.shape != fractional.shape or not np.isfinite(velocities).all():
                raise ValueError("velocities must be a finite Nx3 array")
            object.__setattr__(self, "velocities", velocities)

        forces_present = self.forces is not None
        if forces_present != (self.force_unit is not None):
            raise ValueError("force unit must accompany force values")
        if forces_present:
            forces = _array(self.forces)
            if forces.shape != fractional.shape or not np.isfinite(forces).all():
                raise ValueError("forces must be a finite Nx3 array")
            object.__setattr__(self, "forces", forces)


@dataclass(frozen=True, slots=True)
class PoscarDocument:
    raw_source: str
    source_name: str | None = None
    header: VaspHeader | None = None
    positions: tuple[VaspAtomRow, ...] = ()
    velocity_mode: CoordinateMode | None = None
    velocities: tuple[VaspAtomRow, ...] = ()
    trailing_start: int | None = None

    def render_preserved(self) -> str:
        return self.raw_source


@dataclass(frozen=True, slots=True)
class XdatcarDocument:
    raw_source: str
    source_name: str | None = None
    header: VaspHeader | None = None
    frames: tuple[VaspFrameSpan, ...] = ()


@dataclass(frozen=True, slots=True)
class OutcarDocument:
    raw_source: str
    source_name: str | None = None
    frames: tuple[VaspFrameSpan, ...] = ()
    species_labels: tuple[str, ...] | None = None
    counts: tuple[int, ...] = ()
    lattices: tuple[np.ndarray, ...] = field(default=(), compare=False)

    def __post_init__(self) -> None:
        if len(self.lattices) != len(self.frames):
            raise ValueError("OUTCAR lattice count must match frame count")
        lattices = tuple(_array(value) for value in self.lattices)
        if any(value.shape != (3, 3) for value in lattices):
            raise ValueError("OUTCAR lattices must be 3x3 arrays")
        object.__setattr__(self, "lattices", lattices)


@dataclass(frozen=True, slots=True)
class VasprunDocument:
    raw_source: str
    source_name: str | None = None
    frames: tuple[VaspFrameSpan, ...] = ()
    species_labels: tuple[str, ...] = ()


__all__ = [
    "CoordinateMode",
    "OutcarDocument",
    "PoscarDocument",
    "VaspAtomRow",
    "VaspFrameSpan",
    "VaspHeader",
    "VaspScale",
    "VaspSnapshot",
    "VasprunDocument",
    "XdatcarDocument",
]
