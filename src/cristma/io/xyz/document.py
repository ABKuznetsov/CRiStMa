"""Immutable source and selected-frame records for XYZ/extXYZ."""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal, Mapping

import numpy as np

from cristma.structure import SourceReference


XyzPropertyKind = Literal["S", "I", "R", "L"]


@dataclass(frozen=True, slots=True)
class XyzPropertySpec:
    """One extXYZ per-atom property declaration."""

    name: str
    kind: XyzPropertyKind
    width: int

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("XYZ property name must not be empty")
        if self.kind not in {"S", "I", "R", "L"}:
            raise ValueError(f"unknown XYZ property type: {self.kind!r}")
        if isinstance(self.width, bool) or not isinstance(self.width, int) or self.width <= 0:
            raise ValueError("XYZ property width must be a positive integer")


@dataclass(frozen=True, slots=True)
class XyzFrameSpan:
    """Character offsets and declared size of one complete XYZ frame."""

    index: int
    atom_count: int
    start_offset: int
    end_offset: int
    comment_start_offset: int
    comment_end_offset: int
    atom_rows_start_offset: int

    def __post_init__(self) -> None:
        if isinstance(self.index, bool) or not isinstance(self.index, int) or self.index < 0:
            raise ValueError("XYZ frame index must be a non-negative integer")
        if (
            isinstance(self.atom_count, bool)
            or not isinstance(self.atom_count, int)
            or self.atom_count < 0
        ):
            raise ValueError("XYZ atom count must be a non-negative integer")
        offsets = (
            self.start_offset,
            self.comment_start_offset,
            self.comment_end_offset,
            self.atom_rows_start_offset,
            self.end_offset,
        )
        if any(value < 0 for value in offsets) or tuple(sorted(offsets)) != offsets:
            raise ValueError("XYZ frame offsets must be ordered")


@dataclass(frozen=True, slots=True)
class XyzDocument:
    """Loss-preserving XYZ source plus a cheap complete-frame index."""

    raw_source: str
    source_name: str | None = None
    frames: tuple[XyzFrameSpan, ...] = ()

    def __post_init__(self) -> None:
        previous_end = 0
        for expected_index, frame in enumerate(self.frames):
            if frame.index != expected_index or frame.start_offset < previous_end:
                raise ValueError("XYZ frame spans must be ordered and consecutively indexed")
            if frame.end_offset > len(self.raw_source):
                raise ValueError("XYZ frame span exceeds source length")
            previous_end = frame.end_offset

    def render_preserved(self) -> str:
        return self.raw_source


def _immutable_array(value: object) -> np.ndarray:
    array = np.array(value, copy=True)
    array.flags.writeable = False
    return array


@dataclass(frozen=True, slots=True)
class XyzFrame:
    """One parsed frame ready for canonical scientific mapping."""

    name: str
    atom_count: int
    comment: str
    metadata: Mapping[str, object]
    schema: tuple[XyzPropertySpec, ...]
    columns: Mapping[str, np.ndarray]
    lattice: np.ndarray | None
    pbc: tuple[bool, bool, bool] | None
    source: SourceReference

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("XYZ frame name must not be empty")
        if (
            isinstance(self.atom_count, bool)
            or not isinstance(self.atom_count, int)
            or self.atom_count < 0
        ):
            raise ValueError("XYZ atom count must be a non-negative integer")
        columns = {}
        for name, values in self.columns.items():
            array = _immutable_array(values)
            if array.ndim == 0 or array.shape[0] != self.atom_count:
                raise ValueError(f"XYZ property {name!r} row count does not match atom count")
            columns[name] = array
        object.__setattr__(self, "columns", MappingProxyType(columns))
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))
        if self.lattice is not None:
            lattice = _immutable_array(self.lattice)
            if lattice.shape != (3, 3) or not np.isfinite(lattice).all():
                raise ValueError("XYZ lattice must be a finite 3x3 array")
            if abs(float(np.linalg.det(lattice))) <= 1e-15:
                raise ValueError("XYZ lattice must be non-singular")
            object.__setattr__(self, "lattice", lattice)
        if self.pbc is not None and (
            len(self.pbc) != 3 or any(not isinstance(value, bool) for value in self.pbc)
        ):
            raise ValueError("XYZ pbc must contain three booleans")


__all__ = [
    "XyzDocument",
    "XyzFrame",
    "XyzFrameSpan",
    "XyzPropertyKind",
    "XyzPropertySpec",
]
