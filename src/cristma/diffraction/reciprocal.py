"""Reciprocal metric calculations and bounded integer ellipsoid search."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterator

import numpy as np

from cristma.core import UnitCell

from .models import MillerIndex


def _positive_finite(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a real number")
    converted = float(value)
    if not math.isfinite(converted) or converted <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return converted


@dataclass(frozen=True, slots=True)
class ReciprocalMetric:
    """Reciprocal basis and metric without a ``2*pi`` factor."""

    basis: np.ndarray
    metric: np.ndarray

    def __post_init__(self) -> None:
        for value, name in ((self.basis, "basis"), (self.metric, "metric")):
            if not isinstance(value, np.ndarray) or value.shape != (3, 3):
                raise ValueError(f"reciprocal {name} must be a 3x3 array")
            if not np.all(np.isfinite(value)):
                raise ValueError(f"reciprocal {name} must be finite")
        if not np.allclose(self.metric, self.metric.T, rtol=0.0, atol=1e-14):
            raise ValueError("reciprocal metric must be symmetric")
        if np.min(np.linalg.eigvalsh(self.metric)) <= 0:
            raise ValueError("reciprocal metric must be positive definite")
        basis = np.array(self.basis, dtype=float, copy=True)
        metric = np.array(self.metric, dtype=float, copy=True)
        basis.flags.writeable = False
        metric.flags.writeable = False
        object.__setattr__(self, "basis", basis)
        object.__setattr__(self, "metric", metric)

    @classmethod
    def from_cell(cls, cell: UnitCell) -> ReciprocalMetric:
        if not isinstance(cell, UnitCell):
            raise TypeError("cell must be UnitCell")
        direct_basis = cell.matrix
        reciprocal_basis = np.linalg.inv(direct_basis).T
        reciprocal_metric = np.linalg.inv(cell.metric)
        return cls(reciprocal_basis, reciprocal_metric)

    def norm_squared(self, hkl: MillerIndex) -> float:
        if not isinstance(hkl, MillerIndex):
            raise TypeError("hkl must be MillerIndex")
        if hkl.is_zero:
            raise ValueError("zero Miller index has no reciprocal norm or spacing")
        vector = np.asarray(hkl.as_tuple(), dtype=float)
        value = float(vector @ self.metric @ vector)
        if not math.isfinite(value) or value <= 0:
            raise ValueError("reciprocal metric produced a non-positive norm")
        return value

    def norm(self, hkl: MillerIndex) -> float:
        return math.sqrt(self.norm_squared(hkl))

    def d_spacing(self, hkl: MillerIndex) -> float:
        return 1.0 / self.norm(hkl)


@dataclass(frozen=True, slots=True)
class IntegerEllipsoidResult:
    indices: tuple[MillerIndex, ...]
    integer_points_tested: int
    reflections_within_d_min: int
    complete: bool


def _ellipsoid_leaves(
    metric: np.ndarray,
    limit: float,
) -> Iterator[tuple[int, int, int]]:
    upper = np.linalg.cholesky(metric).T
    values = [0, 0, 0]

    def visit(axis: int, spent: float) -> Iterator[tuple[int, int, int]]:
        if axis < 0:
            yield (values[0], values[1], values[2])
            return
        remaining = max(0.0, limit - spent)
        offset = math.fsum(upper[axis, j] * values[j] for j in range(axis + 1, 3))
        diagonal = float(upper[axis, axis])
        radius = math.sqrt(remaining) / diagonal
        center = -offset / diagonal
        lower = math.ceil(math.nextafter(center - radius, -math.inf))
        upper_bound = math.floor(math.nextafter(center + radius, math.inf))
        for candidate in range(lower, upper_bound + 1):
            row_value = diagonal * candidate + offset
            next_spent = spent + row_value * row_value
            if next_spent <= limit:
                values[axis] = candidate
                yield from visit(axis - 1, next_spent)

    yield from visit(2, 0.0)


def enumerate_integer_ellipsoid(
    reciprocal: ReciprocalMetric,
    d_min: float,
    max_candidates: int | None,
    boundary_tolerance: float,
) -> IntegerEllipsoidResult:
    """Enumerate integer Miller indices inside the physical ``d_min`` ellipsoid."""

    if not isinstance(reciprocal, ReciprocalMetric):
        raise TypeError("reciprocal must be ReciprocalMetric")
    d_min = _positive_finite(d_min, "d_min")
    boundary_tolerance = _positive_finite(boundary_tolerance, "boundary_tolerance")
    if max_candidates is not None and (
        isinstance(max_candidates, bool)
        or not isinstance(max_candidates, int)
        or max_candidates <= 0
    ):
        raise ValueError("max_candidates must be a positive integer or None")

    exact_limit = 1.0 / (d_min * d_min)
    tolerant_limit = exact_limit * (1.0 + boundary_tolerance)
    accepted: list[MillerIndex] = []
    tested = 0
    complete = True
    for values in _ellipsoid_leaves(reciprocal.metric, tolerant_limit):
        if max_candidates is not None and tested >= max_candidates:
            complete = False
            break
        tested += 1
        index = MillerIndex(*values)
        if index.is_zero:
            continue
        vector = np.asarray(values, dtype=float)
        quadratic = float(vector @ reciprocal.metric @ vector)
        if quadratic <= tolerant_limit:
            accepted.append(index)

    indices = tuple(sorted(set(accepted)))
    return IntegerEllipsoidResult(indices, tested, len(indices), complete)


__all__ = ["IntegerEllipsoidResult", "ReciprocalMetric", "enumerate_integer_ellipsoid"]
