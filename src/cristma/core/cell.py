"""Validated crystallographic unit cells."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .values import MeasuredValue


def _required(value: MeasuredValue, name: str) -> float:
    if value.value is None or not math.isfinite(value.value):
        raise ValueError(f"{name} must be a finite reported value")
    return float(value.value)


@dataclass(frozen=True, slots=True)
class UnitCell:
    """Six reported cell parameters with derived Cartesian metric."""

    a: MeasuredValue
    b: MeasuredValue
    c: MeasuredValue
    alpha: MeasuredValue
    beta: MeasuredValue
    gamma: MeasuredValue

    def __post_init__(self) -> None:
        edges = tuple(_required(value, "cell edge") for value in (self.a, self.b, self.c))
        if any(edge <= 0 for edge in edges):
            raise ValueError("cell edge must be positive")

        angles = tuple(
            _required(value, "cell angle")
            for value in (self.alpha, self.beta, self.gamma)
        )
        if any(not 0 < angle < 180 for angle in angles):
            raise ValueError("cell angle must lie strictly between 0 and 180 degrees")

        if self.volume <= 0:
            raise ValueError("cell metric must have positive volume")

    @classmethod
    def cubic(cls, edge: MeasuredValue) -> UnitCell:
        right_angle = MeasuredValue(90.0, None, "90")
        return cls(edge, edge, edge, right_angle, right_angle, right_angle)

    @property
    def matrix(self) -> np.ndarray:
        """Return row-wise Cartesian basis vectors in angstrom."""

        a = _required(self.a, "cell edge")
        b = _required(self.b, "cell edge")
        c = _required(self.c, "cell edge")
        alpha, beta, gamma = (
            math.radians(_required(value, "cell angle"))
            for value in (self.alpha, self.beta, self.gamma)
        )
        cos_alpha = math.cos(alpha)
        cos_beta = math.cos(beta)
        cos_gamma = math.cos(gamma)
        sin_gamma = math.sin(gamma)
        if abs(sin_gamma) < 1e-15:
            raise ValueError("cell gamma produces a singular metric")

        c_x = c * cos_beta
        c_y = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
        c_z_squared = c * c - c_x * c_x - c_y * c_y
        if c_z_squared <= 0:
            raise ValueError("cell parameters produce a non-positive metric")

        matrix = np.array(
            [
                [a, 0.0, 0.0],
                [b * cos_gamma, b * sin_gamma, 0.0],
                [c_x, c_y, math.sqrt(c_z_squared)],
            ],
            dtype=float,
        )
        matrix.flags.writeable = False
        return matrix

    @property
    def metric(self) -> np.ndarray:
        basis = self.matrix
        return basis @ basis.T

    @property
    def volume(self) -> float:
        return float(abs(np.linalg.det(self.matrix)))
