"""Exact VASP scale and coordinate conversion rules."""

from __future__ import annotations

import math

import numpy as np

from .document import VaspScale


def _matrix(value: object, *, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if result.ndim != 2 or result.shape[1:] != (3,) or not np.isfinite(result).all():
        raise ValueError(f"{name} must be a finite Nx3 array")
    return result


def _component_factors(scale: VaspScale, raw_lattice: np.ndarray) -> np.ndarray:
    if len(scale.values) == 3:
        return np.asarray(scale.values, dtype=float)
    reported = scale.values[0]
    if reported > 0:
        return np.full(3, reported, dtype=float)
    determinant = abs(float(np.linalg.det(raw_lattice)))
    if determinant <= 1e-15:
        raise ValueError("negative-volume scaling requires a non-singular raw lattice")
    factor = math.cbrt(abs(reported) / determinant)
    return np.full(3, factor, dtype=float)


def scaled_lattice(scale: VaspScale, raw_lattice: object) -> np.ndarray:
    """Return row-wise lattice vectors in angstrom."""

    raw = _matrix(raw_lattice, name="raw lattice")
    if raw.shape != (3, 3):
        raise ValueError("raw lattice must be a 3x3 array")
    result = np.array(raw * _component_factors(scale, raw), copy=True)
    if abs(float(np.linalg.det(result))) <= 1e-15:
        raise ValueError("scaled lattice is singular")
    result.flags.writeable = False
    return result


def scaled_cartesian(
    scale: VaspScale,
    rows: object,
    raw_lattice: object,
) -> np.ndarray:
    """Apply VASP position scaling to reported Cartesian rows."""

    raw_cell = _matrix(raw_lattice, name="raw lattice")
    if raw_cell.shape != (3, 3):
        raise ValueError("raw lattice must be a 3x3 array")
    result = np.array(
        _matrix(rows, name="Cartesian coordinates")
        * _component_factors(scale, raw_cell),
        copy=True,
    )
    result.flags.writeable = False
    return result


def fractional_from_cartesian(lattice: object, rows: object) -> np.ndarray:
    """Convert row-wise Cartesian positions to fractional coordinates."""

    cell = _matrix(lattice, name="lattice")
    if cell.shape != (3, 3) or abs(float(np.linalg.det(cell))) <= 1e-15:
        raise ValueError("lattice must be a non-singular 3x3 array")
    result = np.array(_matrix(rows, name="Cartesian coordinates") @ np.linalg.inv(cell), copy=True)
    result.flags.writeable = False
    return result


__all__ = ["fractional_from_cartesian", "scaled_cartesian", "scaled_lattice"]
