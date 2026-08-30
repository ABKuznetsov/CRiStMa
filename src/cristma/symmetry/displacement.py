"""Symmetry transformations for atomic displacement parameters."""

from __future__ import annotations

import math

import numpy as np

from cristma.core.values import MeasuredValue
from cristma.structure.crystal import DisplacementParameters

from .affine import AffineOperation


class SymmetryConsistencyError(ValueError):
    """Equivalent symmetry images imply incompatible scientific state."""


def _derived_value(value: float, uncertainty: float | None, source: MeasuredValue) -> MeasuredValue:
    return MeasuredValue(
        value=float(value),
        uncertainty=None if uncertainty is None else float(uncertainty),
        raw=None,
        unit=source.unit,
        missing=source.missing,
    )


def _tensor_values(
    displacement: DisplacementParameters,
) -> tuple[np.ndarray, tuple[float | None, ...], tuple[MeasuredValue, ...]]:
    if displacement.tensor is None or len(displacement.tensor) != 3:
        raise ValueError("anisotropic displacement requires a 3 by 3 tensor")
    rows = displacement.tensor
    if any(len(row) != 3 for row in rows):
        raise ValueError("anisotropic displacement requires a 3 by 3 tensor")
    values = np.array([[item.value for item in row] for row in rows], dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("anisotropic displacement tensor must be finite")
    packed = (rows[0][0], rows[1][1], rows[2][2], rows[0][1], rows[0][2], rows[1][2])
    uncertainties = tuple(
        None if item.uncertainty is None else float(item.uncertainty)
        for item in packed
    )
    return values, uncertainties, packed


def _tensor_linear_map(rotation: np.ndarray) -> np.ndarray:
    basis = []
    for row, column in ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)):
        matrix = np.zeros((3, 3), dtype=float)
        matrix[row, column] = 1.0
        matrix[column, row] = 1.0
        basis.append(matrix)
    columns = []
    for matrix in basis:
        transformed = rotation @ matrix @ rotation.T
        columns.append(
            (transformed[0, 0], transformed[1, 1], transformed[2, 2],
             transformed[0, 1], transformed[0, 2], transformed[1, 2])
        )
    return np.asarray(columns, dtype=float).T


def transform_displacement(
    displacement: DisplacementParameters | None,
    operation: AffineOperation,
) -> DisplacementParameters | None:
    """Transform one displacement model into a symmetry-image orientation."""

    if displacement is None:
        return None
    if displacement.kind in {"U_iso", "B_iso"}:
        return displacement
    if displacement.kind != "U_aniso":
        raise ValueError(f"unsupported displacement kind: {displacement.kind!r}")

    values, uncertainties, packed = _tensor_values(displacement)
    rotation = np.asarray(operation.rotation, dtype=float)
    transformed = rotation @ values @ rotation.T
    transformed = 0.5 * (transformed + transformed.T)

    linear = _tensor_linear_map(rotation)
    propagated: list[float | None] = []
    for coefficients in linear:
        contributors = [
            (coefficient, uncertainty)
            for coefficient, uncertainty in zip(coefficients, uncertainties, strict=True)
            if not math.isclose(coefficient, 0.0, abs_tol=1e-15)
        ]
        if any(uncertainty is None for _, uncertainty in contributors):
            propagated.append(None)
        else:
            propagated.append(
                math.sqrt(
                    math.fsum(
                        (coefficient * float(uncertainty)) ** 2
                        for coefficient, uncertainty in contributors
                    )
                )
            )

    uncertainty_matrix: list[list[float | None]] = [[None] * 3 for _ in range(3)]
    for value, (row, column) in zip(
        propagated,
        ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)),
        strict=True,
    ):
        uncertainty_matrix[row][column] = value
        uncertainty_matrix[column][row] = value

    source = packed[0]
    tensor = tuple(
        tuple(_derived_value(transformed[row, column], uncertainty_matrix[row][column], source)
              for column in range(3))
        for row in range(3)
    )
    return DisplacementParameters(
        kind="U_aniso",
        tensor=tensor,
        reported_kind=displacement.reported_kind,
    )


def displacements_close(
    left: DisplacementParameters | None,
    right: DisplacementParameters | None,
    tolerance: float,
) -> bool:
    """Compare equivalent-image displacement values without source formatting."""

    if left is None or right is None:
        return left is right
    if left.kind != right.kind:
        return False
    if left.kind in {"U_iso", "B_iso"}:
        if left.isotropic is None or right.isotropic is None:
            return left.isotropic is right.isotropic
        if left.isotropic.value is None or right.isotropic.value is None:
            return left.isotropic.value is right.isotropic.value
        return math.isclose(left.isotropic.value, right.isotropic.value, rel_tol=0.0, abs_tol=tolerance)
    if left.kind == "U_aniso":
        left_values, _, _ = _tensor_values(left)
        right_values, _, _ = _tensor_values(right)
        return bool(np.allclose(left_values, right_values, rtol=0.0, atol=tolerance))
    return False


__all__ = ["SymmetryConsistencyError", "displacements_close", "transform_displacement"]
