"""Exact affine coordinate constraints for Wyckoff positions."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import re

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .space_group import SpaceGroupKey


FractionVector = tuple[Fraction, Fraction, Fraction]
FractionMatrix = tuple[FractionVector, FractionVector, FractionVector]
_TERM = re.compile(r"[+-](?:(?:\d+(?:/\d+)?)?[xyz]|\d+(?:/\d+)?)")


def _parse_component(expression: str) -> tuple[FractionVector, Fraction]:
    compact = "".join(expression.split()).lower()
    if not compact:
        raise ValueError(f"invalid affine coordinate expression: {expression!r}")
    signed = compact if compact.startswith(("+", "-")) else "+" + compact
    terms = _TERM.findall(signed)
    if not terms or "".join(terms) != signed:
        raise ValueError(f"invalid affine coordinate expression: {expression!r}")

    coefficients = {axis: Fraction(0) for axis in "xyz"}
    translation = Fraction(0)
    for term in terms:
        sign = Fraction(-1) if term[0] == "-" else Fraction(1)
        body = term[1:]
        if body[-1:] in coefficients:
            variable = body[-1]
            raw_coefficient = body[:-1]
            coefficient = Fraction(1) if not raw_coefficient else Fraction(raw_coefficient)
            coefficients[variable] += sign * coefficient
        else:
            translation += sign * Fraction(body)
    return tuple(coefficients[axis] for axis in "xyz"), translation


def _matrix_rank_over_rationals(matrix: FractionMatrix) -> int:
    rows = [list(row) for row in matrix]
    rank = 0
    column_count = len(rows[0]) if rows else 0
    for column in range(column_count):
        pivot = next(
            (index for index in range(rank, len(rows)) if rows[index][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][column]
        rows[rank] = [value / pivot_value for value in rows[rank]]
        for index, row in enumerate(rows):
            if index == rank or not row[column]:
                continue
            factor = row[column]
            rows[index] = [
                value - factor * pivot_item
                for value, pivot_item in zip(row, rows[rank], strict=True)
            ]
        rank += 1
        if rank == len(rows):
            break
    return rank


@dataclass(frozen=True, slots=True)
class AffineCoordinateMap:
    """Map free fractional parameters to one Wyckoff representative."""

    parameter_matrix: FractionMatrix
    translation: FractionVector
    source: str | None = None

    def __post_init__(self) -> None:
        if len(self.parameter_matrix) != 3 or any(
            len(row) != 3 for row in self.parameter_matrix
        ):
            raise ValueError("parameter matrix must be 3 by 3")
        if len(self.translation) != 3:
            raise ValueError("translation must contain three fractions")
        if any(not isinstance(value, Fraction) for row in self.parameter_matrix for value in row):
            raise TypeError("parameter matrix entries must be Fraction values")
        if any(not isinstance(value, Fraction) for value in self.translation):
            raise TypeError("translation entries must be Fraction values")

    @classmethod
    def from_xyz(cls, text: str) -> "AffineCoordinateMap":
        source = text.strip()
        expression = source[1:-1] if source.startswith("(") and source.endswith(")") else source
        components = expression.split(",")
        if len(components) != 3:
            raise ValueError(f"affine coordinate map requires three components: {text!r}")
        parsed = tuple(_parse_component(component) for component in components)
        return cls(
            parameter_matrix=tuple(item[0] for item in parsed),
            translation=tuple(item[1] for item in parsed),
            source=source,
        )

    @property
    def degrees_of_freedom(self) -> int:
        return _matrix_rank_over_rationals(self.parameter_matrix)

    def evaluate(self, parameters: FractionVector) -> FractionVector:
        if len(parameters) != 3 or any(not isinstance(value, Fraction) for value in parameters):
            raise TypeError("parameters must contain three Fraction values")
        return tuple(
            sum(
                (coefficient * value for coefficient, value in zip(row, parameters, strict=True)),
                start=translation,
            )
            for row, translation in zip(
                self.parameter_matrix,
                self.translation,
                strict=True,
            )
        )


@dataclass(frozen=True, slots=True)
class WyckoffPosition:
    """One setting-specific Wyckoff position from the reference catalog."""

    space_group_key: SpaceGroupKey
    letter: str
    multiplicity: int
    site_symmetry_symbol: str
    representatives: tuple[AffineCoordinateMap, ...]

    def __post_init__(self) -> None:
        if len(self.letter) != 1 or self.letter not in "abcdefghijklmnopqrstuvwxyz":
            raise ValueError("Wyckoff letter must be one lower-case ASCII letter")
        if isinstance(self.multiplicity, bool) or self.multiplicity <= 0:
            raise ValueError("Wyckoff multiplicity must be positive")
        if len(self.representatives) != self.multiplicity:
            raise ValueError("representative count must equal Wyckoff multiplicity")
        if not self.site_symmetry_symbol:
            raise ValueError("site-symmetry symbol must not be empty")

    @property
    def degrees_of_freedom(self) -> int:
        return self.representatives[0].degrees_of_freedom


__all__ = ["AffineCoordinateMap", "WyckoffPosition"]
