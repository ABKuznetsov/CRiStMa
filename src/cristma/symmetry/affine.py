"""Exact rational affine operations used by crystallographic symmetry."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
import re


Vector3 = tuple[Fraction, Fraction, Fraction]
Matrix3 = tuple[Vector3, Vector3, Vector3]

_TERM = re.compile(r"[+-](?:[xyz]|\d+/\d+|\d+(?:\.\d*)?|\.\d+)")


@dataclass(frozen=True, slots=True)
class AffineOperation:
    """An exact fractional-coordinate operation ``x' = R x + t``."""

    rotation: Matrix3
    translation: Vector3
    source: str | None = None
    id: str | None = None

    def normalized(self) -> AffineOperation:
        """Return an equivalent operation with translation in ``[0, 1)``."""

        return AffineOperation(
            rotation=self.rotation,
            translation=tuple(value % 1 for value in self.translation),
            source=self.source,
            id=self.id,
        )

    def apply_fractional(
        self,
        coordinates: tuple[float, float, float],
    ) -> tuple[float, float, float]:
        """Apply the operation and wrap the result into the reference cell."""

        transformed = []
        for row, offset in zip(self.rotation, self.translation, strict=True):
            value = math.fsum(
                float(coefficient) * coordinate
                for coefficient, coordinate in zip(row, coordinates, strict=True)
            ) + float(offset)
            wrapped = value % 1.0
            transformed.append(0.0 if math.isclose(wrapped, 1.0, abs_tol=1e-15) else wrapped)
        return tuple(transformed)


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _format_component(row: Vector3, offset: Fraction) -> str:
    terms: list[str] = []
    for coefficient, variable in zip(row, "xyz", strict=True):
        if coefficient == 0:
            continue
        if coefficient not in {Fraction(-1), Fraction(1)}:
            raise ValueError("symmetry formatting requires unit coordinate coefficients")
        terms.append(
            f"-{variable}"
            if coefficient == -1
            else ("+" if terms else "") + variable
        )
    normalized = offset % 1
    if normalized:
        terms.append(("+" if terms else "") + _fraction_text(normalized))
    return "".join(terms) or "0"


def format_xyz_operation(operation: AffineOperation) -> str:
    """Render an exact affine operation as an ``x,y,z`` triplet."""

    if operation.source:
        return operation.source
    return ",".join(
        _format_component(row, offset)
        for row, offset in zip(
            operation.rotation,
            operation.translation,
            strict=True,
        )
    )


def _parse_component(expression: str) -> tuple[Vector3, Fraction]:
    compact = "".join(expression.split()).lower()
    if not compact:
        raise ValueError(f"Invalid symmetry expression: {expression!r}")
    signed = compact if compact[0] in "+-" else "+" + compact
    terms = _TERM.findall(signed)
    if not terms or "".join(terms) != signed:
        raise ValueError(f"Invalid symmetry expression: {expression!r}")

    coefficients = {"x": Fraction(0), "y": Fraction(0), "z": Fraction(0)}
    translation = Fraction(0)
    for term in terms:
        sign = Fraction(-1) if term[0] == "-" else Fraction(1)
        body = term[1:]
        if body in coefficients:
            coefficients[body] += sign
            if abs(coefficients[body]) > 1:
                raise ValueError(f"Invalid symmetry expression: {expression!r}")
        else:
            try:
                translation += sign * Fraction(body)
            except (ValueError, ZeroDivisionError) as exc:
                raise ValueError(f"Invalid symmetry expression: {expression!r}") from exc

    return (
        (coefficients["x"], coefficients["y"], coefficients["z"]),
        translation,
    )


def parse_xyz_operation(
    text: str,
    operation_id: str | None = None,
) -> AffineOperation:
    """Parse a linear ``x,y,z`` triplet without evaluating arbitrary code."""

    components = text.split(",")
    if len(components) != 3:
        raise ValueError(f"Invalid symmetry operation: {text!r}")
    parsed = [_parse_component(component) for component in components]
    return AffineOperation(
        rotation=tuple(item[0] for item in parsed),
        translation=tuple(item[1] for item in parsed),
        source=text,
        id=operation_id,
    )
