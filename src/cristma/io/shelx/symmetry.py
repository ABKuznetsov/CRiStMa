"""Exact SHELX lattice and explicit-operation construction."""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
import re

from cristma.symmetry.affine import AffineOperation, parse_xyz_operation


_ZERO = Fraction(0)
_HALF = Fraction(1, 2)
_THIRD = Fraction(1, 3)
_TWO_THIRDS = Fraction(2, 3)
_IDENTITY = (
    (Fraction(1), _ZERO, _ZERO),
    (_ZERO, Fraction(1), _ZERO),
    (_ZERO, _ZERO, Fraction(1)),
)
_CENTERING = {
    1: ((_ZERO, _ZERO, _ZERO),),
    2: ((_ZERO, _ZERO, _ZERO), (_HALF, _HALF, _HALF)),
    3: (
        (_ZERO, _ZERO, _ZERO),
        (_TWO_THIRDS, _THIRD, _THIRD),
        (_THIRD, _TWO_THIRDS, _TWO_THIRDS),
    ),
    4: (
        (_ZERO, _ZERO, _ZERO),
        (_ZERO, _HALF, _HALF),
        (_HALF, _ZERO, _HALF),
        (_HALF, _HALF, _ZERO),
    ),
    5: ((_ZERO, _ZERO, _ZERO), (_ZERO, _HALF, _HALF)),
    6: ((_ZERO, _ZERO, _ZERO), (_HALF, _ZERO, _HALF)),
    7: ((_ZERO, _ZERO, _ZERO), (_HALF, _HALF, _ZERO)),
}
_DECIMAL = re.compile(r"(?<![A-Za-z0-9_/])(\d*\.\d+)(?![A-Za-z0-9_/])")


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _component_text(row: tuple[Fraction, Fraction, Fraction], offset: Fraction) -> str:
    pieces: list[str] = []
    for coefficient, variable in zip(row, "xyz", strict=True):
        if coefficient == 0:
            continue
        if coefficient not in {Fraction(-1), Fraction(1)}:
            raise ValueError("SHELX symmetry requires unit coordinate coefficients")
        if coefficient == -1:
            pieces.append(f"-{variable}")
        else:
            pieces.append(("+" if pieces else "") + variable)
    normalized = offset % 1
    if normalized:
        pieces.append(("+" if pieces else "") + _fraction_text(normalized))
    return "".join(pieces) or "0"


def _source(operation: AffineOperation) -> str:
    return ",".join(
        _component_text(row, offset)
        for row, offset in zip(operation.rotation, operation.translation, strict=True)
    )


def parse_shelx_symmetry(text: str) -> AffineOperation:
    """Parse SHELX symmetry, accepting exact decimal translations."""

    rational = _DECIMAL.sub(lambda match: _fraction_text(Fraction(match.group(1))), text)
    return parse_xyz_operation(rational)


def _inverted(operation: AffineOperation) -> AffineOperation:
    return AffineOperation(
        rotation=tuple(tuple(-value for value in row) for row in operation.rotation),
        translation=tuple(-value for value in operation.translation),
    ).normalized()


def _centered(
    operation: AffineOperation,
    translation: tuple[Fraction, Fraction, Fraction],
) -> AffineOperation:
    return AffineOperation(
        rotation=operation.rotation,
        translation=tuple(
            (left + right) % 1
            for left, right in zip(operation.translation, translation, strict=True)
        ),
    )


def build_shelx_operations(
    latt: int,
    reported: tuple[AffineOperation, ...],
) -> tuple[AffineOperation, ...]:
    """Build the full exact operation set implied by LATT and SYMM."""

    lattice = abs(latt)
    if lattice not in _CENTERING:
        raise ValueError(f"unknown SHELX LATT code: {latt}")
    identity = AffineOperation(_IDENTITY, (_ZERO, _ZERO, _ZERO))
    primitive = (identity, *reported)
    if latt > 0:
        primitive = (*primitive, *(_inverted(operation) for operation in primitive))

    unique: list[AffineOperation] = []
    seen = set()
    for operation in primitive:
        for centering in _CENTERING[lattice]:
            candidate = _centered(operation.normalized(), centering)
            key = (candidate.rotation, candidate.translation)
            if key in seen:
                continue
            seen.add(key)
            unique.append(candidate)
    return tuple(
        replace(operation, source=_source(operation), id=f"shelx:op:{index}")
        for index, operation in enumerate(unique, start=1)
    )


__all__ = ["build_shelx_operations", "parse_shelx_symmetry"]
