from __future__ import annotations

from fractions import Fraction

import pytest

from cristma.core import MeasuredValue, UnitCell
from cristma.crystallography import (
    DirectBasisConvention,
    SymmetryContext,
    SymmetryContextInvariantError,
    SymmetrySourceKind,
)
from cristma.symmetry import AffineOperation


def _value(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def _orthogonal_cell(a: float = 5.0, b: float = 5.0, c: float = 7.0) -> UnitCell:
    return UnitCell(
        _value(a),
        _value(b),
        _value(c),
        _value(90.0),
        _value(90.0),
        _value(90.0),
    )


def _operation(
    rotation: tuple[tuple[int, int, int], ...],
    translation: tuple[Fraction | int, Fraction | int, Fraction | int] = (0, 0, 0),
    *,
    operation_id: str | None = None,
) -> AffineOperation:
    return AffineOperation(
        rotation=tuple(
            tuple(Fraction(value) for value in row)
            for row in rotation
        ),
        translation=tuple(Fraction(value) for value in translation),
        id=operation_id,
    )


IDENTITY = _operation(((1, 0, 0), (0, 1, 0), (0, 0, 1)), operation_id="reported-7")
INVERSION = _operation(((-1, 0, 0), (0, -1, 0), (0, 0, -1)), operation_id="reported-2")


def test_explicit_operation_order_does_not_change_context_identity() -> None:
    cell = _orthogonal_cell()

    first = SymmetryContext.from_operations((IDENTITY, INVERSION), cell)
    second = SymmetryContext.from_operations((INVERSION, IDENTITY), cell)

    assert first.operation_keys == second.operation_keys
    assert first.symmetry_action_fingerprint == second.symmetry_action_fingerprint
    assert first.operations == second.operations


def test_context_identity_ignores_reported_operation_ids() -> None:
    cell = _orthogonal_cell()
    renamed_identity = _operation(
        ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        operation_id="anything",
    )
    renamed_inversion = _operation(
        ((-1, 0, 0), (0, -1, 0), (0, 0, -1)),
        operation_id="else",
    )

    first = SymmetryContext.from_operations((IDENTITY, INVERSION), cell)
    second = SymmetryContext.from_operations((renamed_inversion, renamed_identity), cell)

    assert first.symmetry_action_fingerprint == second.symmetry_action_fingerprint


def test_valid_explicit_operations_are_complete_without_setting_id() -> None:
    context = SymmetryContext.from_operations((IDENTITY, INVERSION), _orthogonal_cell())

    assert context.source_kind is SymmetrySourceKind.VALID_EXPLICIT_OPERATIONS
    assert context.setting_id is None
    assert context.basis_convention is DirectBasisConvention.FRACTIONAL_DIRECT


def test_non_closed_operation_set_is_rejected() -> None:
    quarter_turn = _operation(((0, -1, 0), (1, 0, 0), (0, 0, 1)))

    with pytest.raises(SymmetryContextInvariantError) as caught:
        SymmetryContext.from_operations((IDENTITY, quarter_turn), _orthogonal_cell())

    assert caught.value.code == "symmetry.context.group_invalid"


def test_rotation_incompatible_with_cell_metric_is_rejected() -> None:
    swap_xy = _operation(((0, 1, 0), (1, 0, 0), (0, 0, 1)))

    with pytest.raises(SymmetryContextInvariantError) as caught:
        SymmetryContext.from_operations((IDENTITY, swap_xy), _orthogonal_cell(4.0, 5.0))

    assert caught.value.code == "symmetry.context.metric_incompatible"


def test_duplicate_normalized_operations_are_rejected() -> None:
    translated_identity = _operation(
        ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        (1, 0, 0),
    )

    with pytest.raises(SymmetryContextInvariantError) as caught:
        SymmetryContext.from_operations((IDENTITY, translated_identity), _orthogonal_cell())

    assert caught.value.code == "symmetry.context.duplicate_operation"
