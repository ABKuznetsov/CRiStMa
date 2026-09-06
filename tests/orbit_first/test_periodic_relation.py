from __future__ import annotations

from fractions import Fraction
from itertools import product

import pytest

from cristma.core import MeasuredValue, UnitCell
from cristma.crystallography import (
    PeriodicSymmetryRelation,
    SpaceGroupCatalog,
    SymmetryContext,
    canonical_operation_key,
    identity_relation,
)
from cristma.symmetry import AffineOperation


def _value(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def _cell() -> UnitCell:
    return UnitCell.cubic(_value(5.0))


def _operation(
    rotation: tuple[tuple[int, int, int], ...],
    translation: tuple[Fraction | int, Fraction | int, Fraction | int] = (0, 0, 0),
) -> AffineOperation:
    return AffineOperation(
        tuple(tuple(Fraction(value) for value in row) for row in rotation),
        tuple(Fraction(value) for value in translation),
    )


IDENTITY_OPERATION = _operation(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
TWO_FOLD_WITH_TRANSLATION = _operation(
    ((-1, 0, 0), (0, -1, 0), (0, 0, 1)),
    (Fraction(1, 2), 0, Fraction(1, 2)),
)


@pytest.fixture
def context() -> SymmetryContext:
    return SymmetryContext.from_operations(
        (IDENTITY_OPERATION, TWO_FOLD_WITH_TRANSLATION),
        _cell(),
    )


def test_relation_composition_rotates_lattice_translation_and_keeps_carry(
    context: SymmetryContext,
) -> None:
    operation_key = canonical_operation_key(TWO_FOLD_WITH_TRANSLATION)
    left = PeriodicSymmetryRelation(operation_key, (1, 0, 0))
    right = PeriodicSymmetryRelation(operation_key, (0, 1, 0))

    composed = left.compose(right, context)

    assert composed == PeriodicSymmetryRelation(
        canonical_operation_key(IDENTITY_OPERATION),
        (1, -1, 1),
    )


def test_relation_inverse_is_two_sided_exact_identity(context: SymmetryContext) -> None:
    relation = PeriodicSymmetryRelation(
        canonical_operation_key(TWO_FOLD_WITH_TRANSLATION),
        (1, -2, 3),
    )

    inverse = relation.inverse(context)

    assert inverse == PeriodicSymmetryRelation(relation.operation_key, (1, -2, -4))
    assert relation.compose(inverse, context) == identity_relation(context)
    assert inverse.compose(relation, context) == identity_relation(context)


def test_apply_fractional_keeps_the_full_unwrapped_periodic_action(
    context: SymmetryContext,
) -> None:
    relation = PeriodicSymmetryRelation(
        canonical_operation_key(TWO_FOLD_WITH_TRANSLATION),
        (1, -2, 3),
    )

    assert relation.apply_fractional(
        (Fraction(1, 4), Fraction(1, 3), Fraction(2, 5)),
        context,
    ) == (Fraction(5, 4), Fraction(-7, 3), Fraction(39, 10))


def test_relation_rejects_unknown_operation_key(context: SymmetryContext) -> None:
    relation = PeriodicSymmetryRelation("operation:not-in-this-context", (0, 0, 0))

    with pytest.raises(KeyError):
        relation.normalize(context)


def test_group_laws_hold_for_all_operation_pairs_and_small_lattice_vectors(
    context: SymmetryContext,
) -> None:
    translations = ((0, 0, 0), (1, -1, 2), (-2, 1, 0))
    relations = tuple(
        PeriodicSymmetryRelation(operation_key, translation)
        for operation_key, translation in product(context.operation_keys, translations)
    )
    identity = identity_relation(context)

    for relation in relations:
        assert relation.compose(identity, context) == relation
        assert identity.compose(relation, context) == relation
        assert relation.compose(relation.inverse(context), context) == identity
        assert relation.inverse(context).compose(relation, context) == identity

    for left, middle, right in product(relations, repeat=3):
        assert left.compose(middle, context).compose(right, context) == left.compose(
            middle.compose(right, context),
            context,
        )


@pytest.mark.parametrize("setting_id", (2, 525))
def test_every_operation_pair_closes_as_a_periodic_relation(setting_id: int) -> None:
    setting = SpaceGroupCatalog.default().by_setting(setting_id)
    context = SymmetryContext.from_setting(setting, _cell())
    relations = tuple(
        PeriodicSymmetryRelation(operation_key, (0, 0, 0))
        for operation_key in context.operation_keys
    )

    for left, right in product(relations, repeat=2):
        composed = left.compose(right, context)
        assert composed.operation_key in context.operation_keys
        assert all(type(value) is int for value in composed.lattice_translation)

    for relation in relations:
        assert relation.compose(relation.inverse(context), context) == identity_relation(context)
