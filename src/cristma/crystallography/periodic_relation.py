"""Exact affine-periodic relations over a validated symmetry context."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .symmetry_context import (
    SymmetryContext,
    _compose,
    _IDENTITY_ROTATION,
    _matrix_vector,
    _ZERO_TRANSLATION,
    canonical_operation_key,
)


LatticeTranslation = tuple[int, int, int]
ExactFractionalPosition = tuple[Fraction, Fraction, Fraction]


def _integer_vector(values: tuple[object, ...]) -> LatticeTranslation:
    if len(values) != 3:
        raise ValueError("lattice translation must contain three integers")
    converted: list[int] = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError("lattice translation must contain only integers")
        converted.append(value)
    return tuple(converted)


def _exact_integer(value: Fraction, *, name: str) -> int:
    if value.denominator != 1:
        raise ArithmeticError(f"{name} is not an exact lattice translation")
    return value.numerator


@dataclass(frozen=True, slots=True, order=True)
class PeriodicSymmetryRelation:
    """An action ``x -> R x + t + n`` in the full periodic affine group."""

    operation_key: str
    lattice_translation: LatticeTranslation

    def __post_init__(self) -> None:
        if not isinstance(self.operation_key, str) or not self.operation_key:
            raise ValueError("operation_key must be a nonempty string")
        object.__setattr__(
            self,
            "lattice_translation",
            _integer_vector(tuple(self.lattice_translation)),
        )

    def normalize(self, context: SymmetryContext) -> "PeriodicSymmetryRelation":
        """Validate this already-canonical relation against ``context``."""

        context.operation_by_key(self.operation_key)
        return self

    def compose(
        self,
        other: "PeriodicSymmetryRelation",
        context: SymmetryContext,
    ) -> "PeriodicSymmetryRelation":
        """Return ``self`` after ``other`` using exact affine composition."""

        return compose_periodic_relations(self, other, context)

    def inverse(self, context: SymmetryContext) -> "PeriodicSymmetryRelation":
        """Return the exact two-sided inverse within ``context``."""

        return invert_periodic_relation(self, context)

    def apply_fractional(
        self,
        coordinates: tuple[Fraction | int, Fraction | int, Fraction | int],
        context: SymmetryContext,
    ) -> ExactFractionalPosition:
        """Apply the full relation without wrapping into the reference cell."""

        operation = context.operation_by_key(self.operation_key)
        exact_coordinates = tuple(Fraction(value) for value in coordinates)
        rotated = _matrix_vector(operation.rotation, exact_coordinates)
        return tuple(
            rotated[index]
            + operation.translation[index]
            + self.lattice_translation[index]
            for index in range(3)
        )


def identity_relation(context: SymmetryContext) -> PeriodicSymmetryRelation:
    """Return the exact identity element for ``context``."""

    identity_key = canonical_operation_key(
        context.operations[
            next(
                index
                for index, operation in enumerate(context.operations)
                if operation.rotation == _IDENTITY_ROTATION
                and operation.translation == _ZERO_TRANSLATION
            )
        ]
    )
    return PeriodicSymmetryRelation(identity_key, (0, 0, 0))


def compose_periodic_relations(
    left: PeriodicSymmetryRelation,
    right: PeriodicSymmetryRelation,
    context: SymmetryContext,
) -> PeriodicSymmetryRelation:
    """Compose ``left`` after ``right`` and retain the exact lattice carry."""

    left_operation = context.operation_by_key(left.operation_key)
    right_operation = context.operation_by_key(right.operation_key)
    product_operation = _compose(left_operation, right_operation)
    product_key = canonical_operation_key(product_operation)
    context.operation_by_key(product_key)

    raw_translation = tuple(
        value + left_operation.translation[index]
        for index, value in enumerate(
            _matrix_vector(left_operation.rotation, right_operation.translation)
        )
    )
    carry = tuple(
        _exact_integer(
            raw_translation[index] - product_operation.translation[index],
            name="composition carry",
        )
        for index in range(3)
    )
    rotated_right_lattice = _matrix_vector(
        left_operation.rotation,
        tuple(Fraction(value) for value in right.lattice_translation),
    )
    lattice_translation = tuple(
        _exact_integer(rotated_right_lattice[index], name="rotated lattice translation")
        + left.lattice_translation[index]
        + carry[index]
        for index in range(3)
    )
    return PeriodicSymmetryRelation(product_key, lattice_translation)


def invert_periodic_relation(
    relation: PeriodicSymmetryRelation,
    context: SymmetryContext,
) -> PeriodicSymmetryRelation:
    """Invert a relation while preserving the normalization lattice carry."""

    operation = context.operation_by_key(relation.operation_key)
    identity = identity_relation(context)
    inverse_operation = next(
        candidate
        for candidate in context.operations
        if canonical_operation_key(_compose(operation, candidate)) == identity.operation_key
        and canonical_operation_key(_compose(candidate, operation)) == identity.operation_key
    )
    total_translation = tuple(
        operation.translation[index] + relation.lattice_translation[index]
        for index in range(3)
    )
    raw_inverse_translation = tuple(
        -value
        for value in _matrix_vector(inverse_operation.rotation, total_translation)
    )
    inverse_lattice = tuple(
        _exact_integer(
            raw_inverse_translation[index] - inverse_operation.translation[index],
            name="inverse carry",
        )
        for index in range(3)
    )
    return PeriodicSymmetryRelation(
        canonical_operation_key(inverse_operation),
        inverse_lattice,
    )


__all__ = [
    "ExactFractionalPosition",
    "LatticeTranslation",
    "PeriodicSymmetryRelation",
    "compose_periodic_relations",
    "identity_relation",
    "invert_periodic_relation",
]
