"""Exact affine-periodic relations over a validated symmetry context."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .symmetry_context import (
    SymmetryContext,
    _matrix_vector,
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

    return PeriodicSymmetryRelation(context.identity_operation_key, (0, 0, 0))


def compose_periodic_relations(
    left: PeriodicSymmetryRelation,
    right: PeriodicSymmetryRelation,
    context: SymmetryContext,
) -> PeriodicSymmetryRelation:
    """Compose ``left`` after ``right`` and retain the exact lattice carry."""

    left_operation = context.operation_by_key(left.operation_key)
    product_key, carry = context.compose_operation_keys(
        left.operation_key,
        right.operation_key,
    )
    rotated_right_lattice = tuple(
        sum(
            int(left_operation.rotation[row][column])
            * right.lattice_translation[column]
            for column in range(3)
        )
        for row in range(3)
    )
    lattice_translation = tuple(
        rotated_right_lattice[index]
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

    inverse_key = context.inverse_operation_key(relation.operation_key)
    inverse_operation = context.operation_by_key(inverse_key)
    _, carry = context.compose_operation_keys(
        inverse_key,
        relation.operation_key,
    )
    rotated_lattice = tuple(
        sum(
            int(inverse_operation.rotation[row][column])
            * relation.lattice_translation[column]
            for column in range(3)
        )
        for row in range(3)
    )
    inverse_lattice = tuple(
        -rotated_lattice[index] - carry[index]
        for index in range(3)
    )
    return PeriodicSymmetryRelation(inverse_key, inverse_lattice)


__all__ = [
    "ExactFractionalPosition",
    "LatticeTranslation",
    "PeriodicSymmetryRelation",
    "compose_periodic_relations",
    "identity_relation",
    "invert_periodic_relation",
]
