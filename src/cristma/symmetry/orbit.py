"""Symmetry-derived sites with traceable asymmetric-unit provenance."""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Literal

from cristma.structure.crystal import IndependentSite
from cristma.structure.identity import ExpandedAtomRef, ExpandedSite

from .affine import AffineOperation


SymmetryProvenance = Literal[
    "reported",
    "derived",
    "identity_fallback",
    "unreported_identity",
]


@dataclass(frozen=True, slots=True)
class SpaceGroupDefinition:
    """Reported or derived space-group identity and exact operations."""

    operations: tuple[AffineOperation, ...]
    provenance: SymmetryProvenance
    number: int | None = None
    hm_symbol: str | None = None
    hall_symbol: str | None = None
    setting: str | None = None
    origin_choice: str | None = None

    def __post_init__(self) -> None:
        if not self.operations:
            raise ValueError("space group must contain at least one operation")
        if self.provenance not in {
            "reported",
            "derived",
            "identity_fallback",
            "unreported_identity",
        }:
            raise ValueError(f"unknown symmetry provenance: {self.provenance!r}")


def _raw_coordinates(
    operation: AffineOperation,
    coordinates: tuple[float, float, float],
) -> tuple[float, float, float]:
    return tuple(
        math.fsum(
            float(coefficient) * coordinate
            for coefficient, coordinate in zip(row, coordinates, strict=True)
        ) + float(offset)
        for row, offset in zip(
            operation.rotation,
            operation.translation,
            strict=True,
        )
    )


def _wrap_with_translation(
    raw: tuple[float, float, float],
    tolerance: float,
) -> tuple[tuple[float, float, float], tuple[int, int, int]]:
    wrapped = []
    translations = []
    for value in raw:
        nearest = round(value)
        normalized_value = float(nearest) if math.isclose(value, nearest, abs_tol=tolerance) else value
        translation = math.floor(normalized_value)
        coordinate = normalized_value - translation
        if math.isclose(coordinate, 1.0, abs_tol=tolerance):
            coordinate = 0.0
            translation += 1
        wrapped.append(0.0 if math.isclose(coordinate, 0.0, abs_tol=tolerance) else coordinate)
        translations.append(int(translation))
    return tuple(wrapped), tuple(translations)


def _periodically_equal(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
    tolerance: float,
) -> bool:
    return all(
        abs((a - b + 0.5) % 1.0 - 0.5) <= tolerance
        for a, b in zip(left, right, strict=True)
    )


def expand_orbit(
    site: IndependentSite,
    operations: tuple[AffineOperation, ...],
    tolerance: float = 1e-8,
    *,
    structure_id: str | None = None,
) -> tuple[ExpandedSite, ...]:
    """Expand one independent site and merge equivalent special positions."""

    if tolerance <= 0:
        raise ValueError("orbit tolerance must be positive")
    coordinates = tuple(float(value.value) for value in site.fractional)
    expanded: list[ExpandedSite] = []

    for index, operation in enumerate(operations, start=1):
        operation_id = operation.id or f"operation:{index}"
        fractional, translation = _wrap_with_translation(
            _raw_coordinates(operation, coordinates),
            tolerance,
        )
        for existing_index, existing in enumerate(expanded):
            if _periodically_equal(existing.fractional, fractional, tolerance):
                expanded[existing_index] = replace(
                    existing,
                    equivalent_operation_ids=(
                        *existing.equivalent_operation_ids,
                        operation_id,
                    ),
                )
                break
        else:
            expanded.append(
                ExpandedAtomRef(
                    id=(
                        f"expanded:{structure_id or 'unassigned'}:{site.id}:"
                        f"{operation_id}:{','.join(map(str, translation))}"
                    ),
                    structure_id=structure_id,
                    fractional=fractional,
                    source_site_id=site.id,
                    representative_operation_id=operation_id,
                    equivalent_operation_ids=(operation_id,),
                    cell_translation=translation,
                )
            )

    return tuple(expanded)
