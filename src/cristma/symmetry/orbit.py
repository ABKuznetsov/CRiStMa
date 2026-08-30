"""Symmetry-derived sites with traceable asymmetric-unit provenance."""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Literal

import numpy as np

from cristma.core.cell import UnitCell
from cristma.structure.crystal import CrystalStructure, IndependentSite
from cristma.structure.identity import ExpandedAtom, SymmetryImageProvenance
from cristma.structure.properties import AtomicPropertyTable
from cristma.structure.view import AtomicView

from .affine import AffineOperation
from .displacement import SymmetryConsistencyError, displacements_close, transform_displacement


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
        normalization_translation = -math.floor(normalized_value)
        coordinate = normalized_value + normalization_translation
        if math.isclose(coordinate, 1.0, abs_tol=tolerance):
            coordinate = 0.0
            normalization_translation -= 1
        wrapped.append(0.0 if math.isclose(coordinate, 0.0, abs_tol=tolerance) else coordinate)
        translations.append(int(normalization_translation))
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
    cell: UnitCell,
    structure_id: str | None = None,
) -> tuple[ExpandedAtom, ...]:
    """Expand one independent site and merge equivalent special positions."""

    if tolerance <= 0:
        raise ValueError("orbit tolerance must be positive")
    coordinates = tuple(float(value.value) for value in site.fractional)
    expanded: list[ExpandedAtom] = []

    for index, operation in enumerate(operations, start=1):
        operation_id = operation.id or f"operation:{index}"
        fractional, translation = _wrap_with_translation(
            _raw_coordinates(operation, coordinates),
            tolerance,
        )
        image = SymmetryImageProvenance(operation_id, translation)
        displacement = transform_displacement(site.displacement, operation)
        for existing_index, existing in enumerate(expanded):
            if _periodically_equal(existing.fractional, fractional, tolerance):
                if not displacements_close(existing.displacement, displacement, tolerance):
                    raise SymmetryConsistencyError(
                        f"site {site.label!r} has inconsistent anisotropic displacement "
                        "at symmetry-equivalent images"
                    )
                expanded[existing_index] = replace(
                    existing,
                    equivalent_images=(*existing.equivalent_images, image),
                )
                break
        else:
            cartesian = tuple(float(value) for value in np.asarray(fractional) @ cell.matrix)
            position_key = tuple(round(value / tolerance) for value in fractional)
            expanded.append(
                ExpandedAtom(
                    id=(
                        f"expanded:{structure_id or 'unassigned'}:{site.id}:"
                        f"{','.join(map(str, position_key))}"
                    ),
                    structure_id=structure_id,
                    source_site_id=site.id,
                    fractional=fractional,
                    cartesian=cartesian,
                    components=site.components,
                    displacement=displacement,
                    representative_image=image,
                    equivalent_images=(image,),
                )
            )

    return tuple(expanded)


def expand_structure(
    crystal: CrystalStructure,
    tolerance: float = 1e-8,
) -> AtomicView[ExpandedAtom]:
    """Expand all independent sites into one finite reference-cell atomic view."""

    if crystal.space_group is None:
        raise ValueError("crystal symmetry is required for structure expansion")
    atoms = tuple(
        atom
        for site in crystal.sites
        for atom in expand_orbit(
            site,
            crystal.space_group.operations,
            tolerance,
            cell=crystal.cell,
            structure_id=crystal.id,
        )
    )
    ids = tuple(atom.id for atom in atoms)
    if len(set(ids)) != len(ids):
        raise ValueError("symmetry expansion produced duplicate expanded atom IDs")
    return AtomicView(
        atoms=atoms,
        cell=crystal.cell,
        periodic=crystal.periodic,
        properties=AtomicPropertyTable(len(atoms)),
    )
