"""Calculated crystallographic orbits and site stabilizers."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from string import ascii_letters
from typing import Literal

import numpy as np

from cristma.core.cell import UnitCell
from cristma.diagnostics import Diagnostic, Severity
from cristma.structure.crystal import IndependentSite
from cristma.structure.identity import ExpandedAtom, SymmetryImageProvenance
from cristma.symmetry.orbit import expand_orbit

from .space_group import SpaceGroupSetting
from .wyckoff import AffineCoordinateMap, WyckoffPosition


@dataclass(frozen=True, slots=True)
class SiteSymmetry:
    """Symmetry operations that leave one site fixed modulo the lattice."""

    symbol: str | None
    stabilizer_operations: tuple[SymmetryImageProvenance, ...]

    def __post_init__(self) -> None:
        if self.symbol is not None and not self.symbol.strip():
            raise ValueError("site-symmetry symbol must not be empty")
        if not self.stabilizer_operations:
            raise ValueError("site stabilizer must contain at least one operation")

    @property
    def order(self) -> int:
        return len(self.stabilizer_operations)


@dataclass(frozen=True, slots=True)
class CrystallographicOrbit:
    """Finite orbit of one independent site in the reference cell."""

    representative: IndependentSite
    equivalent_sites: tuple[ExpandedAtom, ...]
    multiplicity: int
    stabilizer: tuple[SymmetryImageProvenance, ...]
    site_symmetry: SiteSymmetry

    def __post_init__(self) -> None:
        if self.multiplicity != len(self.equivalent_sites):
            raise ValueError("orbit multiplicity must equal its number of sites")
        if self.stabilizer != self.site_symmetry.stabilizer_operations:
            raise ValueError("orbit and site-symmetry stabilizers disagree")

    @property
    def calculated_multiplicity(self) -> int:
        """Explicit alias used when comparing against source-reported values."""

        return self.multiplicity


WyckoffStatus = Literal["matched", "unresolved", "ambiguous"]


@dataclass(frozen=True, slots=True)
class WyckoffAssignment:
    """Catalog match and source-validation diagnostics for one orbit."""

    position: WyckoffPosition | None
    calculated_multiplicity: int
    status: WyckoffStatus
    site_symmetry: SiteSymmetry
    diagnostics: tuple[Diagnostic, ...] = ()

    def __post_init__(self) -> None:
        if self.calculated_multiplicity <= 0:
            raise ValueError("calculated multiplicity must be positive")
        if self.status not in {"matched", "unresolved", "ambiguous"}:
            raise ValueError(f"unknown Wyckoff assignment status: {self.status!r}")
        if (self.position is None) == (self.status == "matched"):
            raise ValueError("matched status and Wyckoff position disagree")


def _periodically_equal(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
    tolerance: float,
) -> bool:
    return all(
        abs((a - b + 0.5) % 1.0 - 0.5) <= tolerance
        for a, b in zip(left, right, strict=True)
    )


def _as_float_map(
    constraint: AffineCoordinateMap,
) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.asarray(constraint.parameter_matrix, dtype=float),
        np.asarray(constraint.translation, dtype=float),
    )


def _fit_periodic_affine_map(
    constraint: AffineCoordinateMap,
    coordinate: tuple[float, float, float],
    tolerance: float,
) -> np.ndarray | None:
    matrix, translation = _as_float_map(constraint)
    observed = np.asarray(coordinate, dtype=float)
    for lattice_shift in product(range(-2, 3), repeat=3):
        target = observed + np.asarray(lattice_shift, dtype=float) - translation
        parameters, _, _, _ = np.linalg.lstsq(matrix, target, rcond=None)
        calculated = matrix @ parameters + translation
        if _periodically_equal(tuple(calculated), coordinate, tolerance):
            return parameters
    return None


def _unique_periodic_positions(
    positions: tuple[tuple[float, float, float], ...],
    tolerance: float,
) -> tuple[tuple[float, float, float], ...]:
    unique: list[tuple[float, float, float]] = []
    for position in positions:
        wrapped = tuple(float(value % 1.0) for value in position)
        if not any(_periodically_equal(wrapped, item, tolerance) for item in unique):
            unique.append(wrapped)
    return tuple(unique)


def _same_periodic_position_set(
    left: tuple[tuple[float, float, float], ...],
    right: tuple[tuple[float, float, float], ...],
    tolerance: float,
) -> bool:
    if len(left) != len(right):
        return False
    unmatched = list(right)
    for position in left:
        match = next(
            (
                index
                for index, candidate in enumerate(unmatched)
                if _periodically_equal(position, candidate, tolerance)
            ),
            None,
        )
        if match is None:
            return False
        unmatched.pop(match)
    return not unmatched


def _position_matches_orbit(
    position: WyckoffPosition,
    orbit: CrystallographicOrbit,
    tolerance: float,
) -> bool:
    if position.multiplicity != orbit.multiplicity:
        return False
    observed = tuple(float(value.value) % 1.0 for value in orbit.representative.fractional)
    actual_orbit = tuple(atom.fractional for atom in orbit.equivalent_sites)
    for representative in position.coordinate_constraints:
        parameters = _fit_periodic_affine_map(representative, observed, tolerance)
        if parameters is None:
            continue
        generated = tuple(
            tuple(
                float(value)
                for value in (matrix @ parameters + translation)
            )
            for constraint in position.coordinate_constraints
            for matrix, translation in (_as_float_map(constraint),)
        )
        unique = _unique_periodic_positions(generated, tolerance)
        if len(unique) != position.multiplicity:
            continue
        if _same_periodic_position_set(unique, actual_orbit, tolerance):
            return True
    return False


def build_orbit(
    site: IndependentSite,
    setting: SpaceGroupSetting,
    *,
    cell: UnitCell,
    tolerance: float = 1e-6,
    structure_id: str | None = None,
) -> CrystallographicOrbit:
    """Calculate an orbit and stabilizer from a site and one group setting."""

    if tolerance <= 0:
        raise ValueError("orbit tolerance must be positive")
    equivalent_sites = expand_orbit(
        site,
        setting.symmetry_operations,
        tolerance,
        cell=cell,
        structure_id=structure_id,
    )
    source = tuple(float(value.value) % 1.0 for value in site.fractional)
    representative = next(
        (
            atom
            for atom in equivalent_sites
            if _periodically_equal(atom.fractional, source, tolerance)
        ),
        None,
    )
    if representative is None:
        raise ValueError("orbit does not contain its representative site")

    stabilizer = representative.equivalent_images
    multiplicity = len(equivalent_sites)
    if len(setting.symmetry_operations) != multiplicity * len(stabilizer):
        raise ValueError("inconsistent orbit and stabilizer sizes")

    site_symmetry = SiteSymmetry(
        symbol=None,
        stabilizer_operations=stabilizer,
    )
    return CrystallographicOrbit(
        representative=site,
        equivalent_sites=equivalent_sites,
        multiplicity=multiplicity,
        stabilizer=stabilizer,
        site_symmetry=site_symmetry,
    )


def assign_wyckoff(
    orbit: CrystallographicOrbit,
    setting: SpaceGroupSetting,
    *,
    tolerance: float = 1e-6,
) -> WyckoffAssignment:
    """Match a calculated orbit to one setting-specific Wyckoff position."""

    if tolerance <= 0:
        raise ValueError("Wyckoff matching tolerance must be positive")
    matches = tuple(
        position
        for position in setting.wyckoff_positions
        if _position_matches_orbit(position, orbit, tolerance)
    )
    diagnostics: list[Diagnostic] = []
    if len(matches) == 1:
        position = matches[0]
        status: WyckoffStatus = "matched"
        site_symmetry = SiteSymmetry(
            symbol=position.site_symmetry_symbol,
            stabilizer_operations=orbit.stabilizer,
        )
    elif not matches:
        position = None
        status = "unresolved"
        site_symmetry = orbit.site_symmetry
        diagnostics.append(
            Diagnostic(
                Severity.WARNING,
                "crystallography.orbit.wyckoff_unresolved",
                "Calculated orbit does not match a catalog Wyckoff position.",
            )
        )
    else:
        position = None
        status = "ambiguous"
        site_symmetry = orbit.site_symmetry
        diagnostics.append(
            Diagnostic(
                Severity.WARNING,
                "crystallography.orbit.wyckoff_ambiguous",
                "Calculated orbit matches more than one catalog Wyckoff position.",
            )
        )

    site = orbit.representative
    if (
        site.reported_multiplicity is not None
        and site.reported_multiplicity != orbit.multiplicity
    ):
        diagnostics.append(
            Diagnostic(
                Severity.WARNING,
                "crystallography.orbit.reported_multiplicity_mismatch",
                f"Reported multiplicity {site.reported_multiplicity} differs from "
                f"calculated multiplicity {orbit.multiplicity}.",
            )
        )
    reported_letter = None
    if site.wyckoff is not None:
        compact = site.wyckoff.strip()
        if compact and compact[-1] in ascii_letters:
            reported_letter = compact[-1]
    reported_position = next(
        (
            candidate
            for candidate in setting.wyckoff_positions
            if candidate.letter == reported_letter
        ),
        None,
    )
    if (
        reported_position is not None
        and reported_position.multiplicity != orbit.multiplicity
    ):
        diagnostics.append(
            Diagnostic(
                Severity.WARNING,
                "crystallography.orbit.wyckoff_multiplicity_mismatch",
                f"Reported Wyckoff position {reported_position.letter!r} has catalog "
                f"multiplicity {reported_position.multiplicity}, but the calculated "
                f"orbit has multiplicity {orbit.multiplicity}.",
            )
        )
    if position is not None and reported_letter is not None and reported_letter != position.letter:
        diagnostics.append(
            Diagnostic(
                Severity.WARNING,
                "crystallography.orbit.reported_wyckoff_mismatch",
                f"Reported Wyckoff letter {site.wyckoff!r} differs from "
                f"calculated letter {position.letter!r}.",
            )
        )

    return WyckoffAssignment(
        position=position,
        calculated_multiplicity=orbit.multiplicity,
        status=status,
        site_symmetry=site_symmetry,
        diagnostics=tuple(diagnostics),
    )


__all__ = [
    "CrystallographicOrbit",
    "SiteSymmetry",
    "WyckoffAssignment",
    "assign_wyckoff",
    "build_orbit",
]
