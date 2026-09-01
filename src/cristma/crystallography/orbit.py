"""Calculated crystallographic orbits and site stabilizers."""

from __future__ import annotations

from dataclasses import dataclass

from cristma.core.cell import UnitCell
from cristma.structure.crystal import IndependentSite
from cristma.structure.identity import ExpandedAtom, SymmetryImageProvenance
from cristma.symmetry.orbit import expand_orbit

from .space_group import SpaceGroupSetting


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


def _periodically_equal(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
    tolerance: float,
) -> bool:
    return all(
        abs((a - b + 0.5) % 1.0 - 0.5) <= tolerance
        for a, b in zip(left, right, strict=True)
    )


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


__all__ = ["CrystallographicOrbit", "SiteSymmetry", "build_orbit"]
