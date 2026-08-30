"""Canonical periodic crystal structure model."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from types import MappingProxyType
from typing import TYPE_CHECKING, Mapping

from cristma.chemistry.species import ChemicalSpecies, as_species
from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue

from .identity import ExpandedAtomRef, StructureProvenance

if TYPE_CHECKING:
    from cristma.symmetry.orbit import SpaceGroupDefinition


def _frozen_metadata(values: Mapping[str, object]) -> Mapping[str, object]:
    return MappingProxyType(dict(values))


@dataclass(frozen=True, slots=True)
class SiteComponent:
    """One chemical component occupying a crystallographic position."""

    species: ChemicalSpecies | str
    occupancy: MeasuredValue
    oxidation_state: MeasuredValue | None = None
    metadata: Mapping[str, object] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "species", as_species(self.species))
        occupancy = self.occupancy.value
        if occupancy is None or not math.isfinite(occupancy) or occupancy < 0:
            raise ValueError("site component occupancy must be finite and non-negative")
        object.__setattr__(self, "metadata", _frozen_metadata(self.metadata))

    @property
    def element(self) -> str | None:
        return self.species.element


@dataclass(frozen=True, slots=True)
class DisplacementParameters:
    """Reported isotropic or anisotropic atomic displacement data."""

    kind: str
    isotropic: MeasuredValue | None = None
    tensor: tuple[tuple[MeasuredValue, MeasuredValue, MeasuredValue], ...] | None = None
    reported_kind: str | None = None


@dataclass(frozen=True, slots=True)
class IndependentSite:
    """A single independent crystallographic position."""

    id: str
    label: str
    components: tuple[SiteComponent, ...]
    fractional: tuple[MeasuredValue, MeasuredValue, MeasuredValue]
    wyckoff: str | None = None
    reported_multiplicity: int | None = None
    calculated_multiplicity: int | None = None
    disorder_assembly: str | None = None
    disorder_group: str | None = None
    displacement: DisplacementParameters | None = None
    metadata: Mapping[str, object] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        if not self.id or not self.label:
            raise ValueError("site id and label must not be empty")
        if not self.components:
            raise ValueError("site must contain at least one component")
        occupancy = math.fsum(
            float(component.occupancy.value)
            for component in self.components
            if component.occupancy.value is not None
        )
        if occupancy > 1.0 + 1e-6:
            raise ValueError(f"site occupancy exceeds one: {occupancy}")
        for coordinate in self.fractional:
            if coordinate.value is None or not math.isfinite(coordinate.value):
                raise ValueError("fractional coordinate must be finite")
        object.__setattr__(self, "metadata", _frozen_metadata(self.metadata))


@dataclass(frozen=True, slots=True)
class CrystalStructure:
    """A periodic structure whose independent sites are primary coordinates."""

    name: str
    cell: UnitCell
    sites: tuple[IndependentSite, ...]
    id: str | None = None
    space_group: SpaceGroupDefinition | None = None
    formula: str | None = None
    periodic: tuple[bool, bool, bool] = (True, True, True)
    provenance: StructureProvenance = field(default_factory=StructureProvenance)
    metadata: Mapping[str, object] = field(default_factory=dict, compare=False)
    expanded_sites: tuple[ExpandedAtomRef, ...] | None = field(default=None, compare=False, repr=False)

    def __post_init__(self) -> None:
        if not any(self.periodic):
            raise ValueError("crystal structure must be periodic along at least one axis")
        object.__setattr__(self, "metadata", _frozen_metadata(self.metadata))

    @classmethod
    def explicit(
        cls,
        name: str,
        cell: UnitCell,
        sites: tuple[IndependentSite, ...],
        **kwargs: object,
    ) -> CrystalStructure:
        """Build a structure whose source reports sites but no symmetry."""

        from cristma.symmetry.affine import parse_xyz_operation
        from cristma.symmetry.orbit import SpaceGroupDefinition

        space_group = SpaceGroupDefinition(
            operations=(parse_xyz_operation("x,y,z", operation_id="op:1"),),
            provenance="unreported_identity",
        )
        return cls(name, cell, sites, space_group=space_group, **kwargs)


Crystal = CrystalStructure


__all__ = ["Crystal", "CrystalStructure", "DisplacementParameters", "IndependentSite", "SiteComponent"]
