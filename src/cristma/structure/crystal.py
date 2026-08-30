"""Canonical periodic crystal structure model."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from types import MappingProxyType
from typing import TYPE_CHECKING, Mapping

from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue

from .identity import StructureProvenance
from .occupation import SiteComponent

if TYPE_CHECKING:
    from cristma.symmetry.orbit import SpaceGroupDefinition
    from .view import AtomicView


def _frozen_metadata(values: Mapping[str, object]) -> Mapping[str, object]:
    return MappingProxyType(dict(values))


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
        if self.total_occupancy > 1.0 + 1e-12:
            raise ValueError(f"site occupancy exceeds one: {self.total_occupancy}")
        for coordinate in self.fractional:
            if coordinate.value is None or not math.isfinite(coordinate.value):
                raise ValueError("fractional coordinate must be finite")
        object.__setattr__(self, "metadata", _frozen_metadata(self.metadata))

    @property
    def total_occupancy(self) -> float:
        """Return the total chemical occupation of this geometric position."""

        return math.fsum(float(component.occupancy.value) for component in self.components)

    @property
    def vacancy_fraction(self) -> float:
        """Return the unoccupied fraction without creating a vacancy atom."""

        return max(0.0, 1.0 - self.total_occupancy)


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

    def atomic_view(self, *, expanded: bool = True) -> AtomicView:
        """Expose independent or symmetry-expanded sites as numerical rows."""

        if not expanded:
            raise ValueError("independent sites are not an AtomicView; request expanded=True")
        from cristma.symmetry.orbit import expand_structure

        return expand_structure(self)


Crystal = CrystalStructure


__all__ = ["Crystal", "CrystalStructure", "DisplacementParameters", "IndependentSite", "SiteComponent"]
