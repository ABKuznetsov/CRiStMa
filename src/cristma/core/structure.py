"""Canonical asymmetric-unit crystal model."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from types import MappingProxyType
from typing import Mapping

from .cell import UnitCell
from .values import MeasuredValue


def _frozen_metadata(values: Mapping[str, object]) -> Mapping[str, object]:
    return MappingProxyType(dict(values))


@dataclass(frozen=True, slots=True)
class SiteComponent:
    """One chemical component occupying a crystallographic position."""

    element: str
    occupancy: MeasuredValue
    oxidation_state: MeasuredValue | None = None
    metadata: Mapping[str, object] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        if not self.element:
            raise ValueError("site component element must not be empty")
        occupancy = self.occupancy.value
        if occupancy is None or not math.isfinite(occupancy) or occupancy < 0:
            raise ValueError("site component occupancy must be finite and non-negative")
        object.__setattr__(self, "metadata", _frozen_metadata(self.metadata))


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
class Crystal:
    """A crystal whose independent sites are the only primary coordinates."""

    name: str
    cell: UnitCell
    sites: tuple[IndependentSite, ...]
    space_group: object | None = None
    formula: str | None = None
    metadata: Mapping[str, object] = field(default_factory=dict, compare=False)
    expanded_sites: tuple[object, ...] | None = field(
        default=None,
        compare=False,
        repr=False,
    )

    def __post_init__(self) -> None:
        object.__setattr__(self, "metadata", _frozen_metadata(self.metadata))
