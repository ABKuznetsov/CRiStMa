"""Chemical occupation models shared by structural position types."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from types import MappingProxyType
from typing import Mapping

from cristma.chemistry.species import ChemicalSpecies, as_species
from cristma.core.values import MeasuredValue


@dataclass(frozen=True, slots=True)
class SiteComponent:
    """One chemical component occupying a geometric position."""

    species: ChemicalSpecies | str
    occupancy: MeasuredValue
    oxidation_state: MeasuredValue | None = None
    metadata: Mapping[str, object] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "species", as_species(self.species))
        occupancy = self.occupancy.value
        if occupancy is None or not math.isfinite(occupancy) or not 0 <= occupancy <= 1:
            raise ValueError("site component occupancy must lie between zero and one")
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    @property
    def element(self) -> str | None:
        """Return the underlying element symbol when known."""

        return self.species.element


__all__ = ["SiteComponent"]
