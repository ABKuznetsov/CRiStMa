"""Small immutable facade over the default CRiStMa reference catalogs."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from .chemical_reference import ChemicalReference, load_chemical_reference
from .elements import ElementCatalog
from .radii import CovalentRadii


@dataclass(frozen=True, slots=True)
class ReferenceData:
    elements: ElementCatalog
    covalent_radii: CovalentRadii
    chemical: ChemicalReference

    @classmethod
    @lru_cache(maxsize=1)
    def default(cls) -> "ReferenceData":
        return cls(
            elements=ElementCatalog.default(),
            covalent_radii=CovalentRadii.default(),
            chemical=load_chemical_reference(),
        )


__all__ = ["ReferenceData"]
