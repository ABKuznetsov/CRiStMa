"""Curated covalent radii used by structure-search tools."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from types import MappingProxyType
from typing import Mapping

from cristma.chemistry.elements import normalize_element


_CRAFT_COVALENT_RADII = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "K": 2.03, "Ca": 1.76, "Ti": 1.60, "Cr": 1.39, "Mn": 1.39,
    "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Br": 1.20,
    "Sr": 1.95, "Zr": 1.75, "Mo": 1.54, "Ag": 1.45, "I": 1.39, "Ba": 2.15,
    "W": 1.62, "Pt": 1.36, "Au": 1.36, "Pb": 1.46, "U": 1.96,
}


@dataclass(frozen=True, slots=True)
class CovalentRadiusRecord:
    symbol: str
    value: float
    unit: str = "angstrom"
    dataset_id: str = "cristma.covalent_radii.craft"
    dataset_version: str = "1"


class CovalentRadii:
    """Exact lookup with no guessed fallback for missing elements."""

    def __init__(self, records: Mapping[str, CovalentRadiusRecord]) -> None:
        self._records = MappingProxyType(dict(records))

    @classmethod
    @lru_cache(maxsize=1)
    def default(cls) -> "CovalentRadii":
        return cls({
            symbol: CovalentRadiusRecord(symbol, value)
            for symbol, value in _CRAFT_COVALENT_RADII.items()
        })

    def find(self, symbol: str) -> CovalentRadiusRecord:
        normalized = normalize_element(symbol)
        try:
            return self._records[normalized]
        except KeyError as exc:
            raise KeyError(f"No covalent radius for {normalized}") from exc


__all__ = ["CovalentRadii", "CovalentRadiusRecord"]
