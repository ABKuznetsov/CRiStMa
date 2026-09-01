"""Small dependency-free element reference catalog."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from functools import lru_cache

from cristma.chemistry.elements import (
    ELEMENT_SYMBOLS_BY_ATOMIC_NUMBER,
    normalize_element,
)


class ElementCategory(StrEnum):
    """Broad composition-level chemical category."""

    METAL = "metal"
    METALLOID = "metalloid"
    NONMETAL = "nonmetal"
    NOBLE_GAS = "noble_gas"


_METALLOIDS = frozenset({"B", "Si", "Ge", "As", "Sb", "Te"})
_NOBLE_GASES = frozenset({"He", "Ne", "Ar", "Kr", "Xe", "Rn", "Og"})
_NONMETALS = frozenset({"H", "C", "N", "O", "F", "P", "S", "Se", "Cl", "Br", "I", "At", "Ts"})


@dataclass(frozen=True, slots=True)
class ElementRecord:
    """Element identity plus the broad category used by Chemistry."""

    symbol: str
    atomic_number: int
    category: ElementCategory
    dataset_id: str = "cristma.elements"
    dataset_version: str = "1"

    @property
    def is_metal(self) -> bool:
        return self.category is ElementCategory.METAL


class ElementCatalog:
    """Immutable lookup of all IUPAC element identities."""

    def __init__(self, records: tuple[ElementRecord, ...]) -> None:
        self._records = records
        self._by_symbol = {record.symbol: record for record in records}

    @classmethod
    @lru_cache(maxsize=1)
    def default(cls) -> "ElementCatalog":
        records = []
        for atomic_number, symbol in enumerate(ELEMENT_SYMBOLS_BY_ATOMIC_NUMBER, start=1):
            if symbol in _METALLOIDS:
                category = ElementCategory.METALLOID
            elif symbol in _NOBLE_GASES:
                category = ElementCategory.NOBLE_GAS
            elif symbol in _NONMETALS:
                category = ElementCategory.NONMETAL
            else:
                category = ElementCategory.METAL
            records.append(ElementRecord(symbol, atomic_number, category))
        return cls(tuple(records))

    def by_symbol(self, symbol: str) -> ElementRecord:
        return self._by_symbol[normalize_element(symbol)]

    @property
    def records(self) -> tuple[ElementRecord, ...]:
        return self._records


__all__ = ["ElementCatalog", "ElementCategory", "ElementRecord"]
