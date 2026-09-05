"""Typed chemical species used by CrIStMa's canonical structure model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from .elements import normalize_element


@runtime_checkable
class ChemicalSpecies(Protocol):
    """Common runtime contract for known and explicitly unknown species."""

    @property
    def label(self) -> str: ...

    @property
    def element(self) -> str | None: ...

    def require_element(self) -> str:
        """Return the element symbol or fail when no element is known."""


@dataclass(frozen=True, slots=True)
class ElementSpecies:
    """A chemical element without isotope or charge qualification."""

    symbol: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "symbol", normalize_element(self.symbol))

    @property
    def label(self) -> str:
        return self.symbol

    @property
    def element(self) -> str:
        return self.symbol

    def require_element(self) -> str:
        return self.symbol


@dataclass(frozen=True, slots=True)
class IsotopeSpecies:
    """An isotope with a known parent element."""

    element_symbol: str
    mass_number: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "element_symbol", normalize_element(self.element_symbol))
        if isinstance(self.mass_number, bool) or not isinstance(self.mass_number, int):
            raise TypeError("mass_number must be an integer")
        if self.mass_number <= 0:
            raise ValueError("mass_number must be positive")

    @property
    def label(self) -> str:
        return f"{self.mass_number}{self.element_symbol}"

    @property
    def element(self) -> str:
        return self.element_symbol

    def require_element(self) -> str:
        return self.element_symbol


@dataclass(frozen=True, slots=True)
class ChargedSpecies:
    """An element annotated with an integral charge or oxidation state."""

    element_symbol: str
    charge: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "element_symbol", normalize_element(self.element_symbol))
        if isinstance(self.charge, bool) or not isinstance(self.charge, int):
            raise TypeError("charge must be an integer")

    @property
    def label(self) -> str:
        if self.charge == 0:
            return self.element_symbol
        sign = "+" if self.charge > 0 else "-"
        magnitude = abs(self.charge)
        return f"{self.element_symbol}{magnitude}{sign}"

    @property
    def element(self) -> str:
        return self.element_symbol

    def require_element(self) -> str:
        return self.element_symbol


@dataclass(frozen=True, slots=True)
class UnknownSpecies:
    """A source species that cannot yet be assigned to a known element."""

    id: str
    source_label: str | None = None

    def __post_init__(self) -> None:
        if not self.id.strip():
            raise ValueError("unknown species id must not be empty")

    @property
    def label(self) -> str:
        return self.source_label or self.id

    @property
    def element(self) -> None:
        return None

    def require_element(self) -> str:
        raise ValueError(f"Species {self.label!r} has no known element")


def as_species(value: str | ChemicalSpecies) -> ChemicalSpecies:
    """Coerce a legacy element string to a typed chemical species."""

    if isinstance(value, str):
        return ElementSpecies(value)
    if isinstance(value, ChemicalSpecies):
        return value
    raise TypeError(f"Expected a chemical species or element string, got {type(value).__name__}")


__all__ = [
    "ChargedSpecies",
    "ChemicalSpecies",
    "ElementSpecies",
    "IsotopeSpecies",
    "UnknownSpecies",
    "as_species",
]
