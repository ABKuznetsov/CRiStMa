"""Coordinate-free chemical composition."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from types import MappingProxyType
from typing import Mapping

from .elements import normalize_element


def _formula_order(symbols: tuple[str, ...]) -> tuple[str, ...]:
    if "C" in symbols:
        return tuple(
            symbol
            for symbol in ("C", "H", *sorted(set(symbols) - {"C", "H"}))
            if symbol in symbols
        )
    return tuple(sorted(symbols))


def _format_amount(value: float) -> str:
    if math.isclose(value, 1.0, abs_tol=1e-12):
        return ""
    nearest = round(value)
    return str(nearest) if math.isclose(value, nearest, abs_tol=1e-12) else f"{value:g}"


def _is_identity_only(structure: object) -> bool:
    space_group = getattr(structure, "space_group", None)
    operations = getattr(space_group, "operations", ())
    if len(operations) != 1:
        return False
    operation = operations[0].normalized()
    return operation.rotation == (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    ) and operation.translation == (Fraction(0), Fraction(0), Fraction(0))


@dataclass(frozen=True, slots=True)
class Composition:
    """Immutable amounts indexed by normalized element symbol."""

    _amounts: Mapping[str, float]
    normalization_basis: str = "reported_occupancy"

    def __post_init__(self) -> None:
        if not self._amounts:
            raise ValueError("composition must contain occupied elements")
        normalized: dict[str, float] = {}
        for raw_symbol, raw_amount in self._amounts.items():
            symbol = normalize_element(raw_symbol)
            amount = float(raw_amount)
            if not math.isfinite(amount) or amount <= 0:
                raise ValueError("composition amounts must be positive and finite")
            normalized[symbol] = normalized.get(symbol, 0.0) + amount
        object.__setattr__(self, "_amounts", MappingProxyType(dict(sorted(normalized.items()))))

    @classmethod
    def from_mapping(cls, values: Mapping[str, float]) -> "Composition":
        return cls(values)

    @classmethod
    def from_structure(cls, structure: object) -> "Composition":
        independent_sites = getattr(structure, "sites", None)
        if independent_sites is not None:
            positions = independent_sites
            identity_only = _is_identity_only(structure)
            use_multiplicity = True
        else:
            positions = getattr(structure, "atoms", None)
            if positions is None:
                raise TypeError("structure must expose sites or atoms")
            identity_only = True
            use_multiplicity = False

        amounts: dict[str, float] = {}
        for position in positions:
            if use_multiplicity:
                multiplicity = getattr(position, "calculated_multiplicity", None)
                if multiplicity is None:
                    if not identity_only:
                        raise ValueError(
                            f"site {getattr(position, 'label', '?')!r} requires calculated multiplicity"
                        )
                    multiplicity = 1
            else:
                multiplicity = 1
            for component in position.components:
                occupancy = float(component.occupancy.value)
                if occupancy <= 0:
                    continue
                symbol = component.species.require_element()
                amounts[symbol] = amounts.get(symbol, 0.0) + occupancy * multiplicity
        return cls(amounts, normalization_basis="site_occupancy_multiplicity")

    @property
    def elements(self) -> tuple[str, ...]:
        return tuple(self._amounts)

    @property
    def normalized_formula(self) -> str:
        return "".join(
            f"{symbol}{_format_amount(self._amounts[symbol])}"
            for symbol in _formula_order(self.elements)
        )

    def amount(self, symbol: str) -> float:
        return self._amounts[normalize_element(symbol)]

    def as_dict(self) -> dict[str, float]:
        return dict(self._amounts)


__all__ = ["Composition"]
