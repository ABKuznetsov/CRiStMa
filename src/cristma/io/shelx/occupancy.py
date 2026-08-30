"""Symbolic SHELX free-variable occupancy expressions."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True, slots=True)
class ShelxOccupancyExpression:
    """A reported SHELX site-occupation factor and its FVAR dependency."""

    raw: str
    free_variable_index: int | None
    multiplier: float
    complement: bool

    @classmethod
    def parse(cls, token: str) -> ShelxOccupancyExpression:
        """Decode fixed, FVAR, and complement forms without evaluating them."""

        raw = str(token).strip()
        try:
            reported = float(raw)
        except ValueError as error:
            raise ValueError(f"invalid SHELX occupancy expression: {raw!r}") from error
        if not math.isfinite(reported):
            raise ValueError(f"invalid SHELX occupancy expression: {raw!r}")

        magnitude = abs(reported)
        control = int(magnitude // 10.0)
        multiplier = magnitude - 10.0 * control if control else magnitude
        if magnitude > 15.0:
            if control < 2 or not 0.0 <= multiplier < 5.0:
                raise ValueError(f"invalid SHELX occupancy expression: {raw!r}")
            return cls(raw, control, multiplier, reported < 0.0)
        return cls(raw, None, multiplier, False)

    def evaluate(self, free_variables: tuple[float, ...]) -> float:
        """Evaluate physical occupancy against the complete FVAR instruction."""

        if self.free_variable_index is None:
            occupancy = self.multiplier
        else:
            offset = self.free_variable_index - 1
            if offset >= len(free_variables):
                raise ValueError(
                    f"SHELX occupancy references absent FVAR {self.free_variable_index}"
                )
            free_value = float(free_variables[offset])
            occupancy = self.multiplier * (
                1.0 - free_value if self.complement else free_value
            )
        if not math.isfinite(occupancy) or not 0.0 <= occupancy <= 1.0:
            raise ValueError(
                f"SHELX occupancy evaluates outside the physical range: {occupancy!r}"
            )
        return occupancy


__all__ = ["ShelxOccupancyExpression"]
