"""Explicit immutable policies for coordination-shell resolution."""

from __future__ import annotations

from dataclasses import dataclass, replace
import math


@dataclass(frozen=True, slots=True)
class ShellResolutionPolicy:
    candidate_rho_max: float
    distance_group_tolerance: float
    minimum_shell_gap: float
    ambiguity_tolerance: float
    search_rho_max: float

    def __post_init__(self) -> None:
        for name, value in self.get_config().items():
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"{name} must be a positive finite number")
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be a positive finite number")
        if self.search_rho_max <= self.candidate_rho_max:
            raise ValueError("search_rho_max must exceed candidate_rho_max")

    def get_config(self) -> dict[str, float]:
        return {
            "candidate_rho_max": self.candidate_rho_max,
            "distance_group_tolerance": self.distance_group_tolerance,
            "minimum_shell_gap": self.minimum_shell_gap,
            "ambiguity_tolerance": self.ambiguity_tolerance,
            "search_rho_max": self.search_rho_max,
        }

    def clone(self, **changes: float) -> ShellResolutionPolicy:
        return replace(self, **changes)


__all__ = ["ShellResolutionPolicy"]
