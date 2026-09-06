"""Optional Shannon-radius evidence for already identified ionic contacts."""

from __future__ import annotations

from dataclasses import dataclass, replace
import math

from cristma.reference_data import ShannonRadiusRecord

from .models import EvidenceStatus


@dataclass(frozen=True, slots=True)
class ShannonDistanceCheck:
    status: EvidenceStatus
    observed_distance: float
    radius_sum: float
    minimum_distance: float
    distance_ratio: float
    minimum_ratio: float
    excludes_contact: bool = False


@dataclass(frozen=True, slots=True)
class ShannonDistanceValidator:
    """Report overlap evidence without changing contact selection."""

    minimum_ratio: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.minimum_ratio) or not 0 < self.minimum_ratio <= 1:
            raise ValueError("minimum_ratio must lie in (0, 1]")

    def get_config(self) -> dict[str, float]:
        return {"minimum_ratio": self.minimum_ratio}

    def clone(self, **changes: float) -> "ShannonDistanceValidator":
        return replace(self, **changes)

    def evaluate(
        self,
        *,
        distance: float,
        first: ShannonRadiusRecord,
        second: ShannonRadiusRecord,
    ) -> ShannonDistanceCheck:
        if not math.isfinite(distance) or distance <= 0:
            raise ValueError("distance must be positive and finite")
        radius_sum = first.ionic_radius + second.ionic_radius
        minimum_distance = self.minimum_ratio * radius_sum
        status = (
            EvidenceStatus.CONTRADICTORY
            if distance < minimum_distance
            else EvidenceStatus.SUPPORTIVE
        )
        return ShannonDistanceCheck(
            status=status,
            observed_distance=distance,
            radius_sum=radius_sum,
            minimum_distance=minimum_distance,
            distance_ratio=distance / radius_sum,
            minimum_ratio=self.minimum_ratio,
        )


__all__ = ["ShannonDistanceCheck", "ShannonDistanceValidator"]
