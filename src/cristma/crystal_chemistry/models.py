"""Small immutable value types shared by orbit-first crystal chemistry."""
from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
import math

from cristma.chemistry import GrammarOperation, InteractionLayer, InteractionPriority
from cristma.chemistry.species import ChemicalSpecies


class ResolutionStatus(StrEnum):
    RESOLVED = "resolved"
    AMBIGUOUS = "ambiguous"
    INCOMPLETE = "incomplete"
    NOT_APPLICABLE = "not_applicable"


class EvidenceStatus(StrEnum):
    SUPPORTIVE = "supportive"
    NEUTRAL = "neutral"
    CONTRADICTORY = "contradictory"
    NOT_AVAILABLE = "not_available"
    NOT_APPLICABLE = "not_applicable"


@dataclass(frozen=True, slots=True)
class SecondaryEvidence:
    method: str
    status: EvidenceStatus
    message: str

    def __post_init__(self) -> None:
        if not self.method or not self.message:
            raise ValueError("secondary evidence requires method and message")


@dataclass(frozen=True, slots=True)
class ComponentPairInterpretation:
    first_species: ChemicalSpecies
    second_species: ChemicalSpecies
    first_occupancy: float
    second_occupancy: float
    radius_sum: float
    normalized_distance: float
    occupancy_weight: float
    interaction_type: GrammarOperation
    interaction_layer: InteractionLayer
    grammar_priority: InteractionPriority
    centre_elements: tuple[str, ...]
    ligand_elements: tuple[str, ...]

    def __post_init__(self) -> None:
        values = (self.first_occupancy, self.second_occupancy, self.radius_sum,
                  self.normalized_distance, self.occupancy_weight)
        if not all(math.isfinite(value) for value in values):
            raise ValueError("component interpretation values must be finite")
        if not 0 <= self.first_occupancy <= 1 or not 0 <= self.second_occupancy <= 1:
            raise ValueError("component occupancies must lie between zero and one")
        if self.radius_sum <= 0 or self.normalized_distance <= 0:
            raise ValueError("radius sum and normalized distance must be positive")
        if not math.isclose(self.occupancy_weight,
                            self.first_occupancy * self.second_occupancy, abs_tol=1e-12):
            raise ValueError("occupancy weight must equal component occupancy product")
        if not self.centre_elements or not self.ligand_elements:
            raise ValueError("component interpretation requires centre and ligand views")

    @property
    def species_symbols(self) -> tuple[str | None, str | None]:
        return self.first_species.element, self.second_species.element


__all__ = ["ComponentPairInterpretation", "EvidenceStatus", "ResolutionStatus", "SecondaryEvidence"]
