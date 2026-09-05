"""Immutable results shared by inorganic crystal-chemistry tools."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
import math

from cristma.chemistry import GrammarOperation, InteractionLayer, InteractionPriority
from cristma.chemistry.species import ChemicalSpecies
from cristma.crystallography import GeometricContact
from cristma.diagnostics import Diagnostic


class ResolutionStatus(StrEnum):
    RESOLVED = "resolved"
    AMBIGUOUS = "ambiguous"
    INCOMPLETE = "incomplete"
    NOT_APPLICABLE = "not_applicable"


class ContactClassification(StrEnum):
    PRIMARY = "primary"
    SECONDARY = "secondary"


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
        values = (
            self.first_occupancy,
            self.second_occupancy,
            self.radius_sum,
            self.normalized_distance,
            self.occupancy_weight,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("component interpretation values must be finite")
        if not 0 <= self.first_occupancy <= 1 or not 0 <= self.second_occupancy <= 1:
            raise ValueError("component occupancies must lie between zero and one")
        if self.radius_sum <= 0 or self.normalized_distance <= 0:
            raise ValueError("radius sum and normalized distance must be positive")
        expected_weight = self.first_occupancy * self.second_occupancy
        if not math.isclose(self.occupancy_weight, expected_weight, abs_tol=1e-12):
            raise ValueError("occupancy weight must equal the component occupancy product")
        if not self.centre_elements or not self.ligand_elements:
            raise ValueError("component interpretation requires centre and ligand views")

    @property
    def species_symbols(self) -> tuple[str | None, str | None]:
        return (self.first_species.element, self.second_species.element)


@dataclass(frozen=True, slots=True)
class ShellAlternative:
    boundary_group: int
    geometric_CN: int
    relative_gap: float
    internal_spread: float
    strong_contacts_outside: bool


@dataclass(frozen=True, slots=True)
class ResolvedContact:
    geometric_contact: GeometricContact
    interaction_type: GrammarOperation
    interaction_layer: InteractionLayer
    grammar_priority: InteractionPriority
    contact_classification: ContactClassification
    component_interpretations: tuple[ComponentPairInterpretation, ...]
    normalized_distance_min: float
    normalized_distance_max: float
    neighbor_total_occupancy: float
    evidence: tuple[SecondaryEvidence, ...]
    provenance: tuple[tuple[str, object], ...]
    contact_orbit_id: str = ""

    def __post_init__(self) -> None:
        if not self.component_interpretations:
            raise ValueError("resolved contact requires component interpretations")
        distances = tuple(
            item.normalized_distance for item in self.component_interpretations
        )
        if not math.isclose(self.normalized_distance_min, min(distances), abs_tol=1e-12) or not math.isclose(
            self.normalized_distance_max, max(distances), abs_tol=1e-12
        ):
            raise ValueError("normalized-distance bounds must cover every interpretation")
        if any(item.interaction_type is not self.interaction_type for item in self.component_interpretations):
            raise ValueError("component interpretations must share the contact interaction")
        if any(item.interaction_layer is not self.interaction_layer for item in self.component_interpretations):
            raise ValueError("component interpretations must share the interaction layer")
        if any(item.grammar_priority is not self.grammar_priority for item in self.component_interpretations):
            raise ValueError("component interpretations must share the grammar priority")
        if not math.isfinite(self.neighbor_total_occupancy) or not 0 <= self.neighbor_total_occupancy <= 1:
            raise ValueError("neighbor total occupancy must lie between zero and one")
        if not self.contact_orbit_id:
            object.__setattr__(
                self,
                "contact_orbit_id",
                f"contact-orbit:singleton:{self.geometric_contact.contact_id}",
            )


@dataclass(frozen=True, slots=True)
class ResolvedContactOrbit:
    """One symmetry-equivalent family of resolved periodic contacts."""

    contact_orbit_id: str
    representative_contact_id: str
    contacts: tuple[ResolvedContact, ...]

    def __post_init__(self) -> None:
        if not self.contact_orbit_id or not self.representative_contact_id:
            raise ValueError("contact orbit and representative IDs must not be empty")
        if not self.contacts:
            raise ValueError("resolved contact orbit must not be empty")
        contact_ids = tuple(item.geometric_contact.contact_id for item in self.contacts)
        if len(set(contact_ids)) != len(contact_ids):
            raise ValueError("contact orbit members must have unique contact IDs")
        if self.representative_contact_id != contact_ids[0]:
            raise ValueError("representative contact must be the first orbit member")
        if any(item.contact_orbit_id != self.contact_orbit_id for item in self.contacts):
            raise ValueError("contact orbit member has a different orbit ID")


@dataclass(frozen=True, slots=True)
class CoordinationShell:
    source_site_id: str
    center_atom_id: str
    contacts: tuple[ResolvedContact, ...]
    geometric_CN: int
    mean_occupied_neighbors: float
    status: ResolutionStatus
    alternatives: tuple[ShellAlternative, ...] = ()
    evidence: tuple[SecondaryEvidence, ...] = ()
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()

    def __post_init__(self) -> None:
        if not self.source_site_id or not self.center_atom_id:
            raise ValueError("coordination shell identities must not be empty")
        if self.geometric_CN < 0 or not math.isfinite(self.mean_occupied_neighbors):
            raise ValueError("coordination counts must be non-negative and finite")
        if self.status is ResolutionStatus.RESOLVED and self.geometric_CN != len(self.contacts):
            raise ValueError("resolved shell geometric CN must equal its contact count")

    @classmethod
    def resolved(
        cls,
        source_site_id: str,
        center_atom_id: str,
        contacts: tuple[ResolvedContact, ...],
    ) -> CoordinationShell:
        return cls(
            source_site_id=source_site_id,
            center_atom_id=center_atom_id,
            contacts=contacts,
            geometric_CN=len(contacts),
            mean_occupied_neighbors=math.fsum(
                item.neighbor_total_occupancy for item in contacts
            ),
            status=ResolutionStatus.RESOLVED,
        )

    @property
    def diagnostic_codes(self) -> tuple[str, ...]:
        return tuple(item.code for item in self.diagnostics)


@dataclass(frozen=True, slots=True)
class CrystalChemistryResolution:
    contacts: tuple[ResolvedContact, ...]
    coordination_shells: tuple[CoordinationShell, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()
    status: ResolutionStatus = field(
        default=ResolutionStatus.NOT_APPLICABLE,
        kw_only=True,
    )
    contact_orbits: tuple[ResolvedContactOrbit, ...] = field(
        default=(),
        kw_only=True,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.status, ResolutionStatus):
            raise TypeError("crystal-chemistry resolution status must be ResolutionStatus")
        if self.status is ResolutionStatus.NOT_APPLICABLE and (
            self.contacts or self.coordination_shells
        ):
            shell_statuses = {item.status for item in self.coordination_shells}
            derived = (
                ResolutionStatus.INCOMPLETE
                if ResolutionStatus.INCOMPLETE in shell_statuses
                else ResolutionStatus.AMBIGUOUS
                if ResolutionStatus.AMBIGUOUS in shell_statuses
                else ResolutionStatus.RESOLVED
            )
            object.__setattr__(self, "status", derived)
        if not self.contact_orbits:
            return
        orbit_ids = tuple(item.contact_orbit_id for item in self.contact_orbits)
        if len(set(orbit_ids)) != len(orbit_ids):
            raise ValueError(
                "contact orbit IDs must be unique and every contact must belong "
                "to exactly one orbit"
            )
        expected = tuple(item.geometric_contact.contact_id for item in self.contacts)
        observed = tuple(
            item.geometric_contact.contact_id
            for orbit in self.contact_orbits
            for item in orbit.contacts
        )
        if len(observed) != len(set(observed)) or set(observed) != set(expected):
            raise ValueError("every resolved contact must belong to exactly one orbit")

    @property
    def diagnostic_codes(self) -> tuple[str, ...]:
        return tuple(item.code for item in self.diagnostics)


__all__ = [
    "ComponentPairInterpretation",
    "ContactClassification",
    "CoordinationShell",
    "CrystalChemistryResolution",
    "EvidenceStatus",
    "ResolutionStatus",
    "ResolvedContact",
    "ResolvedContactOrbit",
    "SecondaryEvidence",
    "ShellAlternative",
]
