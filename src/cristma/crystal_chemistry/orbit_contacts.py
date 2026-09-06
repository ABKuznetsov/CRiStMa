"""Chemical interpretations calculated once per geometric pair orbit."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from enum import StrEnum
import math

from cristma.chemistry import (
    CandidateInteraction,
    CompositionGrammar,
    GrammarOperation,
    InteractionLayer,
    InteractionPriority,
)
from cristma.crystallography import PairTableStatus, SymmetryContactOrbit, SymmetryPairTable
from cristma.crystallography.symmetry_context import _digest
from cristma.diagnostics import Diagnostic, Severity
from cristma.reference_data import ReferenceData
from cristma.structure import CrystalStructure, SiteComponent

from .contacts import (
    ComponentPairInterpretation,
    EvidenceStatus,
    ResolutionStatus,
    SecondaryEvidence,
)
from .policy import ShellResolutionPolicy


class OrientationMode(StrEnum):
    UNDIRECTED = "undirected"
    ENDPOINT_ROLES = "endpoint_roles"


class EndpointRole(StrEnum):
    CENTER = "center"
    LIGAND = "ligand"
    NETWORK = "network"


EndpointRoles = tuple[EndpointRole, EndpointRole]

_ORIENTED_OPERATIONS = frozenset(
    {
        GrammarOperation.CENTRE_LIGAND_SHELL,
        GrammarOperation.INTERSTITIAL_COORDINATION,
        GrammarOperation.MIXED_ANION_COORDINATION,
    }
)


@dataclass(frozen=True, slots=True)
class ContactInterpretation:
    interpretation_id: str
    interaction_context_id: str
    interaction_type: GrammarOperation
    interaction_layer: InteractionLayer
    grammar_priority: InteractionPriority
    orientation_mode: OrientationMode
    endpoint_roles: EndpointRoles
    component_pair_interpretations: tuple[ComponentPairInterpretation, ...]
    normalized_distance_range: tuple[float, float] | None
    status: ResolutionStatus
    evidence: tuple[SecondaryEvidence, ...]

    def __post_init__(self) -> None:
        if not self.interpretation_id or not self.interaction_context_id:
            raise ValueError("contact interpretation identities must not be empty")
        if len(self.endpoint_roles) != 2:
            raise ValueError("contact interpretation requires two endpoint roles")
        if self.orientation_mode is OrientationMode.UNDIRECTED and self.endpoint_roles != (
            EndpointRole.NETWORK,
            EndpointRole.NETWORK,
        ):
            raise ValueError("undirected interpretations must use network endpoint roles")
        values = tuple(
            item.normalized_distance for item in self.component_pair_interpretations
        )
        if values:
            if self.normalized_distance_range is None:
                raise ValueError("known normalized distances require a range")
            expected = (min(values), max(values))
            if not all(
                math.isclose(observed, wanted, abs_tol=1e-12)
                for observed, wanted in zip(self.normalized_distance_range, expected, strict=True)
            ):
                raise ValueError("normalized-distance range does not cover component records")
        elif self.normalized_distance_range is not None:
            raise ValueError("an empty incomplete interpretation cannot report a distance range")
        if not values and self.status is not ResolutionStatus.INCOMPLETE:
            raise ValueError("an interpretation without radius evidence must be incomplete")


@dataclass(frozen=True, slots=True)
class ResolvedContactOrbit:
    resolved_contact_orbit_id: str
    geometry_orbit_id: str
    interpretations: tuple[ContactInterpretation, ...]
    status: ResolutionStatus
    diagnostics: tuple[Diagnostic, ...]
    provenance: tuple[tuple[str, object], ...]

    def __post_init__(self) -> None:
        if not self.resolved_contact_orbit_id or not self.geometry_orbit_id:
            raise ValueError("resolved and geometry orbit IDs must not be empty")
        if not self.interpretations:
            raise ValueError("resolved contact orbit requires an interpretation")
        ids = tuple(item.interpretation_id for item in self.interpretations)
        if len(set(ids)) != len(ids):
            raise ValueError("contact interpretation IDs must be unique")


@dataclass(frozen=True, slots=True)
class ContactOrbitResolution:
    contact_orbits: tuple[ResolvedContactOrbit, ...]
    status: ResolutionStatus
    diagnostics: tuple[Diagnostic, ...]
    provenance: tuple[tuple[str, object], ...]

    def __post_init__(self) -> None:
        ids = tuple(item.resolved_contact_orbit_id for item in self.contact_orbits)
        if len(set(ids)) != len(ids):
            raise ValueError("resolved contact orbit IDs must be unique")


def _matching_orientation(
    first: SiteComponent,
    second: SiteComponent,
    request: CandidateInteraction,
) -> bool:
    if first.element is None or second.element is None:
        return False
    return (
        first.element in request.first_elements
        and second.element in request.second_elements
    ) or (
        second.element in request.first_elements
        and first.element in request.second_elements
    )


def _endpoint_roles(
    first: SiteComponent,
    second: SiteComponent,
    request: CandidateInteraction,
) -> tuple[tuple[OrientationMode, EndpointRoles], ...]:
    if request.operation not in _ORIENTED_OPERATIONS:
        return ((OrientationMode.UNDIRECTED, (EndpointRole.NETWORK, EndpointRole.NETWORK)),)
    roles: list[tuple[OrientationMode, EndpointRoles]] = []
    if first.element in request.centre_elements and second.element in request.ligand_elements:
        roles.append((OrientationMode.ENDPOINT_ROLES, (EndpointRole.CENTER, EndpointRole.LIGAND)))
    if second.element in request.centre_elements and first.element in request.ligand_elements:
        roles.append((OrientationMode.ENDPOINT_ROLES, (EndpointRole.LIGAND, EndpointRole.CENTER)))
    return tuple(roles)


def _component_record(
    first: SiteComponent,
    second: SiteComponent,
    request: CandidateInteraction,
    distance: float,
    reference: ReferenceData,
) -> tuple[ComponentPairInterpretation | None, tuple[str, ...]]:
    radii: list[float] = []
    missing: list[str] = []
    for component in (first, second):
        try:
            radii.append(reference.covalent_radii.find(component.element).value)
        except (KeyError, TypeError):
            missing.append(component.element or "unknown")
    if missing:
        return None, tuple(sorted(set(missing)))
    radius_sum = radii[0] + radii[1]
    rho = distance / radius_sum
    return (
        ComponentPairInterpretation(
            first_species=first.species,
            second_species=second.species,
            first_occupancy=float(first.occupancy.value),
            second_occupancy=float(second.occupancy.value),
            radius_sum=radius_sum,
            normalized_distance=rho,
            occupancy_weight=float(first.occupancy.value) * float(second.occupancy.value),
            interaction_type=request.operation,
            interaction_layer=request.layer,
            grammar_priority=request.priority,
            centre_elements=request.centre_elements,
            ligand_elements=request.ligand_elements,
        ),
        (),
    )


class ContactOrbitResolver:
    """Resolve chemistry without materializing symmetry-expanded contacts."""

    def __init__(
        self,
        policy: ShellResolutionPolicy,
        reference: ReferenceData | None = None,
    ) -> None:
        if not isinstance(policy, ShellResolutionPolicy):
            raise TypeError("policy must be ShellResolutionPolicy")
        self.policy = policy
        self.reference = reference or ReferenceData.default()

    def _interpret_orbit(
        self,
        orbit: SymmetryContactOrbit,
        first_components: tuple[SiteComponent, ...],
        second_components: tuple[SiteComponent, ...],
        grammar: CompositionGrammar,
    ) -> tuple[tuple[ContactInterpretation, ...], tuple[Diagnostic, ...]]:
        output: list[ContactInterpretation] = []
        diagnostics: list[Diagnostic] = []
        for request in grammar.candidate_interactions:
            records: dict[
                tuple[OrientationMode, EndpointRoles],
                list[ComponentPairInterpretation],
            ] = defaultdict(list)
            missing_by_orientation: dict[
                tuple[OrientationMode, EndpointRoles],
                set[str],
            ] = defaultdict(set)
            order: list[tuple[OrientationMode, EndpointRoles]] = []
            for first in first_components:
                for second in second_components:
                    if not _matching_orientation(first, second, request):
                        continue
                    orientations = _endpoint_roles(first, second, request)
                    for orientation in orientations:
                        if orientation not in order:
                            order.append(orientation)
                        record, missing = _component_record(
                            first,
                            second,
                            request,
                            orbit.representative_distance,
                            self.reference,
                        )
                        if missing:
                            missing_by_orientation[orientation].update(missing)
                        elif record is not None and record.normalized_distance <= self.policy.search_rho_max:
                            records[orientation].append(record)

            for orientation in order:
                component_records = tuple(
                    sorted(
                        set(records[orientation]),
                        key=lambda item: (
                            item.first_species.label,
                            item.second_species.label,
                            item.first_occupancy,
                            item.second_occupancy,
                            item.normalized_distance,
                        ),
                    )
                )
                missing = tuple(sorted(missing_by_orientation[orientation]))
                if not component_records and not missing:
                    continue
                mode, roles = orientation
                status = ResolutionStatus.INCOMPLETE if missing else ResolutionStatus.RESOLVED
                evidence: list[SecondaryEvidence] = []
                if component_records:
                    evidence.append(
                        SecondaryEvidence(
                            "covalent_radii",
                            EvidenceStatus.SUPPORTIVE,
                            "normalized distances use tabulated covalent radii",
                        )
                    )
                for symbol in missing:
                    diagnostics.append(
                        Diagnostic(
                            Severity.WARNING,
                            "crystal_chemistry.contact.radius_missing",
                            f"No covalent radius for {symbol} in a {request.operation.value} component pair",
                        )
                    )
                    evidence.append(
                        SecondaryEvidence(
                            "covalent_radii",
                            EvidenceStatus.NOT_AVAILABLE,
                            f"covalent radius is unavailable for {symbol}",
                        )
                    )
                normalized_range = (
                    (
                        min(item.normalized_distance for item in component_records),
                        max(item.normalized_distance for item in component_records),
                    )
                    if component_records
                    else None
                )
                context_id = "contact-interaction-context:" + _digest(
                    {
                        "operation": request.operation.value,
                        "layer": request.layer.value,
                        "priority": request.priority.value,
                        "first_elements": request.first_elements,
                        "second_elements": request.second_elements,
                        "centre_elements": request.centre_elements,
                        "ligand_elements": request.ligand_elements,
                    }
                )
                interpretation_id = "contact-interpretation:" + _digest(
                    {
                        "geometry_orbit_id": orbit.geometry_orbit_id,
                        "interaction_context_id": context_id,
                        "orientation_mode": mode.value,
                        "endpoint_roles": tuple(role.value for role in roles),
                        "centre_elements": request.centre_elements,
                        "ligand_elements": request.ligand_elements,
                    }
                )
                output.append(
                    ContactInterpretation(
                        interpretation_id,
                        context_id,
                        request.operation,
                        request.layer,
                        request.priority,
                        mode,
                        roles,
                        component_records,
                        normalized_range,
                        status,
                        tuple(evidence),
                    )
                )
        return tuple(output), tuple(dict.fromkeys(diagnostics))

    def resolve(
        self,
        pair_table: SymmetryPairTable,
        structure: CrystalStructure,
        grammar: CompositionGrammar,
    ) -> ContactOrbitResolution:
        if not isinstance(pair_table, SymmetryPairTable):
            raise TypeError("pair_table must be SymmetryPairTable")
        if not isinstance(structure, CrystalStructure):
            raise TypeError("structure must be CrystalStructure")
        if not isinstance(grammar, CompositionGrammar):
            raise TypeError("grammar must be CompositionGrammar")
        sites = {site.id: site for site in structure.sites}
        output: list[ResolvedContactOrbit] = []
        diagnostics: list[Diagnostic] = list(grammar.diagnostics)
        for orbit in pair_table.contact_orbits:
            try:
                first = sites[orbit.first_independent_site_id]
                second = sites[orbit.second_independent_site_id]
            except KeyError as exc:
                raise ValueError("pair table references a site outside the structure") from exc
            interpretations, orbit_diagnostics = self._interpret_orbit(
                orbit,
                first.components,
                second.components,
                grammar,
            )
            diagnostics.extend(orbit_diagnostics)
            if not interpretations:
                continue
            oriented_role_sets: dict[tuple[object, ...], set[EndpointRoles]] = defaultdict(set)
            for item in interpretations:
                if item.orientation_mode is OrientationMode.ENDPOINT_ROLES:
                    oriented_role_sets[
                        (item.interaction_type, item.interaction_layer, item.grammar_priority)
                    ].add(item.endpoint_roles)
            ambiguous = any(len(roles) > 1 for roles in oriented_role_sets.values())
            orbit_status = (
                ResolutionStatus.INCOMPLETE
                if pair_table.status is PairTableStatus.INCOMPLETE
                or any(item.status is ResolutionStatus.INCOMPLETE for item in interpretations)
                else ResolutionStatus.AMBIGUOUS
                if ambiguous
                else ResolutionStatus.RESOLVED
            )
            resolved_id = "resolved-contact-orbit:" + _digest(
                {
                    "geometry_orbit_id": orbit.geometry_orbit_id,
                    "interpretation_ids": tuple(
                        item.interpretation_id for item in interpretations
                    ),
                }
            )
            output.append(
                ResolvedContactOrbit(
                    resolved_id,
                    orbit.geometry_orbit_id,
                    interpretations,
                    orbit_status,
                    orbit_diagnostics,
                    (
                        ("method", "cristma.orbit_contact_resolution:1"),
                        ("reference_version", grammar.reference_version),
                    ),
                )
            )
        aggregate = (
            ResolutionStatus.INCOMPLETE
            if pair_table.status is PairTableStatus.INCOMPLETE
            or any(item.status is ResolutionStatus.INCOMPLETE for item in output)
            else ResolutionStatus.AMBIGUOUS
            if any(item.status is ResolutionStatus.AMBIGUOUS for item in output)
            else ResolutionStatus.RESOLVED
            if output
            else ResolutionStatus.NOT_APPLICABLE
        )
        return ContactOrbitResolution(
            tuple(output),
            aggregate,
            tuple(dict.fromkeys(diagnostics)),
            (
                ("method", "cristma.orbit_contact_resolution:1"),
                ("grammar_method", f"{grammar.method_id}:{grammar.method_version}"),
                ("reference_version", grammar.reference_version),
            ),
        )


__all__ = [
    "ContactInterpretation",
    "ContactOrbitResolution",
    "ContactOrbitResolver",
    "EndpointRole",
    "EndpointRoles",
    "OrientationMode",
    "ResolvedContactOrbit",
]
