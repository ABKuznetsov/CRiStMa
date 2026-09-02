"""Inorganic contact interpretation and coordination-shell resolution."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import statistics
from types import MappingProxyType

from cristma.chemistry import CandidateInteraction, CompositionGrammar, InteractionPriority
from cristma.chemistry.grammar import GrammarOperation
from cristma.crystallography import GeometricContact, geometric_contacts
from cristma.diagnostics import Diagnostic, Severity
from cristma.geometry import NeighborFinder
from cristma.reference_data import ReferenceData
from cristma.structure import AtomicView, CrystalStructure, SiteComponent

from .contacts import (
    ComponentPairInterpretation,
    ContactClassification,
    CoordinationShell,
    CrystalChemistryResolution,
    EvidenceStatus,
    ResolutionStatus,
    SecondaryEvidence,
    ShellAlternative,
    ResolvedContact,
)
from .policy import ShellResolutionPolicy


@dataclass(frozen=True, slots=True)
class InteractionScope:
    operation: GrammarOperation
    priority: InteractionPriority
    centre_elements: tuple[str, ...]
    ligand_elements: tuple[str, ...]

    @classmethod
    def from_request(cls, request: CandidateInteraction) -> InteractionScope:
        return cls(
            request.operation, request.priority,
            request.centre_elements, request.ligand_elements,
        )


@dataclass(frozen=True, slots=True)
class InterpretationOutcome:
    interpretations: tuple[ComponentPairInterpretation, ...]
    diagnostics: tuple[Diagnostic, ...]
    incomplete_scopes: tuple[InteractionScope, ...] = ()

    @property
    def incomplete_interactions(self) -> tuple[GrammarOperation, ...]:
        return tuple(
            sorted({item.operation for item in self.incomplete_scopes}, key=lambda item: item.value)
        )


@dataclass(frozen=True, slots=True)
class BoundaryDecision:
    status: ResolutionStatus
    selected: ShellAlternative | None
    alternatives: tuple[ShellAlternative, ...]
    evidence: tuple[SecondaryEvidence, ...]
    diagnostics: tuple[Diagnostic, ...]


def _group_distances(
    values: tuple[float, ...],
    tolerance: float,
) -> tuple[tuple[float, ...], ...]:
    ordered = tuple(sorted(values))
    if not ordered:
        return ()
    groups: list[list[float]] = [[ordered[0]]]
    for value in ordered[1:]:
        if value - groups[-1][-1] <= tolerance + 1e-12:
            groups[-1].append(value)
        else:
            groups.append([value])
    return tuple(tuple(group) for group in groups)


def _candidate_boundaries(
    groups: tuple[tuple[float, ...], ...],
    policy: ShellResolutionPolicy,
) -> tuple[ShellAlternative, ...]:
    candidates: list[ShellAlternative] = []
    inside: list[float] = []
    for index, group in enumerate(groups[:-1]):
        inside.extend(group)
        rho_last_inside = group[-1]
        rho_first_outside = groups[index + 1][0]
        if rho_last_inside > policy.candidate_rho_max:
            continue
        relative_gap = (rho_first_outside - rho_last_inside) / rho_last_inside
        median = statistics.median(inside)
        internal_spread = (max(inside) - min(inside)) / median
        candidates.append(ShellAlternative(
            boundary_group=index,
            geometric_CN=len(inside),
            relative_gap=relative_gap,
            internal_spread=internal_spread,
            strong_contacts_outside=relative_gap < policy.minimum_shell_gap,
        ))
    return tuple(candidates)


def _secondary_evidence() -> tuple[SecondaryEvidence, ...]:
    return (
        SecondaryEvidence("bvs", EvidenceStatus.NOT_AVAILABLE, "BVS analyzer not supplied"),
        SecondaryEvidence(
            "coordination_geometry",
            EvidenceStatus.NOT_APPLICABLE,
            "geometry validation follows shell resolution",
        ),
    )


def _resolve_rho_values(
    values: tuple[float, ...],
    policy: ShellResolutionPolicy,
) -> BoundaryDecision:
    if not values or any(not math.isfinite(value) or value <= 0 for value in values):
        raise ValueError("normalized distances must be positive and finite")
    groups = _group_distances(values, policy.distance_group_tolerance)
    evidence = _secondary_evidence()
    if len(groups) < 2:
        diagnostic = Diagnostic(
            Severity.WARNING,
            "crystal_chemistry.shell.candidates_insufficient",
            "At least two observed distance groups are required",
        )
        return BoundaryDecision(ResolutionStatus.INCOMPLETE, None, (), evidence, (diagnostic,))

    candidates = _candidate_boundaries(groups, policy)
    significant = tuple(
        item for item in candidates
        if item.relative_gap >= policy.minimum_shell_gap
    )
    if not significant:
        diagnostic = Diagnostic(
            Severity.WARNING,
            "crystal_chemistry.shell.search_boundary_not_observed",
            "No significant shell boundary was followed by an observed outer group",
        )
        return BoundaryDecision(
            ResolutionStatus.INCOMPLETE, None, candidates, evidence, (diagnostic,)
        )

    best_gap = max(item.relative_gap for item in significant)
    survivors = tuple(
        item for item in significant
        if best_gap - item.relative_gap <= policy.ambiguity_tolerance
    )
    best_spread = min(item.internal_spread for item in survivors)
    survivors = tuple(
        item for item in survivors
        if item.internal_spread - best_spread <= policy.ambiguity_tolerance
    )
    without_strong_outside = tuple(
        item for item in survivors if not item.strong_contacts_outside
    )
    if without_strong_outside:
        survivors = without_strong_outside
    if len(survivors) == 1:
        return BoundaryDecision(
            ResolutionStatus.RESOLVED, survivors[0], significant, evidence, ()
        )
    diagnostic = Diagnostic(
        Severity.WARNING,
        "crystal_chemistry.shell.boundary_ambiguous",
        "Lexicographic criteria do not distinguish shell boundaries",
    )
    return BoundaryDecision(
        ResolutionStatus.AMBIGUOUS, None, survivors, evidence, (diagnostic,)
    )


def _matching_interactions(
    first: SiteComponent,
    second: SiteComponent,
    grammar: CompositionGrammar,
) -> tuple[CandidateInteraction, ...]:
    first_symbol, second_symbol = first.element, second.element
    if first_symbol is None or second_symbol is None:
        return ()
    matches: list[CandidateInteraction] = []
    for request in grammar.candidate_interactions:
        forward = (
            first_symbol in request.first_elements
            and second_symbol in request.second_elements
        )
        reverse = (
            second_symbol in request.first_elements
            and first_symbol in request.second_elements
        )
        if forward or reverse:
            matches.append(request)
    return tuple(matches)


def _allowed_radius_sums(
    grammar: CompositionGrammar,
    reference: ReferenceData,
) -> tuple[float, ...]:
    values: list[float] = []
    for request in grammar.candidate_interactions:
        for first in request.first_elements:
            for second in request.second_elements:
                try:
                    values.append(
                        reference.covalent_radii.find(first).value
                        + reference.covalent_radii.find(second).value
                    )
                except KeyError:
                    continue
    return tuple(values)


def derive_search_cutoff(
    grammar: CompositionGrammar,
    reference: ReferenceData,
    policy: ShellResolutionPolicy,
) -> float:
    sums = _allowed_radius_sums(grammar, reference)
    if not sums:
        raise ValueError("grammar has no component pairs with known covalent radii")
    return max(sums) * policy.search_rho_max


def _interpret_contact(
    contact: GeometricContact,
    first_components: tuple[SiteComponent, ...],
    second_components: tuple[SiteComponent, ...],
    grammar: CompositionGrammar,
    reference: ReferenceData,
    policy: ShellResolutionPolicy,
) -> InterpretationOutcome:
    records: list[ComponentPairInterpretation] = []
    diagnostics: list[Diagnostic] = []
    incomplete: set[InteractionScope] = set()
    reported_missing: set[tuple[str, GrammarOperation]] = set()
    for first in first_components:
        for second in second_components:
            requests = _matching_interactions(first, second, grammar)
            if not requests:
                continue
            radii: list[float] = []
            missing_symbols: list[str] = []
            for component in (first, second):
                symbol = component.element
                try:
                    radii.append(reference.covalent_radii.find(symbol).value)
                except (KeyError, TypeError):
                    missing_symbols.append(symbol or "unknown")
            if missing_symbols:
                for request in requests:
                    incomplete.add(InteractionScope.from_request(request))
                    for missing_symbol in missing_symbols:
                        key = (missing_symbol, request.operation)
                        if key not in reported_missing:
                            diagnostics.append(Diagnostic(
                                Severity.WARNING,
                                "crystal_chemistry.contact.radius_missing",
                                f"No covalent radius for {missing_symbol} in a "
                                f"{request.operation.value} component pair",
                            ))
                            reported_missing.add(key)
                continue
            radius_sum = radii[0] + radii[1]
            rho = contact.distance / radius_sum
            if rho > policy.search_rho_max:
                continue
            for request in requests:
                records.append(ComponentPairInterpretation(
                    first_species=first.species,
                    second_species=second.species,
                    first_occupancy=float(first.occupancy.value),
                    second_occupancy=float(second.occupancy.value),
                    radius_sum=radius_sum,
                    normalized_distance=rho,
                    occupancy_weight=(
                        float(first.occupancy.value) * float(second.occupancy.value)
                    ),
                    interaction_type=request.operation,
                    grammar_priority=request.priority,
                    centre_elements=request.centre_elements,
                    ligand_elements=request.ligand_elements,
                ))
    return InterpretationOutcome(
        tuple(records),
        tuple(diagnostics),
        tuple(sorted(
            incomplete,
            key=lambda item: (
                item.operation.value, item.priority.value,
                item.centre_elements, item.ligand_elements,
            ),
        )),
    )


_SHELL_OPERATIONS = frozenset({
    GrammarOperation.CENTRE_LIGAND_SHELL,
    GrammarOperation.INTERSTITIAL_COORDINATION,
    GrammarOperation.MIXED_ANION_COORDINATION,
})


def _component_symbols(atom: object) -> frozenset[str]:
    return frozenset(
        component.element
        for component in getattr(atom, "components")
        if component.element is not None
    )


def _make_resolved_contact(
    contact: GeometricContact,
    interpretations: tuple[ComponentPairInterpretation, ...],
    classification: ContactClassification,
    neighbor_occupancy: float,
) -> ResolvedContact:
    operation = interpretations[0].interaction_type
    priority = interpretations[0].grammar_priority
    return ResolvedContact(
        geometric_contact=contact,
        interaction_type=operation,
        grammar_priority=priority,
        contact_classification=classification,
        component_interpretations=interpretations,
        normalized_distance_min=min(item.normalized_distance for item in interpretations),
        normalized_distance_max=max(item.normalized_distance for item in interpretations),
        neighbor_total_occupancy=neighbor_occupancy,
        evidence=_secondary_evidence(),
        provenance=(("method", "cristma.contact_interpretation:1"),),
    )


@dataclass(frozen=True, slots=True)
class CoordinationShellResolver:
    policy: ShellResolutionPolicy
    reference: ReferenceData = field(default_factory=ReferenceData.default)

    def resolve(
        self,
        structure: CrystalStructure,
        grammar: CompositionGrammar,
    ) -> CrystalChemistryResolution:
        view = structure.atomic_view()
        try:
            cutoff = derive_search_cutoff(grammar, self.reference, self.policy)
        except ValueError:
            diagnostic = Diagnostic(
                Severity.WARNING,
                "crystal_chemistry.contact.radius_missing",
                "No complete covalent-radius pair is available for geometric search",
            )
            shells = tuple(
                CoordinationShell(
                    getattr(atom, "source_site_id", atom.id),
                    atom.id,
                    (), 0, 0.0, ResolutionStatus.INCOMPLETE,
                    diagnostics=(diagnostic,),
                )
                for request in grammar.candidate_interactions
                if request.operation in _SHELL_OPERATIONS
                for atom in view.atoms
                if _component_symbols(atom) & set(request.centre_elements)
            )
            provenance = (
                ("policy", MappingProxyType(self.policy.get_config())),
                ("search_cutoff_angstrom", None),
                ("grammar_method", f"{grammar.method_id}:{grammar.method_version}"),
                ("reference_version", grammar.reference_version),
                ("resolver_method", "cristma.coordination_shell_resolver:1"),
                ("structure_id", structure.id),
            )
            return CrystalChemistryResolution((), shells, (diagnostic,), provenance)
        graph = NeighborFinder(cutoff=cutoff).find(view)
        geometric = geometric_contacts(view, graph)
        atoms = {atom.id: atom for atom in view.atoms}
        interpreted: dict[str, InterpretationOutcome] = {}
        diagnostics: list[Diagnostic] = list(graph.diagnostics)
        for contact in geometric:
            outcome = _interpret_contact(
                contact,
                atoms[contact.first_atom_id].components,
                atoms[contact.second_atom_id].components,
                grammar,
                self.reference,
                self.policy,
            )
            interpreted[contact.contact_id] = outcome
            diagnostics.extend(outcome.diagnostics)

        contacts: list[ResolvedContact] = []
        shells: list[CoordinationShell] = []
        for request in grammar.candidate_interactions:
            if request.operation not in _SHELL_OPERATIONS:
                contacts.extend(self._network_contacts(geometric, interpreted, atoms, request))
                continue
            request_shells, request_contacts, request_diagnostics = self._shells_for_request(
                view, geometric, interpreted, atoms, request
            )
            shells.extend(request_shells)
            contacts.extend(request_contacts)
            diagnostics.extend(request_diagnostics)

        maximum_rho = max(
            (
                item.normalized_distance
                for outcome in interpreted.values()
                for item in outcome.interpretations
            ),
            default=0.0,
        )
        provenance = (
            ("policy", MappingProxyType(self.policy.get_config())),
            ("search_cutoff_angstrom", cutoff),
            ("maximum_observed_rho", maximum_rho),
            ("grammar_method", f"{grammar.method_id}:{grammar.method_version}"),
            ("reference_version", grammar.reference_version),
            ("resolver_method", "cristma.coordination_shell_resolver:1"),
            ("structure_id", structure.id),
        )
        return CrystalChemistryResolution(
            tuple(contacts), tuple(shells), tuple(diagnostics), provenance
        )

    def _matching_records(
        self,
        outcome: InterpretationOutcome,
        request: CandidateInteraction,
    ) -> tuple[ComponentPairInterpretation, ...]:
        return tuple(
            item for item in outcome.interpretations
            if item.interaction_type is request.operation
            and item.grammar_priority is request.priority
            and item.centre_elements == request.centre_elements
            and item.ligand_elements == request.ligand_elements
        )

    def _network_contacts(self, geometric, interpreted, atoms, request):
        results: list[ResolvedContact] = []
        for contact in geometric:
            records = self._matching_records(interpreted[contact.contact_id], request)
            records = tuple(
                item for item in records
                if item.normalized_distance <= self.policy.candidate_rho_max
            )
            if records:
                occupancy = sum(
                    float(item.occupancy.value)
                    for item in atoms[contact.second_atom_id].components
                )
                results.append(_make_resolved_contact(
                    contact, records, ContactClassification.PRIMARY, occupancy
                ))
        return results

    def _shells_for_request(self, view, geometric, interpreted, atoms, request):
        orbit_rows: dict[str, list[tuple[object, list[tuple[float, GeometricContact, tuple[ComponentPairInterpretation, ...], float]]]]] = {}
        for center in view.atoms:
            if not (_component_symbols(center) & set(request.centre_elements)):
                continue
            rows = []
            for contact in geometric:
                if center.id not in {contact.first_atom_id, contact.second_atom_id}:
                    continue
                other_id = (
                    contact.second_atom_id
                    if center.id == contact.first_atom_id
                    else contact.first_atom_id
                )
                other = atoms[other_id]
                if not (_component_symbols(other) & set(request.ligand_elements)):
                    continue
                records = self._matching_records(interpreted[contact.contact_id], request)
                if not records:
                    continue
                rho = statistics.median(item.normalized_distance for item in records)
                occupancy = sum(float(item.occupancy.value) for item in other.components)
                rows.append((rho, contact, records, occupancy))
            source_site_id = getattr(center, "source_site_id", center.id)
            orbit_rows.setdefault(source_site_id, []).append((center, rows))

        shells: list[CoordinationShell] = []
        resolved_contacts: list[ResolvedContact] = []
        diagnostics: list[Diagnostic] = []
        for source_site_id, center_rows in orbit_rows.items():
            request_scope = InteractionScope.from_request(request)
            center_ids = {center.id for center, _ in center_rows}
            missing_outcomes = tuple(
                interpreted[contact.contact_id]
                for contact in geometric
                if center_ids & {contact.first_atom_id, contact.second_atom_id}
                and request_scope in interpreted[contact.contact_id].incomplete_scopes
            )
            if missing_outcomes:
                shell_diagnostics = tuple(dict.fromkeys(
                    diagnostic
                    for outcome in missing_outcomes
                    for diagnostic in outcome.diagnostics
                ))
                for center, _ in center_rows:
                    shells.append(CoordinationShell(
                        source_site_id, center.id, (), 0, 0.0,
                        ResolutionStatus.INCOMPLETE,
                        diagnostics=shell_diagnostics,
                    ))
                continue
            signatures = {
                tuple((round(row[0], 10), row[1].first_source_site_id, row[1].second_source_site_id)
                      for row in sorted(rows, key=lambda item: item[0]))
                for _, rows in center_rows
            }
            if len(signatures) != 1:
                diagnostic = Diagnostic(
                    Severity.WARNING,
                    "crystal_chemistry.shell.symmetry_inconsistent",
                    f"Equivalent centres for {source_site_id} have different contact signatures",
                )
                diagnostics.append(diagnostic)
                for center, _ in center_rows:
                    shells.append(CoordinationShell(
                        source_site_id, center.id, (), 0, 0.0,
                        ResolutionStatus.INCOMPLETE, diagnostics=(diagnostic,),
                    ))
                continue
            template_rows = sorted(center_rows[0][1], key=lambda item: item[0])
            decision = _resolve_rho_values(tuple(row[0] for row in template_rows), self.policy)
            lower_decision = _resolve_rho_values(
                tuple(min(item.normalized_distance for item in row[2]) for row in template_rows),
                self.policy,
            )
            upper_decision = _resolve_rho_values(
                tuple(max(item.normalized_distance for item in row[2]) for row in template_rows),
                self.policy,
            )
            lower_cn = lower_decision.selected.geometric_CN if lower_decision.selected else None
            upper_cn = upper_decision.selected.geometric_CN if upper_decision.selected else None
            if (lower_decision.status, lower_cn) != (upper_decision.status, upper_cn):
                diagnostic = Diagnostic(
                    Severity.WARNING,
                    "crystal_chemistry.shell.mixed_occupancy_disagreement",
                    "Valid component interpretations imply different shell boundaries",
                )
                alternatives = tuple(dict.fromkeys(
                    lower_decision.alternatives + upper_decision.alternatives
                ))
                decision = BoundaryDecision(
                    ResolutionStatus.AMBIGUOUS,
                    None,
                    alternatives,
                    _secondary_evidence(),
                    (diagnostic,),
                )
            for center, rows in center_rows:
                ordered = sorted(rows, key=lambda item: item[0])
                selected_count = decision.selected.geometric_CN if decision.selected else 0
                center_contacts: list[ResolvedContact] = []
                for index, (_, contact, records, occupancy) in enumerate(ordered):
                    classification = (
                        ContactClassification.PRIMARY
                        if index < selected_count else ContactClassification.SECONDARY
                    )
                    resolved = _make_resolved_contact(contact, records, classification, occupancy)
                    if any(
                        item.normalized_distance <= self.policy.candidate_rho_max
                        for item in records
                    ):
                        resolved_contacts.append(resolved)
                    if index < selected_count:
                        center_contacts.append(resolved)
                shells.append(CoordinationShell(
                    source_site_id=source_site_id,
                    center_atom_id=center.id,
                    contacts=tuple(center_contacts),
                    geometric_CN=len(center_contacts),
                    mean_occupied_neighbors=math.fsum(
                        item.neighbor_total_occupancy for item in center_contacts
                    ),
                    status=decision.status,
                    alternatives=decision.alternatives,
                    evidence=decision.evidence,
                    diagnostics=decision.diagnostics,
                    provenance=(("interaction", request.operation.value),),
                ))
                diagnostics.extend(decision.diagnostics)
        return shells, resolved_contacts, diagnostics


__all__ = [
    "BoundaryDecision",
    "CoordinationShellResolver",
    "InterpretationOutcome",
    "derive_search_cutoff",
]
