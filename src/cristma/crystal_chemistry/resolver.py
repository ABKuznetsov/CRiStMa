"""Inorganic contact interpretation and coordination-shell resolution."""

from __future__ import annotations

from dataclasses import dataclass
import math
import statistics

from cristma.chemistry import CandidateInteraction, CompositionGrammar
from cristma.chemistry.grammar import GrammarOperation
from cristma.crystallography import GeometricContact
from cristma.diagnostics import Diagnostic, Severity
from cristma.reference_data import ReferenceData
from cristma.structure import SiteComponent

from .contacts import (
    ComponentPairInterpretation,
    EvidenceStatus,
    ResolutionStatus,
    SecondaryEvidence,
    ShellAlternative,
)
from .policy import ShellResolutionPolicy


@dataclass(frozen=True, slots=True)
class InterpretationOutcome:
    interpretations: tuple[ComponentPairInterpretation, ...]
    diagnostics: tuple[Diagnostic, ...]
    incomplete_interactions: tuple[GrammarOperation, ...] = ()


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
    return max(sums) * policy.candidate_rho_max


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
    incomplete: set[GrammarOperation] = set()
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
                    incomplete.add(request.operation)
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
            if rho > policy.candidate_rho_max:
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
                ))
    return InterpretationOutcome(
        tuple(records),
        tuple(diagnostics),
        tuple(sorted(incomplete, key=lambda item: item.value)),
    )


__all__ = ["InterpretationOutcome", "derive_search_cutoff"]
