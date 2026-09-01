"""Inorganic contact interpretation and coordination-shell resolution."""

from __future__ import annotations

from dataclasses import dataclass

from cristma.chemistry import CandidateInteraction, CompositionGrammar
from cristma.chemistry.grammar import GrammarOperation
from cristma.crystallography import GeometricContact
from cristma.diagnostics import Diagnostic, Severity
from cristma.reference_data import ReferenceData
from cristma.structure import SiteComponent

from .contacts import ComponentPairInterpretation
from .policy import ShellResolutionPolicy


@dataclass(frozen=True, slots=True)
class InterpretationOutcome:
    interpretations: tuple[ComponentPairInterpretation, ...]
    diagnostics: tuple[Diagnostic, ...]
    incomplete_interactions: tuple[GrammarOperation, ...] = ()


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
