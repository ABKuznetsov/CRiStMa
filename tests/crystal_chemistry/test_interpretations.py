from __future__ import annotations

import pytest

from cristma.chemistry import (
    CandidateInteraction,
    CompositionGrammar,
    DecompositionMode,
    GrammarOperation,
    InteractionLayer,
    InteractionPriority,
)
from cristma.chemistry.evidence import ChemicalEvidence
from cristma.core.values import MeasuredValue
from cristma.crystallography import GeometricContact
from cristma.crystal_chemistry import ShellResolutionPolicy
from cristma.crystal_chemistry.resolver import _interpret_contact, derive_search_cutoff
from cristma.reference_data import ReferenceData
from cristma.structure import SiteComponent


POLICY = ShellResolutionPolicy(1.6, 0.01, 0.08, 0.01)
EVIDENCE = (ChemicalEvidence("test.grammar", "analytic grammar"),)


def component(symbol: str, occupancy: float) -> SiteComponent:
    return SiteComponent(symbol, MeasuredValue(occupancy, None, str(occupancy)))


def interaction(
    first: tuple[str, ...],
    second: tuple[str, ...],
    operation: GrammarOperation,
) -> CandidateInteraction:
    return CandidateInteraction(
        first_elements=first,
        second_elements=second,
        operation=operation,
        layer=InteractionLayer.PRIMARY_COORDINATION,
        priority=InteractionPriority.PRIMARY,
        centre_elements=first,
        ligand_elements=second,
        evidence=EVIDENCE,
    )


def grammar(*items: CandidateInteraction) -> CompositionGrammar:
    return CompositionGrammar(
        mode=DecompositionMode.CATION_ANION_SUBSYSTEM,
        candidate_interactions=items,
        confidence=1.0,
        evidence=EVIDENCE,
        diagnostics=(),
        reference_version="test",
    )


def contact(distance: float = 2.4) -> GeometricContact:
    return GeometricContact(
        "contact:M|X|finite", "M", "X", None, distance, (distance, 0.0, 0.0),
        "site:M", "site:X", "neighbor_graph",
    )


def test_same_pair_retains_every_matching_interaction_context() -> None:
    requests = grammar(
        interaction(("Ca", "Sr"), ("F", "O"), GrammarOperation.CENTRE_LIGAND_SHELL),
        interaction(("Ca", "Sr"), ("F", "O"), GrammarOperation.INTERSTITIAL_COORDINATION),
    )

    outcome = _interpret_contact(
        contact(),
        (component("Ca", 0.7), component("Sr", 0.3)),
        (component("O", 0.8), component("F", 0.2)),
        requests,
        ReferenceData.default(),
        POLICY,
    )

    assert len(outcome.interpretations) == 8
    assert {item.interaction_type for item in outcome.interpretations} == {
        GrammarOperation.CENTRE_LIGAND_SHELL,
        GrammarOperation.INTERSTITIAL_COORDINATION,
    }


def test_missing_primary_radius_is_reported_not_silently_skipped() -> None:
    requests = grammar(
        interaction(("Xe",), ("O",), GrammarOperation.CENTRE_LIGAND_SHELL)
    )

    outcome = _interpret_contact(
        contact(), (component("Xe", 1.0),), (component("O", 1.0),),
        requests, ReferenceData.default(), POLICY,
    )

    assert outcome.interpretations == ()
    assert {item.code for item in outcome.diagnostics} == {
        "crystal_chemistry.contact.radius_missing"
    }
    assert outcome.incomplete_interactions == (
        GrammarOperation.CENTRE_LIGAND_SHELL,
    )


def test_missing_second_component_radius_names_the_missing_element() -> None:
    requests = grammar(
        interaction(("O",), ("Xe",), GrammarOperation.CENTRE_LIGAND_SHELL)
    )

    outcome = _interpret_contact(
        contact(), (component("O", 1.0),), (component("Xe", 1.0),),
        requests, ReferenceData.default(), POLICY,
    )

    assert "Xe" in outcome.diagnostics[0].message


def test_search_cutoff_uses_largest_available_allowed_radius_sum() -> None:
    requests = grammar(
        interaction(("Ca",), ("O",), GrammarOperation.CENTRE_LIGAND_SHELL),
        interaction(("Xe",), ("O",), GrammarOperation.INTERSTITIAL_COORDINATION),
    )

    cutoff = derive_search_cutoff(requests, ReferenceData.default(), POLICY)

    assert cutoff == pytest.approx((1.76 + 0.66) * 1.6)
