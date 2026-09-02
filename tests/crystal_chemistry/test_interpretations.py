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


POLICY = ShellResolutionPolicy(1.6, 0.01, 0.08, 0.01, 2.0)
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
        layer=InteractionLayer.COORDINATION,
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


def test_component_interpretation_preserves_interaction_layer() -> None:
    request = CandidateInteraction(
        first_elements=("S",),
        second_elements=("S",),
        operation=GrammarOperation.INTRA_SUBSYSTEM_BONDS,
        layer=InteractionLayer.INTRA_SUBSYSTEM,
        priority=InteractionPriority.ALLOWED,
        centre_elements=("S",),
        ligand_elements=("S",),
        evidence=EVIDENCE,
    )

    outcome = _interpret_contact(
        contact(2.0),
        (component("S", 1.0),),
        (component("S", 1.0),),
        grammar(request),
        ReferenceData.default(),
        POLICY,
    )

    assert outcome.interpretations[0].interaction_layer is InteractionLayer.INTRA_SUBSYSTEM


def test_same_operation_retains_distinct_centre_views() -> None:
    forward = interaction(
        ("Ca",), ("O",), GrammarOperation.CENTRE_LIGAND_SHELL
    )
    reverse = CandidateInteraction(
        first_elements=("Ca",), second_elements=("O",),
        operation=GrammarOperation.CENTRE_LIGAND_SHELL,
        layer=InteractionLayer.COORDINATION,
        priority=InteractionPriority.ALLOWED,
        centre_elements=("O",), ligand_elements=("Ca",), evidence=EVIDENCE,
    )

    outcome = _interpret_contact(
        contact(), (component("Ca", 1.0),), (component("O", 1.0),),
        grammar(forward, reverse), ReferenceData.default(), POLICY,
    )

    assert {(item.centre_elements, item.ligand_elements) for item in outcome.interpretations} == {
        (("Ca",), ("O",)),
        (("O",), ("Ca",)),
    }


def test_missing_primary_radius_is_reported_not_silently_skipped() -> None:
    requests = grammar(
        interaction(("Cf",), ("O",), GrammarOperation.CENTRE_LIGAND_SHELL)
    )

    outcome = _interpret_contact(
        contact(), (component("Cf", 1.0),), (component("O", 1.0),),
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
        interaction(("O",), ("Cf",), GrammarOperation.CENTRE_LIGAND_SHELL)
    )

    outcome = _interpret_contact(
        contact(), (component("O", 1.0),), (component("Cf", 1.0),),
        requests, ReferenceData.default(), POLICY,
    )

    assert "Cf" in outcome.diagnostics[0].message


def test_search_cutoff_uses_largest_available_allowed_radius_sum() -> None:
    requests = grammar(
        interaction(("Ca",), ("O",), GrammarOperation.CENTRE_LIGAND_SHELL),
        interaction(("Xe",), ("O",), GrammarOperation.INTERSTITIAL_COORDINATION),
    )

    cutoff = derive_search_cutoff(requests, ReferenceData.default(), POLICY)

    assert cutoff == pytest.approx((1.76 + 0.66) * 2.0)
