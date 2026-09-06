from __future__ import annotations

from dataclasses import fields
from fractions import Fraction

import pytest

from cristma.crystal_chemistry import (
    ContactAnalysisResult,
    ContactAnalyzer,
    ResolutionStatus,
    ShellResolutionPolicy,
    aggregate_contact_analysis_status,
)
from cristma.chemistry import (
    CandidateInteraction,
    ChemicalEvidence,
    CompositionGrammar,
    DecompositionMode,
    GrammarOperation,
    InteractionLayer,
    InteractionPriority,
)
from cristma.core import MeasuredValue, UnitCell
from cristma.crystallography import PairTableStatus, SymmetryContext
from cristma.structure import CrystalStructure, IndependentSite, SiteComponent
from cristma.symmetry import AffineOperation


@pytest.mark.parametrize(
    ("pair_status", "shell_statuses", "expected"),
    (
        (PairTableStatus.INCOMPLETE, (), ResolutionStatus.INCOMPLETE),
        (PairTableStatus.COMPLETE, (ResolutionStatus.INCOMPLETE,), ResolutionStatus.INCOMPLETE),
        (PairTableStatus.COMPLETE, (ResolutionStatus.AMBIGUOUS,), ResolutionStatus.AMBIGUOUS),
        (PairTableStatus.COMPLETE, (ResolutionStatus.RESOLVED,), ResolutionStatus.RESOLVED),
        (PairTableStatus.COMPLETE, (), ResolutionStatus.NOT_APPLICABLE),
    ),
)
def test_contact_analysis_status_order(pair_status, shell_statuses, expected) -> None:
    assert aggregate_contact_analysis_status(pair_status, shell_statuses) is expected


def test_contact_analysis_result_does_not_store_expanded_contacts() -> None:
    names = {item.name for item in fields(ContactAnalysisResult)}
    assert "contacts" not in names
    assert "expanded_contacts" not in names


def test_network_contacts_remain_available_when_coordination_is_not_applicable() -> None:
    value = lambda number: MeasuredValue(number, None, str(number))
    cell = UnitCell.cubic(value(10.0))
    structure = CrystalStructure(
        "network-only",
        cell,
        (
            IndependentSite(
                "C", "C", (SiteComponent("C", value(1.0)),),
                (value(0.0), value(0.0), value(0.0)),
            ),
            IndependentSite(
                "N", "N", (SiteComponent("N", value(1.0)),),
                (value(0.14), value(0.0), value(0.0)),
            ),
        ),
    )
    identity = AffineOperation(
        (
            (Fraction(1), Fraction(0), Fraction(0)),
            (Fraction(0), Fraction(1), Fraction(0)),
            (Fraction(0), Fraction(0), Fraction(1)),
        ),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    context = SymmetryContext.from_operations((identity,), cell)
    grammar = CompositionGrammar(
        DecompositionMode.COVALENT_NETWORK,
        (
            CandidateInteraction(
                ("C",), ("N",), GrammarOperation.COVALENT_NETWORK,
                InteractionLayer.STRUCTURAL, InteractionPriority.PRIMARY,
                ("C",), ("N",), (ChemicalEvidence("fixture", "network"),),
            ),
        ),
        1.0,
        (ChemicalEvidence("fixture", "grammar"),),
        (),
        "fixture",
    )

    result = ContactAnalyzer(
        ShellResolutionPolicy(1.6, 0.01, 0.08, 0.01, 2.0)
    ).analyze(structure, context, grammar)

    assert result.status is ResolutionStatus.NOT_APPLICABLE
    assert result.contact_orbits
    assert result.contact_incidence_orbits
    assert result.coordination_shell_orbits == ()


def test_empty_grammar_returns_not_applicable_without_inventing_search_cutoff() -> None:
    value = lambda number: MeasuredValue(number, None, str(number))
    cell = UnitCell.cubic(value(10.0))
    structure = CrystalStructure(
        "no-applicable-chemistry",
        cell,
        (
            IndependentSite(
                "X", "X", (SiteComponent("C", value(1.0)),),
                (value(0.0), value(0.0), value(0.0)),
            ),
        ),
    )
    identity = AffineOperation(
        (
            (Fraction(1), Fraction(0), Fraction(0)),
            (Fraction(0), Fraction(1), Fraction(0)),
            (Fraction(0), Fraction(0), Fraction(1)),
        ),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    grammar = CompositionGrammar(
        DecompositionMode.COVALENT_NETWORK,
        (), 0.0, (ChemicalEvidence("fixture", "no request"),), (), "fixture",
    )

    result = ContactAnalyzer(
        ShellResolutionPolicy(1.6, 0.01, 0.08, 0.01, 2.0)
    ).analyze(
        structure,
        SymmetryContext.from_operations((identity,), cell),
        grammar,
    )

    assert result.status is ResolutionStatus.NOT_APPLICABLE
    assert result.pair_table.status is PairTableStatus.COMPLETE
    assert result.pair_table.cutoff == 0.0
    assert result.contact_orbits == ()
