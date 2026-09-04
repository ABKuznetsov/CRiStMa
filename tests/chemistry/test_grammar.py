from __future__ import annotations

import pytest

from cristma.chemistry import (
    ChemistryAnalyzer,
    Composition,
    GrammarOperation,
    InteractionLayer,
    InteractionPriority,
)


def signatures(formula: dict[str, float]) -> set[tuple[tuple[str, ...], tuple[str, ...], str, str]]:
    grammar = ChemistryAnalyzer().analyze(Composition.from_mapping(formula)).grammar
    return {
        (
            item.first_elements,
            item.second_elements,
            item.operation.value,
            item.priority.value,
        )
        for item in grammar.candidate_interactions
    }


@pytest.mark.parametrize(
    ("formula", "pairs"),
    (
        ({"Fe": 1}, {(('Fe',), ('Fe',), "metallic_coordination", "primary")}),
        ({"Si": 1}, {(('Si',), ('Si',), "covalent_network", "primary")}),
        ({"Fe": 1, "Al": 1}, {(('Al',), ('Fe',), "metallic_coordination", "primary")}),
        ({"Ca": 1, "O": 1}, {(('Ca',), ('O',), "interstitial_coordination", "allowed")}),
        ({"Na": 1, "Cl": 1}, {(('Cl',), ('Na',), "centre_ligand_shell", "primary")}),
        (
            {"Fe": 1, "S": 2},
            {
                (('Fe',), ('S',), "centre_ligand_shell", "primary"),
                (('S',), ('S',), "intra_subsystem_bonds", "allowed"),
            },
        ),
    ),
)
def test_grammar_returns_concrete_search_pairs(formula, pairs) -> None:
    assert signatures(formula) == pairs


def test_organic_grammar_searches_intramolecular_covalent_contacts() -> None:
    grammar = ChemistryAnalyzer().analyze(
        Composition.from_mapping({"C": 2, "H": 6, "O": 1})
    ).grammar

    assert grammar.candidate_interactions[0].operation is GrammarOperation.COVALENT_NETWORK
    assert grammar.candidate_interactions[0].first_elements == ("C", "H", "O")
    assert grammar.candidate_interactions[0].priority is InteractionPriority.PRIMARY


def test_metal_organic_grammar_separates_organic_and_metal_donor_searches() -> None:
    grammar = ChemistryAnalyzer().analyze(
        Composition.from_mapping({"Zn": 1, "C": 2, "H": 4, "N": 2})
    ).grammar

    assert {item.operation for item in grammar.candidate_interactions} == {
        GrammarOperation.COVALENT_NETWORK,
        GrammarOperation.CENTRE_LIGAND_SHELL,
    }
    metal_donor = next(
        item
        for item in grammar.candidate_interactions
        if item.operation is GrammarOperation.CENTRE_LIGAND_SHELL
    )
    assert metal_donor.centre_elements == ("Zn",)
    assert metal_donor.ligand_elements == ("N",)


def test_directed_coordination_operations_are_distinct() -> None:
    assert GrammarOperation.INTERSTITIAL_COORDINATION.value == "interstitial_coordination"
    assert GrammarOperation.MIXED_ANION_COORDINATION.value == "mixed_anion_coordination"


def test_mixed_anion_grammar_keeps_both_anions_and_excludes_metal_centres() -> None:
    grammar = ChemistryAnalyzer().analyze(
        Composition.from_mapping({"La": 1, "O": 1, "F": 1})
    ).grammar

    assert len(grammar.candidate_interactions) == 1
    interaction = grammar.candidate_interactions[0]
    assert interaction.centre_elements == ("La",)
    assert interaction.ligand_elements == ("F", "O")
    assert interaction.operation is GrammarOperation.MIXED_ANION_COORDINATION


def test_oxide_template_separates_structural_former_from_interstitial_cation() -> None:
    grammar = ChemistryAnalyzer().analyze(
        Composition.from_mapping({"Ca": 1, "Mo": 1, "O": 4})
    ).grammar

    assert {
        (
            item.centre_elements,
            item.ligand_elements,
            item.operation,
            item.layer,
            item.priority,
        )
        for item in grammar.candidate_interactions
    } == {
        (
            ("Mo",),
            ("O",),
            GrammarOperation.CENTRE_LIGAND_SHELL,
            "structural",
            InteractionPriority.PRIMARY,
        ),
        (
            ("Ca",),
            ("O",),
            GrammarOperation.INTERSTITIAL_COORDINATION,
            "interstitial",
            InteractionPriority.ALLOWED,
        ),
    }


@pytest.mark.parametrize(
    ("formula", "expected"),
    (
        (
            {"Li": 1, "B": 3, "O": 5},
            {
                ("B", "O", "structural"),
                ("Li", "O", "interstitial"),
            },
        ),
        (
            {"Fe": 1, "S": 2},
            {
                ("Fe", "S", "coordination"),
                ("S", "S", "intra_subsystem"),
            },
        ),
    ),
)
def test_reference_templates_assign_structural_interaction_layers(
    formula: dict[str, float],
    expected: set[tuple[str, str, str]],
) -> None:
    grammar = ChemistryAnalyzer().analyze(Composition.from_mapping(formula)).grammar

    assert {
        (item.centre_elements[0], item.ligand_elements[0], item.layer.value)
        for item in grammar.candidate_interactions
    } == expected
