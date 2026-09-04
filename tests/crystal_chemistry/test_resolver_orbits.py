from __future__ import annotations

import math

import pytest

from cristma.chemistry import (
    CandidateInteraction,
    ChemistryAnalyzer,
    Composition,
    CompositionGrammar,
    DecompositionMode,
    GrammarOperation,
    InteractionLayer,
    InteractionPriority,
)
from cristma.chemistry.evidence import ChemicalEvidence
from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
from cristma.crystallography import GeometricContact, geometric_contacts
from cristma.crystal_chemistry import (
    ComponentPairInterpretation,
    ContactClassification,
    CoordinationShellResolver,
    ResolutionStatus,
    ShellResolutionPolicy,
)
from cristma.crystal_chemistry.resolver import InterpretationOutcome, _interpret_contact
from cristma.geometry import NeighborFinder
from cristma.reference_data import ReferenceData
from cristma.structure import CrystalStructure, IndependentSite, SiteComponent
from cristma.symmetry import SpaceGroupDefinition, parse_xyz_operation


POLICY = ShellResolutionPolicy(1.6, 0.01, 0.08, 0.01, 2.0)


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def site(
    label: str,
    symbol: str,
    fractional: tuple[float, float, float],
    occupancy: float = 1.0,
) -> IndependentSite:
    return IndependentSite(
        id=f"site:{label}",
        label=label,
        components=(SiteComponent(symbol, number(occupancy)),),
        fractional=tuple(number(value) for value in fractional),
        calculated_multiplicity=2,
    )


def equivalent_tetrahedra(*, ligand_occupancy: float = 1.0) -> CrystalStructure:
    center = (0.25, 0.25, 0.25)
    delta = 2.0 / (20.0 * math.sqrt(3.0))
    ligand_rows = tuple(
        tuple(center[axis] + delta * signs[axis] for axis in range(3))
        for signs in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    )
    outer = (center[0] + 3.0 / 20.0, center[1], center[2])
    sites = (site("Ca", "Ca", center),) + tuple(
        site(f"O{index}", "O", row, ligand_occupancy)
        for index, row in enumerate((*ligand_rows, outer), start=1)
    )
    symmetry = SpaceGroupDefinition(
        operations=(
            parse_xyz_operation("x,y,z", operation_id="op:1"),
            parse_xyz_operation("-x,-y,-z", operation_id="op:2"),
        ),
        provenance="derived",
    )
    return CrystalStructure(
        "equivalent tetrahedra",
        UnitCell.cubic(number(20.0)),
        sites,
        id="structure:equivalent-tetrahedra",
        space_group=symmetry,
    )


def mixed_boundary_structure() -> CrystalStructure:
    structure = equivalent_tetrahedra()
    center = (0.25, 0.25, 0.25)
    mixed = IndependentSite(
        id="site:mixed", label="mixed",
        components=(
            SiteComponent("O", number(0.5)),
            SiteComponent("F", number(0.5)),
        ),
        fractional=(number(center[0] + 2.3 / 20.0), number(center[1]), number(center[2])),
        calculated_multiplicity=2,
    )
    outer = site("outer", "O", (center[0] + 2.7 / 20.0, center[1], center[2]))
    return CrystalStructure(
        "mixed boundary", structure.cell, structure.sites[:-1] + (mixed, outer),
        id="structure:mixed-boundary", space_group=structure.space_group,
    )


def rocksalt_naf() -> CrystalStructure:
    symmetry = SpaceGroupDefinition(
        operations=tuple(
            parse_xyz_operation(expression, operation_id=f"op:{index}")
            for index, expression in enumerate((
                "x,y,z",
                "x,y+1/2,z+1/2",
                "x+1/2,y,z+1/2",
                "x+1/2,y+1/2,z",
            ), start=1)
        ),
        provenance="derived",
    )
    return CrystalStructure(
        "NaF",
        UnitCell.cubic(number(4.634)),
        (
            site("Na", "Na", (0.0, 0.0, 0.0)),
            site("F", "F", (0.5, 0.5, 0.5)),
        ),
        id="structure:naf",
        space_group=symmetry,
    )
def resolve(structure: CrystalStructure):
    composition = Composition.from_structure(structure)
    grammar = ChemistryAnalyzer().analyze(composition).grammar
    return CoordinationShellResolver(POLICY).resolve(structure, grammar)


def test_resolver_projects_one_orbit_decision_to_every_equivalent_center() -> None:
    structure = equivalent_tetrahedra()

    result = resolve(structure)

    shells = [
        shell for shell in result.coordination_shells
        if shell.source_site_id == "site:Ca"
    ]
    expected_centers = [
        atom for atom in structure.atomic_view().atoms
        if atom.source_site_id == "site:Ca"
    ]
    assert len(shells) == len(expected_centers) == 2
    assert {shell.status for shell in shells} == {ResolutionStatus.RESOLVED}
    assert {shell.geometric_CN for shell in shells} == {4}


def test_search_horizon_observes_outer_group_for_single_distance_shell() -> None:
    result = resolve(rocksalt_naf())

    shells = [
        shell for shell in result.coordination_shells
        if shell.source_site_id == "site:Na"
    ]

    assert shells
    assert {shell.status for shell in shells} == {ResolutionStatus.RESOLVED}
    assert {shell.geometric_CN for shell in shells} == {6}


def test_vacancy_changes_mean_occupancy_not_geometric_cn() -> None:
    result = resolve(equivalent_tetrahedra(ligand_occupancy=0.75))

    shells = [
        shell for shell in result.coordination_shells
        if shell.source_site_id == "site:Ca"
    ]
    assert {shell.geometric_CN for shell in shells} == {4}
    assert all(shell.mean_occupied_neighbors == pytest.approx(3.0) for shell in shells)


def test_result_records_reproducible_method_provenance() -> None:
    result = resolve(equivalent_tetrahedra())

    provenance = dict(result.provenance)
    assert provenance["policy"] == POLICY.get_config()
    assert provenance["search_cutoff_angstrom"] == pytest.approx(4.84)
    assert provenance["grammar_method"] == "cristma.composition_grammar:2"
    assert provenance["resolver_method"] == "cristma.coordination_shell_resolver:1"
    assert provenance["structure_id"] == "structure:equivalent-tetrahedra"


def test_missing_contact_in_one_equivalent_environment_is_incomplete() -> None:
    structure = equivalent_tetrahedra()
    view = structure.atomic_view()
    chemistry = ChemistryAnalyzer().analyze(Composition.from_structure(structure))
    request = chemistry.grammar.candidate_interactions[0]
    graph = NeighborFinder(derive_cutoff := (1.76 + 0.66) * POLICY.candidate_rho_max).find(view)
    geometric = list(geometric_contacts(view, graph))
    centers = [atom for atom in view.atoms if atom.source_site_id == "site:Ca"]
    removed = next(
        item for item in geometric
        if centers[1].id in {item.first_atom_id, item.second_atom_id}
        and item.distance < 2.1
    )
    geometric.remove(removed)
    atoms = {atom.id: atom for atom in view.atoms}
    interpreted = {
        item.contact_id: _interpret_contact(
            item,
            atoms[item.first_atom_id].components,
            atoms[item.second_atom_id].components,
            chemistry.grammar,
            ReferenceData.default(),
            POLICY,
        )
        for item in geometric
    }

    shells, _, diagnostics = CoordinationShellResolver(POLICY)._shells_for_request(
        view, tuple(geometric), interpreted, atoms, request
    )

    assert derive_cutoff > 0
    assert all(shell.status is ResolutionStatus.INCOMPLETE for shell in shells)
    assert "crystal_chemistry.shell.symmetry_inconsistent" in {
        item.code for item in diagnostics
    }


def test_mixed_components_with_different_boundaries_remain_ambiguous() -> None:
    structure = mixed_boundary_structure()
    evidence = (ChemicalEvidence("test", "mixed ligand grammar"),)
    request = CandidateInteraction(
        first_elements=("Ca",), second_elements=("F", "O"),
        operation=GrammarOperation.CENTRE_LIGAND_SHELL,
        layer=InteractionLayer.COORDINATION,
        priority=InteractionPriority.PRIMARY,
        centre_elements=("Ca",), ligand_elements=("F", "O"), evidence=evidence,
    )
    grammar = CompositionGrammar(
        DecompositionMode.CATION_ANION_SUBSYSTEM, (request,), 1.0,
        evidence, (), "test",
    )

    result = CoordinationShellResolver(POLICY).resolve(structure, grammar)

    shells = [item for item in result.coordination_shells if item.source_site_id == "site:Ca"]
    assert all(item.status is ResolutionStatus.AMBIGUOUS for item in shells)
    assert all(
        "crystal_chemistry.shell.mixed_occupancy_disagreement" in item.diagnostic_codes
        for item in shells
    )


def test_primary_scope_with_no_available_radius_is_incomplete() -> None:
    structure = CrystalStructure.explicit(
        "californium oxide",
        UnitCell.cubic(number(10.0)),
        (site("Cf", "Cf", (0.5, 0.5, 0.5)), site("O", "O", (0.7, 0.5, 0.5))),
        id="structure:missing-radius",
    )
    evidence = (ChemicalEvidence("test", "missing radius grammar"),)
    request = CandidateInteraction(
        ("Cf",), ("O",), GrammarOperation.CENTRE_LIGAND_SHELL,
        InteractionLayer.COORDINATION, InteractionPriority.PRIMARY,
        ("Cf",), ("O",), evidence,
    )
    grammar = CompositionGrammar(
        DecompositionMode.CATION_ANION_SUBSYSTEM, (request,), 1.0,
        evidence, (), "test",
    )

    result = CoordinationShellResolver(POLICY).resolve(structure, grammar)

    assert result.coordination_shells[0].status is ResolutionStatus.INCOMPLETE
    assert "crystal_chemistry.contact.radius_missing" in result.diagnostic_codes


def test_partly_missing_primary_mixed_site_radius_blocks_resolution() -> None:
    base = equivalent_tetrahedra()
    mixed_center = IndependentSite(
        id="site:Ca", label="CaCf",
        components=(
            SiteComponent("Ca", number(0.5)),
            SiteComponent("Cf", number(0.5)),
        ),
        fractional=base.sites[0].fractional,
        calculated_multiplicity=2,
    )
    structure = CrystalStructure(
        "mixed centre", base.cell, (mixed_center,) + base.sites[1:],
        id="structure:partly-missing-radius", space_group=base.space_group,
    )
    evidence = (ChemicalEvidence("test", "mixed centre grammar"),)
    request = CandidateInteraction(
        ("Ca", "Cf"), ("O",), GrammarOperation.CENTRE_LIGAND_SHELL,
        InteractionLayer.COORDINATION, InteractionPriority.PRIMARY,
        ("Ca", "Cf"), ("O",), evidence,
    )
    grammar = CompositionGrammar(
        DecompositionMode.CATION_ANION_SUBSYSTEM, (request,), 1.0,
        evidence, (), "test",
    )

    result = CoordinationShellResolver(POLICY).resolve(structure, grammar)

    shells = [item for item in result.coordination_shells if item.source_site_id == "site:Ca"]
    assert all(item.status is ResolutionStatus.INCOMPLETE for item in shells)
    assert "crystal_chemistry.contact.radius_missing" in result.diagnostic_codes


def test_network_boundaries_are_resolved_per_crystallographic_site_pair() -> None:
    structure = CrystalStructure.explicit(
        "inequivalent sulfur pairs",
        UnitCell.cubic(number(20.0)),
        (
            site("S1", "S", (0.1, 0.1, 0.1)),
            site("S2", "S", (0.6, 0.6, 0.6)),
        ),
        id="structure:inequivalent-sulfur-pairs",
    )
    atoms = {atom.id: atom for atom in structure.atomic_view().atoms}
    atoms_by_site = {atom.source_site_id: atom for atom in atoms.values()}
    evidence = (ChemicalEvidence("test", "site-pair-scoped network"),)
    request = CandidateInteraction(
        ("S",), ("S",), GrammarOperation.INTRA_SUBSYSTEM_BONDS,
        InteractionLayer.INTRA_SUBSYSTEM, InteractionPriority.PRIMARY,
        ("S",), ("S",), evidence,
    )
    rows = (
        ("s1-near", "site:S1", 1.00, (1, 0, 0)),
        ("s1-outer", "site:S1", 1.30, (2, 0, 0)),
        ("s2-near-a", "site:S2", 1.15, (0, 1, 0)),
        ("s2-near-b", "site:S2", 1.16, (0, 0, 1)),
        ("s2-outer", "site:S2", 1.50, (0, 2, 0)),
    )
    geometric = []
    interpreted = {}
    for contact_id, source_site_id, rho, translation in rows:
        atom = atoms_by_site[source_site_id]
        contact = GeometricContact(
            contact_id,
            atom.id,
            atom.id,
            translation,
            2.0 * rho,
            (2.0 * rho, 0.0, 0.0),
            source_site_id,
            source_site_id,
            "test",
        )
        record = ComponentPairInterpretation(
            atom.components[0].species,
            atom.components[0].species,
            1.0,
            1.0,
            2.0,
            rho,
            1.0,
            request.operation,
            request.layer,
            request.priority,
            request.centre_elements,
            request.ligand_elements,
        )
        geometric.append(contact)
        interpreted[contact_id] = InterpretationOutcome((record,), ())

    contacts, diagnostics = CoordinationShellResolver(POLICY)._network_contacts(
        tuple(geometric), interpreted, atoms, request
    )

    classifications = {
        item.geometric_contact.contact_id: item.contact_classification
        for item in contacts
    }
    assert diagnostics == ()
    assert classifications == {
        "s1-near": ContactClassification.PRIMARY,
        "s1-outer": ContactClassification.SECONDARY,
        "s2-near-a": ContactClassification.PRIMARY,
        "s2-near-b": ContactClassification.PRIMARY,
        "s2-outer": ContactClassification.SECONDARY,
    }
