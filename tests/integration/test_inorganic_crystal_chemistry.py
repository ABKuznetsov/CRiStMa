from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pytest

import cristma
from cristma.chemistry import (
    ChemistryAnalyzer,
    Composition,
    GrammarOperation,
    InteractionLayer,
)
from cristma.crystal_chemistry import (
    ContactClassification,
    CoordinationShellResolver,
    CrystalChemistryResolution,
    PolyhedronBuilder,
    ResolutionStatus,
    ShellResolutionPolicy,
    StructuralGraphBuilder,
    StructuralUnitBuilder,
    StructuralUnitKind,
)
from cristma.structure import CrystalStructure


FIXTURES = Path(__file__).parents[1] / "fixtures" / "crystal_chemistry"
ACCEPTANCE_POLICY = ShellResolutionPolicy(1.60, 0.01, 0.08, 0.01, 2.0)


@dataclass(frozen=True)
class Calculation:
    structure: CrystalStructure
    result: CrystalChemistryResolution


def calculate(filename: str) -> Calculation:
    structure = cristma.read(FIXTURES / filename).structures[0]
    chemistry = ChemistryAnalyzer().analyze(Composition.from_structure(structure))
    result = CoordinationShellResolver(ACCEPTANCE_POLICY).resolve(
        structure, chemistry.grammar
    )
    return Calculation(structure, result)


def shell_elements(calculation: Calculation):
    atoms = {atom.id: atom for atom in calculation.structure.atomic_view().atoms}
    return tuple(
        (atoms[shell.center_atom_id].components[0].element, shell)
        for shell in calculation.result.coordination_shells
    )


def build_unit_graph(filename: str):
    calculation = calculate(filename)
    view = calculation.structure.atomic_view()
    polyhedra = tuple(
        built.polyhedron
        for shell in calculation.result.coordination_shells
        if shell.status is ResolutionStatus.RESOLVED
        for built in (PolyhedronBuilder().build(shell, view),)
        if built.polyhedron is not None
    )
    units = StructuralUnitBuilder().build(calculation.result, polyhedra).units
    return StructuralGraphBuilder().build(units, calculation.result.contacts)


@pytest.mark.parametrize(
    ("filename", "expected"),
    (
        (
            "CaMoO4_9009632.cif",
            {
                ("Ca-O", GrammarOperation.INTERSTITIAL_COORDINATION, InteractionLayer.INTERSTITIAL),
                ("Mo-O", GrammarOperation.CENTRE_LIGAND_SHELL, InteractionLayer.STRUCTURAL),
            },
        ),
        (
            "LiB3O5_3000122.cif",
            {
                ("B-O", GrammarOperation.CENTRE_LIGAND_SHELL, InteractionLayer.STRUCTURAL),
                ("Li-O", GrammarOperation.INTERSTITIAL_COORDINATION, InteractionLayer.INTERSTITIAL),
            },
        ),
        (
            "FeS2_9000594.cif",
            {
                ("Fe-S", GrammarOperation.CENTRE_LIGAND_SHELL, InteractionLayer.COORDINATION),
                ("S-S", GrammarOperation.INTRA_SUBSYSTEM_BONDS, InteractionLayer.INTRA_SUBSYSTEM),
            },
        ),
    ),
)
def test_resolved_contacts_preserve_reference_interaction_roles(
    filename: str,
    expected: set[tuple[str, GrammarOperation, InteractionLayer]],
) -> None:
    calculation = calculate(filename)
    atoms = {atom.id: atom for atom in calculation.structure.atomic_view().atoms}

    actual = {
        (
            "-".join(sorted((
                atoms[contact.geometric_contact.first_atom_id].components[0].element,
                atoms[contact.geometric_contact.second_atom_id].components[0].element,
            ))),
            contact.interaction_type,
            contact.interaction_layer,
        )
        for contact in calculation.result.contacts
    }

    assert actual == expected


@pytest.mark.parametrize(
    ("filename", "expected_layers"),
    (
        (
            "CaMoO4_9009632.cif",
            {InteractionLayer.STRUCTURAL, InteractionLayer.INTERSTITIAL},
        ),
        (
            "LiB3O5_3000122.cif",
            {InteractionLayer.STRUCTURAL, InteractionLayer.INTERSTITIAL},
        ),
        (
            "FeS2_9000594.cif",
            {InteractionLayer.COORDINATION, InteractionLayer.INTRA_SUBSYSTEM},
        ),
    ),
)
def test_structural_graph_preserves_contact_semantics(
    filename: str,
    expected_layers: set[InteractionLayer],
) -> None:
    graph = build_unit_graph(filename)

    assert any(unit.kind is StructuralUnitKind.POLYHEDRON for unit in graph.units)
    actual_layers = {
        layer
        for connection in graph.connections
        for layer in connection.interaction_layers
    }
    assert actual_layers >= expected_layers


def test_alpha_si3n4_rounded_special_position_is_not_duplicated() -> None:
    calculation = calculate("Si3N4_9013139.cif")
    composition = Composition.from_structure(calculation.structure)

    assert composition.as_dict() == {"N": 16.0, "Si": 12.0}
    shells = [shell for element, shell in shell_elements(calculation) if element == "Si"]
    assert shells
    assert {shell.status for shell in shells} == {ResolutionStatus.RESOLVED}
    assert {shell.geometric_CN for shell in shells} == {4}


@pytest.mark.parametrize(
    ("filename", "operation"),
    (
        ("SiC_9008856.cif", GrammarOperation.COVALENT_NETWORK),
        ("FeS2_9000594.cif", GrammarOperation.INTRA_SUBSYSTEM_BONDS),
        ("NiAl_B2_analytic.cif", GrammarOperation.METALLIC_COORDINATION),
    ),
)
def test_network_materials_return_contacts_without_forced_polyhedra(
    filename: str,
    operation: GrammarOperation,
) -> None:
    calculation = calculate(filename)

    assert any(
        contact.interaction_type is operation
        for contact in calculation.result.contacts
    )
    if filename in {"SiC_9008856.cif", "NiAl_B2_analytic.cif"}:
        assert calculation.result.coordination_shells == ()


@pytest.mark.parametrize(
    ("filename", "center", "coordination"),
    (
        ("NaF_9007457.cif", "Na", 6),
        ("Si3N4_9013139.cif", "Si", 4),
        ("FeS2_9000594.cif", "Fe", 6),
        ("CaMoO4_9009632.cif", "Mo", 4),
    ),
)
def test_established_shells_are_discovered_not_supplied(
    filename: str,
    center: str,
    coordination: int,
) -> None:
    calculation = calculate(filename)
    shells = [
        shell for element, shell in shell_elements(calculation)
        if element == center and shell.status is ResolutionStatus.RESOLVED
    ]

    assert shells
    assert {shell.geometric_CN for shell in shells} == {coordination}
    assert all(
        PolyhedronBuilder().build(shell, calculation.structure.atomic_view()).polyhedron
        is not None
        for shell in shells
    )


def test_calcium_diazenide_keeps_ca_n_shell_and_n_n_contact() -> None:
    calculation = calculate("CaN2_analytic.cif")
    calcium_shells = [
        shell for element, shell in shell_elements(calculation)
        if element == "Ca" and shell.status is ResolutionStatus.RESOLVED
    ]

    assert calcium_shells
    assert any(
        contact.interaction_type is GrammarOperation.INTRA_SUBSYSTEM_BONDS
        and contact.geometric_contact.distance == pytest.approx(1.202, abs=1e-4)
        for contact in calculation.result.contacts
    )


def test_lithium_triborate_resolves_both_bo3_and_bo4_units() -> None:
    calculation = calculate("LiB3O5_3000122.cif")
    boron_shells = [
        shell for element, shell in shell_elements(calculation)
        if element == "B" and shell.status is ResolutionStatus.RESOLVED
    ]

    assert {shell.geometric_CN for shell in boron_shells} == {3, 4}


def test_bismuth_telluride_retains_contacts_beyond_primary_shell() -> None:
    calculation = calculate("Bi2Te3_9011962.cif")

    assert any(
        contact.contact_classification is ContactClassification.SECONDARY
        for contact in calculation.result.contacts
    )


def test_sodium_phosphide_resolves_each_distinct_na_environment() -> None:
    calculation = calculate("Na3P_1010294.cif")
    sodium_shells = [
        shell for element, shell in shell_elements(calculation)
        if element == "Na" and shell.status is ResolutionStatus.RESOLVED
    ]

    assert sodium_shells
    assert {shell.geometric_CN for shell in sodium_shells} == {3, 9}


def test_anorthite_never_forces_an_unresolved_tetrahedron() -> None:
    calculation = calculate("anorthite_9000361.cif")
    view = calculation.structure.atomic_view()
    framework_shells = [
        shell for element, shell in shell_elements(calculation)
        if element in {"Al", "Si"}
    ]

    assert framework_shells
    for shell in framework_shells:
        built = PolyhedronBuilder().build(shell, view)
        if shell.status is ResolutionStatus.RESOLVED:
            assert shell.geometric_CN == 4
            assert built.polyhedron is not None
        else:
            assert shell.status in {ResolutionStatus.AMBIGUOUS, ResolutionStatus.INCOMPLETE}
            assert built.polyhedron is None
