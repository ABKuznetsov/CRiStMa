from __future__ import annotations

from fractions import Fraction

import cristma.crystal_chemistry as crystal_chemistry
from cristma.chemistry import (
    CandidateInteraction, ChemicalEvidence, CompositionGrammar, DecompositionMode,
    GrammarOperation, InteractionLayer, InteractionPriority,
)
from cristma.core import MeasuredValue, UnitCell
from cristma.crystal_chemistry import (
    ContactAnalyzer,
    PeriodicConnectivityAnalyzer,
    PolyhedronOrbitBuilder,
    ReferenceCell,
    RingFinder,
    ShellResolutionPolicy,
    ShellRole,
    StructuralBlockFinder,
    StructuralGraphBuilder,
    StructuralRepresentationBuilder,
    StructuralSelectionPolicy,
)
from cristma.crystallography import SymmetryContext
from cristma.structure import CrystalStructure, IndependentSite, SiteComponent
from cristma.symmetry import AffineOperation


def _value(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def _result():
    cell = UnitCell.cubic(_value(10.0))
    sites = (
        IndependentSite("Ca", "Ca", (SiteComponent("Ca", _value(1.0)),),
                        (_value(0.0), _value(0.0), _value(0.0))),
        IndependentSite("O", "O", (SiteComponent("O", _value(1.0)),),
                        (_value(0.2), _value(0.0), _value(0.0))),
    )
    structure = CrystalStructure("pair", cell, sites)
    identity = AffineOperation(
        ((Fraction(1), Fraction(0), Fraction(0)),
         (Fraction(0), Fraction(1), Fraction(0)),
         (Fraction(0), Fraction(0), Fraction(1))),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    context = SymmetryContext.from_operations((identity,), cell)
    request = CandidateInteraction(
        ("Ca",), ("O",), GrammarOperation.CENTRE_LIGAND_SHELL,
        InteractionLayer.STRUCTURAL, InteractionPriority.PRIMARY,
        ("Ca",), ("O",), (ChemicalEvidence("fixture", "pair"),),
    )
    grammar = CompositionGrammar(
        DecompositionMode.STRUCTURAL_ANION_SUBSYSTEM, (request,), 1.0,
        (ChemicalEvidence("fixture", "grammar"),), (), "fixture",
    )
    return ContactAnalyzer(ShellResolutionPolicy(1.6, 0.01, 0.08, 0.01, 2.0)).analyze(
        structure, context, grammar
    )


def test_public_contact_route_materializes_reference_cell_from_orbits() -> None:
    result = _result()

    assert result.contact_orbits
    assert result.contacts == result.materialize_contacts(ReferenceCell())
    assert len(result.contacts) == sum(
        geometry.multiplicity_in_reference_cell
        for geometry in result.pair_table.contact_orbits
        if geometry.geometry_orbit_id in {item.geometry_orbit_id for item in result.contact_orbits}
    )
    assert all(item.resolved_contact_orbit_id for item in result.contacts)


def test_legacy_expanded_first_symbols_are_not_public() -> None:
    assert not hasattr(crystal_chemistry, "CoordinationShellResolver")
    assert not hasattr(crystal_chemistry, "CrystalChemistryResolution")
    assert not hasattr(crystal_chemistry, "ContactClassification")
    assert not hasattr(crystal_chemistry, "LegacyResolvedContactOrbit")
    assert not hasattr(crystal_chemistry, "StructuralUnit")
    assert not hasattr(crystal_chemistry, "StructuralConnection")
    assert not hasattr(crystal_chemistry, "PeriodicUnitRef")
    assert not hasattr(crystal_chemistry, "StructuralRing")
    assert not hasattr(crystal_chemistry, "StructuralBlockOrbit")


def test_materialization_is_a_dead_end_for_all_scientific_results() -> None:
    result = _result()
    policy = StructuralSelectionPolicy(
        frozenset((InteractionLayer.STRUCTURAL,)),
        frozenset((ShellRole.PRIMARY, ShellRole.SECONDARY)),
    )

    def calculate():
        polyhedra = PolyhedronOrbitBuilder().build(result)
        graph = StructuralGraphBuilder().build(result, polyhedra)
        representation = StructuralRepresentationBuilder(policy).build(graph)
        connectivity = PeriodicConnectivityAnalyzer().analyze(representation)
        blocks = StructuralBlockFinder().find(representation, connectivity)
        rings = RingFinder().find(representation, blocks)
        return (
            result.coordination_shell_orbits,
            polyhedra,
            graph,
            connectivity,
            rings,
        )

    before = calculate()
    _ = result.contacts
    after = calculate()

    assert after == before
