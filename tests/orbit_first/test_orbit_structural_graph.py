from __future__ import annotations

from dataclasses import fields
from fractions import Fraction

import pytest

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
from cristma.crystal_chemistry import (
    ContactAnalyzer,
    PeriodicConnectivityAnalyzer,
    PolyhedronOrbitBuilder,
    ShellResolutionPolicy,
    StructuralGraphBuilder,
    StructuralRepresentationBuilder,
    StructuralSelectionPolicy,
    StructuralUnitOrbit,
    StructuralUnitKind,
    StructuralConnectionOrbit,
    StructuralConnectionKind,
    StructuralRepresentation,
    ShellRole,
    StructuralBlockFinder,
    integer_translation_lattice_basis,
)
from cristma.crystallography import SymmetryContext
from cristma.crystallography import PeriodicSymmetryRelation
from cristma.structure import CrystalStructure, IndependentSite, PeriodicAtomRef, SiteComponent
from cristma.symmetry import AffineOperation
from cristma.crystal_chemistry.structural_units import _ordered_planar_face


def _value(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


IDENTITY = AffineOperation(
    (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    ),
    (Fraction(0), Fraction(0), Fraction(0)),
)


def _cell(a: float, b: float, c: float) -> UnitCell:
    right = _value(90.0)
    return UnitCell(_value(a), _value(b), _value(c), right, right, right)


def _site(site_id: str, fractional: tuple[float, float, float], symbol: str) -> IndependentSite:
    return IndependentSite(
        site_id,
        site_id,
        (SiteComponent(symbol, _value(1.0)),),
        tuple(_value(item) for item in fractional),
    )


def _rank_fixture(rank: int):
    if rank == 0:
        cell = _cell(10.0, 10.0, 10.0)
        sites = (
            _site("C", (0.0, 0.0, 0.0), "C"),
            _site("N", (0.14, 0.0, 0.0), "N"),
        )
        first, second = ("C",), ("N",)
    else:
        edges = (2.0, 10.0, 10.0) if rank == 1 else (2.0, 2.0, 10.0) if rank == 2 else (2.0, 2.0, 2.0)
        cell = _cell(*edges)
        sites = (_site("C", (0.0, 0.0, 0.0), "C"),)
        first = second = ("C",)
    structure = CrystalStructure(f"rank-{rank}", cell, sites)
    context = SymmetryContext.from_operations((IDENTITY,), cell)
    request = CandidateInteraction(
        first, second, GrammarOperation.COVALENT_NETWORK,
        InteractionLayer.STRUCTURAL, InteractionPriority.PRIMARY,
        first, second, (ChemicalEvidence("fixture", "network"),),
    )
    grammar = CompositionGrammar(
        DecompositionMode.COVALENT_NETWORK,
        (request,), 1.0, (ChemicalEvidence("fixture", "grammar"),), (), "fixture",
    )
    contact_result = ContactAnalyzer(
        ShellResolutionPolicy(1.6, 0.01, 0.08, 0.01, 2.0)
    ).analyze(structure, context, grammar)
    polyhedra = PolyhedronOrbitBuilder().build(contact_result)
    graph = StructuralGraphBuilder().build(contact_result, polyhedra)
    representation = StructuralRepresentationBuilder(
        StructuralSelectionPolicy(
            frozenset({InteractionLayer.STRUCTURAL}),
            frozenset({ShellRole.PRIMARY}),
        )
    ).build(graph)
    return graph, PeriodicConnectivityAnalyzer().analyze(representation)


def test_structural_graph_edges_reference_contact_orbits_only() -> None:
    graph, _ = _rank_fixture(1)

    assert graph.connection_orbits
    assert all(item.source_resolved_contact_orbit_ids for item in graph.connection_orbits)
    assert all(not hasattr(item, "contact_id") for item in graph.connection_orbits)
    assert all(not hasattr(item, "source_contact_ids") for item in graph.connection_orbits)
    assert all(item.periodic_relation.operation_key for item in graph.connection_orbits)


def test_one_physical_contact_orbit_creates_only_one_connection_orbit() -> None:
    graph, _ = _rank_fixture(1)
    sources = tuple(
        source
        for connection in graph.connection_orbits
        for source in connection.source_resolved_contact_orbit_ids
    )

    assert len(sources) == len(set(sources))


def test_structural_units_do_not_store_expanded_contact_identity() -> None:
    names = {item.name for item in fields(StructuralUnitOrbit)}
    assert "units" not in names
    assert "source_contact_ids" not in names
    assert "source_resolved_contact_orbit_ids" in names


def test_integer_translation_basis_is_primitive() -> None:
    assert integer_translation_lattice_basis(((2, 0, 0), (3, 0, 0))) == ((1, 0, 0),)


def test_planar_unit_face_follows_polygon_boundary() -> None:
    coordinates = ((1.0, 1.0, 0.0), (-1.0, -1.0, 0.0), (1.0, -1.0, 0.0), (-1.0, 1.0, 0.0))
    face = _ordered_planar_face(coordinates)
    squared_edges = {
        sum((coordinates[face[(index + 1) % 4]][axis] - coordinates[face[index]][axis]) ** 2 for axis in range(3))
        for index in range(4)
    }

    assert squared_edges == {4.0}


def _affine_self_loop_rank(operation: AffineOperation) -> tuple[int, tuple[tuple[int, int, int], ...]]:
    context = SymmetryContext.from_operations((IDENTITY, operation), _cell(2.0, 2.0, 2.0))
    unit = StructuralUnitOrbit(
        "unit", StructuralUnitKind.ATOM, "A", None, None,
        (PeriodicAtomRef("image", (0, 0, 0)),), (), 1, (("fixture", "affine"),),
    )
    relation = PeriodicSymmetryRelation(
        next(key for key in context.operation_keys if key != context.identity_operation_key),
        (0, 0, 0),
    )
    edge = StructuralConnectionOrbit(
        "edge", "unit", "unit", relation, StructuralConnectionKind.DIRECT_CONTACT,
        (), (InteractionLayer.STRUCTURAL,), (), ("contact",), ("interpretation",), 1,
    )
    policy = StructuralSelectionPolicy(
        frozenset({InteractionLayer.STRUCTURAL}), frozenset({ShellRole.PRIMARY})
    )
    representation = StructuralRepresentation(
        "representation", (unit,), (edge,), policy, _symmetry_context=context
    )
    component = PeriodicConnectivityAnalyzer().analyze(representation).components[0]
    return component.rank, component.periodic_generators


def test_finite_inversion_cycle_has_rank_zero() -> None:
    inversion = AffineOperation(
        ((Fraction(-1), Fraction(0), Fraction(0)),
         (Fraction(0), Fraction(-1), Fraction(0)),
         (Fraction(0), Fraction(0), Fraction(-1))),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    assert _affine_self_loop_rank(inversion)[0] == 0


def test_screw_cycle_contributes_its_intrinsic_translation() -> None:
    screw = AffineOperation(
        ((Fraction(-1), Fraction(0), Fraction(0)),
         (Fraction(0), Fraction(-1), Fraction(0)),
         (Fraction(0), Fraction(0), Fraction(1))),
        (Fraction(0), Fraction(0), Fraction(1, 2)),
    )
    assert _affine_self_loop_rank(screw) == (1, ((0, 0, 1),))


def test_structural_blocks_keep_quotient_graph_identity() -> None:
    graph, connectivity = _rank_fixture(1)
    representation = StructuralRepresentationBuilder(
        StructuralSelectionPolicy(
            frozenset({InteractionLayer.STRUCTURAL}),
            frozenset({ShellRole.PRIMARY}),
        )
    ).build(graph)
    result = StructuralBlockFinder().find(representation, connectivity)

    assert len(result.blocks) == 1
    block = result.blocks[0]
    assert block.rank == 1
    assert block.unit_orbit_ids == tuple(item.unit_orbit_id for item in representation.unit_orbits)
    assert not hasattr(block, "atom_refs")
    assert not hasattr(block, "unit_ids")


@pytest.mark.parametrize("expected_rank", (0, 1, 2, 3))
def test_periodic_rank_comes_from_exact_relation_translations(expected_rank: int) -> None:
    _, connectivity = _rank_fixture(expected_rank)

    assert len(connectivity.components) == 1
    assert connectivity.components[0].rank == expected_rank
