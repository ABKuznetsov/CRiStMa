from __future__ import annotations

from fractions import Fraction

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
    ContactIncidenceBuilder,
    ContactOrbitResolver,
    ShellResolutionPolicy,
)
from cristma.crystallography import AsymmetricUnitMapper, SymmetryContext, SymmetryPairFinder
from cristma.structure import CrystalStructure, IndependentSite, SiteComponent
from cristma.symmetry import AffineOperation


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
INVERSION = AffineOperation(
    (
        (Fraction(-1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(-1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(-1)),
    ),
    (Fraction(0), Fraction(0), Fraction(0)),
)


def _site(site_id: str, x: float, components) -> IndependentSite:
    return IndependentSite(
        site_id,
        site_id,
        tuple(SiteComponent(symbol, _value(occupancy)) for symbol, occupancy in components),
        (_value(x), _value(0.0), _value(0.0)),
    )


def _request(
    first: tuple[str, ...],
    second: tuple[str, ...],
    operation: GrammarOperation,
    *,
    centers: tuple[str, ...],
    ligands: tuple[str, ...],
) -> CandidateInteraction:
    return CandidateInteraction(
        first,
        second,
        operation,
        InteractionLayer.COORDINATION,
        InteractionPriority.PRIMARY,
        centers,
        ligands,
        (ChemicalEvidence("fixture", "fixture interaction"),),
    )


def _analyze(cell: UnitCell, sites, request, policy=None, operations=(IDENTITY,)):
    structure = CrystalStructure("incidence", cell, tuple(sites))
    context = SymmetryContext.from_operations(tuple(operations), cell)
    mapping = AsymmetricUnitMapper().build(structure, context)
    table = SymmetryPairFinder(cutoff=1.2).find(structure, context, mapping)
    grammar = CompositionGrammar(
        DecompositionMode.COVALENT_NETWORK,
        (request,),
        1.0,
        (ChemicalEvidence("fixture", "fixture grammar"),),
        (),
        "fixture",
    )
    selected_policy = policy or ShellResolutionPolicy(1.6, 0.01, 0.08, 0.01, 2.0)
    contacts = ContactOrbitResolver(selected_policy).resolve(table, structure, grammar)
    incidences = ContactIncidenceBuilder().build(
        table,
        contacts.contact_orbits,
        structure,
        mapping,
        context,
    )
    return incidences


def test_one_pair_orbit_gives_two_local_incidences_in_simple_chain() -> None:
    cell = UnitCell(
        _value(1.0), _value(10.0), _value(10.0),
        _value(90.0), _value(90.0), _value(90.0),
    )
    incidences = _analyze(
        cell,
        (_site("A", 0.0, (("C", 1.0),)),),
        _request(
            ("C",), ("C",), GrammarOperation.COVALENT_NETWORK,
            centers=("C",), ligands=("C",),
        ),
    )

    assert len(incidences) == 1
    assert incidences[0].incidence_multiplicity_per_center == 2


def test_effective_occupancy_uses_only_participating_component() -> None:
    cell = UnitCell.cubic(_value(10.0))
    incidences = _analyze(
        cell,
        (
            _site("M", 0.0, (("Ca", 1.0),)),
            _site("X", 0.1, (("O", 0.6), ("F", 0.4))),
        ),
        _request(
            ("Ca",), ("O",), GrammarOperation.CENTRE_LIGAND_SHELL,
            centers=("Ca",), ligands=("O",),
        ),
    )

    assert len(incidences) == 1
    assert incidences[0].center_independent_site_id == "M"
    assert incidences[0].ligand_independent_site_id == "X"
    assert incidences[0].effective_neighbor_occupancy == 0.6


def test_incidence_identity_does_not_include_shell_policy() -> None:
    cell = UnitCell.cubic(_value(10.0))
    sites = (
        _site("M", 0.0, (("Ca", 1.0),)),
        _site("X", 0.1, (("O", 1.0),)),
    )
    request = _request(
        ("Ca",), ("O",), GrammarOperation.CENTRE_LIGAND_SHELL,
        centers=("Ca",), ligands=("O",),
    )

    first = _analyze(cell, sites, request, ShellResolutionPolicy(1.5, 0.01, 0.08, 0.01, 2.0))
    second = _analyze(cell, sites, request, ShellResolutionPolicy(1.7, 0.02, 0.10, 0.02, 2.2))

    assert tuple(item.incidence_orbit_id for item in first) == tuple(
        item.incidence_orbit_id for item in second
    )


def test_undirected_pair_with_distinct_sites_has_one_incidence_per_center() -> None:
    cell = UnitCell.cubic(_value(10.0))
    incidences = _analyze(
        cell,
        (
            _site("A", 0.0, (("C", 1.0),)),
            _site("B", 0.1, (("N", 1.0),)),
        ),
        _request(
            ("C",), ("N",), GrammarOperation.COVALENT_NETWORK,
            centers=("C",), ligands=("N",),
        ),
    )

    assert {(item.center_independent_site_id, item.ligand_independent_site_id) for item in incidences} == {
        ("A", "B"),
        ("B", "A"),
    }


def test_special_position_incidence_is_invariant_to_operation_order() -> None:
    cell = UnitCell.cubic(_value(10.0))
    sites = (
        _site("M", 0.0, (("Ca", 1.0),)),
        _site("X", 0.1, (("O", 1.0),)),
    )
    request = _request(
        ("Ca",), ("O",), GrammarOperation.CENTRE_LIGAND_SHELL,
        centers=("Ca",), ligands=("O",),
    )

    first = _analyze(cell, sites, request, operations=(IDENTITY, INVERSION))
    second = _analyze(cell, sites, request, operations=(INVERSION, IDENTITY))

    assert first == second
    assert len(first) == 1
    assert first[0].incidence_multiplicity_per_center == 2
