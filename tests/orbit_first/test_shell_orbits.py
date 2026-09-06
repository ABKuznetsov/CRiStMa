from __future__ import annotations

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
    ContactIncidenceBuilder,
    ContactOrbitResolver,
    CoordinationShellOrbitResolver,
    ResolutionStatus,
    ShellResolutionPolicy,
)
from cristma.crystallography import AsymmetricUnitMapper, SymmetryContext, SymmetryPairFinder
from cristma.structure import CrystalStructure, IndependentSite, SiteComponent
from cristma.symmetry import AffineOperation


def _value(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def _site(site_id: str, x: float, symbol: str, occupancy: float = 1.0) -> IndependentSite:
    return IndependentSite(
        site_id,
        site_id,
        (SiteComponent(symbol, _value(occupancy)),),
        (_value(x), _value(0.0), _value(0.0)),
    )


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


def _shells(
    policy: ShellResolutionPolicy,
    positions: tuple[float, float, float] = (0.10, 0.105, 0.20),
):
    cell = UnitCell.cubic(_value(20.0))
    structure = CrystalStructure(
        "weighted-shell",
        cell,
        (
            _site("M", 0.0, "Ca"),
            _site("X1", positions[0], "O", 0.6),
            _site("X2", positions[1], "O"),
            _site("X3", positions[2], "O"),
        ),
    )
    context = SymmetryContext.from_operations((IDENTITY, INVERSION), cell)
    mapping = AsymmetricUnitMapper().build(structure, context)
    table = SymmetryPairFinder(cutoff=4.5).find(structure, context, mapping)
    request = CandidateInteraction(
        ("Ca",),
        ("O",),
        GrammarOperation.CENTRE_LIGAND_SHELL,
        InteractionLayer.COORDINATION,
        InteractionPriority.PRIMARY,
        ("Ca",),
        ("O",),
        (ChemicalEvidence("fixture", "fixture interaction"),),
    )
    grammar = CompositionGrammar(
        DecompositionMode.COVALENT_NETWORK,
        (request,),
        1.0,
        (ChemicalEvidence("fixture", "fixture grammar"),),
        (),
        "fixture",
    )
    contact_resolution = ContactOrbitResolver(policy).resolve(table, structure, grammar)
    incidences = ContactIncidenceBuilder().build(
        table,
        contact_resolution.contact_orbits,
        structure,
        mapping,
        context,
    )
    shells = CoordinationShellOrbitResolver(policy).resolve(
        table,
        contact_resolution.contact_orbits,
        incidences,
    )
    return shells, incidences


def test_shell_cn_and_occupancy_use_incidence_weights() -> None:
    shells, _ = _shells(ShellResolutionPolicy(1.6, 0.01, 0.20, 0.01, 2.0))

    shell = shells[0]
    assert shell.status is ResolutionStatus.RESOLVED
    assert shell.selected is not None
    assert shell.selected.geometric_CN == 4
    assert shell.selected.mean_occupied_neighbors == pytest.approx(3.2)


def test_ambiguous_shell_keeps_alternatives_without_selecting_one() -> None:
    shells, incidences = _shells(
        ShellResolutionPolicy(1.6, 0.01, 0.20, 1.0, 2.0),
        (0.10, 0.15, 0.20),
    )

    shell = shells[0]
    assert shell.status is ResolutionStatus.AMBIGUOUS
    assert shell.selected_alternative is None
    assert len(shell.alternatives) >= 2
    middle = next(
        item.incidence_orbit_id
        for item in incidences
        if item.ligand_independent_site_id == "X2"
    )
    assert any(middle in item.primary_incidence_ids for item in shell.alternatives)
    assert any(middle in item.secondary_incidence_ids for item in shell.alternatives)


def test_network_only_contacts_do_not_create_coordination_shells() -> None:
    cell = UnitCell.cubic(_value(10.0))
    structure = CrystalStructure(
        "network",
        cell,
        (_site("C", 0.0, "C"), _site("N", 0.14, "N")),
    )
    context = SymmetryContext.from_operations((IDENTITY,), cell)
    mapping = AsymmetricUnitMapper().build(structure, context)
    table = SymmetryPairFinder(cutoff=2.0).find(structure, context, mapping)
    request = CandidateInteraction(
        ("C",), ("N",), GrammarOperation.COVALENT_NETWORK,
        InteractionLayer.STRUCTURAL, InteractionPriority.PRIMARY,
        ("C",), ("N",), (ChemicalEvidence("fixture", "network"),),
    )
    grammar = CompositionGrammar(
        DecompositionMode.COVALENT_NETWORK,
        (request,), 1.0, (ChemicalEvidence("fixture", "grammar"),), (), "fixture",
    )
    policy = ShellResolutionPolicy(1.6, 0.01, 0.08, 0.01, 2.0)
    contacts = ContactOrbitResolver(policy).resolve(table, structure, grammar)
    incidences = ContactIncidenceBuilder().build(
        table, contacts.contact_orbits, structure, mapping, context,
    )

    assert CoordinationShellOrbitResolver(policy).resolve(
        table, contacts.contact_orbits, incidences,
    ) == ()
