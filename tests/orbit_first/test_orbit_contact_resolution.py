from __future__ import annotations

from dataclasses import fields
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
    ContactOrbitResolver,
    EndpointRole,
    OrientationMode,
    ResolutionStatus,
    ShellResolutionPolicy,
)
from cristma.crystallography import AsymmetricUnitMapper, SymmetryContext, SymmetryPairFinder
from cristma.reference_data import CovalentRadii, CovalentRadiusRecord, ReferenceData
from cristma.structure import CrystalStructure, IndependentSite, SiteComponent
from cristma.symmetry import AffineOperation


def _value(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def _site(
    site_id: str,
    x: float,
    components: tuple[tuple[str, float], ...],
) -> IndependentSite:
    return IndependentSite(
        site_id,
        site_id,
        tuple(SiteComponent(symbol, _value(occupancy)) for symbol, occupancy in components),
        (_value(x), _value(0.0), _value(0.0)),
    )


def _pair_table(first: IndependentSite, second: IndependentSite):
    cell = UnitCell.cubic(_value(10.0))
    structure = CrystalStructure("chemistry", cell, (first, second))
    identity = AffineOperation(
        (
            (Fraction(1), Fraction(0), Fraction(0)),
            (Fraction(0), Fraction(1), Fraction(0)),
            (Fraction(0), Fraction(0), Fraction(1)),
        ),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    context = SymmetryContext.from_operations((identity,), cell)
    mapping = AsymmetricUnitMapper().build(structure, context)
    table = SymmetryPairFinder(cutoff=2.1).find(structure, context, mapping)
    return structure, table


def _request(operation: GrammarOperation) -> CandidateInteraction:
    evidence = (ChemicalEvidence("fixture", "fixture interaction"),)
    return CandidateInteraction(
        ("Ca", "Sr"),
        ("O",),
        operation,
        (
            InteractionLayer.COORDINATION
            if operation is GrammarOperation.CENTRE_LIGAND_SHELL
            else InteractionLayer.STRUCTURAL
        ),
        InteractionPriority.PRIMARY,
        ("Ca", "Sr"),
        ("O",),
        evidence,
    )


def _grammar(*requests: CandidateInteraction) -> CompositionGrammar:
    return CompositionGrammar(
        DecompositionMode.COVALENT_NETWORK,
        tuple(requests),
        1.0,
        (ChemicalEvidence("fixture", "fixture grammar"),),
        (),
        "fixture",
    )


POLICY = ShellResolutionPolicy(1.60, 0.01, 0.08, 0.01, 2.0)


def test_mixed_pair_preserves_every_supported_interpretation() -> None:
    structure, pair_table = _pair_table(
        _site("A", 0.0, (("Ca", 0.7), ("Sr", 0.3))),
        _site("B", 0.2, (("O", 1.0),)),
    )
    grammar = _grammar(
        _request(GrammarOperation.CENTRE_LIGAND_SHELL),
        _request(GrammarOperation.COVALENT_NETWORK),
    )

    resolved = ContactOrbitResolver(POLICY).resolve(pair_table, structure, grammar)
    orbit = resolved.contact_orbits[0]

    assert tuple(item.interaction_type for item in orbit.interpretations) == (
        GrammarOperation.CENTRE_LIGAND_SHELL,
        GrammarOperation.COVALENT_NETWORK,
    )
    assert tuple(len(item.component_pair_interpretations) for item in orbit.interpretations) == (2, 2)
    assert orbit.interpretations[0].orientation_mode is OrientationMode.ENDPOINT_ROLES
    assert orbit.interpretations[0].endpoint_roles == (
        EndpointRole.CENTER,
        EndpointRole.LIGAND,
    )
    assert orbit.interpretations[1].orientation_mode is OrientationMode.UNDIRECTED
    assert resolved.status is ResolutionStatus.RESOLVED


def test_missing_radius_marks_only_affected_interpretation_incomplete() -> None:
    structure, pair_table = _pair_table(
        _site("A", 0.0, (("Sr", 1.0),)),
        _site("B", 0.2, (("O", 1.0),)),
    )
    defaults = ReferenceData.default()
    reference = ReferenceData(
        defaults.elements,
        CovalentRadii((CovalentRadiusRecord("O", 0.66),)),
        defaults.shannon_radii,
        defaults.chemical,
    )

    resolved = ContactOrbitResolver(POLICY, reference).resolve(
        pair_table,
        structure,
        _grammar(_request(GrammarOperation.CENTRE_LIGAND_SHELL)),
    )

    assert resolved.contact_orbits[0].status is ResolutionStatus.INCOMPLETE
    assert resolved.contact_orbits[0].interpretations[0].status is ResolutionStatus.INCOMPLETE
    assert resolved.contact_orbits[0].interpretations[0].normalized_distance_range is None
    assert "crystal_chemistry.contact.radius_missing" in {
        diagnostic.code for diagnostic in resolved.diagnostics
    }


def test_orbit_result_contains_no_materialized_contact_or_expanded_atom_state() -> None:
    structure, pair_table = _pair_table(
        _site("A", 0.0, (("Ca", 1.0),)),
        _site("B", 0.2, (("O", 1.0),)),
    )

    resolved = ContactOrbitResolver(POLICY).resolve(
        pair_table,
        structure,
        _grammar(_request(GrammarOperation.CENTRE_LIGAND_SHELL)),
    )

    orbit_field_names = {item.name for item in fields(resolved.contact_orbits[0])}
    assert "contacts" not in orbit_field_names
    assert "contact_id" not in orbit_field_names
    assert "expanded_atom_id" not in orbit_field_names
