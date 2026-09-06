from __future__ import annotations

from fractions import Fraction

from cristma.core import MeasuredValue, UnitCell
from cristma.crystallography import (
    AsymmetricUnitMapper,
    PeriodicSymmetryRelation,
    SymmetryContext,
    SymmetryPairCandidate,
    canonical_operation_key,
    canonical_pair_relation,
)
from cristma.structure import CrystalStructure, IndependentSite, SiteComponent
from cristma.symmetry import AffineOperation


def _value(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def _operation(rotation) -> AffineOperation:
    return AffineOperation(
        tuple(tuple(Fraction(value) for value in row) for row in rotation),
        (Fraction(0), Fraction(0), Fraction(0)),
    )


IDENTITY = _operation(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
INVERSION = _operation(((-1, 0, 0), (0, -1, 0), (0, 0, -1)))


def _site(site_id: str, x: float) -> IndependentSite:
    return IndependentSite(
        site_id,
        site_id,
        (SiteComponent("C", _value(1.0)),),
        (_value(x), _value(0.0), _value(0.0)),
    )


def test_special_position_descriptions_collapse_under_endpoint_stabilizer() -> None:
    cell = UnitCell.cubic(_value(10.0))
    structure = CrystalStructure("special", cell, (_site("A", 0.0), _site("B", 0.1)))
    context = SymmetryContext.from_operations((IDENTITY, INVERSION), cell)
    mapping = AsymmetricUnitMapper().build(structure, context)
    forward = SymmetryPairCandidate(
        "A",
        "B",
        PeriodicSymmetryRelation(canonical_operation_key(IDENTITY), (0, 0, 0)),
        1.0,
        (1.0, 0.0, 0.0),
    )
    symmetry_related = SymmetryPairCandidate(
        "A",
        "B",
        PeriodicSymmetryRelation(canonical_operation_key(INVERSION), (0, 0, 0)),
        1.0,
        (-1.0, 0.0, 0.0),
    )

    descriptors = {
        canonical_pair_relation(
            candidate,
            mapping.by_site_id["A"],
            mapping.by_site_id["B"],
            context,
        )
        for candidate in (forward, symmetry_related)
    }

    assert len(descriptors) == 1
    descriptor = descriptors.pop()
    assert descriptor[0:2] == ("A", "B")


def test_undirected_self_relation_and_inverse_have_one_identity() -> None:
    cell = UnitCell(
        _value(1.0), _value(10.0), _value(10.0),
        _value(90.0), _value(90.0), _value(90.0),
    )
    structure = CrystalStructure("chain", cell, (_site("A", 0.0),))
    context = SymmetryContext.from_operations((IDENTITY,), cell)
    mapping = AsymmetricUnitMapper().build(structure, context)
    operation_key = canonical_operation_key(IDENTITY)

    descriptors = {
        canonical_pair_relation(
            SymmetryPairCandidate(
                "A", "A", PeriodicSymmetryRelation(operation_key, translation),
                1.0, (float(translation[0]), 0.0, 0.0),
            ),
            mapping.by_site_id["A"],
            mapping.by_site_id["A"],
            context,
        )
        for translation in ((1, 0, 0), (-1, 0, 0))
    }

    assert len(descriptors) == 1
    assert descriptors.pop()[2].lattice_translation == (-1, 0, 0)
