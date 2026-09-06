from __future__ import annotations

from fractions import Fraction

from cristma.core import MeasuredValue, UnitCell
from cristma.crystallography import (
    AsymmetricUnitMapper,
    PairTableStatus,
    SymmetryContext,
    SymmetryPairFinder,
)
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


def _site(site_id: str, x: float = 0.0) -> IndependentSite:
    return IndependentSite(
        site_id,
        site_id,
        (SiteComponent("C", _value(1.0)),),
        (_value(x), _value(0.0), _value(0.0)),
    )


def _chain(edge: float):
    cell = UnitCell(
        _value(edge), _value(10.0), _value(10.0),
        _value(90.0), _value(90.0), _value(90.0),
    )
    structure = CrystalStructure("chain", cell, (_site("A"),))
    context = SymmetryContext.from_operations((IDENTITY,), cell)
    mapping = AsymmetricUnitMapper().build(structure, context)
    return SymmetryPairFinder(cutoff=1.2).find(structure, context, mapping)


def test_geometry_id_ignores_metric_only_cell_change() -> None:
    first = _chain(1.0)
    second = _chain(1.05)

    assert tuple(item.geometry_orbit_id for item in first.contact_orbits) == tuple(
        item.geometry_orbit_id for item in second.contact_orbits
    )
    assert tuple(item.representative_distance for item in first.contact_orbits) != tuple(
        item.representative_distance for item in second.contact_orbits
    )


def test_pair_multiplicity_is_not_chain_coordination_number() -> None:
    table = _chain(1.0)

    assert table.status is PairTableStatus.COMPLETE
    assert len(table.contact_orbits) == 1
    assert table.contact_orbits[0].multiplicity_in_reference_cell == 1


def test_special_position_pair_orbit_has_exact_reference_cell_multiplicity() -> None:
    cell = UnitCell.cubic(_value(10.0))
    structure = CrystalStructure("special", cell, (_site("A", 0.0), _site("B", 0.1)))
    context = SymmetryContext.from_operations((IDENTITY, INVERSION), cell)
    mapping = AsymmetricUnitMapper().build(structure, context)

    table = SymmetryPairFinder(cutoff=1.1).find(structure, context, mapping)

    assert len(table.contact_orbits) == 1
    assert table.contact_orbits[0].multiplicity_in_reference_cell == 2


def test_search_limit_propagates_to_table_and_every_returned_orbit() -> None:
    cell = UnitCell.cubic(_value(1.0))
    structure = CrystalStructure("limited", cell, (_site("A"),))
    context = SymmetryContext.from_operations((IDENTITY,), cell)
    mapping = AsymmetricUnitMapper().build(structure, context)

    table = SymmetryPairFinder(cutoff=2.1, max_candidates=10).find(
        structure, context, mapping
    )

    assert table.status is PairTableStatus.INCOMPLETE
    assert all(
        orbit.status is PairTableStatus.INCOMPLETE for orbit in table.contact_orbits
    )
