from __future__ import annotations

import pytest

from cristma.chemistry import Composition
from cristma.chemistry.species import UnknownSpecies
from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
from cristma.structure import CrystalStructure, IndependentSite, MolecularAtom, MolecularStructure, SiteComponent
from cristma.symmetry.affine import parse_xyz_operation
from cristma.symmetry.orbit import SpaceGroupDefinition


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def component(symbol: str, occupancy: float = 1.0) -> SiteComponent:
    return SiteComponent(symbol, number(occupancy))


def test_composition_normalizes_symbols_and_formula() -> None:
    composition = Composition.from_mapping({"ca": 1, "O": 1})

    assert composition.elements == ("Ca", "O")
    assert composition.normalized_formula == "CaO"
    assert composition.amount("ca") == pytest.approx(1.0)


def test_molecular_structure_counts_each_position_once() -> None:
    molecule = MolecularStructure(
        "water",
        (
            MolecularAtom("O1", "O1", (component("O"),), (0.0, 0.0, 0.0)),
            MolecularAtom("H1", "H1", (component("H"),), (1.0, 0.0, 0.0)),
            MolecularAtom("H2", "H2", (component("H"),), (-1.0, 0.0, 0.0)),
        ),
    )

    composition = Composition.from_structure(molecule)

    assert composition.as_dict() == {"H": 2.0, "O": 1.0}
    assert composition.normalized_formula == "H2O"


def test_independent_sites_use_occupancy_times_calculated_multiplicity() -> None:
    site = IndependentSite(
        id="site:Ca1",
        label="Ca1",
        components=(component("Ca", 0.5),),
        fractional=(number(0), number(0), number(0)),
        calculated_multiplicity=4,
    )
    crystal = CrystalStructure(
        "calcium",
        UnitCell.cubic(number(4)),
        (site,),
        space_group=SpaceGroupDefinition(
            operations=(
                parse_xyz_operation("x,y,z", operation_id="op:1"),
                parse_xyz_operation("-x,-y,-z", operation_id="op:2"),
            ),
            provenance="reported",
        ),
    )

    assert Composition.from_structure(crystal).amount("Ca") == pytest.approx(2.0)


def test_explicit_identity_structure_uses_multiplicity_one() -> None:
    site = IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(component("Si"),),
        fractional=(number(0), number(0), number(0)),
    )
    crystal = CrystalStructure.explicit("silicon", UnitCell.cubic(number(5.4)), (site,))

    assert Composition.from_structure(crystal).amount("Si") == pytest.approx(1.0)


def test_nonidentity_structure_without_calculated_multiplicity_is_rejected() -> None:
    site = IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(component("Si"),),
        fractional=(number(0.1), number(0.2), number(0.3)),
    )
    crystal = CrystalStructure(
        "silicon",
        UnitCell.cubic(number(5.4)),
        (site,),
        space_group=SpaceGroupDefinition(
            operations=(
                parse_xyz_operation("x,y,z", operation_id="op:1"),
                parse_xyz_operation("-x,-y,-z", operation_id="op:2"),
            ),
            provenance="reported",
        ),
    )

    with pytest.raises(ValueError, match="calculated multiplicity"):
        Composition.from_structure(crystal)


def test_unknown_occupied_species_is_rejected() -> None:
    atom = MolecularAtom(
        "X1",
        "X1",
        (SiteComponent(UnknownSpecies("unknown:X1"), number(1)),),
        (0.0, 0.0, 0.0),
    )

    with pytest.raises(ValueError, match="no known element"):
        Composition.from_structure(MolecularStructure("unknown", (atom,)))
