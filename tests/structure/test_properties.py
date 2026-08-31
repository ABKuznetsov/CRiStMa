import numpy as np
import pytest

from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
from cristma.structure import (
    AtomicProperty,
    AtomicPropertyTable,
    CrystalStructure,
    IndependentSite,
    MolecularAtom,
    MolecularStructure,
    SiteComponent,
)


def test_property_values_are_immutable_and_typed() -> None:
    prop = AtomicProperty("magnetic_moment", np.array([1.0, -1.0]), unit="mu_B")
    table = AtomicPropertyTable(2, (prop,))

    assert table["magnetic_moment"].unit == "mu_B"
    with pytest.raises(ValueError):
        table["magnetic_moment"].values[0] = 0


def test_property_length_must_match_atoms() -> None:
    with pytest.raises(ValueError, match="leading dimension"):
        AtomicPropertyTable(2, (AtomicProperty("charge", np.array([0.0])),))


def test_crystal_property_rows_must_match_independent_sites() -> None:
    value = MeasuredValue(1.0, None, "1")
    site = IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(SiteComponent("Si", value),),
        fractional=(
            MeasuredValue(0.0, None, "0"),
            MeasuredValue(0.0, None, "0"),
            MeasuredValue(0.0, None, "0"),
        ),
    )
    table = AtomicPropertyTable(
        1,
        (AtomicProperty("selective_dynamics", np.array([[True, False, True]])),),
    )

    crystal = CrystalStructure.explicit(
        "demo",
        UnitCell.cubic(MeasuredValue(4.0, None, "4")),
        (site,),
        properties=table,
    )

    assert crystal.properties is table

    with pytest.raises(ValueError, match="independent sites"):
        CrystalStructure.explicit(
            "invalid",
            crystal.cell,
            (site,),
            properties=AtomicPropertyTable(0),
        )


def test_molecular_properties_reach_atomic_view() -> None:
    value = MeasuredValue(1.0, None, "1")
    atoms = (
        MolecularAtom("atom:O", "O", (SiteComponent("O", value),), (0.0, 0.0, 0.0)),
        MolecularAtom("atom:H1", "H1", (SiteComponent("H", value),), (1.0, 0.0, 0.0)),
        MolecularAtom("atom:H2", "H2", (SiteComponent("H", value),), (-1.0, 0.0, 0.0)),
    )
    table = AtomicPropertyTable(
        3,
        (AtomicProperty("charge", np.array([-0.8, 0.4, 0.4])),),
    )

    molecule = MolecularStructure("water", atoms, properties=table)

    assert molecule.properties is table
    assert molecule.atomic_view().properties is table


def test_molecular_property_rows_must_match_atoms() -> None:
    value = MeasuredValue(1.0, None, "1")
    atom = MolecularAtom("atom:H", "H", (SiteComponent("H", value),), (0.0, 0.0, 0.0))

    with pytest.raises(ValueError, match="molecular atoms"):
        MolecularStructure("invalid", (atom,), properties=AtomicPropertyTable(0))
