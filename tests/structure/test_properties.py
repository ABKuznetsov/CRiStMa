import numpy as np
import pytest

from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
from cristma.structure import (
    AtomicProperty,
    AtomicPropertyTable,
    CrystalStructure,
    IndependentSite,
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
