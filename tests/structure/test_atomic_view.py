import numpy as np

from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
from cristma.structure import (
    CrystalStructure,
    IndependentSite,
    MolecularAtom,
    MolecularStructure,
    SiteComponent,
)


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def test_molecular_atomic_view_has_cartesian_coordinates_and_no_cell() -> None:
    molecule = MolecularStructure(
        "water",
        atoms=(MolecularAtom("atom:O", "O", "O", (1.0, 2.0, 3.0)),),
    )

    view = molecule.atomic_view()

    assert np.array_equal(view.cartesian, [[1.0, 2.0, 3.0]])
    assert view.fractional is None
    assert view.cell is None


def test_mixed_crystal_site_stays_one_geometric_row() -> None:
    site = IndependentSite(
        id="site:M1",
        label="M1",
        components=(
            SiteComponent("Fe", number(0.75)),
            SiteComponent("Mg", number(0.25)),
        ),
        fractional=(number(0.5), number(0.0), number(0.0)),
    )
    crystal = CrystalStructure.explicit(
        "mixed",
        UnitCell.cubic(number(4.0)),
        (site,),
        id="structure:mixed",
    )

    view = crystal.atomic_view()

    assert view.ids == ("expanded:structure:mixed:site:M1:op:1:0,0,0",)
    assert view.source_site_ids == ("site:M1",)
    assert np.array_equal(view.cartesian, [[2.0, 0.0, 0.0]])
    assert len(view.properties["site_components"].values[0]) == 2
