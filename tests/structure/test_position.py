import numpy as np
import pytest

from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
from cristma.structure import (
    AtomicPropertyTable,
    AtomicView,
    CrystalStructure,
    ExpandedAtom,
    IndependentSite,
    MolecularAtom,
    MolecularStructure,
    SiteComponent,
)
from cristma.symmetry import expand_structure


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def component(element: str) -> SiteComponent:
    return SiteComponent(element, number(1.0))


def one_site_p1_crystal() -> CrystalStructure:
    site = IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(component("Si"),),
        fractional=(number(0.25), number(0.5), number(0.75)),
    )
    return CrystalStructure.explicit(
        "silicon",
        UnitCell.cubic(number(4.0)),
        (site,),
        id="structure:si",
    )


def test_crystal_view_retains_expanded_atom_objects() -> None:
    crystal = one_site_p1_crystal()

    view = expand_structure(crystal)

    assert isinstance(view.atoms[0], ExpandedAtom)
    assert view.atoms[0].source_site_id == crystal.sites[0].id
    assert view.cell is crystal.cell
    assert view.cell_matrix.flags.writeable is False
    with pytest.raises(ValueError):
        view.cartesian[0, 0] = 10.0


def test_molecular_view_retains_molecular_atom_objects() -> None:
    atom = MolecularAtom(
        id="atom:C1",
        label="C1",
        components=(component("C"),),
        cartesian=(1.0, 2.0, 3.0),
    )

    view = MolecularStructure("fragment", atoms=(atom,)).atomic_view()

    assert view.atoms == (atom,)
    assert np.array_equal(view.cartesian, [[1.0, 2.0, 3.0]])
    assert view.fractional is None
    assert view.cell is None
    assert view.periodic == (False, False, False)


def test_periodic_view_requires_cell_and_fractional_coordinates() -> None:
    atom = MolecularAtom("atom:C1", "C1", (component("C"),), (0.0, 0.0, 0.0))

    with pytest.raises(ValueError, match="periodic"):
        AtomicView(
            atoms=(atom,),
            cell=None,
            periodic=(True, False, False),
            properties=AtomicPropertyTable(1),
        )
