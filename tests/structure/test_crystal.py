import pytest

from cristma.chemistry.species import ElementSpecies
from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
from cristma.structure import Crystal, CrystalStructure, IndependentSite, SiteComponent
from cristma.symmetry.affine import parse_xyz_operation
from cristma.symmetry.orbit import expand_orbit


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def test_component_occupancy_above_one_is_invalid() -> None:
    with pytest.raises(ValueError, match="occupancy must lie between zero and one"):
        SiteComponent("Ca", number(1.01))


def test_site_reports_total_occupancy_and_vacancy_fraction() -> None:
    site = IndependentSite(
        id="site:M1",
        label="M1",
        components=(
            SiteComponent("Ca", number(0.6)),
            SiteComponent("Sr", number(0.2)),
        ),
        fractional=(number(0), number(0), number(0)),
    )

    assert site.total_occupancy == pytest.approx(0.8)
    assert site.vacancy_fraction == pytest.approx(0.2)


def test_public_crystal_name_and_compatibility_alias() -> None:
    site = IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(SiteComponent("Si", number(1)),),
        fractional=(number(0), number(0), number(0)),
    )
    structure = CrystalStructure("Si", UnitCell.cubic(number(5.43)), (site,))

    assert isinstance(structure, Crystal)
    assert structure.periodic == (True, True, True)
    assert structure.sites[0].components[0].species == ElementSpecies("Si")
    assert structure.sites[0].components[0].element == "Si"


def test_identity_for_explicit_dft_sites_is_not_reported_p1() -> None:
    structure = CrystalStructure.explicit(
        name="simulation",
        cell=UnitCell.cubic(number(4)),
        sites=(),
    )

    assert structure.space_group.provenance == "unreported_identity"
    assert structure.space_group.number is None


def test_expanded_atom_identity_resolves_to_independent_site() -> None:
    site = IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(SiteComponent("Si", number(1)),),
        fractional=(number(0), number(0), number(0)),
    )
    atom = expand_orbit(
        site,
        (parse_xyz_operation("x,y,z", operation_id="op:1"),),
        structure_id="structure:Si",
    )[0]

    assert atom.structure_id == "structure:Si"
    assert atom.source_site_id == "site:Si1"
    assert atom.independent_site_id == "site:Si1"
    assert atom.id == "expanded:structure:Si:site:Si1:op:1:0,0,0"
