from __future__ import annotations

from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
from cristma.crystallography import SpaceGroupCatalog, build_orbit
from cristma.structure import IndependentSite, SiteComponent


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value=value, uncertainty=None, raw=str(value))


def site_at(x: float, y: float, z: float) -> IndependentSite:
    return IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(SiteComponent("Si", number(1.0)),),
        fractional=(number(x), number(y), number(z)),
    )


def cubic_cell() -> UnitCell:
    return UnitCell.cubic(number(4.0))


def test_p_minus_one_general_position_has_multiplicity_two() -> None:
    setting = SpaceGroupCatalog.default().by_setting(2)

    orbit = build_orbit(site_at(0.1, 0.2, 0.3), setting, cell=cubic_cell())

    assert orbit.multiplicity == 2
    assert orbit.calculated_multiplicity == 2
    assert len(orbit.stabilizer) == 1
    assert orbit.site_symmetry.order == 1
    assert orbit.site_symmetry.symbol is None


def test_p_minus_one_origin_is_special_position() -> None:
    setting = SpaceGroupCatalog.default().by_setting(2)

    orbit = build_orbit(site_at(0.0, 0.0, 0.0), setting, cell=cubic_cell())

    assert orbit.multiplicity == 1
    assert len(orbit.stabilizer) == 2
    assert orbit.site_symmetry.order == 2
    assert {image.operation_id for image in orbit.stabilizer} == {
        "hall:2:op:1",
        "hall:2:op:2",
    }


def test_orbit_keeps_calculated_facts_separate_from_reported_cif_values() -> None:
    setting = SpaceGroupCatalog.default().by_setting(2)
    reported_site = IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(SiteComponent("Si", number(1.0)),),
        fractional=(number(0.1), number(0.2), number(0.3)),
        wyckoff="a",
        reported_multiplicity=1,
    )

    orbit = build_orbit(reported_site, setting, cell=cubic_cell())

    assert orbit.multiplicity == 2
    assert orbit.representative.reported_multiplicity == 1
    assert not hasattr(orbit, "diagnostics")
