from __future__ import annotations

from dataclasses import replace

from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
import pytest

from cristma.crystallography import SpaceGroupCatalog, assign_wyckoff, build_orbit
from cristma.structure import IndependentSite, SiteComponent


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value=value, uncertainty=None, raw=str(value))


def site_at(
    x: float,
    y: float,
    z: float,
    *,
    wyckoff: str | None = None,
    reported_multiplicity: int | None = None,
) -> IndependentSite:
    return IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(SiteComponent("Si", number(1.0)),),
        fractional=(number(x), number(y), number(z)),
        wyckoff=wyckoff,
        reported_multiplicity=reported_multiplicity,
    )


def cubic_cell() -> UnitCell:
    return UnitCell.cubic(number(4.0))


def tetragonal_cell() -> UnitCell:
    right = number(90.0)
    return UnitCell(number(4.0), number(4.0), number(6.0), right, right, right)


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


@pytest.mark.parametrize(
    ("fractional", "letter", "multiplicity", "site_symmetry"),
    [
        ((0.0, 0.0, 0.0), "a", 2, "-4.."),
        ((0.0, 0.0, 0.5), "b", 2, "-4.."),
        ((0.0, 0.5, 0.23), "c", 2, "2.mm"),
        ((0.0, 0.0, 0.23), "d", 4, "2.."),
        ((0.17, 0.67, 0.23), "e", 4, "..m"),
        ((0.17, 0.29, 0.23), "f", 8, "1"),
    ],
)
def test_p421m_wyckoff_positions_are_identified(
    fractional: tuple[float, float, float],
    letter: str,
    multiplicity: int,
    site_symmetry: str,
) -> None:
    setting = SpaceGroupCatalog.default().by_setting(390)
    orbit = build_orbit(site_at(*fractional), setting, cell=tetragonal_cell())

    assignment = assign_wyckoff(orbit, setting)

    assert assignment.position is not None
    assert assignment.position.letter == letter
    assert assignment.calculated_multiplicity == multiplicity
    assert assignment.site_symmetry.symbol == site_symmetry
    assert assignment.status == "matched"


def test_rounded_special_coordinate_respects_explicit_tolerance() -> None:
    setting = SpaceGroupCatalog.default().by_setting(390)
    orbit = build_orbit(
        site_at(0.0, 0.5, 0.2300004),
        setting,
        cell=tetragonal_cell(),
    )

    assignment = assign_wyckoff(orbit, setting, tolerance=1e-6)

    assert assignment.position is not None
    assert assignment.position.letter == "c"


def test_reported_wyckoff_and_multiplicity_mismatches_are_diagnostics() -> None:
    setting = SpaceGroupCatalog.default().by_setting(390)
    orbit = build_orbit(
        site_at(0.0, 0.0, 0.0, wyckoff="b", reported_multiplicity=4),
        setting,
        cell=tetragonal_cell(),
    )

    assignment = assign_wyckoff(orbit, setting)
    codes = {diagnostic.code for diagnostic in assignment.diagnostics}

    assert assignment.position is not None
    assert assignment.position.letter == "a"
    assert "crystallography.orbit.reported_wyckoff_mismatch" in codes
    assert "crystallography.orbit.reported_multiplicity_mismatch" in codes


def test_reported_wyckoff_catalog_multiplicity_is_checked_independently() -> None:
    setting = SpaceGroupCatalog.default().by_setting(390)
    orbit = build_orbit(
        site_at(0.0, 0.0, 0.0, wyckoff="d"),
        setting,
        cell=tetragonal_cell(),
    )

    assignment = assign_wyckoff(orbit, setting)

    assert "crystallography.orbit.wyckoff_multiplicity_mismatch" in {
        diagnostic.code for diagnostic in assignment.diagnostics
    }


def test_unresolved_and_ambiguous_assignments_are_explicit() -> None:
    setting = SpaceGroupCatalog.default().by_setting(390)
    orbit = build_orbit(site_at(0.0, 0.0, 0.0), setting, cell=tetragonal_cell())

    unresolved = assign_wyckoff(
        orbit,
        replace(setting, wyckoff_positions=()),
    )
    duplicate = replace(setting.wyckoff_positions[-1], letter="q")
    ambiguous = assign_wyckoff(
        orbit,
        replace(
            setting,
            wyckoff_positions=(*setting.wyckoff_positions, duplicate),
        ),
    )

    assert unresolved.status == "unresolved"
    assert unresolved.position is None
    assert ambiguous.status == "ambiguous"
    assert ambiguous.position is None
