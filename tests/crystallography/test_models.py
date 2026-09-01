from __future__ import annotations

from fractions import Fraction

import pytest

from cristma.crystallography import (
    AffineCoordinateMap,
    SpaceGroupSetting,
    WyckoffPosition,
)
from cristma.symmetry import SpaceGroupDefinition, parse_xyz_operation


def test_space_group_setting_builds_existing_definition() -> None:
    identity = parse_xyz_operation("x,y,z", operation_id="hall:1:op:1")
    setting = SpaceGroupSetting(
        setting_id=1,
        number=1,
        hall_symbol="P 1",
        choice="",
        hm_short="P 1",
        hm_full="P 1",
        point_group="1",
        centering="P",
        crystal_system="triclinic",
        symmetry_operations=(identity,),
        wyckoff_positions=(),
    )

    definition = setting.definition(provenance="derived")

    assert isinstance(definition, SpaceGroupDefinition)
    assert definition.number == 1
    assert definition.hall_symbol == "P 1"
    assert definition.operations == (identity,)


def test_affine_coordinate_map_evaluates_exact_parameterization() -> None:
    representative = AffineCoordinateMap.from_xyz("x,x+1/2,-z")

    assert representative.evaluate(
        (Fraction(1, 4), Fraction(0), Fraction(1, 3))
    ) == (Fraction(1, 4), Fraction(3, 4), Fraction(-1, 3))
    assert representative.degrees_of_freedom == 2


def test_affine_coordinate_map_accepts_integer_coefficients() -> None:
    representative = AffineCoordinateMap.from_xyz("2x-y,0,z")

    assert representative.evaluate(
        (Fraction(1, 3), Fraction(1, 6), Fraction(1, 4))
    ) == (Fraction(1, 2), Fraction(0), Fraction(1, 4))


def test_wyckoff_position_rejects_wrong_representative_count() -> None:
    with pytest.raises(ValueError, match="multiplicity"):
        WyckoffPosition(
            setting_id=390,
            letter="a",
            multiplicity=2,
            site_symmetry_symbol="-4..",
            coordinate_constraints=(AffineCoordinateMap.from_xyz("0,0,0"),),
        )


def test_wyckoff_position_accepts_spglib_overflow_letter() -> None:
    position = WyckoffPosition(
        setting_id=227,
        letter="A",
        multiplicity=8,
        site_symmetry_symbol="1",
        coordinate_constraints=tuple(
            AffineCoordinateMap.from_xyz("x,y,z") for _ in range(8)
        ),
    )

    assert position.letter == "A"


def test_space_group_setting_rejects_dataset_number_outside_catalog() -> None:
    with pytest.raises(ValueError, match="Hall number"):
        SpaceGroupSetting(
            setting_id=531,
            number=1,
            hall_symbol="invalid",
            choice="",
            hm_short="P 1",
            hm_full="P 1",
            point_group="1",
            centering="P",
            crystal_system="triclinic",
            symmetry_operations=(
                parse_xyz_operation("x,y,z", operation_id="invalid"),
            ),
            wyckoff_positions=(),
        )
