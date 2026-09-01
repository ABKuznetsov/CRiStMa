from __future__ import annotations

from fractions import Fraction

import pytest

from cristma.crystallography import (
    AffineCoordinateMap,
    SpaceGroupKey,
    SpaceGroupRecord,
    WyckoffPosition,
)
from cristma.symmetry import SpaceGroupDefinition, parse_xyz_operation


def test_space_group_record_builds_existing_definition() -> None:
    identity = parse_xyz_operation("x,y,z", operation_id="hall:1:op:1")
    record = SpaceGroupRecord(
        key=SpaceGroupKey(1, "P 1", ""),
        number=1,
        hm_short="P 1",
        hm_full="P 1",
        point_group="1",
        centering="P",
        crystal_system="triclinic",
        operations=(identity,),
        wyckoff_positions=(),
    )

    definition = record.definition(provenance="derived")

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
            space_group_key=SpaceGroupKey(390, "P -4 2ab", ""),
            letter="a",
            multiplicity=2,
            site_symmetry_symbol="-4..",
            representatives=(AffineCoordinateMap.from_xyz("0,0,0"),),
        )


def test_space_group_key_rejects_dataset_number_outside_catalog() -> None:
    with pytest.raises(ValueError, match="Hall number"):
        SpaceGroupKey(531, "invalid", "")
