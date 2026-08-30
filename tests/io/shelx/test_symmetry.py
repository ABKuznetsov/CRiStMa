from fractions import Fraction

import pytest

from cristma.io.shelx.symmetry import build_shelx_operations
from cristma.symmetry.affine import parse_xyz_operation


@pytest.mark.parametrize(
    ("latt", "translations"),
    [
        (-1, {(0, 0, 0)}),
        (-2, {(0, 0, 0), (Fraction(1, 2),) * 3}),
        (
            -3,
            {
                (0, 0, 0),
                (Fraction(2, 3), Fraction(1, 3), Fraction(1, 3)),
                (Fraction(1, 3), Fraction(2, 3), Fraction(2, 3)),
            },
        ),
        (
            -4,
            {
                (0, 0, 0),
                (0, Fraction(1, 2), Fraction(1, 2)),
                (Fraction(1, 2), 0, Fraction(1, 2)),
                (Fraction(1, 2), Fraction(1, 2), 0),
            },
        ),
        (-5, {(0, 0, 0), (0, Fraction(1, 2), Fraction(1, 2))}),
        (-6, {(0, 0, 0), (Fraction(1, 2), 0, Fraction(1, 2))}),
        (-7, {(0, 0, 0), (Fraction(1, 2), Fraction(1, 2), 0)}),
    ],
)
def test_latt_generates_exact_centering_translations(latt: int, translations: set) -> None:
    operations = build_shelx_operations(latt, ())

    assert {operation.translation for operation in operations} == translations
    assert operations[0].translation == (0, 0, 0)


def test_positive_latt_adds_inversion_and_negative_does_not() -> None:
    non_centrosymmetric = build_shelx_operations(-1, ())
    centrosymmetric = build_shelx_operations(1, ())

    assert len(non_centrosymmetric) == 1
    assert len(centrosymmetric) == 2
    assert centrosymmetric[1].rotation == (
        (Fraction(-1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(-1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(-1)),
    )


def test_explicit_operations_are_combined_and_deduplicated_exactly() -> None:
    operation = parse_xyz_operation("-x+1/2,y+1/2,-z+1/2")

    operations = build_shelx_operations(1, (operation, operation))

    assert len(operations) == 4
    assert [item.id for item in operations] == [
        "shelx:op:1",
        "shelx:op:2",
        "shelx:op:3",
        "shelx:op:4",
    ]


def test_p21n_example_produces_four_operations() -> None:
    reported = parse_xyz_operation("1/2-x,1/2+y,1/2-z")

    operations = build_shelx_operations(1, (reported,))

    assert len(operations) == 4
    assert len({(item.rotation, item.translation) for item in operations}) == 4
    assert all(item.source for item in operations)
    assert all(
        parse_xyz_operation(item.source).normalized().rotation == item.rotation
        and parse_xyz_operation(item.source).normalized().translation == item.translation
        for item in operations
    )


@pytest.mark.parametrize("latt", [0, 8, -8])
def test_unknown_latt_code_is_rejected(latt: int) -> None:
    with pytest.raises(ValueError, match="SHELX LATT"):
        build_shelx_operations(latt, ())
