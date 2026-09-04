from fractions import Fraction

import pytest

from cristma.symmetry.affine import AffineOperation, format_xyz_operation, parse_xyz_operation


def test_parses_rational_affine_operation_exactly():
    operation = parse_xyz_operation("-x+1/2, y+1/3, z")

    assert operation.rotation[0] == (
        Fraction(-1),
        Fraction(0),
        Fraction(0),
    )
    assert operation.translation == (
        Fraction(1, 2),
        Fraction(1, 3),
        Fraction(0),
    )
    assert operation.apply_fractional((0.1, 0.2, 0.3)) == (
        0.4,
        0.5333333333333333,
        0.3,
    )


def test_operation_normalizes_integer_translation():
    operation = parse_xyz_operation("x+1,y-1,z")

    assert operation.normalized().translation == (Fraction(0),) * 3


def test_parser_rejects_non_linear_expression_without_eval():
    with pytest.raises(ValueError, match="Invalid symmetry expression"):
        parse_xyz_operation("x*y,y,z")


def test_format_xyz_operation_renders_exact_operation_without_source_text() -> None:
    parsed = parse_xyz_operation("-x+1/2,y+1/3,z")
    operation = AffineOperation(parsed.rotation, parsed.translation)

    assert format_xyz_operation(operation) == "-x+1/2,y+1/3,z"
