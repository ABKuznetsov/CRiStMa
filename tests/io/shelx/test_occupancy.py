import math

import pytest

from cristma.io.shelx.occupancy import ShelxOccupancyExpression


@pytest.mark.parametrize(
    ("token", "expected"),
    [("11.00000", 1.0), ("10.50000", 0.5), ("1.00000", 1.0)],
)
def test_fixed_occupancy_codes_evaluate_without_fvar(token: str, expected: float) -> None:
    expression = ShelxOccupancyExpression.parse(token)

    assert expression.free_variable_index is None
    assert expression.evaluate(()) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("token", "expected", "complement"),
    [
        ("21.00000", 0.7, False),
        ("20.50000", 0.35, False),
        ("-21.00000", 0.3, True),
        ("-20.50000", 0.15, True),
    ],
)
def test_fvar_codes_retain_dependency_and_evaluate(
    token: str,
    expected: float,
    complement: bool,
) -> None:
    expression = ShelxOccupancyExpression.parse(token)

    assert expression.free_variable_index == 2
    assert expression.multiplier in {0.5, 1.0}
    assert expression.complement is complement
    assert expression.evaluate((0.55, 0.7)) == pytest.approx(expected)


def test_missing_referenced_fvar_is_an_error() -> None:
    expression = ShelxOccupancyExpression.parse("31")

    with pytest.raises(ValueError, match="FVAR 3"):
        expression.evaluate((0.55, 0.7))


@pytest.mark.parametrize("token", ["nan", "inf", "25", "-25", "nonsense"])
def test_invalid_occupancy_expression_is_rejected(token: str) -> None:
    with pytest.raises(ValueError, match="SHELX occupancy"):
        ShelxOccupancyExpression.parse(token)


@pytest.mark.parametrize("fvar", [-0.1, 1.1, math.nan])
def test_unphysical_evaluated_occupancy_is_rejected(fvar: float) -> None:
    expression = ShelxOccupancyExpression.parse("21")

    with pytest.raises(ValueError, match="physical range"):
        expression.evaluate((0.55, fvar))
