import pytest

from cristma.chemistry import element_from_atomic_number


def test_atomic_number_maps_to_iupac_symbol() -> None:
    assert element_from_atomic_number(1) == "H"
    assert element_from_atomic_number(14) == "Si"
    assert element_from_atomic_number(118) == "Og"


@pytest.mark.parametrize("value", [0, 119, -1])
def test_out_of_range_atomic_number_is_rejected(value: int) -> None:
    with pytest.raises(ValueError, match="between 1 and 118"):
        element_from_atomic_number(value)


def test_boolean_is_not_an_atomic_number() -> None:
    with pytest.raises(TypeError, match="integer"):
        element_from_atomic_number(True)
