import numpy as np
import pytest

from cristma.structure import AtomicProperty, AtomicPropertyTable


def test_property_values_are_immutable_and_typed() -> None:
    prop = AtomicProperty("magnetic_moment", np.array([1.0, -1.0]), unit="mu_B")
    table = AtomicPropertyTable(2, (prop,))

    assert table["magnetic_moment"].unit == "mu_B"
    with pytest.raises(ValueError):
        table["magnetic_moment"].values[0] = 0


def test_property_length_must_match_atoms() -> None:
    with pytest.raises(ValueError, match="leading dimension"):
        AtomicPropertyTable(2, (AtomicProperty("charge", np.array([0.0])),))
