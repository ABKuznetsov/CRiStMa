from __future__ import annotations

import pytest

from cristma.reference_data import CovalentRadii


def test_covalent_radius_exposes_value_units_and_provenance() -> None:
    radius = CovalentRadii.default().find("O")

    assert radius.value == pytest.approx(0.66)
    assert radius.unit == "angstrom"
    assert radius.dataset_id == "cristma.covalent_radii.craft"


def test_covalent_radius_normalizes_element_symbol() -> None:
    assert CovalentRadii.default().find("fe").value == pytest.approx(1.32)


def test_covalent_radius_never_guesses_a_missing_value() -> None:
    with pytest.raises(KeyError, match="No covalent radius"):
        CovalentRadii.default().find("He")
