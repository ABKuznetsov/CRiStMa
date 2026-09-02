from __future__ import annotations

import pytest

from cristma.reference_data import CovalentRadii


def test_covalent_radius_exposes_value_units_and_provenance() -> None:
    radius = CovalentRadii.default().find("O")

    assert radius.value == pytest.approx(0.66)
    assert radius.unit == "angstrom"
    assert radius.dataset_id == "cristma.covalent_radii.cordero_2008"


def test_covalent_radius_normalizes_element_symbol() -> None:
    # Use the largest published variant when bonding/spin is unknown so the
    # geometric candidate search cannot miss a physically possible contact.
    assert CovalentRadii.default().find("fe").value == pytest.approx(1.52)


def test_covalent_radius_never_guesses_a_missing_value() -> None:
    with pytest.raises(KeyError, match="No covalent radius"):
        CovalentRadii.default().find("Cf")


def test_covalent_catalog_covers_heavy_chalcogenide_pair() -> None:
    radii = CovalentRadii.default()

    assert radii.find("Bi").value == pytest.approx(1.48)
    assert radii.find("Te").value == pytest.approx(1.38)


def test_covalent_catalog_contains_every_cordero_element_without_fill_values() -> None:
    radii = CovalentRadii.default()

    assert len(radii.records) == 101
    assert len({record.symbol for record in radii.records}) == 96
    assert radii.find("H").value == pytest.approx(0.31)
    assert radii.find("He").value == pytest.approx(0.28)
    assert radii.find("Cm").value == pytest.approx(1.69)


def test_covalent_catalog_preserves_published_variants() -> None:
    radii = CovalentRadii.default()

    assert [item.value for item in radii.find_variants("C")] == pytest.approx(
        [0.76, 0.73, 0.69]
    )
    assert [item.value for item in radii.find_variants("Fe")] == pytest.approx(
        [1.32, 1.52]
    )
