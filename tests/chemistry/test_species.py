import pytest

from cristma.chemistry.species import (
    ChargedSpecies,
    ElementSpecies,
    IsotopeSpecies,
    UnknownSpecies,
    as_species,
)


def test_string_becomes_normalized_element_species() -> None:
    species = as_species("si")

    assert species == ElementSpecies("Si")
    assert species.element == "Si"


def test_isotope_and_charge_keep_element_identity() -> None:
    assert IsotopeSpecies("C", 13).element == "C"
    assert ChargedSpecies("Fe", 3).label == "Fe3+"


def test_unknown_species_is_explicit_and_not_an_element() -> None:
    species = UnknownSpecies("species:1", source_label="type 1")

    assert species.element is None
    with pytest.raises(ValueError, match="known element"):
        species.require_element()
