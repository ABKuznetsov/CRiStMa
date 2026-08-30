import pytest

from cristma.io.shelx.parser import parse_shelx
from cristma.io.shelx.sfac import extract_sfac_entries


def test_element_list_sfac_maps_every_label_without_external_chemistry_library() -> None:
    document = parse_shelx("SFAC C H N O S AG\nEND\n").document

    entries = extract_sfac_entries(document.records)

    assert [entry.source_label for entry in entries] == ["C", "H", "N", "O", "S", "AG"]
    assert [entry.species.require_element() for entry in entries] == [
        "C",
        "H",
        "N",
        "O",
        "S",
        "Ag",
    ]


def test_one_scattering_factor_record_per_element_maps_one_entry_each() -> None:
    source = (
        "SFAC C 2.31 20.84 1.02 10.21 1.59 0.57 0.87 51.65 0.22 0 0 0\n"
        "SFAC O 3.05 13.28 2.29 5.70 1.55 0.32 0.87 32.91 0.25 0 0 0\n"
        "END\n"
    )

    entries = extract_sfac_entries(parse_shelx(source).document.records)

    assert [entry.species.require_element() for entry in entries] == ["C", "O"]
    assert all(entry.coefficients for entry in entries)


def test_ionic_and_special_sfac_labels_retain_source_and_resolve_element() -> None:
    entries = extract_sfac_entries(
        parse_shelx("SFAC Ca2+ $C O\nEND\n").document.records
    )

    assert [entry.source_label for entry in entries] == ["Ca2+", "$C", "O"]
    assert [entry.species.require_element() for entry in entries] == ["Ca", "C", "O"]


def test_invalid_sfac_label_has_record_span() -> None:
    document = parse_shelx("SFAC NotAnElement\nEND\n").document

    with pytest.raises(ValueError, match="NotAnElement.*line 1"):
        extract_sfac_entries(document.records)
