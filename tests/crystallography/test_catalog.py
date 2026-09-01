from __future__ import annotations

import pytest

from cristma.crystallography import SpaceGroupCatalog


def test_default_catalog_has_all_hall_settings() -> None:
    catalog = SpaceGroupCatalog.default()

    assert len(catalog) == 530
    assert {setting.number for setting in catalog.settings} == set(range(1, 231))
    assert catalog.by_setting(390).number == 113
    assert catalog.by_setting(390).hall_symbol == "P -4 2ab"


def test_number_lookup_preserves_setting_ambiguity() -> None:
    records = SpaceGroupCatalog.default().by_number(5)

    assert len(records) == 9
    assert len({setting.setting_id for setting in records}) == len(records)


def test_hall_symbol_lookup_normalizes_whitespace_and_case() -> None:
    setting = SpaceGroupCatalog.default().by_hall("  p   -4   2AB ")

    assert setting.setting_id == 390


def test_duplicated_hall_symbol_requires_explicit_hall_number() -> None:
    with pytest.raises(LookupError, match="ambiguous Hall symbol"):
        SpaceGroupCatalog.default().by_hall("C 2 2 -1ac")


def test_default_catalog_is_cached() -> None:
    assert SpaceGroupCatalog.default() is SpaceGroupCatalog.default()


def test_catalog_rejects_unknown_lookup_values() -> None:
    catalog = SpaceGroupCatalog.default()

    with pytest.raises(KeyError):
        catalog.by_setting(531)
    with pytest.raises(KeyError):
        catalog.by_hall("not a Hall symbol")


def test_catalog_exposes_wyckoff_positions_without_new_scientific_state() -> None:
    catalog = SpaceGroupCatalog.default()

    assert catalog.wyckoff_positions(390) is catalog.by_setting(390).wyckoff_positions
    assert [position.letter for position in catalog.wyckoff_positions(390)] == [
        "f",
        "e",
        "d",
        "c",
        "b",
        "a",
    ]
