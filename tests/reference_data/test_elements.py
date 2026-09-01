from __future__ import annotations

import pytest

from cristma.reference_data import ElementCatalog, ElementCategory


def test_element_catalog_identifies_metal_metalloid_and_nonmetal() -> None:
    catalog = ElementCatalog.default()

    assert catalog.by_symbol("Fe").category is ElementCategory.METAL
    assert catalog.by_symbol("Fe").is_metal
    assert catalog.by_symbol("si").category is ElementCategory.METALLOID
    assert not catalog.by_symbol("O").is_metal


def test_element_catalog_normalizes_symbol_and_preserves_atomic_number() -> None:
    record = ElementCatalog.default().by_symbol("fe")

    assert record.symbol == "Fe"
    assert record.atomic_number == 26
    assert record.dataset_id == "cristma.elements"


def test_element_catalog_rejects_unknown_symbol() -> None:
    with pytest.raises(ValueError, match="Unknown chemical element"):
        ElementCatalog.default().by_symbol("Xx")
