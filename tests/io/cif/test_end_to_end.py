from pathlib import Path

import pytest

import cristma
from cristma.core.structure import Crystal
from cristma.io.cif.document import CifDocument


FIXTURES = Path("tests/fixtures/cif")


def test_public_read_maps_real_inorganic_cif():
    result = cristma.read(FIXTURES / "lithium_triborate.cif")

    assert result.ok
    assert result.structures
    assert {
        component.element
        for site in result.structures[0].sites
        for component in site.components
    } >= {"Li", "B", "O"}
    assert result.source_info.format == "cif"


def test_public_read_text_probes_cif_and_maps_structure():
    source = (FIXTURES / "lithium_triborate.cif").read_text(encoding="utf-8")

    result = cristma.read_text(source, source_name="memory.cif")

    assert result.ok
    assert result.structures[0].formula == "B3 Li O5"
    assert result.source_info.name == "memory.cif"
    assert result.source_info.format == "cif"


def test_public_preserve_round_trip_keeps_mixed_occupancy_fixture(tmp_path):
    source = FIXTURES / "mixed_occupancy_positions.cif"
    result = cristma.read(source)
    target = tmp_path / "copy.cif"

    cristma.write(result.document, target, mode="preserve")

    assert target.read_bytes() == source.read_bytes()


def test_public_canonical_write_can_be_read_back(tmp_path):
    result = cristma.read(FIXTURES / "lithium_triborate.cif")
    original = result.structures[0]
    target = tmp_path / "canonical.cif"

    cristma.write(original, target, mode="canonical")
    restored = cristma.read(target).structures[0]

    assert [site.label for site in restored.sites] == [site.label for site in original.sites]
    assert restored.cell.a.value == original.cell.a.value


@pytest.mark.parametrize(
    ("value_type", "mode"),
    ((CifDocument, "canonical"), (Crystal, "preserve")),
)
def test_public_write_rejects_incompatible_mode(value_type, mode, tmp_path):
    result = cristma.read(FIXTURES / "lithium_triborate.cif")
    value = result.document if value_type is CifDocument else result.structures[0]

    with pytest.raises(ValueError, match="mode"):
        cristma.write(value, tmp_path / "invalid.cif", mode=mode)
