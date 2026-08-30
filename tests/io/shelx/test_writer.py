from pathlib import Path

import pytest

import cristma
from cristma.io.shelx import ShelxWriteOptions
from cristma.io.shelx.writer import write_crystal_shelx


SOURCE = (
    "TITL demo in P 1\n"
    "CELL 0.71073 10 11 12 90 91 92\n"
    "LATT -1\n"
    "SFAC C O\n"
    "UNIT 1 0.75\n"
    "FVAR 0.55 0.75\n"
    "C1 1 0.1 0.2 0.3 11 0.05\n"
    "O1 2 0.2 0.3 0.4 21 0.01 0.02 0.03 0.004 0.005 0.006\n"
    "HKLF 4\n"
    "END\n"
)


def test_canonical_writer_requires_measurement_wavelength() -> None:
    crystal = cristma.read_text(SOURCE, format="shelx").structures[0]

    with pytest.raises(ValueError, match="wavelength"):
        write_crystal_shelx(crystal)


def test_canonical_writer_round_trips_scientific_structure_and_fvar() -> None:
    original = cristma.read_text(SOURCE, format="shelx").structures[0]

    text = write_crystal_shelx(
        original,
        options=ShelxWriteOptions(wavelength=0.71073),
    )
    result = cristma.read_text(text, format="shelx")
    restored = result.structures[0]

    assert result.ok
    assert "CELL 0.71073 10 11 12 90 91 92" in text
    assert "FVAR 0.55 0.75" in text
    assert "O1 2 0.2 0.3 0.4 21" in text
    assert "HKLF 4\nEND\n" in text
    assert restored.cell == original.cell
    assert [site.label for site in restored.sites] == ["C1", "O1"]
    assert [site.total_occupancy for site in restored.sites] == pytest.approx([1.0, 0.75])
    assert restored.sites[1].displacement.tensor == original.sites[1].displacement.tensor
    assert [
        (operation.rotation, operation.translation)
        for operation in restored.space_group.operations
    ] == [
        (operation.rotation, operation.translation)
        for operation in original.space_group.operations
    ]


def test_canonical_writer_computes_unit_from_expanded_contents() -> None:
    crystal = cristma.read_text(
        SOURCE.replace("LATT -1", "LATT 1"),
        format="shelx",
    ).structures[0]

    text = write_crystal_shelx(
        crystal,
        options=ShelxWriteOptions(wavelength=0.71073),
    )

    assert "UNIT 2 1.5" in text


def test_top_level_write_dispatches_document_and_explicit_shelx_crystal(
    tmp_path: Path,
) -> None:
    result = cristma.read_text(SOURCE, format="shelx")
    preserved = tmp_path / "copy.res"
    canonical = tmp_path / "canonical.ins"

    cristma.write(result.document, preserved, mode="preserve")
    cristma.write(
        result.structures[0],
        canonical,
        format="shelx",
        mode="canonical",
        options=ShelxWriteOptions(wavelength=0.71073),
    )

    assert preserved.read_text(encoding="utf-8") == SOURCE
    assert cristma.read(canonical).ok


def test_cif_default_write_remains_backward_compatible(tmp_path: Path) -> None:
    crystal = cristma.read_text(SOURCE, format="shelx").structures[0]
    path = tmp_path / "default.cif"

    cristma.write(crystal, path)

    assert path.read_text(encoding="utf-8").startswith("data_demo")
