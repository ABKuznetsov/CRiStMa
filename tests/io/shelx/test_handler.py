from pathlib import Path

import cristma
from cristma.io.shelx.document import ShelxDocument


SOURCE = (
    "TITL demo in P 1\n"
    "CELL 0.71073 10 10 10 90 90 90\n"
    "LATT -1\n"
    "SFAC C\n"
    "C1 1 0.1 0.2 0.3 11 0.05\n"
    "END\n"
)


def test_read_text_auto_detects_shelx_and_returns_canonical_structure() -> None:
    result = cristma.read_text(SOURCE, source_name="memory.res")

    assert result.ok
    assert isinstance(result.document, ShelxDocument)
    assert result.source_info.format == "shelx"
    assert result.structures[0].name == "demo"
    assert result.structures[0].sites[0].label == "C1"


def test_explicit_res_and_ins_aliases_select_shelx() -> None:
    assert cristma.read_text(SOURCE, format="res").source_info.format == "shelx"
    assert cristma.read_text(SOURCE, format="ins").source_info.format == "shelx"


def test_read_path_auto_detects_res_without_application_dispatch(tmp_path: Path) -> None:
    path = tmp_path / "model.res"
    path.write_text(SOURCE, encoding="utf-8")

    result = cristma.read(path)

    assert result.ok
    assert result.source_info.name == str(path)
    assert result.source_info.format == "shelx"


def test_invalid_shelx_keeps_document_and_mapping_diagnostics() -> None:
    result = cristma.read_text("TITL invalid\nEND\n", format="shelx")

    assert isinstance(result.document, ShelxDocument)
    assert not result.ok
    assert not result.structures
    assert "shelx.map.cell_missing" in [item.code for item in result.diagnostics]
