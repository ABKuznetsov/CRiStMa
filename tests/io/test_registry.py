from pathlib import Path

import pytest

from cristma.io.formats import FormatCapabilities, FormatDescriptor
from cristma.io.registry import FormatRegistry
from cristma.io.result import ReadResult


class StubHandler:
    name = "stub"
    suffixes = (".stub",)

    def probe(self, source: str) -> float:
        return 1.0 if source.startswith("STUB") else 0.0

    def read_text(self, source: str, source_name: str | None = None):
        return (source_name, source)


def test_registry_uses_content_probe_without_suffix(tmp_path: Path):
    path = tmp_path / "unknown.data"
    path.write_text("STUB value", encoding="utf-8")
    registry = FormatRegistry((StubHandler(),))

    assert registry.read(path) == (str(path), "STUB value")


def test_registry_rejects_equal_best_probe_scores(tmp_path: Path):
    path = tmp_path / "unknown.data"
    path.write_text("STUB value", encoding="utf-8")
    first = StubHandler()
    second = StubHandler()
    second.name = "other"
    registry = FormatRegistry((first, second))

    with pytest.raises(ValueError, match="Ambiguous structure format.*other.*stub"):
        registry.read(path)


class ResultHandler(StubHandler):
    def read_text(self, source: str, source_name: str | None = None):
        return ReadResult(document=source)


def test_registry_reports_latin1_decoding_fallback(tmp_path: Path):
    path = tmp_path / "sample.stub"
    path.write_bytes(b"STUB caf\xe9\r\n")
    registry = FormatRegistry((ResultHandler(),))

    result = registry.read(path)

    assert result.source_info.encoding == "latin-1"
    assert result.source_info.newline == "\r\n"
    assert "io.encoding_fallback" in {item.code for item in result.diagnostics}


def test_special_basename_beats_uninformative_suffix_without_loading_handler():
    loaded = []
    descriptor = FormatDescriptor(
        name="vasp",
        aliases=("poscar",),
        suffixes=(),
        basenames=("POSCAR", "CONTCAR"),
        probe=lambda source: 0.0,
        factory=lambda: loaded.append(True),
        capabilities=FormatCapabilities(text=True),
    )
    registry = FormatRegistry((descriptor,))

    selected = registry.select("title", basename="POSCAR")

    assert selected.name == "vasp"
    assert loaded == []


def test_handler_factory_runs_only_when_selected_reader_is_used():
    loaded = []
    handler = StubHandler()
    descriptor = FormatDescriptor(
        name="stub",
        aliases=(),
        suffixes=(".stub",),
        basenames=(),
        probe=handler.probe,
        factory=lambda: loaded.append(handler) or handler,
        capabilities=FormatCapabilities(text=True),
    )
    registry = FormatRegistry((descriptor,))

    assert loaded == []
    assert registry.read_text("STUB value") == (None, "STUB value")
    assert loaded == [handler]
