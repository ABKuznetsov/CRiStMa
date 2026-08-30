from pathlib import Path

import pytest

from cristma.io.cif.mapper import map_cif_structures
from cristma.io.cif.parser import parse_cif
from cristma.io.result import ReadResult


FIXTURE_DIR = Path(__file__).parents[2] / "fixtures" / "cif"


@pytest.fixture
def read_fixture():
    def read(name: str) -> ReadResult:
        parsed = parse_cif(
            (FIXTURE_DIR / name).read_text(encoding="utf-8"),
            source_name=name,
        )
        structures, mapped = map_cif_structures(parsed.document)
        return ReadResult(
            parsed.document,
            structures,
            parsed.diagnostics + mapped,
            parsed.source_info,
        )

    return read
