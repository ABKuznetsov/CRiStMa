from __future__ import annotations

import json
from pathlib import Path

import pytest
import spglib

from tools.compile_spglib_crystallography import compile_catalog


SPGLIB_COMMIT = "12355c77fb7c505a55f52cae36341d73b781a065"
FIXTURES = Path(__file__).parents[1] / "fixtures" / "spglib"
SPG = FIXTURES / "spg_minimal.csv"
WYCKOFF = FIXTURES / "Wyckoff_minimal.csv"


def test_compiler_emits_canonical_catalogs(tmp_path: Path) -> None:
    space_groups, wyckoffs = compile_catalog(
        SPG,
        WYCKOFF,
        tmp_path,
        upstream_commit=SPGLIB_COMMIT,
        compiled_date="2026-09-01",
    )

    groups = json.loads(space_groups.read_text(encoding="utf-8"))
    positions = json.loads(wyckoffs.read_text(encoding="utf-8"))

    assert [record["hall_number"] for record in groups["records"]] == [1, 2, 390]
    assert groups["records"][2]["hall_symbol"] == "P -4 2ab"
    assert len(groups["records"][2]["operations"]) == 8
    assert positions["records"]["390"][0]["letter"] == "f"
    assert positions["records"]["390"][-1]["letter"] == "a"
    assert positions["records"]["390"][-1]["multiplicity"] == 2


def test_compiler_output_is_byte_deterministic(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    kwargs = {
        "upstream_commit": SPGLIB_COMMIT,
        "compiled_date": "2026-09-01",
    }

    compile_catalog(SPG, WYCKOFF, first, **kwargs)
    compile_catalog(SPG, WYCKOFF, second, **kwargs)

    assert (first / "space_groups.json").read_bytes() == (
        second / "space_groups.json"
    ).read_bytes()
    assert (first / "wyckoff_positions.json").read_bytes() == (
        second / "wyckoff_positions.json"
    ).read_bytes()


def test_compiler_rejects_wrong_spglib_version(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(spglib, "__version__", "2.6.0")

    with pytest.raises(RuntimeError, match="spglib 2.7.0 is required"):
        compile_catalog(
            SPG,
            WYCKOFF,
            tmp_path,
            upstream_commit=SPGLIB_COMMIT,
            compiled_date="2026-09-01",
        )
