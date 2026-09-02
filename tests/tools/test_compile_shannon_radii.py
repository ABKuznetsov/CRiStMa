from __future__ import annotations

import json
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


_SCRIPT = Path(__file__).parents[2] / "tools" / "compile_shannon_radii.py"
_SPEC = spec_from_file_location("compile_shannon_radii", _SCRIPT)
assert _SPEC is not None and _SPEC.loader is not None
_MODULE = module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)
compile_csv = _MODULE.compile_csv


def test_compiler_inherits_blank_charge_within_the_same_element(tmp_path) -> None:
    source = tmp_path / "Shannon_Radii.csv"
    source.write_text(
        "Element,Charge,Coordination,Spin State,Crystal Radius,Ionic Radius\n"
        "Eu,3,VI,,1.087,0.947\n"
        "Eu,,VII,,1.15,1.01\n",
        encoding="utf-8",
    )
    destination = tmp_path / "shannon_radii.json"

    compile_csv(source, destination)

    records = json.loads(destination.read_text(encoding="utf-8"))["records"]
    assert [item["oxidation_state"] for item in records] == [3, 3]
