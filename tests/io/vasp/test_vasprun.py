import numpy as np
import pytest

from cristma.io.vasp import mapper
from cristma.io.vasp.vasprun import parse_vasprun


def _calculation(x: float) -> str:
    return f"""<calculation>
  <structure><crystal><varray name="basis"><v>2 0 0</v><v>0 2 0</v><v>0 0 2</v></varray></crystal>
    <varray name="positions"><v>{x} 0 0</v><v>0.5 0.5 0.5</v></varray></structure>
  <varray name="forces"><v>1 0 0</v><v>-1 0 0</v></varray>
</calculation>"""


THREE_STEP_XML = """<?xml version="1.0"?>
<modeling>
<atominfo><array name="atoms"><set>
  <rc><c>Na</c><c>1</c></rc><rc><c>Cl</c><c>2</c></rc>
</set></array></atominfo>
""" + _calculation(0.0) + _calculation(0.1) + _calculation(0.2) + "</modeling>"


def test_vasprun_indexes_ionic_steps_and_maps_only_on_access(monkeypatch) -> None:
    calls = []
    original = mapper.map_vasp_snapshot

    def record(snapshot):
        calls.append(snapshot)
        return original(snapshot)

    monkeypatch.setattr(mapper, "map_vasp_snapshot", record)
    result = parse_vasprun(THREE_STEP_XML, "vasprun.xml")

    assert result.ok
    assert len(result.structures) == 3
    assert calls == []
    assert result.structures.references[-1].role == "final"
    final = result.structures.final
    assert len(calls) == 1
    assert final.properties["force"].values.shape == (2, 3)
    assert result.structures.final is final
    assert len(calls) == 1


def test_truncated_xml_keeps_only_complete_calculations() -> None:
    source = THREE_STEP_XML.replace(_calculation(0.1) + _calculation(0.2) + "</modeling>", "<calculation><structure>")

    result = parse_vasprun(source, "vasprun.xml")

    assert len(result.structures) == 1
    assert any(item.code == "vasp.vasprun.xml_incomplete" for item in result.diagnostics)


def test_default_namespace_is_accepted() -> None:
    source = THREE_STEP_XML.replace("<modeling>", '<modeling xmlns="urn:vasp">', 1)

    result = parse_vasprun(source, "vasprun.xml")

    assert len(result.structures) == 3
    assert [site.components[0].element for site in result.structures.final.sites] == ["Na", "Cl"]


def test_initial_only_single_point_structure_is_available() -> None:
    source = """<modeling><atominfo><array name="atoms"><set><rc><c>Si</c></rc></set></array></atominfo>
<structure name="initialpos"><crystal><varray name="basis"><v>1 0 0</v><v>0 1 0</v><v>0 0 1</v></varray></crystal>
<varray name="positions"><v>0 0 0</v></varray></structure></modeling>"""

    result = parse_vasprun(source, "vasprun.xml")

    assert result.ok
    assert len(result.structures) == 1
    assert result.structures.final.sites[0].components[0].element == "Si"


def test_atom_count_mismatch_is_rejected_when_frame_is_loaded() -> None:
    source = THREE_STEP_XML.replace("<v>0.5 0.5 0.5</v>", "", 1)
    result = parse_vasprun(source, "vasprun.xml")

    with pytest.raises(ValueError, match="position count"):
        result.structures[0]


def test_malformed_force_array_is_not_guessed() -> None:
    source = THREE_STEP_XML.replace("<v>1 0 0</v>", "<v>not-a-force</v>", 1)
    result = parse_vasprun(source, "vasprun.xml")

    with pytest.raises(ValueError, match="forces"):
        result.structures[0]


def test_unicode_before_frame_keeps_character_source_offsets() -> None:
    source = THREE_STEP_XML.replace("<modeling>", "<modeling><!-- расчёт -->", 1)
    result = parse_vasprun(source, "vasprun.xml")
    reference = result.structures.references[0]

    fragment = source[reference.source.start_offset : reference.source.end_offset]
    assert fragment.startswith("<calculation>")
    assert fragment.endswith("</calculation>")
    assert np.allclose(result.structures[0].cell.matrix, np.diag([2.0, 2.0, 2.0]), atol=1e-12)
