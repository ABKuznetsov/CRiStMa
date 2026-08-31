from cristma.io.vasp import mapper
from cristma.io.vasp.xdatcar import parse_xdatcar


THREE_FRAME_XDATCAR = """Silicon trajectory
1.0
2 0 0
0 2 0
0 0 2
Si
1
Direct configuration=     1
0.0 0.0 0.0
Direct configuration=     2
0.1 0.0 0.0
Direct configuration=     3
0.2 0.0 0.0
"""


def test_xdatcar_indexes_complete_frames_without_mapping_them(monkeypatch) -> None:
    calls = []
    original = mapper.map_vasp_snapshot

    def record(snapshot):
        calls.append(snapshot)
        return original(snapshot)

    monkeypatch.setattr(mapper, "map_vasp_snapshot", record)

    result = parse_xdatcar(THREE_FRAME_XDATCAR, "XDATCAR")

    assert result.ok
    assert len(result.structures) == 3
    assert calls == []
    assert result.structures[-1].name.endswith("configuration 3")
    assert len(calls) == 1
    assert result.structures.final is result.structures[-1]
    assert len(calls) == 1


def test_incomplete_trailing_configuration_is_diagnostic_only() -> None:
    source = THREE_FRAME_XDATCAR + "Direct configuration= 4\n"

    result = parse_xdatcar(source, "XDATCAR")

    assert len(result.structures) == 3
    assert result.structures.references[-1].role == "final"
    assert any(item.code == "vasp.xdatcar.frame_incomplete" for item in result.diagnostics)


def test_reported_configuration_and_source_span_are_preserved() -> None:
    result = parse_xdatcar(THREE_FRAME_XDATCAR, "run/XDATCAR")
    reference = result.structures.references[1]

    assert reference.metadata["configuration"] == 2
    assert reference.source.source_name == "run/XDATCAR"
    assert reference.source.start_offset < reference.source.end_offset
    assert result.structures[1].provenance.source is reference.source


def test_vasp4_xdatcar_species_are_explicitly_unknown() -> None:
    source = """unknown trajectory
1
1 0 0
0 1 0
0 0 1
1
Direct configuration= 1
0 0 0
"""

    result = parse_xdatcar(source, "XDATCAR")

    assert result.ok
    assert result.structures[0].sites[0].components[0].element is None
    assert any(item.code == "vasp.map.species_unresolved" for item in result.diagnostics)
