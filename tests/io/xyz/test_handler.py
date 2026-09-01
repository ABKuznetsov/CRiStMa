from pathlib import Path

import cristma
from cristma.io.xyz import mapper
from cristma.io.xyz.parser import parse_xyz
from cristma.structure import MolecularStructure


THREE_FRAME_XYZ = """1
frame 1
H 0 0 0
1
frame 2
H 1 0 0
1
frame 3
H 2 0 0
"""


def test_multiframe_xyz_maps_only_requested_frame(monkeypatch) -> None:
    calls = []
    original = mapper.map_xyz_frame

    def record(frame):
        calls.append(frame)
        return original(frame)

    monkeypatch.setattr(mapper, "map_xyz_frame", record)

    result = parse_xyz(THREE_FRAME_XYZ, "trajectory.xyz")

    assert len(result.structures) == 3
    assert calls == []
    assert result.structures.final.name == "frame 3"
    assert len(calls) == 1
    assert result.structures.final is result.structures[-1]
    assert len(calls) == 1


def test_references_preserve_roles_spans_and_changing_schemas() -> None:
    source = (
        "1\none\nH 0 0 0\n"
        "1\nProperties=Z:I:1:pos:R:3:charge:R:1\n1 1 0 0 0.2\n"
    )

    result = parse_xyz(source, "changing.xyz")

    first, final = result.structures.references
    assert first.role == "intermediate"
    assert final.role == "final"
    assert first.source.end_offset <= final.source.start_offset
    assert result.structures[1].properties["charge"].values.tolist() == [0.2]


def test_complete_frames_survive_incomplete_tail() -> None:
    result = parse_xyz(THREE_FRAME_XYZ + "2\nbroken\nH 0 0 0\n", "tail.xyz")

    assert len(result.structures) == 3
    assert any(item.code == "xyz.frame.incomplete" for item in result.diagnostics)


def test_public_read_maps_xyz(tmp_path: Path) -> None:
    path = tmp_path / "molecule.xyz"
    path.write_text("1\nhelium\nHe 0 0 0\n", encoding="utf-8")

    result = cristma.read(path)

    assert result.ok
    assert isinstance(result.structures[0], MolecularStructure)
    assert result.source_info.format == "xyz"

