from pathlib import Path

import pytest

import cristma


FIXTURES = Path(__file__).parents[2] / "fixtures" / "vasp"


@pytest.mark.parametrize(
    ("name", "expected_frames"),
    [("POSCAR", 1), ("XDATCAR", 3), ("OUTCAR", 2), ("vasprun.xml", 2)],
)
def test_reference_vasp_fixture_reads_through_public_api(name: str, expected_frames: int) -> None:
    result = cristma.read(FIXTURES / name)

    assert result.ok
    assert len(result.structures) == expected_frames
    assert result.structures[-1].atomic_view().cartesian.shape == (1, 3)
