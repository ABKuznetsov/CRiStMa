from pathlib import Path

import numpy as np
import pytest

import cristma
from cristma.structure import CrystalStructure, MolecularStructure


FIXTURES = Path(__file__).parents[2] / "fixtures" / "xyz"


@pytest.mark.parametrize(
    ("name", "frames"),
    [
        ("water.xyz", 1),
        ("water.extxyz", 1),
        ("silicon.extxyz", 1),
        ("trajectory.xyz", 3),
    ],
)
def test_xyz_fixture_reads_through_public_api(name: str, frames: int) -> None:
    result = cristma.read(FIXTURES / name)

    assert result.ok
    assert len(result.structures) == frames


def test_plain_and_extended_water_have_equal_geometry() -> None:
    plain = cristma.read(FIXTURES / "water.xyz").structures[0]
    extended = cristma.read(FIXTURES / "water.extxyz").structures[0]

    assert isinstance(plain, MolecularStructure)
    assert isinstance(extended, MolecularStructure)
    assert np.allclose(plain.atomic_view().cartesian, extended.atomic_view().cartesian)
    assert extended.properties["charge"].values.tolist() == [-0.84, 0.42, 0.42]


def test_mixed_trajectory_frames_keep_independent_scientific_types() -> None:
    structures = cristma.read(FIXTURES / "trajectory.xyz").structures

    assert isinstance(structures[0], MolecularStructure)
    assert isinstance(structures[1], MolecularStructure)
    assert isinstance(structures[2], CrystalStructure)
    assert structures[1].properties["step"].values.tolist() == [2]
