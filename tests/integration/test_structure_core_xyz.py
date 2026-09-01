from pathlib import Path

import numpy as np

import cristma
from cristma.geometry import CoordinationAnalyzer, NeighborFinder


XYZ_FIXTURES = Path(__file__).parents[1] / "fixtures" / "xyz"
VASP_FIXTURES = Path(__file__).parents[1] / "fixtures" / "vasp"


def _coordination_numbers(structure) -> tuple[int, ...]:
    view = structure.atomic_view()
    graph = NeighborFinder(cutoff=2.01).find(view)
    result = CoordinationAnalyzer().analyze(view, graph)
    return tuple(environment.coordination_number for environment in result.environments)


def test_extxyz_and_poscar_silicon_have_equal_periodic_coordination() -> None:
    xyz = cristma.read(XYZ_FIXTURES / "silicon.extxyz").structures[0]
    poscar = cristma.read(VASP_FIXTURES / "POSCAR").structures[0]

    assert np.allclose(xyz.cell.matrix, poscar.cell.matrix)
    assert np.allclose(xyz.atomic_view().fractional, poscar.atomic_view().fractional)
    assert _coordination_numbers(xyz) == _coordination_numbers(poscar) == (6,)
