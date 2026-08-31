from pathlib import Path

import numpy as np

import cristma
from cristma.geometry import CoordinationAnalyzer, NeighborFinder


FIXTURES = Path(__file__).parents[1] / "fixtures" / "vasp"
CIF_EQUIVALENT = """data_silicon
_cell_length_a 2
_cell_length_b 2
_cell_length_c 2
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0 0 0
"""


def _coordination_numbers(structure) -> tuple[int, ...]:
    view = structure.atomic_view()
    graph = NeighborFinder(cutoff=2.01).find(view)
    result = CoordinationAnalyzer().analyze(view, graph)
    return tuple(environment.coordination_number for environment in result.environments)


def test_equivalent_cif_and_poscar_have_equal_geometry_and_coordination() -> None:
    cif = cristma.read_text(CIF_EQUIVALENT, format="cif").structures[0]
    poscar = cristma.read(FIXTURES / "POSCAR").structures[0]

    assert np.allclose(cif.cell.matrix, poscar.cell.matrix, atol=1e-12)
    assert np.allclose(cif.atomic_view().fractional, poscar.atomic_view().fractional, atol=1e-12)
    assert _coordination_numbers(cif) == _coordination_numbers(poscar) == (6,)


def test_all_trajectory_families_produce_same_final_cell_and_typed_forces() -> None:
    outcar = cristma.read(FIXTURES / "OUTCAR").structures.final
    xml = cristma.read(FIXTURES / "vasprun.xml").structures.final

    assert np.allclose(outcar.cell.matrix, xml.cell.matrix, atol=1e-12)
    assert outcar.properties["force"].unit == xml.properties["force"].unit == "eV/angstrom"
