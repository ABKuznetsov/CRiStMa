import numpy as np

import cristma
from cristma.geometry import CoordinationAnalyzer, NeighborFinder


SHELX = (
    "TITL salt in P 1\n"
    "CELL 1.5406 5 5 5 90 90 90\n"
    "LATT -1\n"
    "SFAC Na Cl\n"
    "UNIT 1 1\n"
    "NA1 1 0 0 0 11 0.02\n"
    "CL1 2 0.5 0.5 0.5 11 0.02\n"
    "HKLF 4\n"
    "END\n"
)

CIF = """data_salt
_cell_length_a 5
_cell_length_b 5
_cell_length_c 5
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_symop_operation_xyz 'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Na1 Na 0 0 0 1
Cl1 Cl 0.5 0.5 0.5 1
"""


def _signature(crystal: object) -> tuple[object, ...]:
    view = crystal.atomic_view()
    graph = NeighborFinder(cutoff=4.4).find(view)
    coordination = CoordinationAnalyzer().analyze(view, graph)
    return (
        np.asarray(crystal.cell.matrix),
        tuple(atom.components[0].element for atom in view.atoms),
        np.asarray(view.fractional),
        tuple(environment.coordination_number for environment in coordination.environments),
    )


def test_shelx_runs_through_structure_geometry_and_coordination() -> None:
    crystal = cristma.read_text(SHELX, format="shelx").structures[0]

    cell, elements, fractional, coordination = _signature(crystal)

    assert cell.shape == (3, 3)
    assert elements == ("Na", "Cl")
    assert fractional.shape == (2, 3)
    assert coordination == (8, 8)


def test_equivalent_cif_and_shelx_have_equivalent_scientific_geometry() -> None:
    shelx = _signature(cristma.read_text(SHELX, format="shelx").structures[0])
    cif = _signature(cristma.read_text(CIF, format="cif").structures[0])

    assert np.allclose(shelx[0], cif[0])
    assert shelx[1] == cif[1]
    assert np.allclose(shelx[2], cif[2])
    assert shelx[3] == cif[3]
