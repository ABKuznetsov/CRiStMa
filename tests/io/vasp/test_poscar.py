import numpy as np
import pytest

from cristma.io.vasp.mapper import map_vasp_snapshot
from cristma.io.vasp.poscar import parse_poscar, poscar_snapshot


POSCAR = """Silicon oxide
2.0
1 0 0
0 1 0
0 0 1
Si O
1 2
Selective dynamics
Cartesian
0.0 0.0 0.0 T F T Si
0.5 0.5 0.5 F F F O
0.25 0.25 0.25 T T T O
Cartesian
0.01 0.02 0.03
0.04 0.05 0.06
0.07 0.08 0.09
predictor-corrector data stays untouched
"""


def test_vasp5_cartesian_selective_and_velocity_sections_are_preserved() -> None:
    result = parse_poscar(POSCAR, "CONTCAR")
    document = result.document

    assert result.ok
    assert document.header.species_labels == ("Si", "O")
    assert document.header.counts == (1, 2)
    assert document.header.coordinate_mode == "cartesian"
    assert document.positions[0].selective == (True, False, True)
    assert document.velocity_mode == "cartesian"
    assert len(document.velocities) == 3
    assert document.render_preserved() == POSCAR
    assert document.trailing_start is not None


def test_cartesian_positions_scale_but_velocities_do_not() -> None:
    document = parse_poscar(POSCAR, "CONTCAR").document
    snapshot = poscar_snapshot(document)

    assert snapshot.lattice.tolist() == np.diag([2.0, 2.0, 2.0]).tolist()
    assert snapshot.fractional[1].tolist() == [0.5, 0.5, 0.5]
    assert snapshot.velocities[0].tolist() == [0.01, 0.02, 0.03]
    assert snapshot.velocity_mode == "cartesian"
    assert snapshot.velocity_unit == "angstrom/fs"


def test_vasp4_species_remain_explicitly_unknown() -> None:
    source = """mystery material
1
2 0 0
0 2 0
0 0 2
2
Direct
0 0 0
0.5 0.5 0.5
"""

    result = parse_poscar(source, "POSCAR")
    crystal = map_vasp_snapshot(poscar_snapshot(result.document))

    assert result.ok
    assert [site.components[0].element for site in crystal.sites] == [None, None]
    assert any(item.code == "vasp.map.species_unresolved" for item in result.diagnostics)


def test_direct_velocities_keep_reported_direct_convention() -> None:
    source = """demo
1
1 0 0
0 1 0
0 0 1
Si
1
Direct
0 0 0
Direct
0.1 0.2 0.3
"""

    snapshot = poscar_snapshot(parse_poscar(source, "POSCAR").document)

    assert snapshot.velocities.tolist() == [[0.1, 0.2, 0.3]]
    assert snapshot.velocity_mode == "direct"
    assert snapshot.velocity_unit == "direct_lattice_vector/timestep"


def test_incomplete_declared_position_block_is_an_error() -> None:
    source = """broken
1
1 0 0
0 1 0
0 0 1
Si
2
Direct
0 0 0
"""

    result = parse_poscar(source, "POSCAR")

    assert not result.ok
    assert any(item.code == "vasp.poscar.positions_incomplete" for item in result.diagnostics)
    with pytest.raises(ValueError, match="position count"):
        poscar_snapshot(result.document)


def test_invalid_selective_flag_is_not_silently_accepted() -> None:
    source = POSCAR.replace("T F T Si", "T MAYBE T Si")

    result = parse_poscar(source, "POSCAR")

    assert not result.ok
    assert any(item.code == "vasp.poscar.selective_flag_invalid" for item in result.diagnostics)


def test_zero_count_species_creates_no_sites() -> None:
    source = """vacant type
1
1 0 0
0 1 0
0 0 1
Na Cl
0 1
Direct
0.5 0.5 0.5
"""

    crystal = map_vasp_snapshot(poscar_snapshot(parse_poscar(source).document))

    assert len(crystal.sites) == 1
    assert crystal.sites[0].components[0].element == "Cl"
