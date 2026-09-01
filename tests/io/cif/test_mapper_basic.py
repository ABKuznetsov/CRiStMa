import pytest

from cristma.io.cif.mapper import map_cif_structures
from cristma.io.cif.parser import parse_cif
from cristma.symmetry import expand_structure


MINIMAL = """data_si
_cell_length_a 5.43
_cell_length_b 5.43
_cell_length_c 5.43
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt 'P -1'
loop_
_space_group_symop_operation_xyz
'x,y,z'
'-x,-y,-z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Si1 Si 0 0 0 1
"""


def test_maps_asymmetric_site_and_exact_reported_symmetry():
    document = parse_cif(MINIMAL).document

    structures, diagnostics = map_cif_structures(document)

    assert not [item for item in diagnostics if item.severity.value == "error"]
    crystal = structures[0]
    assert crystal.cell.volume == pytest.approx(5.43**3)
    assert crystal.sites[0].components[0].element == "Si"
    assert crystal.space_group.provenance == "reported"
    assert len(crystal.space_group.operations) == 2
    assert len(expand_structure(crystal).atoms) == 1


def test_missing_cell_does_not_fabricate_structure():
    document = parse_cif(
        "data_a\n_atom_site_label Si1\n_atom_site_fract_x 0\n"
    ).document

    structures, diagnostics = map_cif_structures(document)

    assert structures == ()
    assert {item.code for item in diagnostics} >= {"cif.map.cell_missing"}


def test_incomplete_coordinate_row_is_error():
    source = MINIMAL.replace("Si1 Si 0 0 0 1", "Si1 Si 0 ? 0 1")

    structures, diagnostics = map_cif_structures(parse_cif(source).document)

    assert structures == ()
    assert "cif.map.coordinate_missing" in {item.code for item in diagnostics}


def test_missing_operations_are_derived_from_unique_hm_symbol():
    source = MINIMAL.replace(
        "loop_\n_space_group_symop_operation_xyz\n'x,y,z'\n'-x,-y,-z'\n",
        "",
    )

    structures, diagnostics = map_cif_structures(parse_cif(source).document)

    assert structures[0].space_group.provenance == "derived"
    assert "cif.map.symmetry_operations_derived" in {item.code for item in diagnostics}


def test_maps_formula_hall_symbol_and_publication_metadata():
    source = MINIMAL.replace(
        "_space_group_name_H-M_alt 'P -1'",
        """_space_group_name_H-M_alt 'P -1'
_space_group_name_Hall '-P 1'
_space_group_IT_number 2
_chemical_formula_sum 'Si'
_journal_name_full 'Journal of Tests'
_journal_year 2026
_journal_paper_doi '10.1/example'""",
    )

    structures, diagnostics = map_cif_structures(parse_cif(source).document)

    assert not [item for item in diagnostics if item.severity.value == "error"]
    crystal = structures[0]
    assert crystal.formula == "Si"
    assert crystal.space_group.hall_symbol == "-P 1"
    assert crystal.space_group.number == 2
    assert crystal.metadata["journal"] == "Journal of Tests"
    assert crystal.metadata["doi"] == "10.1/example"
