from cristma.io.cif.mapper import map_cif_structures
from cristma.io.cif.parser import parse_cif


def _single_site_cif(occupancy: float) -> str:
    return f"""data_a
_cell_length_a 5
_cell_length_b 5
_cell_length_c 5
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_space_group_symop_operation_xyz
'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Ca1 Ca 0 0 0 {occupancy}
"""


def test_source_occupancy_above_one_rejects_canonical_structure() -> None:
    structures, diagnostics = map_cif_structures(
        parse_cif(_single_site_cif(1.2)).document
    )

    assert not structures
    assert "cif.map.occupancy_out_of_range" in {
        item.code for item in diagnostics
    }


def test_source_negative_occupancy_rejects_canonical_structure() -> None:
    structures, diagnostics = map_cif_structures(
        parse_cif(_single_site_cif(-0.2)).document
    )

    assert not structures
    assert "cif.map.occupancy_out_of_range" in {
        item.code for item in diagnostics
    }


def test_explicit_disorder_total_above_one_rejects_canonical_structure() -> None:
    source = """data_a
_cell_length_a 5
_cell_length_b 5
_cell_length_c 5
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_space_group_symop_operation_xyz
'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_disorder_assembly
_atom_site_disorder_group
Ca1 Ca 0 0 0 0.7 A 1
Sr1 Sr 0 0 0 0.6 A 1
"""

    structures, diagnostics = map_cif_structures(parse_cif(source).document)

    assert not structures
    assert "cif.map.occupancy_total_exceeds_one" in {
        item.code for item in diagnostics
    }


def test_complementary_disorder_rows_form_one_mixed_site(read_fixture):
    result = read_fixture("mixed_disorder.cif")

    assert result.ok
    site = result.structures[0].sites[0]
    assert [(component.element, component.occupancy.value) for component in site.components] == [
        ("La", 0.6),
        ("Zr", 0.4),
    ]
    assert site.disorder_assembly == "A"
    assert site.disorder_group == "1"
    assert site.wyckoff == "1a"


def test_maps_anisotropic_u_tensor_by_atom_label(read_fixture):
    result = read_fixture("anisotropic.cif")

    assert result.ok
    site = result.structures[0].sites[0]
    assert site.displacement.kind == "U_aniso"
    assert site.displacement.tensor[0][1].value == 0.0012
    assert site.displacement.tensor[1][0].value == 0.0012


def test_maps_reported_oxidation_and_checks_orbit_multiplicity(read_fixture):
    result = read_fixture("anisotropic.cif")

    site = result.structures[0].sites[0]
    assert site.components[0].oxidation_state.value == 4
    assert site.reported_multiplicity == 1
    assert site.calculated_multiplicity == 1
    assert "cif.map.multiplicity_mismatch" not in {
        item.code for item in result.diagnostics
    }


def test_coincident_full_sites_are_not_merged():
    source = ("""data_a
_cell_length_a 5
_cell_length_b 5
_cell_length_c 5
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_space_group_symop_operation_xyz
'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
La1 La 0 0 0 1
Zr1 Zr 0 0 0 1
""")

    structures, diagnostics = map_cif_structures(parse_cif(source).document)

    assert len(structures[0].sites) == 2
    assert "cif.map.coincident_sites_unmerged" in {
        item.code for item in diagnostics
    }


def test_reported_multiplicity_mismatch_is_visible(read_fixture):
    source = (
        (read_fixture("anisotropic.cif").document.raw_source)
        .replace("Zr1 Zr 0 0 0 1 4 1a 1", "Zr1 Zr 0 0 0 1 4 2a 2")
    )

    structures, diagnostics = map_cif_structures(parse_cif(source).document)

    assert structures[0].sites[0].reported_multiplicity == 2
    assert structures[0].sites[0].calculated_multiplicity == 1
    assert "cif.map.multiplicity_mismatch" in {item.code for item in diagnostics}


def test_non_positive_anisotropic_tensor_is_retained_with_warning(read_fixture):
    source = read_fixture("anisotropic.cif").document.raw_source.replace(
        "Zr1 0.010 0.011 0.012",
        "Zr1 -0.010 0.011 0.012",
    )

    structures, diagnostics = map_cif_structures(parse_cif(source).document)

    assert structures[0].sites[0].displacement.tensor[0][0].value == -0.01
    assert "cif.map.adp_not_positive_semidefinite" in {
        item.code for item in diagnostics
    }
