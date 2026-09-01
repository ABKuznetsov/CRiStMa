from __future__ import annotations

import cristma


def _p421m_cif(*, operations: str = "", site_tail: str = "") -> str:
    return f"""data_demo
_cell_length_a 5
_cell_length_b 5
_cell_length_c 7
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_Hall 'P -4 2ab'
{operations}
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
{site_tail or 'Si1 Si 0 0 0'}
"""


def test_cif_without_operations_resolves_exact_hall_catalog_entry() -> None:
    result = cristma.read_text(_p421m_cif(), format="cif")

    crystal = result.structures[0]
    codes = {diagnostic.code for diagnostic in result.diagnostics}
    assert crystal.space_group is not None
    assert crystal.space_group.provenance == "derived"
    assert crystal.space_group.number == 113
    assert len(crystal.space_group.operations) == 8
    assert "cif.map.symmetry_operations_derived" in codes
    assert "cif.map.symmetry_operations_missing" not in codes


def test_explicit_operations_remain_authoritative_when_catalog_disagrees() -> None:
    source = _p421m_cif(
        operations="""loop_
_space_group_symop_operation_xyz
'x,y,z'
""",
    )

    result = cristma.read_text(source, format="cif")

    crystal = result.structures[0]
    assert crystal.space_group is not None
    assert crystal.space_group.provenance == "reported"
    assert len(crystal.space_group.operations) == 1
    assert "cif.map.space_group_operations_mismatch" in {
        diagnostic.code for diagnostic in result.diagnostics
    }


def test_ambiguous_hm_symbol_does_not_guess_a_setting() -> None:
    source = _p421m_cif().replace(
        "_space_group_name_Hall 'P -4 2ab'",
        "_space_group_name_H-M_alt 'C 2'",
    )

    result = cristma.read_text(source, format="cif")

    assert result.structures[0].space_group.provenance == "identity_fallback"
    assert "cif.map.space_group_lookup_ambiguous" in {
        diagnostic.code for diagnostic in result.diagnostics
    }


def test_cif_site_is_validated_against_calculated_wyckoff_position() -> None:
    source = _p421m_cif(
        site_tail="""loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_wyckoff_symbol
_atom_site_symmetry_multiplicity
Si1 Si 0 0 0 d 4""",
    ).replace(
        "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nloop_\n",
        "loop_\n",
    )

    result = cristma.read_text(source, format="cif")
    codes = {diagnostic.code for diagnostic in result.diagnostics}

    assert result.structures[0].sites[0].calculated_multiplicity == 2
    assert "crystallography.orbit.reported_wyckoff_mismatch" in codes
    assert "crystallography.orbit.reported_multiplicity_mismatch" in codes
