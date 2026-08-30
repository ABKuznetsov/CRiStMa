from cristma.io.cif.document import replace_scalar
from cristma.io.cif.mapper import map_cif_structures
from cristma.io.cif.parser import parse_cif
from cristma.io.cif.writer import write_cif_document, write_crystal_cif


def test_unchanged_document_round_trips_byte_for_byte():
    source = (
        "data_a\r\n"
        "# keep me\r\n"
        "_local_unknown 'A B'\r\n"
        "_cell_length_a 5.0\r\n"
    )
    document = parse_cif(source).document

    assert write_cif_document(document, mode="preserve") == source


def test_scalar_edit_preserves_unknown_content_and_newlines():
    source = (
        "data_a\r\n"
        "# keep me\r\n"
        "_local_unknown 'A B'\r\n"
        "_cell_length_a 5.0\r\n"
    )
    document = replace_scalar(
        parse_cif(source).document,
        "a",
        "_cell_length_a",
        "5.1(2)",
    )

    rendered = write_cif_document(document, mode="preserve")

    assert "# keep me\r\n_local_unknown 'A B'" in rendered
    assert "_cell_length_a 5.1(2)" in rendered


def test_canonical_writer_emits_parseable_asymmetric_structure():
    source = """data_si
_cell_length_a 5.43
_cell_length_b 5.43
_cell_length_c 5.43
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
Si1 Si 0 0 0 1
"""
    parsed_source = parse_cif(source)
    original, original_diagnostics = map_cif_structures(parsed_source.document)
    assert not [
        item for item in original_diagnostics if item.severity.value == "error"
    ]

    text = write_crystal_cif(original[0], block_name="silicon")
    reparsed = parse_cif(text)
    structures, diagnostics = map_cif_structures(reparsed.document)

    assert reparsed.ok
    assert not [item for item in diagnostics if item.severity.value == "error"]
    assert structures[0].sites[0].id.endswith("Si1:0")
    assert "_space_group_symop_operation_xyz" in text


def test_canonical_writer_preserves_mixed_components(read_fixture):
    original = read_fixture("mixed_disorder.cif").structures[0]

    text = write_crystal_cif(original)
    structures, diagnostics = map_cif_structures(parse_cif(text).document)

    assert not [item for item in diagnostics if item.severity.value == "error"]
    assert [component.element for component in structures[0].sites[0].components] == [
        "La",
        "Zr",
    ]


def test_canonical_writer_preserves_anisotropic_tensor(read_fixture):
    original = read_fixture("anisotropic.cif").structures[0]

    text = write_crystal_cif(original)
    structures, diagnostics = map_cif_structures(parse_cif(text).document)

    assert not [item for item in diagnostics if item.severity.value == "error"]
    assert structures[0].sites[0].displacement.tensor[0][1].value == 0.0012
