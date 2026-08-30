from cristma.io.cif.parser import parse_cif


def test_parser_keeps_multiple_blocks_scalars_and_loops():
    source = """data_first
_cell_length_a 5.0
loop_
_atom_site_label
_atom_site_fract_x
Si1 0
data_second
_audit_note ?
"""

    result = parse_cif(source)

    assert result.ok
    assert [block.name for block in result.document.blocks] == ["first", "second"]
    assert result.document.blocks[0].scalar("_CELL_LENGTH_A").value == "5.0"
    assert result.document.blocks[0].loops[0].rows == (("Si1", "0"),)
    assert result.document.blocks[1].scalar("_audit_note").value == "?"


def test_short_loop_row_is_reported_without_fabricated_cell():
    result = parse_cif("data_a\nloop_\n_a\n_b\n1\n")

    assert not result.ok
    assert result.diagnostics[-1].code == "cif.parse.loop_width"
    assert result.document.blocks[0].loops[0].rows == ()
    assert result.document.blocks[0].loops[0].incomplete_values == ("1",)


def test_parser_retains_unknown_tag_and_comment():
    source = "data_a\n# instrument note\n_local_detector_mode fast\n"

    result = parse_cif(source)

    assert result.document.blocks[0].scalar("_local_detector_mode").value == "fast"
    assert result.document.blocks[0].comments[0].raw == "# instrument note"
    assert result.document.raw_source == source


def test_scalar_without_value_is_reported_and_not_fabricated():
    result = parse_cif("data_a\n_cell_length_a\n_cell_length_b 5\n")

    assert not result.ok
    assert result.document.blocks[0].scalar("_cell_length_a") is None
    assert result.document.blocks[0].scalar("_cell_length_b").value == "5"
    assert "cif.parse.scalar_value_missing" in {item.code for item in result.diagnostics}


def test_data_item_before_first_block_is_reported():
    result = parse_cif("_audit_note orphan\ndata_a\n_tag value\n")

    assert not result.ok
    assert [block.name for block in result.document.blocks] == ["a"]
    assert result.document.blocks[0].scalar("_tag").value == "value"
