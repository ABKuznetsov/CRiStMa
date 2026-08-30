from cristma.io.shelx.mapper import map_shelx_structures
from cristma.io.shelx.parser import parse_shelx


def test_atom_mapping_error_points_to_physical_source_line() -> None:
    source = (
        "CELL 0.71073 10 10 10 90 90 90\nLATT -1\nSFAC C\n"
        "C1 1 bad 0.2 0.3 11 0.05\nEND\n"
    )

    structures, diagnostics = map_shelx_structures(parse_shelx(source).document)
    error = next(item for item in diagnostics if item.code == "shelx.map.atom_invalid")

    assert not structures
    assert error.span.start.line == 4
    assert error.span.end.line == 4
