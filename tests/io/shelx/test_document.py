from dataclasses import replace

import pytest

from cristma.io.shelx.document import ShelxSourceEdit
from cristma.io.shelx.parser import parse_shelx


def test_physical_lines_retain_text_newlines_and_source_spans() -> None:
    source = "TITL demo\r\n\r\nEND\r\n"

    document = parse_shelx(source, source_name="demo.res").document

    assert [(line.text, line.newline) for line in document.physical_lines] == [
        ("TITL demo", "\r\n"),
        ("", "\r\n"),
        ("END", "\r\n"),
    ]
    assert document.physical_lines[0].span.start.offset == 0
    assert document.physical_lines[0].span.end.offset == len("TITL demo")
    assert document.physical_lines[2].span.start.line == 3
    assert document.render_preserved() == source


def test_source_edit_changes_only_its_half_open_span() -> None:
    source = "TITL demo\nCELL 0.71073 10 10 10 90 90 90\nEND\n"
    document = parse_shelx(source).document
    start = source.index("demo")
    edited = replace(
        document,
        edits=(ShelxSourceEdit(start=start, end=start + 4, replacement="changed"),),
    )

    assert edited.render_preserved() == source.replace("demo", "changed", 1)


def test_overlapping_source_edits_are_rejected_when_rendered() -> None:
    document = parse_shelx("TITL demo\nEND\n").document
    edited = replace(
        document,
        edits=(
            ShelxSourceEdit(0, 5, "TITLE"),
            ShelxSourceEdit(4, 8, "overlap"),
        ),
    )

    with pytest.raises(ValueError, match="overlapping SHELX source edits"):
        edited.render_preserved()


@pytest.mark.parametrize("start,end", [(-1, 0), (3, 2), (0, 99)])
def test_invalid_source_edit_bounds_are_rejected(start: int, end: int) -> None:
    document = parse_shelx("END\n").document
    edited = replace(
        document,
        edits=(ShelxSourceEdit(start, end, "x"),),
    )

    with pytest.raises(ValueError, match="invalid SHELX source edit span"):
        edited.render_preserved()
