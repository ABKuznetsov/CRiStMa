from cristma.io.shelx.document import (
    ShelxAtomRecord,
    ShelxBlankRecord,
    ShelxCommentRecord,
    ShelxInstructionRecord,
    ShelxQPeakRecord,
    ShelxUnknownRecord,
)
from cristma.io.shelx.parser import parse_shelx


def test_parser_classifies_records_without_losing_original_text() -> None:
    source = (
        "TITL demo\n"
        "REM keep this\n"
        "\n"
        "C1 1 0.1 0.2 0.3 11 0.05\n"
        "FROB arbitrary payload\n"
        "END\n"
        "Q1 1 0.2 0.3 0.4 11 0.05 1.2\n"
    )

    result = parse_shelx(source, source_name="demo.res")

    assert result.ok
    assert [type(record) for record in result.document.records] == [
        ShelxInstructionRecord,
        ShelxCommentRecord,
        ShelxBlankRecord,
        ShelxAtomRecord,
        ShelxUnknownRecord,
        ShelxInstructionRecord,
        ShelxQPeakRecord,
    ]
    assert result.document.records[0].keyword == "TITL"
    assert result.document.records[1].fields == ("keep", "this")
    assert result.document.records[-1].after_end
    assert result.document.render_preserved() == source


def test_parser_assembles_continuation_into_one_logical_record() -> None:
    source = (
        "C1 1 0.106507 1.433210 0.096736 11.00000 0.05548 =\r\n"
        "       0.03094 0.01587 -0.00145 -0.00054\r\n"
        "END\r\n"
    )

    result = parse_shelx(source)
    atom = result.document.records[0]

    assert isinstance(atom, ShelxAtomRecord)
    assert atom.fields == (
        "1",
        "0.106507",
        "1.433210",
        "0.096736",
        "11.00000",
        "0.05548",
        "0.03094",
        "0.01587",
        "-0.00145",
        "-0.00054",
    )
    assert atom.physical_line_indices == (0, 1)
    assert atom.span.start.line == 1
    assert atom.span.end.line == 2
    assert result.document.render_preserved() == source


def test_inline_comment_is_excluded_from_fields_but_preserved() -> None:
    source = "SYMM -X, Y, -Z ! twofold operation\nEND\n"

    result = parse_shelx(source)
    record = result.document.records[0]

    assert isinstance(record, ShelxInstructionRecord)
    assert record.keyword == "SYMM"
    assert record.fields == ("-X,", "Y,", "-Z")
    assert record.inline_comment == "twofold operation"
    assert result.document.render_preserved() == source
