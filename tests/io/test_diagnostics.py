from cristma.io.diagnostics import (
    Diagnostic,
    Severity,
    SourcePosition,
    SourceSpan,
)
from cristma.io.result import ReadResult, SourceInfo


def test_diagnostic_carries_stable_code_and_source_span():
    span = SourceSpan(
        start=SourcePosition(offset=7, line=2, column=3),
        end=SourcePosition(offset=12, line=2, column=8),
    )
    diagnostic = Diagnostic(Severity.WARNING, "cif.loop.width", "short row", span)
    result = ReadResult(document=None, diagnostics=(diagnostic,))

    assert result.ok
    assert result.diagnostics[0].code == "cif.loop.width"
    assert result.diagnostics[0].span == span


def test_error_diagnostic_makes_result_not_ok():
    diagnostic = Diagnostic(Severity.ERROR, "io.empty", "empty source")

    assert not ReadResult(document=None, diagnostics=(diagnostic,)).ok


def test_read_result_preserves_decoding_and_newline_metadata():
    info = SourceInfo(
        name="sample.cif",
        format="cif",
        encoding="utf-8-sig",
        newline="\r\n",
    )

    assert ReadResult(document=None, source_info=info).source_info == info
