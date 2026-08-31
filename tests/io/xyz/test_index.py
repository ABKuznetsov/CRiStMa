from cristma.io.xyz.index import index_xyz


TWO_FRAME_SOURCE = """2
first
H 0 0 0
H 1 0 0
1
second
He 0 0 0
"""


def test_document_preserves_exact_source_and_indexes_frames() -> None:
    document, diagnostics = index_xyz(TWO_FRAME_SOURCE, "trajectory.xyz")

    assert diagnostics == ()
    assert document.render_preserved() == TWO_FRAME_SOURCE
    assert [frame.atom_count for frame in document.frames] == [2, 1]
    assert document.raw_source[
        document.frames[1].start_offset : document.frames[1].end_offset
    ] == "1\nsecond\nHe 0 0 0\n"


def test_zero_atom_frame_and_trailing_blank_lines_are_valid() -> None:
    source = "0\nempty frame\n\n\n"

    document, diagnostics = index_xyz(source, "empty.xyz")

    assert diagnostics == ()
    assert len(document.frames) == 1
    assert document.frames[0].atom_count == 0


def test_truncated_tail_is_not_indexed() -> None:
    source = TWO_FRAME_SOURCE + "2\ntruncated\nH 0 0 0\n"

    document, diagnostics = index_xyz(source, "run.xyz")

    assert len(document.frames) == 2
    assert any(item.code == "xyz.frame.incomplete" for item in diagnostics)


def test_unicode_and_crlf_offsets_slice_exact_frame() -> None:
    source = "1\r\nрасчёт воды\r\nO 0 0 0\r\n"

    document, diagnostics = index_xyz(source, "модель.xyz")
    span = document.frames[0]

    assert diagnostics == ()
    assert document.raw_source[span.start_offset : span.end_offset] == source
    assert document.raw_source[
        span.comment_start_offset : span.comment_end_offset
    ] == "расчёт воды\r\n"


def test_blank_line_between_frames_is_recovered_with_diagnostic() -> None:
    source = "1\nfirst\nH 0 0 0\n\n1\nsecond\nHe 0 0 0\n"

    document, diagnostics = index_xyz(source)

    assert len(document.frames) == 2
    assert any(item.code == "xyz.frame.blank_between_frames" for item in diagnostics)


def test_invalid_count_stops_indexing() -> None:
    source = TWO_FRAME_SOURCE + "not-a-count\ncomment\n"

    document, diagnostics = index_xyz(source)

    assert len(document.frames) == 2
    assert any(item.code == "xyz.frame.count_invalid" for item in diagnostics)


def test_negative_count_is_invalid() -> None:
    document, diagnostics = index_xyz("-1\ninvalid\n")

    assert len(document.frames) == 0
    assert any(item.code == "xyz.frame.count_invalid" for item in diagnostics)
