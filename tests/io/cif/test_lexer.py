from cristma.io.cif.lexer import lex_cif
from cristma.io.cif.tokens import CifTokenKind


def test_lexer_preserves_quotes_comments_and_spans():
    source = "data_a\n_tag 'two words' # note\n"

    tokens, diagnostics = lex_cif(source)

    assert not diagnostics
    assert [token.kind for token in tokens] == [
        CifTokenKind.DATA,
        CifTokenKind.TAG,
        CifTokenKind.VALUE,
        CifTokenKind.COMMENT,
    ]
    assert tokens[2].value == "two words"
    assert tokens[2].raw == "'two words'"
    assert source[tokens[2].span.start.offset : tokens[2].span.end.offset] == "'two words'"


def test_lexer_reads_semicolon_text_only_from_column_one():
    source = "data_a\n_note\n;line one\nline two\n;\n"

    tokens, diagnostics = lex_cif(source)

    assert not diagnostics
    assert tokens[-1].value == "line one\nline two"
    assert tokens[-1].raw == ";line one\nline two\n;"


def test_reserved_words_are_case_insensitive_but_raw_text_is_preserved():
    tokens, diagnostics = lex_cif("DATA_Test\nLOOP_\n_tag 1\n")

    assert not diagnostics
    assert [token.kind for token in tokens[:2]] == [
        CifTokenKind.DATA,
        CifTokenKind.LOOP,
    ]
    assert tokens[0].value == "Test"
    assert tokens[0].raw == "DATA_Test"


def test_unterminated_quote_reports_location():
    _tokens, diagnostics = lex_cif("data_a\n_tag 'broken\n")

    assert diagnostics[0].code == "cif.lex.unterminated_quote"
    assert diagnostics[0].span.start.line == 2


def test_indented_semicolon_is_an_ordinary_value():
    tokens, diagnostics = lex_cif("data_a\n_tag  ;value\n")

    assert not diagnostics
    assert tokens[-1].kind is CifTokenKind.VALUE
    assert tokens[-1].value == ";value"
