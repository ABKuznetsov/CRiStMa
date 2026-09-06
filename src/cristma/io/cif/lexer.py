"""A source-preserving CIF 1.1 lexical state machine."""

from __future__ import annotations

from dataclasses import dataclass

from cristma.io.diagnostics import (
    Diagnostic,
    Severity,
    SourcePosition,
    SourceSpan,
)

from .tokens import CifToken, CifTokenKind


@dataclass(slots=True)
class _Cursor:
    source: str
    offset: int = 0
    line: int = 1
    column: int = 1

    def position(self) -> SourcePosition:
        return SourcePosition(self.offset, self.line, self.column)

    def advance_to(self, target: int) -> None:
        while self.offset < target:
            character = self.source[self.offset]
            self.offset += 1
            if character == "\n":
                self.line += 1
                self.column = 1
            else:
                self.column += 1


def _kind_and_value(raw: str) -> tuple[CifTokenKind, str]:
    folded = raw.casefold()
    if folded.startswith("data_"):
        return CifTokenKind.DATA, raw[5:]
    if folded == "global_":
        return CifTokenKind.GLOBAL, raw
    if folded == "loop_":
        return CifTokenKind.LOOP, raw
    if folded.startswith("save_"):
        return CifTokenKind.SAVE, raw[5:]
    if folded == "stop_":
        return CifTokenKind.STOP, raw
    if raw.startswith("_"):
        return CifTokenKind.TAG, raw
    return CifTokenKind.VALUE, raw


def _without_final_newline(value: str) -> str:
    if value.endswith("\r\n"):
        return value[:-2]
    if value.endswith(("\n", "\r")):
        return value[:-1]
    return value


def lex_cif(
    source: str,
) -> tuple[tuple[CifToken, ...], tuple[Diagnostic, ...]]:
    """Tokenize CIF source while retaining every token's original span."""

    cursor = _Cursor(source)
    tokens: list[CifToken] = []
    diagnostics: list[Diagnostic] = []
    length = len(source)

    while cursor.offset < length:
        character = source[cursor.offset]
        if character.isspace():
            cursor.advance_to(cursor.offset + 1)
            continue

        start = cursor.position()
        start_offset = cursor.offset

        if character == "#":
            end = cursor.offset
            while end < length and source[end] not in "\r\n":
                end += 1
            cursor.advance_to(end)
            raw = source[start_offset:end]
            tokens.append(
                CifToken(
                    CifTokenKind.COMMENT,
                    raw[1:],
                    raw,
                    SourceSpan(start, cursor.position()),
                )
            )
            continue

        if character == ";" and cursor.column == 1:
            closing = None
            search = cursor.offset + 1
            while search < length:
                newline = source.find("\n", search)
                if newline < 0:
                    break
                candidate = newline + 1
                if candidate < length and source[candidate] == ";":
                    closing = candidate
                    break
                search = candidate
            if closing is None:
                cursor.advance_to(length)
                diagnostics.append(
                    Diagnostic(
                        Severity.ERROR,
                        "cif.lex.unterminated_text",
                        "Unterminated semicolon-delimited CIF text field",
                        SourceSpan(start, cursor.position()),
                    )
                )
                break
            content = _without_final_newline(source[start_offset + 1 : closing])
            cursor.advance_to(closing + 1)
            tokens.append(
                CifToken(
                    CifTokenKind.VALUE,
                    content,
                    source[start_offset : closing + 1],
                    SourceSpan(start, cursor.position()),
                )
            )
            continue

        if character in "'\"":
            quote = character
            end = cursor.offset + 1
            while end < length and source[end] not in "\r\n":
                if source[end] == quote and (
                    end + 1 == length or source[end + 1].isspace()
                ):
                    break
                end += 1
            if end >= length or source[end] != quote:
                while end < length and source[end] not in "\r\n":
                    end += 1
                cursor.advance_to(end)
                diagnostics.append(
                    Diagnostic(
                        Severity.ERROR,
                        "cif.lex.unterminated_quote",
                        "Unterminated quoted CIF value",
                        SourceSpan(start, cursor.position()),
                    )
                )
                continue
            cursor.advance_to(end + 1)
            raw = source[start_offset : end + 1]
            tokens.append(
                CifToken(
                    CifTokenKind.VALUE,
                    raw[1:-1],
                    raw,
                    SourceSpan(start, cursor.position()),
                )
            )
            continue

        end = cursor.offset
        while end < length and not source[end].isspace() and source[end] != "#":
            end += 1
        cursor.advance_to(end)
        raw = source[start_offset:end]
        kind, value = _kind_and_value(raw)
        tokens.append(CifToken(kind, value, raw, SourceSpan(start, cursor.position())))

    return tuple(tokens), tuple(diagnostics)
