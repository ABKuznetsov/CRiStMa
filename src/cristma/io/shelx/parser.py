"""Line-oriented parser for loss-preserving SHELX documents."""

from __future__ import annotations

from cristma.io.diagnostics import SourcePosition, SourceSpan
from cristma.io.result import ReadResult, SourceInfo

from .document import (
    ShelxAtomRecord,
    ShelxBlankRecord,
    ShelxCommentRecord,
    ShelxDocument,
    ShelxInstructionRecord,
    ShelxPhysicalLine,
    ShelxQPeakRecord,
    ShelxRecord,
    ShelxUnknownRecord,
)


_INSTRUCTIONS = frozenset(
    {
        "AFIX", "ANIS", "BASF", "BIND", "BLOC", "BOND", "CELL", "CONF",
        "CONN", "DAMP", "DANG", "DEFS", "DELU", "DFIX", "EADP", "END",
        "EQIV", "EXTI", "EXYZ", "FEND", "FLAT", "FMAP", "FRAG", "FREE",
        "FVAR", "GRID", "HFIX", "HKLF", "HTAB", "LATT", "LIST", "MERG",
        "MORE", "MOVE", "MPLA", "OMIT", "PART", "PLAN", "REM", "RESI",
        "RIGU", "RTAB", "SADI", "SAME", "SFAC", "SIMU", "SIZE", "SPEC",
        "STIR", "SWAT", "SYMM", "TEMP", "TITL", "TWIN", "UNIT", "WGHT",
        "WIGL", "ZERR",
    }
)


def _split_line_ending(raw: str) -> tuple[str, str]:
    if raw.endswith("\r\n"):
        return raw[:-2], "\r\n"
    if raw.endswith(("\n", "\r")):
        return raw[:-1], raw[-1]
    return raw, ""


def _physical_lines(source: str) -> tuple[ShelxPhysicalLine, ...]:
    rows: list[ShelxPhysicalLine] = []
    offset = 0
    for line_number, raw in enumerate(source.splitlines(keepends=True), start=1):
        text, newline = _split_line_ending(raw)
        rows.append(
            ShelxPhysicalLine(
                text=text,
                newline=newline,
                span=SourceSpan(
                    SourcePosition(offset, line_number, 1),
                    SourcePosition(offset + len(text), line_number, len(text) + 1),
                ),
            )
        )
        offset += len(raw)
    return tuple(rows)


def _without_inline_comment(text: str) -> tuple[str, str | None]:
    body, marker, comment = text.partition("!")
    return body, comment.strip() if marker else None


def _continues(text: str) -> bool:
    body, _comment = _without_inline_comment(text)
    return body.rstrip().endswith("=")


def _logical_text(lines: tuple[ShelxPhysicalLine, ...], indices: tuple[int, ...]) -> tuple[str, str | None]:
    parts: list[str] = []
    inline_comment: str | None = None
    for index in indices:
        body, comment = _without_inline_comment(lines[index].text)
        stripped = body.strip()
        if stripped.endswith("="):
            stripped = stripped[:-1].rstrip()
        if stripped:
            parts.append(stripped)
        if comment is not None:
            inline_comment = comment
    return " ".join(parts), inline_comment


def _looks_numeric(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True


def _record_type(keyword: str | None, fields: tuple[str, ...]) -> type[ShelxRecord]:
    if keyword is None:
        return ShelxBlankRecord
    if keyword == "REM":
        return ShelxCommentRecord
    if keyword in _INSTRUCTIONS:
        return ShelxInstructionRecord
    if keyword.startswith("Q") and keyword[1:].isdigit() and fields and _looks_numeric(fields[0]):
        return ShelxQPeakRecord
    if fields and _looks_numeric(fields[0]) and len(fields) >= 6:
        return ShelxAtomRecord
    return ShelxUnknownRecord


def parse_shelx(source: str, source_name: str | None = None) -> ReadResult:
    """Parse physical and logical SHELX records without scientific mapping."""

    lines = _physical_lines(source)
    records: list[ShelxRecord] = []
    after_hklf = False
    after_end = False
    index = 0
    while index < len(lines):
        indices = [index]
        while _continues(lines[indices[-1]].text) and indices[-1] + 1 < len(lines):
            indices.append(indices[-1] + 1)
        logical, inline_comment = _logical_text(lines, tuple(indices))
        tokens = logical.split()
        keyword = tokens[0].upper() if tokens else None
        fields = tuple(tokens[1:])
        record_class = _record_type(keyword, fields)
        span = SourceSpan(lines[indices[0]].span.start, lines[indices[-1]].span.end)
        record = record_class(
            keyword=keyword,
            fields=fields,
            physical_line_indices=tuple(indices),
            span=span,
            inline_comment=inline_comment,
            after_hklf=after_hklf,
            after_end=after_end,
        )
        records.append(record)
        if keyword == "HKLF":
            after_hklf = True
        if keyword == "END":
            after_end = True
        index = indices[-1] + 1

    newline = "\r\n" if "\r\n" in source else "\r" if "\r" in source else "\n"
    document = ShelxDocument(source, lines, tuple(records), source_name=source_name)
    return ReadResult(
        document=document,
        source_info=SourceInfo(
            name=source_name,
            format="shelx",
            encoding="utf-8",
            newline=newline,
        ),
    )


__all__ = ["parse_shelx"]
