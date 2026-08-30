"""Line-oriented parser for loss-preserving SHELX documents."""

from __future__ import annotations

from cristma.core.cell import UnitCell
from cristma.core.values import parse_measured_value
from cristma.io.diagnostics import Diagnostic, Severity, SourcePosition, SourceSpan
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
from .records import (
    ShelxCellInstruction,
    ShelxEndInstruction,
    ShelxFvarInstruction,
    ShelxHklfInstruction,
    ShelxLattInstruction,
    ShelxPartInstruction,
    ShelxResiInstruction,
    ShelxSfacInstruction,
    ShelxSymmInstruction,
    ShelxZerrInstruction,
)
from .symmetry import parse_shelx_symmetry


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


def _typed_instruction(record: ShelxRecord) -> ShelxRecord:
    if not isinstance(record, ShelxInstructionRecord) or record.keyword is None:
        return record
    common = {
        "keyword": record.keyword,
        "fields": record.fields,
        "physical_line_indices": record.physical_line_indices,
        "span": record.span,
        "inline_comment": record.inline_comment,
        "after_hklf": record.after_hklf,
        "after_end": record.after_end,
    }
    fields = record.fields
    if record.keyword == "CELL":
        if len(fields) != 7:
            raise ValueError("CELL requires wavelength and six cell values")
        wavelength = parse_measured_value(fields[0], unit="angstrom")
        cell_values = (
            *(parse_measured_value(value, unit="angstrom") for value in fields[1:4]),
            *(parse_measured_value(value, unit="degree") for value in fields[4:7]),
        )
        return ShelxCellInstruction(**common, wavelength=wavelength, cell=UnitCell(*cell_values))
    if record.keyword == "ZERR":
        if len(fields) != 7:
            raise ValueError("ZERR requires Z and six cell uncertainties")
        uncertainties = (
            *(parse_measured_value(value, unit="angstrom") for value in fields[1:4]),
            *(parse_measured_value(value, unit="degree") for value in fields[4:7]),
        )
        return ShelxZerrInstruction(
            **common,
            formula_units=parse_measured_value(fields[0]),
            cell_uncertainties=tuple(uncertainties),
        )
    if record.keyword == "LATT":
        if len(fields) != 1:
            raise ValueError("LATT requires one integer code")
        return ShelxLattInstruction(**common, code=int(fields[0]))
    if record.keyword == "SYMM":
        return ShelxSymmInstruction(
            **common,
            operation=parse_shelx_symmetry(" ".join(fields)),
        )
    if record.keyword == "SFAC":
        if not fields:
            raise ValueError("SFAC requires at least one entry")
        return ShelxSfacInstruction(**common, entries=fields)
    if record.keyword == "FVAR":
        if not fields:
            raise ValueError("FVAR requires at least one value")
        return ShelxFvarInstruction(
            **common,
            values=tuple(parse_measured_value(value) for value in fields),
        )
    if record.keyword == "PART":
        if not fields:
            raise ValueError("PART requires a part number")
        return ShelxPartInstruction(
            **common,
            part=int(fields[0]),
            occupancy_code=fields[1] if len(fields) > 1 else None,
        )
    if record.keyword == "RESI":
        residue_number = None
        residue_class = None
        if fields:
            try:
                residue_number = int(fields[0])
            except ValueError:
                residue_class = fields[0]
            else:
                residue_class = fields[1] if len(fields) > 1 else None
        return ShelxResiInstruction(
            **common,
            residue_number=residue_number,
            residue_class=residue_class,
        )
    if record.keyword == "HKLF":
        if not fields:
            raise ValueError("HKLF requires a code")
        return ShelxHklfInstruction(**common, code=int(fields[0]))
    if record.keyword == "END":
        return ShelxEndInstruction(**common)
    return record


def parse_shelx(source: str, source_name: str | None = None) -> ReadResult:
    """Parse physical and logical SHELX records without scientific mapping."""

    lines = _physical_lines(source)
    records: list[ShelxRecord] = []
    diagnostics: list[Diagnostic] = []
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
        try:
            record = _typed_instruction(record)
        except (TypeError, ValueError) as error:
            diagnostics.append(
                Diagnostic(
                    Severity.ERROR,
                    f"shelx.parse.invalid_{keyword.casefold()}",
                    str(error),
                    span,
                )
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
        diagnostics=tuple(diagnostics),
        source_info=SourceInfo(
            name=source_name,
            format="shelx",
            encoding="utf-8",
            newline=newline,
        ),
    )


__all__ = ["parse_shelx"]
