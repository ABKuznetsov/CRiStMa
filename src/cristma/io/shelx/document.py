"""Loss-preserving source model for SHELX instruction documents."""

from __future__ import annotations

from dataclasses import dataclass

from cristma.io.diagnostics import SourceSpan


@dataclass(frozen=True, slots=True)
class ShelxPhysicalLine:
    """One physical source line without its retained newline terminator."""

    text: str
    newline: str
    span: SourceSpan


@dataclass(frozen=True, slots=True)
class ShelxSourceEdit:
    """A half-open source replacement applied only during preserve writing."""

    start: int
    end: int
    replacement: str


@dataclass(frozen=True, slots=True)
class ShelxRecord:
    """One logical record assembled from one or more physical lines."""

    keyword: str | None
    fields: tuple[str, ...]
    physical_line_indices: tuple[int, ...]
    span: SourceSpan
    inline_comment: str | None = None
    after_hklf: bool = False
    after_end: bool = False


@dataclass(frozen=True, slots=True)
class ShelxInstructionRecord(ShelxRecord):
    """A recognized SHELX instruction not yet scientifically interpreted."""


@dataclass(frozen=True, slots=True)
class ShelxAtomRecord(ShelxRecord):
    """An atom-like logical record."""


@dataclass(frozen=True, slots=True)
class ShelxQPeakRecord(ShelxRecord):
    """A difference-Fourier peak retained outside canonical chemistry."""


@dataclass(frozen=True, slots=True)
class ShelxCommentRecord(ShelxRecord):
    """A REM comment record."""


@dataclass(frozen=True, slots=True)
class ShelxBlankRecord(ShelxRecord):
    """A blank physical line retained in record order."""


@dataclass(frozen=True, slots=True)
class ShelxUnknownRecord(ShelxRecord):
    """An uninterpreted non-atom record retained without data loss."""


@dataclass(frozen=True, slots=True)
class ShelxDocument:
    """An immutable SHELX document with exact original source text."""

    raw_source: str
    physical_lines: tuple[ShelxPhysicalLine, ...]
    records: tuple[ShelxRecord, ...]
    source_name: str | None = None
    edits: tuple[ShelxSourceEdit, ...] = ()

    def render_preserved(self) -> str:
        """Render source edits while retaining every untouched character."""

        ordered = sorted(self.edits, key=lambda item: (item.start, item.end))
        for edit in ordered:
            if edit.start < 0 or edit.end < edit.start or edit.end > len(self.raw_source):
                raise ValueError("invalid SHELX source edit span")
        for previous, current in zip(ordered, ordered[1:]):
            if current.start < previous.end:
                raise ValueError("overlapping SHELX source edits")
        rendered = self.raw_source
        for edit in reversed(ordered):
            rendered = rendered[: edit.start] + edit.replacement + rendered[edit.end :]
        return rendered


__all__ = [
    "ShelxAtomRecord",
    "ShelxBlankRecord",
    "ShelxCommentRecord",
    "ShelxDocument",
    "ShelxInstructionRecord",
    "ShelxPhysicalLine",
    "ShelxQPeakRecord",
    "ShelxRecord",
    "ShelxSourceEdit",
    "ShelxUnknownRecord",
]
