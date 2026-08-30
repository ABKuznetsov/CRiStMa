"""Machine-readable diagnostics with precise source locations."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Severity(str, Enum):
    """Severity of a parser or scientific mapping diagnostic."""

    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


@dataclass(frozen=True, slots=True)
class SourcePosition:
    """Zero-based byte/character offset and one-based line and column."""

    offset: int
    line: int
    column: int


@dataclass(frozen=True, slots=True)
class SourceSpan:
    """Half-open source range."""

    start: SourcePosition
    end: SourcePosition


@dataclass(frozen=True, slots=True)
class Diagnostic:
    """A stable diagnostic suitable for both APIs and user interfaces."""

    severity: Severity
    code: str
    message: str
    span: SourceSpan | None = None
    recovery: str | None = None
