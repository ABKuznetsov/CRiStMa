"""Toolbox-neutral machine-readable diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Severity(str, Enum):
    """Severity of a scientific or source diagnostic."""

    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


@dataclass(frozen=True, slots=True)
class SourcePosition:
    """Zero-based offset and one-based source line and column."""

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
    """Stable diagnostic suitable for APIs and user interfaces."""

    severity: Severity
    code: str
    message: str
    span: SourceSpan | None = None
    recovery: str | None = None


__all__ = ["Diagnostic", "Severity", "SourcePosition", "SourceSpan"]
