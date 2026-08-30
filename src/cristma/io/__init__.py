"""Structure input and output contracts."""

from .diagnostics import Diagnostic, Severity, SourcePosition, SourceSpan
from .registry import FormatHandler, FormatRegistry
from .result import ReadResult, SourceInfo

__all__ = [
    "Diagnostic",
    "FormatHandler",
    "FormatRegistry",
    "ReadResult",
    "Severity",
    "SourceInfo",
    "SourcePosition",
    "SourceSpan",
]
