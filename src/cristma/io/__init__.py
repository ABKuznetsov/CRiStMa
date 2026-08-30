"""Structure input and output contracts."""

from .diagnostics import Diagnostic, Severity, SourcePosition, SourceSpan
from .result import ReadResult, SourceInfo

__all__ = [
    "Diagnostic",
    "ReadResult",
    "Severity",
    "SourceInfo",
    "SourcePosition",
    "SourceSpan",
]
