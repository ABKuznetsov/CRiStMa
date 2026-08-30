"""Structure input and output contracts."""

from .diagnostics import Diagnostic, Severity, SourcePosition, SourceSpan
from .registry import FormatHandler, FormatRegistry
from .result import ReadResult, SourceInfo
from .source import (
    DecodedSource,
    MappingSourceResolver,
    ResolvedSource,
    SourceResolver,
    decode_bytes,
    decode_source,
)

__all__ = [
    "Diagnostic",
    "DecodedSource",
    "FormatHandler",
    "FormatRegistry",
    "ReadResult",
    "MappingSourceResolver",
    "ResolvedSource",
    "Severity",
    "SourceInfo",
    "SourceResolver",
    "SourcePosition",
    "SourceSpan",
    "decode_bytes",
    "decode_source",
]
