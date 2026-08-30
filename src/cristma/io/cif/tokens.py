"""Token types emitted by the native CIF lexer."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from cristma.io.diagnostics import SourceSpan


class CifTokenKind(str, Enum):
    COMMENT = "comment"
    DATA = "data"
    GLOBAL = "global"
    LOOP = "loop"
    SAVE = "save"
    STOP = "stop"
    TAG = "tag"
    VALUE = "value"


@dataclass(frozen=True, slots=True)
class CifToken:
    kind: CifTokenKind
    value: str
    raw: str
    span: SourceSpan
