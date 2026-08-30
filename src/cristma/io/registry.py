"""Content-aware dispatch for structure document formats."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Protocol, runtime_checkable

from .diagnostics import Diagnostic, Severity
from .result import ReadResult, SourceInfo


@runtime_checkable
class FormatHandler(Protocol):
    """Reader contract registered with :class:`FormatRegistry`."""

    name: str
    suffixes: tuple[str, ...]

    def probe(self, source: str) -> float:
        """Return content confidence from zero to one."""

    def read_text(self, source: str, source_name: str | None = None) -> object:
        """Read an already decoded source string."""


def _newline_style(source: str) -> str:
    if "\r\n" in source:
        return "\r\n"
    if "\r" in source:
        return "\r"
    return "\n"


class FormatRegistry:
    """Select structure readers by explicit name, suffix, and content."""

    def __init__(self, handlers: tuple[FormatHandler, ...] = ()) -> None:
        names = [handler.name.casefold() for handler in handlers]
        if len(names) != len(set(names)):
            raise ValueError("format handler names must be unique")
        self._handlers = tuple(handlers)

    @property
    def handlers(self) -> tuple[FormatHandler, ...]:
        return self._handlers

    def register(self, handler: FormatHandler) -> None:
        if any(existing.name.casefold() == handler.name.casefold() for existing in self._handlers):
            raise ValueError(f"format handler already registered: {handler.name}")
        self._handlers = (*self._handlers, handler)

    def select(
        self,
        source: str,
        *,
        suffix: str = "",
        format: str | None = None,
    ) -> FormatHandler:
        if format is not None:
            for handler in self._handlers:
                if handler.name.casefold() == format.casefold():
                    return handler
            raise ValueError(f"Unknown structure format: {format}")

        normalized_suffix = suffix.casefold()
        scored: list[tuple[float, FormatHandler]] = []
        for handler in self._handlers:
            confidence = float(handler.probe(source))
            if not 0.0 <= confidence <= 1.0:
                raise ValueError(
                    f"format probe {handler.name!r} returned invalid confidence {confidence}"
                )
            suffixes = {item.casefold() for item in handler.suffixes}
            if normalized_suffix and normalized_suffix in suffixes:
                confidence = max(confidence, 0.6)
            if confidence > 0:
                scored.append((confidence, handler))

        if not scored:
            raise ValueError("No registered structure format recognized the source")
        best_score = max(item[0] for item in scored)
        best = [handler for score, handler in scored if score == best_score]
        if len(best) != 1:
            names = ", ".join(sorted(handler.name for handler in best))
            raise ValueError(f"Ambiguous structure format: {names}")
        return best[0]

    def read(
        self,
        path: str | Path,
        *,
        format: str | None = None,
    ) -> object:
        source_path = Path(path)
        raw = source_path.read_bytes()
        decoding_diagnostic: Diagnostic | None = None
        try:
            source = raw.decode("utf-8-sig")
            encoding = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf-8"
        except UnicodeDecodeError:
            source = raw.decode("latin-1")
            encoding = "latin-1"
            decoding_diagnostic = Diagnostic(
                Severity.WARNING,
                "io.encoding_fallback",
                "Source is not valid UTF-8; decoded as Latin-1",
                recovery="latin-1",
            )

        handler = self.select(
            source,
            suffix=source_path.suffix,
            format=format,
        )
        result = handler.read_text(source, source_name=str(source_path))
        if not isinstance(result, ReadResult):
            return result

        diagnostics = result.diagnostics
        if decoding_diagnostic is not None:
            diagnostics = (*diagnostics, decoding_diagnostic)
        return replace(
            result,
            diagnostics=diagnostics,
            source_info=SourceInfo(
                name=str(source_path),
                format=handler.name,
                encoding=encoding,
                newline=_newline_style(source),
            ),
        )

    def read_text(
        self,
        source: str,
        *,
        format: str | None = None,
        source_name: str | None = None,
    ) -> object:
        """Select a handler and read an already decoded source string."""

        suffix = Path(source_name).suffix if source_name is not None else ""
        handler = self.select(source, suffix=suffix, format=format)
        result = handler.read_text(source, source_name=source_name)
        if not isinstance(result, ReadResult):
            return result
        return replace(
            result,
            source_info=SourceInfo(
                name=source_name,
                format=handler.name,
                encoding="utf-8",
                newline=_newline_style(source),
            ),
        )
