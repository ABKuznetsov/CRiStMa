"""Descriptor-first, content-aware dispatch for structure formats."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from threading import RLock

from .formats import FormatDescriptor, FormatHandler, descriptor_for
from .result import ReadResult, SourceInfo
from .source import decode_source


class FormatRegistry:
    """Select descriptors without loading their reader implementations."""

    def __init__(
        self,
        formats: tuple[FormatDescriptor | FormatHandler, ...] = (),
    ) -> None:
        self._descriptors = tuple(
            item if isinstance(item, FormatDescriptor) else descriptor_for(item)
            for item in formats
        )
        self._handler_cache: dict[str, FormatHandler] = {}
        self._lock = RLock()
        self._validate_unique_names()

    def _validate_unique_names(self) -> None:
        claimed: dict[str, str] = {}
        for descriptor in self._descriptors:
            for name in (descriptor.name, *descriptor.aliases):
                folded = name.casefold()
                if folded in claimed:
                    raise ValueError(
                        f"format name or alias {name!r} is already used by {claimed[folded]!r}"
                    )
                claimed[folded] = descriptor.name

    @property
    def descriptors(self) -> tuple[FormatDescriptor, ...]:
        return self._descriptors

    @property
    def handlers(self) -> tuple[FormatDescriptor, ...]:
        """Compatibility view; descriptors deliberately remain unloaded."""

        return self._descriptors

    def register(self, value: FormatDescriptor | FormatHandler) -> None:
        descriptor = value if isinstance(value, FormatDescriptor) else descriptor_for(value)
        self._descriptors = (*self._descriptors, descriptor)
        try:
            self._validate_unique_names()
        except ValueError:
            self._descriptors = self._descriptors[:-1]
            raise

    def select(
        self,
        source: str,
        *,
        suffix: str = "",
        basename: str = "",
        format: str | None = None,
    ) -> FormatDescriptor:
        if format is not None:
            requested = format.casefold()
            for descriptor in self._descriptors:
                if requested in {
                    descriptor.name.casefold(),
                    *(alias.casefold() for alias in descriptor.aliases),
                }:
                    return descriptor
            raise ValueError(f"Unknown structure format: {format}")

        normalized_suffix = suffix.casefold()
        normalized_basename = basename.casefold()
        scored: list[tuple[float, FormatDescriptor]] = []
        for descriptor in self._descriptors:
            confidence = float(descriptor.probe(source))
            if not 0.0 <= confidence <= 1.0:
                raise ValueError(
                    f"format probe {descriptor.name!r} returned invalid confidence {confidence}"
                )
            if normalized_basename and normalized_basename in {
                item.casefold() for item in descriptor.basenames
            }:
                confidence = max(confidence, 0.7)
            if normalized_suffix and normalized_suffix in {
                item.casefold() for item in descriptor.suffixes
            }:
                confidence = max(confidence, 0.6)
            if confidence > 0:
                scored.append((confidence, descriptor))

        if not scored:
            raise ValueError("No registered structure format recognized the source")
        best_score = max(score for score, _descriptor in scored)
        best = [descriptor for score, descriptor in scored if score == best_score]
        if len(best) != 1:
            names = ", ".join(sorted(descriptor.name for descriptor in best))
            raise ValueError(f"Ambiguous structure format: {names}")
        return best[0]

    def _handler(self, descriptor: FormatDescriptor) -> FormatHandler:
        key = descriptor.name.casefold()
        with self._lock:
            handler = self._handler_cache.get(key)
            if handler is None:
                handler = descriptor.factory()
                if not isinstance(handler, FormatHandler):
                    raise TypeError(
                        f"format factory {descriptor.name!r} returned an invalid handler"
                    )
                self._handler_cache[key] = handler
            return handler

    def read(self, path: str | Path, *, format: str | None = None) -> object:
        source_path = Path(path)
        decoded = decode_source(source_path)
        logical_path = Path(decoded.logical_name or source_path.name)
        descriptor = self.select(
            decoded.text,
            suffix=logical_path.suffix,
            basename=logical_path.name,
            format=format,
        )
        result = self._handler(descriptor).read_text(
            decoded.text,
            source_name=str(source_path),
        )
        if not isinstance(result, ReadResult):
            return result
        return replace(
            result,
            diagnostics=(*result.diagnostics, *decoded.diagnostics),
            source_info=SourceInfo(
                name=str(source_path),
                format=descriptor.name,
                encoding=decoded.encoding,
                newline=decoded.newline,
            ),
        )

    def read_text(
        self,
        source: str,
        *,
        format: str | None = None,
        source_name: str | None = None,
    ) -> object:
        source_path = Path(source_name) if source_name is not None else None
        descriptor = self.select(
            source,
            suffix=source_path.suffix if source_path is not None else "",
            basename=source_path.name if source_path is not None else "",
            format=format,
        )
        result = self._handler(descriptor).read_text(source, source_name=source_name)
        if not isinstance(result, ReadResult):
            return result
        newline = "\r\n" if "\r\n" in source else "\r" if "\r" in source else "\n"
        return replace(
            result,
            source_info=SourceInfo(
                name=source_name,
                format=descriptor.name,
                encoding="utf-8",
                newline=newline,
            ),
        )


__all__ = ["FormatHandler", "FormatRegistry"]
