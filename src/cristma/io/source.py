"""Deterministic source decoding and explicit auxiliary-file resolution."""

from __future__ import annotations

import bz2
from dataclasses import dataclass
import gzip
import lzma
from pathlib import Path, PurePosixPath
from types import MappingProxyType
from typing import Literal, Mapping, Protocol, runtime_checkable

from .diagnostics import Diagnostic, Severity


Compression = Literal["gzip", "bzip2", "xz"]


@dataclass(frozen=True, slots=True)
class DecodedSource:
    """Original bytes plus deterministic decoded text and source metadata."""

    raw: bytes
    text: str
    logical_name: str | None
    encoding: str
    newline: str
    compression: Compression | None
    diagnostics: tuple[Diagnostic, ...] = ()


@dataclass(frozen=True, slots=True)
class ResolvedSource:
    """One explicitly supplied auxiliary source."""

    reference: str
    raw: bytes
    source_name: str | None = None


@runtime_checkable
class SourceResolver(Protocol):
    """Capability passed to readers that are allowed to resolve references."""

    def resolve(
        self,
        reference: str,
        *,
        from_source: str | None,
    ) -> ResolvedSource | None: ...


def _newline_style(source: str) -> str:
    if "\r\n" in source:
        return "\r\n"
    if "\r" in source:
        return "\r"
    return "\n"


def _decompress(raw: bytes) -> tuple[bytes, Compression | None, str | None]:
    if raw.startswith(b"\x1f\x8b"):
        return gzip.decompress(raw), "gzip", ".gz"
    if raw.startswith(b"BZh"):
        return bz2.decompress(raw), "bzip2", ".bz2"
    if raw.startswith(b"\xfd7zXZ\x00"):
        return lzma.decompress(raw), "xz", ".xz"
    return raw, None, None


def decode_bytes(raw: bytes, source_name: str | None = None) -> DecodedSource:
    """Decode bytes after magic-based compression detection."""

    original = bytes(raw)
    payload, compression, compression_suffix = _decompress(original)
    diagnostics: tuple[Diagnostic, ...] = ()
    try:
        text = payload.decode("utf-8-sig")
        encoding = "utf-8-sig" if payload.startswith(b"\xef\xbb\xbf") else "utf-8"
    except UnicodeDecodeError:
        text = payload.decode("latin-1")
        encoding = "latin-1"
        diagnostics = (
            Diagnostic(
                Severity.WARNING,
                "io.encoding_fallback",
                "Source is not valid UTF-8; decoded as Latin-1",
                recovery="latin-1",
            ),
        )

    logical_name = source_name
    if (
        logical_name is not None
        and compression_suffix is not None
        and logical_name.casefold().endswith(compression_suffix)
    ):
        logical_name = logical_name[: -len(compression_suffix)]
    return DecodedSource(
        raw=original,
        text=text,
        logical_name=logical_name,
        encoding=encoding,
        newline=_newline_style(text),
        compression=compression,
        diagnostics=diagnostics,
    )


def decode_source(path: str | Path) -> DecodedSource:
    """Read and decode a source path without inspecting neighboring files."""

    source_path = Path(path)
    return decode_bytes(source_path.read_bytes(), str(source_path))


def _safe_reference(reference: str) -> str | None:
    candidate = PurePosixPath(reference)
    if candidate.is_absolute() or ".." in candidate.parts:
        return None
    normalized = str(candidate)
    if normalized in {"", "."}:
        return None
    return normalized


class MappingSourceResolver:
    """Resolve only sources supplied explicitly in an immutable mapping."""

    __slots__ = ("_sources",)

    def __init__(self, sources: Mapping[str, bytes]) -> None:
        normalized: dict[str, bytes] = {}
        for name, raw in sources.items():
            key = _safe_reference(name)
            if key is None:
                raise ValueError(f"unsafe source mapping name: {name!r}")
            normalized[key] = bytes(raw)
        self._sources = MappingProxyType(normalized)

    def resolve(
        self,
        reference: str,
        *,
        from_source: str | None,
    ) -> ResolvedSource | None:
        key = _safe_reference(reference)
        if key is None:
            return None
        candidates = [key]
        if from_source is not None:
            parent = PurePosixPath(from_source).parent
            if str(parent) not in {"", "."}:
                combined = _safe_reference(str(parent / key))
                if combined is not None:
                    candidates.insert(0, combined)
        for candidate in candidates:
            raw = self._sources.get(candidate)
            if raw is not None:
                return ResolvedSource(candidate, raw, source_name=candidate)
        return None


__all__ = [
    "DecodedSource",
    "MappingSourceResolver",
    "ResolvedSource",
    "SourceResolver",
    "decode_bytes",
    "decode_source",
]
