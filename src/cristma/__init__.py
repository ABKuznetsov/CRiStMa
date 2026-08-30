"""Stable public API for the Qt-free CRiStMa scientific library."""

from __future__ import annotations

from pathlib import Path

from .core.structure import Crystal
from .io.cif.document import CifDocument
from .io.cif.handler import CifFormatHandler
from .io.cif.writer import write_cif_document, write_crystal_cif
from .io.registry import FormatRegistry
from .io.result import ReadResult

__version__ = "0.1.0"

_formats = FormatRegistry((CifFormatHandler(),))


def read(path: str | Path, *, format: str | None = None) -> ReadResult:
    """Read a structure file through the native format registry."""

    result = _formats.read(path, format=format)
    if not isinstance(result, ReadResult):
        raise TypeError("structure format handler returned an invalid read result")
    return result


def read_text(
    source: str,
    *,
    format: str | None = None,
    source_name: str | None = None,
) -> ReadResult:
    """Read already decoded structure text through content probing."""

    result = _formats.read_text(
        source,
        format=format,
        source_name=source_name,
    )
    if not isinstance(result, ReadResult):
        raise TypeError("structure format handler returned an invalid read result")
    return result


def write(
    value: CifDocument | Crystal,
    path: str | Path,
    *,
    mode: str | None = None,
) -> None:
    """Write a preserved CIF document or canonical asymmetric-unit crystal."""

    if isinstance(value, CifDocument):
        selected_mode = "preserve" if mode is None else mode
        if selected_mode != "preserve":
            raise ValueError("CifDocument write mode must be 'preserve'")
        rendered = write_cif_document(value, mode=selected_mode)
    elif isinstance(value, Crystal):
        selected_mode = "canonical" if mode is None else mode
        if selected_mode != "canonical":
            raise ValueError("Crystal write mode must be 'canonical'")
        rendered = write_crystal_cif(value)
    else:
        raise TypeError("write accepts only CifDocument or Crystal")
    Path(path).write_bytes(rendered.encode("utf-8"))


__all__ = ["ReadResult", "__version__", "read", "read_text", "write"]
