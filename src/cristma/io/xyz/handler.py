"""Format-registry adapter for native XYZ structure sources."""

from __future__ import annotations

from cristma.io.result import ReadResult

from .probe import probe_xyz


class XyzFormatHandler:
    """Read plain XYZ and extXYZ through the application-neutral API."""

    name = "xyz"
    suffixes = (".xyz", ".extxyz")

    def probe(self, source: str) -> float:
        return probe_xyz(source)

    def read_text(self, source: str, source_name: str | None = None) -> ReadResult:
        from .parser import parse_xyz

        return parse_xyz(source, source_name)


__all__ = ["XyzFormatHandler"]
