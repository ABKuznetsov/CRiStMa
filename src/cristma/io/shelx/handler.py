"""Format-registry adapter for native SHELX documents."""

from __future__ import annotations

from dataclasses import replace

from cristma.io.result import ReadResult

from .mapper import map_shelx_structures
from .parser import parse_shelx
from .probe import probe_shelx


class ShelxFormatHandler:
    """Recognize, parse, and scientifically map RES/INS source text."""

    name = "shelx"
    suffixes = (".res", ".ins")

    def probe(self, source: str) -> float:
        return probe_shelx(source)

    def read_text(
        self,
        source: str,
        source_name: str | None = None,
    ) -> ReadResult:
        parsed = parse_shelx(source, source_name=source_name)
        structures, mapping_diagnostics = map_shelx_structures(parsed.document)
        return replace(
            parsed,
            structures=structures,
            diagnostics=(*parsed.diagnostics, *mapping_diagnostics),
        )


__all__ = ["ShelxFormatHandler"]
