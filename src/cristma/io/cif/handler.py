"""Format-registry adapter for native CIF documents."""

from __future__ import annotations

from dataclasses import replace

from cristma.io.result import ReadResult
from cristma.structure import StructureCollection

from .mapper import map_cif_structures
from .parser import parse_cif
from .probe import probe_cif


class CifFormatHandler:
    """Recognize, parse, and scientifically map CIF 1.1 text."""

    name = "cif"
    suffixes = (".cif",)

    def probe(self, source: str) -> float:
        return probe_cif(source)

    def read_text(
        self,
        source: str,
        source_name: str | None = None,
    ) -> ReadResult:
        parsed = parse_cif(source, source_name=source_name)
        structures, mapping_diagnostics = map_cif_structures(parsed.document)
        return replace(
            parsed,
            structures=StructureCollection.from_structures(structures),
            diagnostics=(*parsed.diagnostics, *mapping_diagnostics),
        )
