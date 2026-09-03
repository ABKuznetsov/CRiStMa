"""Format-registry adapter for native PDB sources."""

from __future__ import annotations

from cristma.io.result import ReadResult

from .probe import probe_pdb


class PdbFormatHandler:
    name = "pdb"
    suffixes = (".pdb",)

    def probe(self, source: str) -> float:
        return probe_pdb(source)

    def read_text(self, source: str, source_name: str | None = None) -> ReadResult:
        from .parser import parse_pdb

        return parse_pdb(source, source_name)


__all__ = ["PdbFormatHandler"]
