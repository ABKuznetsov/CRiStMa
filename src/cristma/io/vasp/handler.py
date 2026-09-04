"""Format-registry adapter for native VASP structure sources."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from cristma.io.result import ReadResult

from .probe import probe_vasp


class VaspFormatHandler:
    """Dispatch VASP structural sources behind one application-neutral API."""

    name = "vasp"
    suffixes = (".vasp", ".xml")

    def probe(self, source: str) -> float:
        return probe_vasp(source)

    def read_text(self, source: str, source_name: str | None = None) -> ReadResult:
        basename = Path(source_name).name.casefold() if source_name else ""
        prefix = source[:100_000]
        if basename == "xdatcar" or "Direct configuration=" in prefix:
            from .xdatcar import parse_xdatcar

            return parse_xdatcar(source, source_name)
        if basename == "outcar" or ("POSITION" in prefix and "TOTAL-FORCE" in prefix):
            from .outcar import parse_outcar

            return parse_outcar(source, source_name)
        if basename == "vasprun.xml" or "<modeling" in prefix:
            from .vasprun import parse_vasprun

            return parse_vasprun(source, source_name)

        from .mapper import map_vasp_snapshot
        from .poscar import parse_poscar, poscar_snapshot

        parsed = parse_poscar(source, source_name)
        if parsed.document.header is None or not parsed.ok:
            return parsed
        try:
            structure = map_vasp_snapshot(poscar_snapshot(parsed.document))
        except ValueError as error:
            from cristma.diagnostics import Diagnostic, Severity

            return replace(
                parsed,
                diagnostics=(
                    *parsed.diagnostics,
                    Diagnostic(Severity.ERROR, "vasp.map.invalid_structure", str(error)),
                ),
            )
        return replace(parsed, structures=(structure,))


__all__ = ["VaspFormatHandler"]
