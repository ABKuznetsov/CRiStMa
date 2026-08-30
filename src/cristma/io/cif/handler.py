"""Format-registry adapter for native CIF documents."""

from __future__ import annotations

from dataclasses import replace

from cristma.io.result import ReadResult

from .lexer import lex_cif
from .mapper import map_cif_structures
from .parser import parse_cif
from .tokens import CifTokenKind


class CifFormatHandler:
    """Recognize, parse, and scientifically map CIF 1.1 text."""

    name = "cif"
    suffixes = (".cif",)

    def probe(self, source: str) -> float:
        tokens, _diagnostics = lex_cif(source)
        significant = tuple(
            token for token in tokens if token.kind is not CifTokenKind.COMMENT
        )
        if significant and significant[0].kind is CifTokenKind.DATA:
            return 1.0
        if any(
            token.kind is CifTokenKind.TAG
            and token.value.casefold().startswith(("_cell_", "_atom_site_"))
            for token in significant
        ):
            return 0.8
        return 0.0

    def read_text(
        self,
        source: str,
        source_name: str | None = None,
    ) -> ReadResult:
        parsed = parse_cif(source, source_name=source_name)
        structures, mapping_diagnostics = map_cif_structures(parsed.document)
        return replace(
            parsed,
            structures=structures,
            diagnostics=(*parsed.diagnostics, *mapping_diagnostics),
        )
