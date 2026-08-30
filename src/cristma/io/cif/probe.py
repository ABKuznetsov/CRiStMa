"""Lightweight CIF content recognition without parser or mapper imports."""

from __future__ import annotations

from .lexer import lex_cif
from .tokens import CifTokenKind


def probe_cif(source: str) -> float:
    tokens, _diagnostics = lex_cif(source)
    significant = tuple(token for token in tokens if token.kind is not CifTokenKind.COMMENT)
    if significant and significant[0].kind is CifTokenKind.DATA:
        return 1.0
    if any(
        token.kind is CifTokenKind.TAG
        and token.value.casefold().startswith(("_cell_", "_atom_site_"))
        for token in significant
    ):
        return 0.8
    return 0.0


__all__ = ["probe_cif"]
