"""Scientific evidence records produced by Chemistry."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class ChemicalEvidence:
    code: str
    message: str
    elements: tuple[str, ...] = ()
    reference_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.code or not self.message:
            raise ValueError("chemical evidence requires code and message")


__all__ = ["ChemicalEvidence"]
