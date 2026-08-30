"""Common result objects returned by structure readers."""

from __future__ import annotations

from dataclasses import dataclass

from cristma.structure import Structure, StructureCollection, StructureSequence

from .diagnostics import Diagnostic, Severity


@dataclass(frozen=True, slots=True)
class SourceInfo:
    """Properties of the decoded source that matter for round trips."""

    name: str | None = None
    format: str | None = None
    encoding: str = "utf-8"
    newline: str = "\n"


@dataclass(frozen=True, slots=True)
class ReadResult:
    """Parsed document, mapped structures, and all emitted diagnostics."""

    document: object | None
    structures: StructureCollection | StructureSequence | tuple[Structure, ...] = ()
    diagnostics: tuple[Diagnostic, ...] = ()
    source_info: SourceInfo | None = None

    def __post_init__(self) -> None:
        if isinstance(self.structures, tuple):
            object.__setattr__(
                self,
                "structures",
                StructureCollection.from_structures(self.structures),
            )

    @property
    def ok(self) -> bool:
        """Return true when the result contains no error diagnostics."""

        return not any(item.severity is Severity.ERROR for item in self.diagnostics)
