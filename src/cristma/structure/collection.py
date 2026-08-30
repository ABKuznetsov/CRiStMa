"""Finite structure collections and indexed lazy frame sequences."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass, field
from threading import RLock
from types import MappingProxyType
from typing import Literal, Mapping, Protocol, overload, runtime_checkable

from .identity import SourceReference
from .molecular import Structure


StructureRole = Literal["model", "primary", "intermediate", "final"]
_ROLES = frozenset({"model", "primary", "intermediate", "final"})


def _validate_role(role: str) -> None:
    if role not in _ROLES:
        raise ValueError(f"unknown structure role: {role!r}")


@dataclass(frozen=True, slots=True)
class StructureEntry:
    """One eagerly available structure and its role in a source."""

    structure: Structure
    role: StructureRole = "model"
    source_index: int | None = None
    source: SourceReference | None = None

    def __post_init__(self) -> None:
        _validate_role(self.role)


@dataclass(frozen=True, slots=True)
class FrameReference:
    """Cheap index record from which a structure can be loaded later."""

    index: int
    role: StructureRole = "intermediate"
    source: SourceReference | None = None
    metadata: Mapping[str, object] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        if isinstance(self.index, bool) or not isinstance(self.index, int) or self.index < 0:
            raise ValueError("frame index must be a non-negative integer")
        _validate_role(self.role)
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))


@runtime_checkable
class StructureSeries(Protocol):
    """Sequence-like contract shared by eager and lazy multi-structure data."""

    def __len__(self) -> int: ...

    def __getitem__(self, index: int) -> Structure: ...

    def __iter__(self) -> Iterator[Structure]: ...

    @property
    def primary(self) -> Structure | None: ...

    @property
    def final(self) -> Structure | None: ...


class StructureCollection(Sequence[Structure]):
    """Immutable finite collection with explicit scientific roles."""

    __slots__ = ("_entries",)

    def __init__(self, entries: tuple[StructureEntry, ...] = ()) -> None:
        primary_count = sum(entry.role == "primary" for entry in entries)
        final_count = sum(entry.role == "final" for entry in entries)
        if primary_count > 1:
            raise ValueError("collection may contain at most one primary structure")
        if final_count > 1:
            raise ValueError("collection may contain at most one final structure")
        self._entries = tuple(entries)

    @classmethod
    def from_structures(
        cls,
        structures: Sequence[Structure],
        *,
        primary_index: int | None = None,
        final_index: int | None = None,
    ) -> StructureCollection:
        values = tuple(structures)
        for name, index in (("primary_index", primary_index), ("final_index", final_index)):
            if index is not None and not 0 <= index < len(values):
                raise IndexError(f"{name} is outside the structure collection")
        if primary_index is not None and primary_index == final_index:
            raise ValueError("one structure cannot have both primary and final roles")
        entries = tuple(
            StructureEntry(
                structure,
                role=(
                    "primary"
                    if index == primary_index
                    else "final"
                    if index == final_index
                    else "model"
                ),
                source_index=index,
            )
            for index, structure in enumerate(values)
        )
        return cls(entries)

    @property
    def entries(self) -> tuple[StructureEntry, ...]:
        return self._entries

    @property
    def primary(self) -> Structure | None:
        return next((entry.structure for entry in self._entries if entry.role == "primary"), None)

    @property
    def final(self) -> Structure | None:
        return next((entry.structure for entry in self._entries if entry.role == "final"), None)

    def __len__(self) -> int:
        return len(self._entries)

    @overload
    def __getitem__(self, index: int) -> Structure: ...

    @overload
    def __getitem__(self, index: slice) -> StructureCollection: ...

    def __getitem__(self, index: int | slice) -> Structure | StructureCollection:
        if isinstance(index, slice):
            return StructureCollection(self._entries[index])
        return self._entries[index].structure

    def __iter__(self) -> Iterator[Structure]:
        return (entry.structure for entry in self._entries)


class StructureSequence(Sequence[Structure]):
    """Indexed lazy structures, caching each successfully loaded frame once."""

    __slots__ = ("_references", "_loader", "_cache", "_lock")

    def __init__(
        self,
        references: Sequence[FrameReference],
        loader: Callable[[FrameReference], Structure],
    ) -> None:
        self._references = tuple(references)
        self._loader = loader
        self._cache: dict[int, Structure] = {}
        self._lock = RLock()

    @property
    def references(self) -> tuple[FrameReference, ...]:
        return self._references

    def _load_position(self, position: int) -> Structure:
        normalized = position if position >= 0 else len(self._references) + position
        reference = self._references[normalized]
        with self._lock:
            if normalized not in self._cache:
                self._cache[normalized] = self._loader(reference)
            return self._cache[normalized]

    @property
    def primary(self) -> Structure | None:
        position = next(
            (index for index, reference in enumerate(self._references) if reference.role == "primary"),
            None,
        )
        return None if position is None else self._load_position(position)

    @property
    def final(self) -> Structure | None:
        position = next(
            (index for index, reference in enumerate(self._references) if reference.role == "final"),
            None,
        )
        return None if position is None else self._load_position(position)

    def __len__(self) -> int:
        return len(self._references)

    @overload
    def __getitem__(self, index: int) -> Structure: ...

    @overload
    def __getitem__(self, index: slice) -> StructureCollection: ...

    def __getitem__(self, index: int | slice) -> Structure | StructureCollection:
        if isinstance(index, slice):
            positions = range(*index.indices(len(self)))
            return StructureCollection(
                tuple(
                    StructureEntry(
                        self._load_position(position),
                        role=self._references[position].role,
                        source_index=self._references[position].index,
                        source=self._references[position].source,
                    )
                    for position in positions
                )
            )
        return self._load_position(index)


__all__ = [
    "FrameReference",
    "StructureCollection",
    "StructureEntry",
    "StructureRole",
    "StructureSequence",
    "StructureSeries",
]
