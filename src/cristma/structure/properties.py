"""Immutable typed per-atom property arrays."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True, slots=True)
class PropertyProvenance:
    """Origin of a property independently of the enclosing structure."""

    source_name: str | None = None
    source_field: str | None = None
    method: str | None = None


def _immutable_array(value: np.ndarray | object) -> np.ndarray:
    array = np.array(value, copy=True)
    array.flags.writeable = False
    return array


@dataclass(frozen=True, slots=True)
class AtomicProperty:
    """One named property whose leading dimension follows the atoms."""

    name: str
    values: np.ndarray
    unit: str | None = None
    missing: np.ndarray | None = None
    source_name: str | None = None
    provenance: PropertyProvenance = field(default_factory=PropertyProvenance)

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("atomic property name must not be empty")
        values = _immutable_array(self.values)
        if values.ndim == 0:
            raise ValueError("atomic property values must have a leading dimension")
        object.__setattr__(self, "values", values)
        if self.missing is not None:
            missing = _immutable_array(np.asarray(self.missing, dtype=bool))
            if missing.shape != values.shape:
                raise ValueError("missing mask shape must match property values")
            object.__setattr__(self, "missing", missing)


class AtomicPropertyTable(Mapping[str, AtomicProperty]):
    """Immutable mapping of property name to arrays for a fixed atom count."""

    __slots__ = ("_atom_count", "_properties")

    def __init__(
        self,
        atom_count: int,
        properties: tuple[AtomicProperty, ...] = (),
    ) -> None:
        if isinstance(atom_count, bool) or not isinstance(atom_count, int) or atom_count < 0:
            raise ValueError("atom_count must be a non-negative integer")
        mapped: dict[str, AtomicProperty] = {}
        for prop in properties:
            if prop.name in mapped:
                raise ValueError(f"duplicate atomic property: {prop.name!r}")
            if prop.values.shape[0] != atom_count:
                raise ValueError(
                    f"property {prop.name!r} leading dimension "
                    f"{prop.values.shape[0]} does not match atom count {atom_count}"
                )
            mapped[prop.name] = prop
        self._atom_count = atom_count
        self._properties = mapped

    @property
    def atom_count(self) -> int:
        return self._atom_count

    def __getitem__(self, name: str) -> AtomicProperty:
        return self._properties[name]

    def __iter__(self) -> Iterator[str]:
        return iter(self._properties)

    def __len__(self) -> int:
        return len(self._properties)


__all__ = ["AtomicProperty", "AtomicPropertyTable", "PropertyProvenance"]
