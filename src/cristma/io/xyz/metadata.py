"""Dependency-free extXYZ comment metadata parsing."""

from __future__ import annotations

from dataclasses import dataclass
import re
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .document import XyzPropertySpec


_INTEGER = re.compile(r"^[+-]?\d+$")
_REAL = re.compile(r"^[+-]?(?:(?:\d+\.\d*)|(?:\.\d+)|(?:\d+))(?:[EeDd][+-]?\d+)?$")
_TRUE = frozenset({"T", "TRUE"})
_FALSE = frozenset({"F", "FALSE"})


@dataclass(frozen=True, slots=True)
class XyzMetadata:
    """Typed structural and arbitrary entries from one XYZ comment line."""

    values: Mapping[str, object]
    schema: tuple[XyzPropertySpec, ...] = ()
    lattice: np.ndarray | None = None
    pbc: tuple[bool, bool, bool] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "values", MappingProxyType(dict(self.values)))
        if self.lattice is not None:
            lattice = np.array(self.lattice, dtype=float, copy=True)
            lattice.flags.writeable = False
            object.__setattr__(self, "lattice", lattice)


def parse_property_schema(value: str) -> tuple[XyzPropertySpec, ...]:
    """Parse the extXYZ colon-separated property triplets."""

    parts = value.split(":")
    if not parts or len(parts) % 3:
        raise ValueError("XYZ Properties must contain complete triplets")
    result: list[XyzPropertySpec] = []
    names: set[str] = set()
    for offset in range(0, len(parts), 3):
        name, kind, width_text = parts[offset : offset + 3]
        if name in names:
            raise ValueError(f"duplicate XYZ property name: {name!r}")
        try:
            width = int(width_text)
        except ValueError as exc:
            raise ValueError("XYZ property width must be an integer") from exc
        result.append(XyzPropertySpec(name, kind, width))
        names.add(name)
    return tuple(result)


def _scan_entries(comment: str) -> list[tuple[str, str, bool]] | None:
    entries: list[tuple[str, str, bool]] = []
    index = 0
    while index < len(comment):
        while index < len(comment) and comment[index].isspace():
            index += 1
        if index == len(comment):
            break
        start = index
        while index < len(comment) and not comment[index].isspace() and comment[index] != "=":
            index += 1
        key = comment[start:index]
        while index < len(comment) and comment[index].isspace():
            index += 1
        if not key or index == len(comment) or comment[index] != "=":
            if not entries:
                return None
            raise ValueError(f"malformed extXYZ metadata near {comment[start:]!r}")
        index += 1
        while index < len(comment) and comment[index].isspace():
            index += 1
        quoted = index < len(comment) and comment[index] == '"'
        if quoted:
            index += 1
            characters: list[str] = []
            while index < len(comment):
                character = comment[index]
                if character == "\\":
                    index += 1
                    if index == len(comment):
                        raise ValueError("unterminated escape in extXYZ metadata")
                    characters.append(comment[index])
                    index += 1
                elif character == '"':
                    index += 1
                    break
                else:
                    characters.append(character)
                    index += 1
            else:
                raise ValueError("unterminated quoted extXYZ metadata value")
            value = "".join(characters)
        else:
            value_start = index
            while index < len(comment) and not comment[index].isspace():
                index += 1
            value = comment[value_start:index]
            if not value:
                raise ValueError(f"missing extXYZ metadata value for {key!r}")
        entries.append((key, value, quoted))
    return entries


def _logical(value: str) -> bool:
    normalized = value.upper()
    if normalized in _TRUE:
        return True
    if normalized in _FALSE:
        return False
    raise ValueError(f"invalid logical value: {value!r}")


def _typed_value(value: str, *, quoted: bool) -> object:
    normalized = value.upper()
    if normalized in _TRUE | _FALSE:
        return _logical(value)
    if _INTEGER.fullmatch(value):
        return int(value)
    if _REAL.fullmatch(value):
        return float(value.replace("D", "E").replace("d", "e"))
    if quoted:
        items = value.split()
        if len(items) > 1 and all(item.upper() in _TRUE | _FALSE for item in items):
            return tuple(_logical(item) for item in items)
        if len(items) > 1 and all(_INTEGER.fullmatch(item) for item in items):
            return tuple(int(item) for item in items)
        if len(items) > 1 and all(_REAL.fullmatch(item) for item in items):
            return tuple(float(item.replace("D", "E").replace("d", "e")) for item in items)
    return value


def parse_xyz_metadata(comment: str) -> XyzMetadata:
    """Parse one plain XYZ or extXYZ comment without guessing periodicity."""

    entries = _scan_entries(comment)
    if entries is None:
        return XyzMetadata({"comment": comment}) if comment else XyzMetadata({})
    seen: set[str] = set()
    arbitrary: dict[str, object] = {}
    schema: tuple[XyzPropertySpec, ...] = ()
    lattice: np.ndarray | None = None
    pbc: tuple[bool, bool, bool] | None = None
    for key, raw_value, quoted in entries:
        if key in seen:
            raise ValueError(f"duplicate extXYZ metadata key: {key!r}")
        seen.add(key)
        if key == "Properties":
            schema = parse_property_schema(raw_value)
        elif key == "Lattice":
            parts = raw_value.split()
            if len(parts) != 9:
                raise ValueError("Lattice must contain nine real values")
            try:
                lattice = np.asarray([float(item) for item in parts], dtype=float).reshape(3, 3)
            except ValueError as exc:
                raise ValueError("Lattice must contain nine real values") from exc
            if not np.isfinite(lattice).all() or abs(float(np.linalg.det(lattice))) <= 1e-15:
                raise ValueError("Lattice must be finite and non-singular")
        elif key == "pbc":
            parts = raw_value.split()
            if len(parts) != 3:
                raise ValueError("pbc must contain three logical values")
            try:
                pbc = tuple(_logical(item) for item in parts)  # type: ignore[assignment]
            except ValueError as exc:
                raise ValueError("pbc must contain three logical values") from exc
        else:
            arbitrary[key] = _typed_value(raw_value, quoted=quoted)
    return XyzMetadata(arbitrary, schema, lattice, pbc)


__all__ = ["XyzMetadata", "parse_property_schema", "parse_xyz_metadata"]
