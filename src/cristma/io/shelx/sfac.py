"""Scientific interpretation of SHELX SFAC element declarations."""

from __future__ import annotations

from dataclasses import dataclass
import re

from cristma.chemistry.elements import normalize_element
from cristma.chemistry.species import ElementSpecies
from cristma.core.values import MeasuredValue, parse_measured_value

from .document import ShelxRecord
from .records import ShelxSfacInstruction


@dataclass(frozen=True, slots=True)
class ShelxScatteringEntry:
    """One indexed SHELX scattering entry with retained source label."""

    source_label: str
    species: ElementSpecies
    coefficients: tuple[MeasuredValue, ...] = ()


def _is_number(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True


def _entry(
    label: str,
    coefficients: tuple[str, ...],
    *,
    line: int,
) -> ShelxScatteringEntry:
    normalized_source = label.lstrip("$")
    match = re.fullmatch(r"([A-Za-z]{1,2})(?:(?:\d+)?[+-])?", normalized_source)
    if match is None:
        raise ValueError(f"invalid SFAC label {label!r} at line {line}")
    try:
        element = normalize_element(match.group(1))
    except ValueError as error:
        raise ValueError(f"invalid SFAC label {label!r} at line {line}") from error
    return ShelxScatteringEntry(
        source_label=label,
        species=ElementSpecies(element),
        coefficients=tuple(parse_measured_value(value) for value in coefficients),
    )


def extract_sfac_entries(records: tuple[ShelxRecord, ...]) -> tuple[ShelxScatteringEntry, ...]:
    """Collect indexed elements from list and coefficient SFAC forms."""

    result: list[ShelxScatteringEntry] = []
    for record in records:
        if not isinstance(record, ShelxSfacInstruction):
            continue
        fields = record.entries
        line = record.span.start.line
        if len(fields) > 1 and any(_is_number(value) for value in fields[1:]):
            result.append(_entry(fields[0], fields[1:], line=line))
            continue
        result.extend(_entry(label, (), line=line) for label in fields)
    return tuple(result)


__all__ = ["ShelxScatteringEntry", "extract_sfac_entries"]
