"""Scientific values that retain reported uncertainty and missing semantics."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import re


class MissingKind(str, Enum):
    """Meaning of a reported or absent scientific value."""

    PRESENT = "present"
    ABSENT = "absent"
    UNKNOWN = "unknown"
    INAPPLICABLE = "inapplicable"


@dataclass(frozen=True, slots=True)
class MeasuredValue:
    """A normalized number together with its original reported form."""

    value: float | None
    uncertainty: float | None
    raw: str | None
    unit: str | None = None
    missing: MissingKind = MissingKind.PRESENT


_NUMBER = re.compile(
    r"^(?P<number>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    r"(?:\((?P<uncertainty>\d+)\))?$"
)


def parse_measured_value(
    token: str | None,
    *,
    unit: str | None = None,
) -> MeasuredValue:
    """Parse a CIF-style number without discarding its reported token."""

    if token is None:
        return MeasuredValue(None, None, None, unit, MissingKind.ABSENT)

    raw = str(token).strip()
    if raw == "?":
        return MeasuredValue(None, None, raw, unit, MissingKind.UNKNOWN)
    if raw == ".":
        return MeasuredValue(None, None, raw, unit, MissingKind.INAPPLICABLE)

    match = _NUMBER.fullmatch(raw)
    if match is None:
        raise ValueError(f"Invalid measured number: {raw!r}")

    number_token = match.group("number")
    value = float(number_token)
    uncertainty_digits = match.group("uncertainty")
    uncertainty = None

    if uncertainty_digits is not None:
        mantissa, exponent_token = re.split(r"[eE]", number_token, maxsplit=1) \
            if "e" in number_token.lower() else (number_token, "0")
        decimals = len(mantissa.partition(".")[2])
        exponent = int(exponent_token)
        uncertainty = int(uncertainty_digits) * 10.0 ** (exponent - decimals)

    return MeasuredValue(value, uncertainty, raw, unit, MissingKind.PRESENT)
