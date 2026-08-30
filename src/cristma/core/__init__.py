"""Canonical scientific data models."""

from .cell import UnitCell
from .values import MeasuredValue, MissingKind, parse_measured_value


def __getattr__(name: str) -> object:
    """Load legacy structure exports without cycling the public namespace."""

    if name in {"Crystal", "DisplacementParameters", "IndependentSite", "SiteComponent"}:
        from . import structure

        return getattr(structure, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    "Crystal",
    "DisplacementParameters",
    "IndependentSite",
    "MeasuredValue",
    "MissingKind",
    "SiteComponent",
    "UnitCell",
    "parse_measured_value",
]
