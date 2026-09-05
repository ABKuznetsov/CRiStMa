"""Stable diagnostics and invariant failures for diffraction calculations."""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType
from typing import Any


SEARCH_LIMIT_REACHED = "diffraction.reflections.search_limit_reached"
INCOMPATIBLE_CELL_AND_SETTING = "diffraction.reflections.incompatible_cell_and_setting"
NON_INTEGRAL_RECIPROCAL_ACTION = "diffraction.reflections.non_integral_reciprocal_action"
ORBIT_METRIC_MISMATCH = "diffraction.reflections.orbit_metric_mismatch"
INCONSISTENT_PHASE_BUCKETS = "diffraction.extinction.inconsistent_phase_buckets"


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return frozenset(_freeze(item) for item in value)
    return value


class DiffractionInvariantError(RuntimeError):
    """A mathematical contradiction in diffraction inputs or implementation."""

    def __init__(self, code: str, message: str, evidence: Mapping[str, Any]) -> None:
        if not code.strip():
            raise ValueError("invariant error code must not be empty")
        if not message.strip():
            raise ValueError("invariant error message must not be empty")
        super().__init__(message)
        self.code = code
        self.evidence = _freeze(evidence)


__all__ = [
    "DiffractionInvariantError",
    "INCOMPATIBLE_CELL_AND_SETTING",
    "INCONSISTENT_PHASE_BUCKETS",
    "NON_INTEGRAL_RECIPROCAL_ACTION",
    "ORBIT_METRIC_MISMATCH",
    "SEARCH_LIMIT_REACHED",
]
