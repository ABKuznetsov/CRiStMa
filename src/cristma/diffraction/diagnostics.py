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
STRUCTURE_FACTOR_SETTING_MISMATCH = "diffraction.structure_factor.setting_mismatch"
STRUCTURE_FACTOR_CELL_MISMATCH = (
    "diffraction.structure_factor.cell_fingerprint_mismatch"
)
STRUCTURE_FACTOR_SYMMETRY_MISMATCH = "diffraction.structure_factor.symmetry_mismatch"
STRUCTURE_FACTOR_UNSUPPORTED_ANISOTROPIC_ADP = (
    "diffraction.structure_factor.unsupported_anisotropic_adp"
)
STRUCTURE_FACTOR_UNSUPPORTED_SPECIES = (
    "diffraction.structure_factor.unsupported_species"
)
STRUCTURE_FACTOR_EXTINCT_NONZERO = (
    "diffraction.structure_factor.extinct_reflection_nonzero"
)
STRUCTURE_FACTOR_ANISOTROPIC_ADP_APPROXIMATED = (
    "diffraction.structure_factor.anisotropic_adp_approximated_by_ueq"
)
POWDER_MISSING_FRIEDEL_MATE = "diffraction.powder.missing_friedel_mate"
POWDER_NONRECIPROCAL_FRIEDEL_LINK = (
    "diffraction.powder.nonreciprocal_friedel_link"
)
POWDER_FRIEDEL_D_SPACING_MISMATCH = (
    "diffraction.powder.friedel_d_spacing_mismatch"
)


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
    "POWDER_FRIEDEL_D_SPACING_MISMATCH",
    "POWDER_MISSING_FRIEDEL_MATE",
    "POWDER_NONRECIPROCAL_FRIEDEL_LINK",
    "SEARCH_LIMIT_REACHED",
    "STRUCTURE_FACTOR_CELL_MISMATCH",
    "STRUCTURE_FACTOR_ANISOTROPIC_ADP_APPROXIMATED",
    "STRUCTURE_FACTOR_EXTINCT_NONZERO",
    "STRUCTURE_FACTOR_SETTING_MISMATCH",
    "STRUCTURE_FACTOR_SYMMETRY_MISMATCH",
    "STRUCTURE_FACTOR_UNSUPPORTED_ANISOTROPIC_ADP",
    "STRUCTURE_FACTOR_UNSUPPORTED_SPECIES",
]
