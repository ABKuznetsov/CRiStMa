"""Catalog-backed structural crystallography tools."""

from .space_group import SpaceGroupSetting
from .wyckoff import AffineCoordinateMap, WyckoffPosition
from .catalog import SpaceGroupCatalog
from .orbit import (
    CrystallographicOrbit,
    SiteSymmetry,
    WyckoffAssignment,
    assign_wyckoff,
    build_orbit,
)
from .local_geometry import GeometricContact, geometric_contacts
from .symmetry_context import (
    DEFAULT_METRIC_TOLERANCE,
    DirectBasisConvention,
    SymmetryContext,
    SymmetryContextInvariantError,
    SymmetrySourceKind,
    canonical_operation_key,
)

__all__ = [
    "AffineCoordinateMap",
    "DEFAULT_METRIC_TOLERANCE",
    "DirectBasisConvention",
    "SpaceGroupSetting",
    "SpaceGroupCatalog",
    "SymmetryContext",
    "SymmetryContextInvariantError",
    "SymmetrySourceKind",
    "CrystallographicOrbit",
    "GeometricContact",
    "SiteSymmetry",
    "WyckoffAssignment",
    "WyckoffPosition",
    "assign_wyckoff",
    "build_orbit",
    "canonical_operation_key",
    "geometric_contacts",
]
