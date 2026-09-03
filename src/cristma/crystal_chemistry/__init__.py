"""Reusable inorganic crystal-chemistry calculations and results."""

from .contacts import (
    ComponentPairInterpretation,
    ContactClassification,
    CoordinationShell,
    CrystalChemistryResolution,
    EvidenceStatus,
    ResolutionStatus,
    ResolvedContact,
    SecondaryEvidence,
    ShellAlternative,
)
from .policy import ShellResolutionPolicy
from .polyhedra import (
    CoordinationPolyhedron,
    PolyhedronBuildResult,
    PolyhedronBuilder,
    polyhedron_face_signature,
)
from .resolver import CoordinationShellResolver
from .shannon_distance import ShannonDistanceCheck, ShannonDistanceValidator
from .structural_units import (
    StructuralUnit,
    StructuralUnitBuildResult,
    StructuralUnitBuilder,
    StructuralUnitKind,
)
from .structural_graph import (
    StructuralConnection,
    StructuralConnectionKind,
    StructuralGraphBuilder,
    StructuralUnitGraph,
)

__all__ = [
    "ComponentPairInterpretation",
    "ContactClassification",
    "CoordinationShell",
    "CoordinationShellResolver",
    "CoordinationPolyhedron",
    "CrystalChemistryResolution",
    "EvidenceStatus",
    "ResolutionStatus",
    "ResolvedContact",
    "PolyhedronBuildResult",
    "PolyhedronBuilder",
    "SecondaryEvidence",
    "ShellAlternative",
    "ShellResolutionPolicy",
    "ShannonDistanceCheck",
    "ShannonDistanceValidator",
    "StructuralUnit",
    "StructuralUnitBuildResult",
    "StructuralUnitBuilder",
    "StructuralUnitKind",
    "StructuralConnection",
    "StructuralConnectionKind",
    "StructuralGraphBuilder",
    "StructuralUnitGraph",
    "polyhedron_face_signature",
]
