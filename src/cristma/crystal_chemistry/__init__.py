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
    "polyhedron_face_signature",
]
