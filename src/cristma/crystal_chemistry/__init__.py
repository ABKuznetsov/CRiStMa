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
from .resolver import CoordinationShellResolver

__all__ = [
    "ComponentPairInterpretation",
    "ContactClassification",
    "CoordinationShell",
    "CoordinationShellResolver",
    "CrystalChemistryResolution",
    "EvidenceStatus",
    "ResolutionStatus",
    "ResolvedContact",
    "SecondaryEvidence",
    "ShellAlternative",
    "ShellResolutionPolicy",
]
