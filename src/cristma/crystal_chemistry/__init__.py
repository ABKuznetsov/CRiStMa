"""Reusable inorganic crystal-chemistry calculations and results."""

from .contacts import (
    ComponentPairInterpretation,
    ContactClassification,
    CoordinationShell,
    CrystalChemistryResolution,
    EvidenceStatus,
    ResolutionStatus,
    ResolvedContact,
    ResolvedContactOrbit,
    SecondaryEvidence,
    ShellAlternative,
)
from .policy import ShellResolutionPolicy
from .polyhedra import (
    CoordinationPolyhedron,
    CoordinationPolyhedronOrbit,
    FaceSignature,
    PolyhedronBuildResult,
    PolyhedronBuilder,
    PolyhedronVertex,
    canonical_face_signature,
    polyhedron_face_signature,
    unique_hull_edges,
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
from .representation import (
    StructuralRepresentation,
    StructuralRepresentationBuilder,
    StructuralSelectionPolicy,
)
from .periodic_connectivity import (
    PeriodicComponent,
    PeriodicConnectivityAnalyzer,
    PeriodicConnectivityResult,
)
from .structural_blocks import (
    StructuralBlock,
    StructuralBlockClassification,
    StructuralBlockFinder,
    StructuralBlockResult,
)
from .ring_finder import RingFinder
from .rings import (
    PeriodicUnitRef,
    RingAnalysisResult,
    RingAnalysisStatus,
    RingSearchPolicy,
    StructuralRing,
    StructuralRingOrbit,
    StructuralRingScope,
)

__all__ = [
    "ComponentPairInterpretation",
    "ContactClassification",
    "CoordinationShell",
    "CoordinationShellResolver",
    "CoordinationPolyhedron",
    "CoordinationPolyhedronOrbit",
    "CrystalChemistryResolution",
    "EvidenceStatus",
    "ResolutionStatus",
    "ResolvedContact",
    "ResolvedContactOrbit",
    "PolyhedronBuildResult",
    "PolyhedronBuilder",
    "PolyhedronVertex",
    "FaceSignature",
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
    "StructuralRepresentation",
    "StructuralRepresentationBuilder",
    "StructuralSelectionPolicy",
    "PeriodicComponent",
    "PeriodicConnectivityAnalyzer",
    "PeriodicConnectivityResult",
    "StructuralBlock",
    "StructuralBlockClassification",
    "StructuralBlockFinder",
    "StructuralBlockResult",
    "PeriodicUnitRef",
    "RingAnalysisResult",
    "RingAnalysisStatus",
    "RingFinder",
    "RingSearchPolicy",
    "StructuralRing",
    "StructuralRingOrbit",
    "StructuralRingScope",
    "polyhedron_face_signature",
    "canonical_face_signature",
    "unique_hull_edges",
]
