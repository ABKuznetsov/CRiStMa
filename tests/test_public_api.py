from cristma import __version__
from cristma.geometry import CoordinationAnalyzer, NeighborFinder, PeriodicNeighborGraph
from cristma.structure import AtomicView, ExpandedAtom, PeriodicAtomRef
from cristma.symmetry import SymmetryImageProvenance, expand_structure


def test_public_package_has_semantic_version():
    assert __version__ == "0.1.0.dev0"


def test_structure_core_subpackages_export_intended_types() -> None:
    assert AtomicView and ExpandedAtom and PeriodicAtomRef
    assert SymmetryImageProvenance and expand_structure
    assert NeighborFinder and PeriodicNeighborGraph and CoordinationAnalyzer
