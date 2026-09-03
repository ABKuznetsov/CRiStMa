from pathlib import Path

import cristma
from cristma import __version__
from cristma.io.shelx import ShelxDocument, ShelxOccupancyExpression, ShelxWriteOptions
from cristma.geometry import CoordinationAnalyzer, NeighborFinder, PeriodicNeighborGraph
from cristma.chemistry import ChemistryAnalyzer, Composition
from cristma.reference_data import ReferenceData
from cristma.structure import AtomicView, ExpandedAtom, PeriodicAtomRef
from cristma.symmetry import SymmetryImageProvenance, expand_structure


def test_public_package_has_semantic_version():
    assert __version__ == "0.1.0.dev0"


def test_structure_core_subpackages_export_intended_types() -> None:
    assert AtomicView and ExpandedAtom and PeriodicAtomRef
    assert SymmetryImageProvenance and expand_structure
    assert NeighborFinder and PeriodicNeighborGraph and CoordinationAnalyzer


def test_chemistry_subpackages_export_intended_tools() -> None:
    assert Composition and ChemistryAnalyzer and ReferenceData


def test_shelx_package_exports_intended_format_controls() -> None:
    assert ShelxDocument and ShelxOccupancyExpression and ShelxWriteOptions


def test_public_read_maps_poscar(tmp_path: Path) -> None:
    path = tmp_path / "POSCAR"
    path.write_text(
        "Silicon\n1\n1 0 0\n0 1 0\n0 0 1\nSi\n1\nDirect\n0 0 0\n",
        encoding="utf-8",
    )

    result = cristma.read(path)

    assert result.ok
    assert result.structures[0].name == "Silicon"


def test_structure_formats_exposes_immutable_read_capabilities() -> None:
    formats = cristma.structure_formats()

    assert isinstance(formats, tuple)
    assert {item.name for item in formats} >= {"cif", "shelx", "vasp", "xyz", "pdb"}
    assert next(item for item in formats if item.name == "pdb").suffixes == (".pdb",)
    vasp = next(item for item in formats if item.name == "vasp")
    assert {"POSCAR", "CONTCAR", "XDATCAR", "OUTCAR", "vasprun.xml"} <= set(
        vasp.basenames
    )
