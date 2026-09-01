import numpy as np
import pytest

from cristma.io.xyz import index_xyz
from cristma.io.xyz.mapper import map_xyz_frame
from cristma.io.xyz.parser import load_xyz_frame
from cristma.structure import CrystalStructure, FrameReference, MolecularStructure, SourceReference


def _map(source: str, name: str = "sample.xyz"):
    document, diagnostics = index_xyz(source, name)
    assert not diagnostics
    span = document.frames[0]
    reference = FrameReference(
        index=0,
        source=SourceReference(name, "xyz", "frame:0", span.start_offset, span.end_offset),
    )
    return map_xyz_frame(load_xyz_frame(document, reference))


def test_plain_xyz_maps_to_molecule() -> None:
    structure = _map("2\nwater\nO 0 0 0\nH 0.96 0 0\n")

    assert isinstance(structure, MolecularStructure)
    assert structure.periodic == (False, False, False)
    assert structure.name == "water"
    assert structure.atoms[0].components[0].element == "O"
    assert structure.atoms[0].id == "atom:1"
    assert structure.provenance.source.source_name == "sample.xyz"


def test_lattice_without_pbc_remains_molecular() -> None:
    structure = _map(
        '1\nLattice="2 0 0 0 2 0 0 0 2" Properties=species:S:1:pos:R:3\nSi 0 0 0\n'
    )

    assert isinstance(structure, MolecularStructure)
    assert structure.metadata["xyz_lattice"] == (
        (2.0, 0.0, 0.0),
        (0.0, 2.0, 0.0),
        (0.0, 0.0, 2.0),
    )


def test_explicit_periodic_extxyz_maps_to_identity_crystal() -> None:
    structure = _map(
        '1\nLattice="0 2 0 -3 0 0 0 0 4" '
        'Properties=species:S:1:pos:R:3:forces:R:3 pbc="T T F"\n'
        'Si -1.5 0.5 3 9 8 7\n'
    )

    assert isinstance(structure, CrystalStructure)
    assert structure.periodic == (True, True, False)
    assert structure.space_group.provenance == "unreported_identity"
    assert np.allclose(structure.sites[0].fractional[0].value, 0.25)
    assert np.allclose(structure.atomic_view().cartesian, [[0.5, 1.5, 3.0]])
    assert structure.properties["forces"].values.tolist() == [[9.0, 8.0, 7.0]]


def test_nonstructural_columns_are_preserved_without_inferred_units() -> None:
    structure = _map(
        "2\nProperties=Z:I:1:pos:R:3:charge:R:1:tag:S:1\n"
        "14 0 0 0 -0.2 core\n8 1 0 0 0.4 shell\n"
    )

    assert [atom.components[0].element for atom in structure.atoms] == ["Si", "O"]
    assert set(structure.properties) == {"charge", "tag"}
    assert structure.properties["charge"].unit is None
    assert structure.properties["charge"].provenance.source_field == "Properties:charge"
    assert structure.properties["tag"].values.tolist() == ["core", "shell"]


def test_unknown_species_is_represented_explicitly() -> None:
    structure = _map("1\nunknown\nXx 0 0 0\n")

    assert structure.atoms[0].components[0].element is None
    assert structure.atoms[0].components[0].species.label == "Xx"


def test_true_pbc_without_lattice_is_rejected() -> None:
    source = '1\nProperties=species:S:1:pos:R:3 pbc="T F F"\nSi 0 0 0\n'
    document, _ = index_xyz(source)
    span = document.frames[0]
    reference = FrameReference(0, source=SourceReference(None, "xyz", "frame:0"))

    with pytest.raises(ValueError, match="Lattice"):
        load_xyz_frame(document, reference)
