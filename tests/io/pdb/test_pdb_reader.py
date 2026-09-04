from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import cristma
from cristma.structure import CrystalStructure, MolecularStructure


DATA = Path(__file__).parent / "data"


def test_pdb_with_cryst1_maps_to_crystal_structure() -> None:
    result = cristma.read(DATA / "crystal.pdb")

    assert result.ok
    structure = result.structures[0]
    assert isinstance(structure, CrystalStructure)
    assert tuple(
        value.value for value in (structure.cell.a, structure.cell.b, structure.cell.c)
    ) == pytest.approx((10.0, 11.0, 12.0))
    assert structure.space_group is not None
    assert structure.space_group.hm_symbol == "P 1"
    assert [site.label for site in structure.sites] == ["CA1", "O1"]
    assert [site.total_occupancy for site in structure.sites] == pytest.approx(
        [0.75, 1.0]
    )
    np.testing.assert_allclose(
        [[value.value for value in site.fractional] for site in structure.sites],
        [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]],
        rtol=0.0,
        atol=1e-12,
    )


def test_pdb_without_cryst1_maps_to_molecular_structure() -> None:
    result = cristma.read(DATA / "molecule.pdb")

    assert result.ok
    structure = result.structures[0]
    assert isinstance(structure, MolecularStructure)
    assert [atom.label for atom in structure.atoms] == ["C1", "O1"]
    assert [
        atom.components[0].species.require_element() for atom in structure.atoms
    ] == ["C", "O"]
    np.testing.assert_allclose(
        [atom.cartesian for atom in structure.atoms],
        [[0.0, 0.0, 0.0], [1.23, 0.0, 0.0]],
        rtol=0.0,
        atol=1e-12,
    )


def test_pdb_atom_name_fallback_respects_fixed_column_alignment() -> None:
    source = (
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00\n"
        "HETATM    2 CA   UNK A   2       1.000   0.000   0.000  1.00 10.00\n"
        "END\n"
    )

    result = cristma.read_text(source, format="pdb", source_name="aligned.pdb")

    assert result.ok
    assert [
        atom.components[0].species.require_element()
        for atom in result.structures[0].atoms
    ] == ["C", "Ca"]


def test_pdb_digit_prefixed_hydrogen_name_does_not_become_heavy_element() -> None:
    source = (
        "ATOM      1 1HG  CYS A   1       0.000   0.000   0.000  1.00 10.00\n"
        "ATOM      2 2HE  GLN A   1       1.000   0.000   0.000  1.00 10.00\n"
        "END\n"
    )

    result = cristma.read_text(source, format="pdb", source_name="hydrogen.pdb")

    assert result.ok
    assert [
        atom.components[0].species.require_element()
        for atom in result.structures[0].atoms
    ] == ["H", "H"]


def test_pdb_models_map_to_separate_canonical_structures() -> None:
    source = (
        "MODEL        1\n"
        "ATOM      1  C1  MOL A   1       0.000   0.000   0.000  1.00 10.00           C\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "ATOM      1  C1  MOL A   1       0.500   0.000   0.000  1.00 10.00           C\n"
        "ENDMDL\nEND\n"
    )

    result = cristma.read_text(source, format="pdb", source_name="trajectory.pdb")

    assert result.ok
    assert len(result.structures) == 2
    assert [item.id for item in result.structures] == ["model:1", "model:2"]
    assert [item.atoms[0].cartesian[0] for item in result.structures] == [0.0, 0.5]


def test_pdb_alternate_locations_remain_distinct_disorder_positions() -> None:
    source = (
        "ATOM      1  C1 AMOL A   1       0.000   0.000   0.000  0.60 10.00           C\n"
        "ATOM      2  C1 BMOL A   1       0.200   0.000   0.000  0.40 10.00           C\n"
        "END\n"
    )

    result = cristma.read_text(source, format="pdb", source_name="disorder.pdb")

    assert result.ok
    atoms = result.structures[0].atoms
    assert len(atoms) == 2
    assert [atom.metadata["alternate_location"] for atom in atoms] == ["A", "B"]
    assert [atom.components[0].occupancy.value for atom in atoms] == [0.6, 0.4]


def test_pdb_invalid_coordinate_is_an_error_without_partial_structure() -> None:
    source = (
        "ATOM      1  C1  MOL A   1       broken   0.000   0.000  1.00 10.00           C\n"
        "END\n"
    )

    result = cristma.read_text(source, format="pdb", source_name="broken.pdb")

    assert not result.ok
    assert not result.structures
    assert {item.code for item in result.diagnostics} == {
        "pdb.parse.value_invalid",
        "pdb.map.atoms_missing",
    }


def test_pdb_unresolved_space_group_is_explicit_identity_fallback() -> None:
    source = (
        "CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 NOT-A-GROUP   1\n"
        "ATOM      1  C1  MOL A   1       0.000   0.000   0.000  1.00 10.00           C\n"
        "END\n"
    )

    result = cristma.read_text(source, format="pdb", source_name="unknown-sg.pdb")

    assert result.ok
    assert result.structures[0].space_group.provenance == "identity_fallback"
    assert {item.code for item in result.diagnostics} == {
        "pdb.map.space_group_unresolved"
    }
