import pytest

from cristma.structure import (
    MolecularAtom,
    MolecularBond,
    MolecularGroup,
    MolecularStructure,
    Structure,
)


def test_molecule_has_no_artificial_periodic_cell() -> None:
    atoms = (
        MolecularAtom("atom:C1", "C1", "C", (0.0, 0.0, 0.0)),
        MolecularAtom("atom:O1", "O1", "O", (1.2, 0.0, 0.0)),
    )
    molecule = MolecularStructure(
        name="CO",
        atoms=atoms,
        bonds=(MolecularBond("bond:1", "atom:C1", "atom:O1", order=2.0),),
    )

    assert isinstance(molecule, Structure)
    assert molecule.cell is None
    assert molecule.periodic == (False, False, False)


def test_bond_rejects_missing_atom_identity() -> None:
    atom = MolecularAtom("atom:C1", "C1", "C", (0.0, 0.0, 0.0))

    with pytest.raises(ValueError, match="unknown atom"):
        MolecularStructure(
            "bad",
            atoms=(atom,),
            bonds=(MolecularBond("bond:1", "atom:C1", "atom:X", 1.0),),
        )


def test_group_references_existing_atoms_without_implying_rigidity() -> None:
    molecule = MolecularStructure(
        "fragment",
        atoms=(MolecularAtom("atom:C1", "C1", "C", (0.0, 0.0, 0.0)),),
        groups=(MolecularGroup("group:1", "fragment", ("atom:C1",)),),
    )

    assert molecule.groups[0].atom_ids == ("atom:C1",)
    assert molecule.groups[0].rigid is False
