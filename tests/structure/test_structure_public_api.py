def test_stable_structure_namespace_is_importable() -> None:
    from cristma.structure import (
        AtomicView,
        CrystalStructure,
        MolecularStructure,
        StructureCollection,
        StructureSequence,
    )

    assert all(
        value is not None
        for value in (
            AtomicView,
            CrystalStructure,
            MolecularStructure,
            StructureCollection,
            StructureSequence,
        )
    )
