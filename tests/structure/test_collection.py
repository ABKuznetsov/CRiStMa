from cristma.core.values import MeasuredValue
from cristma.structure import MolecularAtom, MolecularStructure, SiteComponent, StructureCollection


def molecule(name: str) -> MolecularStructure:
    return MolecularStructure(
        name,
        (MolecularAtom(f"atom:{name}", name, (SiteComponent("H", MeasuredValue(1.0, None, "1")),), (0.0, 0.0, 0.0)),),
    )


def test_collection_is_sequence_with_primary_and_final_roles() -> None:
    first, last = molecule("first"), molecule("last")
    collection = StructureCollection.from_structures(
        (first, last),
        primary_index=0,
        final_index=1,
    )

    assert tuple(collection) == (first, last)
    assert collection.primary is first
    assert collection.final is last
