from cristma.io.result import ReadResult
from cristma.structure import StructureCollection, StructureSequence


def test_result_converts_legacy_structure_tuple_to_collection() -> None:
    result = ReadResult(document=None, structures=())

    assert isinstance(result.structures, StructureCollection)
    assert len(result.structures) == 0


def test_result_keeps_explicit_lazy_sequence() -> None:
    sequence = StructureSequence((), lambda reference: None)

    assert ReadResult(document=None, structures=sequence).structures is sequence
