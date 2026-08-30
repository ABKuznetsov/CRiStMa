import pytest

from cristma.core.values import MeasuredValue

from cristma.structure import (
    FrameReference,
    MolecularAtom,
    MolecularStructure,
    SiteComponent,
    StructureCollection,
    StructureSequence,
)


def molecule(name: str) -> MolecularStructure:
    return MolecularStructure(
        name,
        (MolecularAtom(f"atom:{name}", name, (SiteComponent("H", MeasuredValue(1.0, None, "1")),), (0.0, 0.0, 0.0)),),
    )


def test_sequence_does_not_materialize_frames_until_requested() -> None:
    loaded: list[int] = []
    refs = (FrameReference(0), FrameReference(1, role="final"))

    def load(reference: FrameReference) -> MolecularStructure:
        loaded.append(reference.index)
        return molecule(str(reference.index))

    sequence = StructureSequence(refs, load)

    assert loaded == []
    assert sequence.final.name == "1"
    assert loaded == [1]
    assert sequence.final is sequence.final
    assert loaded == [1]


def test_sequence_slice_materializes_only_selected_frames() -> None:
    loaded: list[int] = []
    refs = tuple(FrameReference(index) for index in range(3))

    def load(reference: FrameReference) -> MolecularStructure:
        loaded.append(reference.index)
        return molecule(str(reference.index))

    selected = StructureSequence(refs, load)[1:]

    assert isinstance(selected, StructureCollection)
    assert [structure.name for structure in selected] == ["1", "2"]
    assert loaded == [1, 2]


def test_failed_frame_load_is_not_cached() -> None:
    attempts = 0

    def load(reference: FrameReference) -> MolecularStructure:
        nonlocal attempts
        attempts += 1
        raise ValueError("broken frame")

    sequence = StructureSequence((FrameReference(0),), load)

    with pytest.raises(ValueError, match="broken frame"):
        sequence[0]
    with pytest.raises(ValueError, match="broken frame"):
        sequence[0]
    assert attempts == 2
