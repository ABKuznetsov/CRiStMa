"""Lazy structural reader for VASP vasprun.xml sources."""

from __future__ import annotations

import xml.etree.ElementTree as ET
import xml.parsers.expat

import numpy as np

from cristma.chemistry import UnknownSpecies
from cristma.diagnostics import Diagnostic, Severity
from cristma.io.result import ReadResult, SourceInfo
from cristma.structure import FrameReference, SourceReference, StructureSequence

from .document import VaspFrameSpan, VaspSnapshot, VasprunDocument


def _local(name: str) -> str:
    return name.rsplit("}", 1)[-1].rsplit(":", 1)[-1]


def _character_boundaries(source: str) -> dict[int, int]:
    boundaries = {0: 0}
    byte_offset = 0
    for char_offset, character in enumerate(source, start=1):
        byte_offset += len(character.encode("utf-8"))
        boundaries[byte_offset] = char_offset
    return boundaries


def _closed_element_end(data: bytes, close_start: int) -> int:
    end = data.find(b">", close_start)
    if end < 0:
        raise ValueError("XML closing tag is incomplete")
    return end + 1


def _extract_species(fragment: str) -> tuple[str, ...]:
    root = ET.fromstring(fragment)
    atoms_array = next(
        (
            element
            for element in root.iter()
            if _local(element.tag) == "array" and element.attrib.get("name") == "atoms"
        ),
        None,
    )
    if atoms_array is None:
        return ()
    species = []
    for row in atoms_array.iter():
        if _local(row.tag) != "rc":
            continue
        columns = [child for child in row if _local(child.tag) == "c"]
        if columns and columns[0].text and columns[0].text.strip():
            species.append(columns[0].text.strip())
    return tuple(species)


def parse_vasprun(source: str, source_name: str | None = None) -> ReadResult:
    """Index complete XML structural frames without mapping their arrays."""

    data = source.encode("utf-8")
    boundaries = _character_boundaries(source)
    calculations: list[tuple[int, int]] = []
    named_structures: list[tuple[int, int, str]] = []
    atominfo_span: tuple[int, int] | None = None
    stack: list[tuple[str, int, dict[str, str]]] = []
    diagnostics: list[Diagnostic] = []
    parser = xml.parsers.expat.ParserCreate()

    def start(name: str, attributes: dict[str, str]) -> None:
        stack.append((_local(name), parser.CurrentByteIndex, dict(attributes)))

    def end(name: str) -> None:
        nonlocal atominfo_span
        local = _local(name)
        opened_local, start_byte, attributes = stack.pop()
        if opened_local != local:
            return
        end_byte = _closed_element_end(data, parser.CurrentByteIndex)
        if local == "calculation":
            calculations.append((start_byte, end_byte))
        elif local == "atominfo":
            atominfo_span = (start_byte, end_byte)
        elif local == "structure" and attributes.get("name"):
            named_structures.append((start_byte, end_byte, attributes["name"]))

    parser.StartElementHandler = start
    parser.EndElementHandler = end
    try:
        parser.Parse(data, True)
    except (xml.parsers.expat.ExpatError, ValueError) as error:
        diagnostics.append(
            Diagnostic(
                Severity.WARNING,
                "vasp.vasprun.xml_incomplete",
                f"XML ended before a complete document: {error}",
            )
        )

    species: tuple[str, ...] = ()
    if atominfo_span is not None:
        start_char = boundaries[atominfo_span[0]]
        end_char = boundaries[atominfo_span[1]]
        try:
            species = _extract_species(source[start_char:end_char])
        except (ET.ParseError, ValueError) as error:
            diagnostics.append(
                Diagnostic(Severity.ERROR, "vasp.vasprun.atominfo_invalid", str(error))
            )
    if not species:
        diagnostics.append(
            Diagnostic(
                Severity.WARNING,
                "vasp.map.species_unresolved",
                "vasprun.xml contains no complete atom species table",
            )
        )

    selected: list[tuple[int, int, int, str]] = []
    if calculations:
        selected = [(*span, index + 1, "calculation") for index, span in enumerate(calculations)]
    else:
        selected = [
            (start_byte, end_byte, index + 1, name)
            for index, (start_byte, end_byte, name) in enumerate(named_structures)
        ]
    frames = tuple(
        VaspFrameSpan(
            index,
            boundaries[start_byte],
            boundaries[end_byte],
            reported_index,
        )
        for index, (start_byte, end_byte, reported_index, _kind) in enumerate(selected)
    )
    document = VasprunDocument(source, source_name, frames, species)
    references = tuple(
        FrameReference(
            frame.index,
            role="final" if index == len(frames) - 1 else "intermediate",
            source=SourceReference(
                source_name,
                "vasp-xml",
                f"{selected[index][3]}:{frame.reported_index}",
                frame.start_offset,
                frame.end_offset,
            ),
            metadata={"xml_kind": selected[index][3], "ionic_step": frame.reported_index},
        )
        for index, frame in enumerate(frames)
    )

    def load(reference: FrameReference):
        from . import mapper

        return mapper.map_vasp_snapshot(load_vasprun_snapshot(document, reference))

    newline = "\r\n" if "\r\n" in source else "\n"
    return ReadResult(
        document=document,
        structures=StructureSequence(references, load),
        diagnostics=tuple(diagnostics),
        source_info=SourceInfo(source_name, "vasp-xml", "utf-8", newline),
    )


def _named_varray(root: ET.Element, name: str) -> ET.Element | None:
    return next(
        (
            element
            for element in root.iter()
            if _local(element.tag) == "varray" and element.attrib.get("name") == name
        ),
        None,
    )


def _numeric_varray(root: ET.Element, name: str, *, required: bool) -> np.ndarray | None:
    varray = _named_varray(root, name)
    if varray is None:
        if required:
            raise ValueError(f"vasprun {name} array is missing")
        return None
    rows = []
    for vector in varray:
        if _local(vector.tag) != "v":
            continue
        try:
            row = tuple(float(value) for value in (vector.text or "").split())
        except ValueError as error:
            raise ValueError(f"vasprun {name} contains malformed numeric row") from error
        if len(row) != 3 or not np.isfinite(row).all():
            raise ValueError(f"vasprun {name} contains malformed numeric row")
        rows.append(row)
    return np.asarray(rows, dtype=float)


def load_vasprun_snapshot(
    document: VasprunDocument,
    reference: FrameReference,
) -> VaspSnapshot:
    """Parse one indexed calculation or named structure fragment."""

    frame = document.frames[reference.index]
    fragment = document.raw_source[frame.start_offset : frame.end_offset]
    try:
        root = ET.fromstring(fragment)
    except ET.ParseError as error:
        raise ValueError(f"indexed vasprun frame is not standalone XML: {error}") from error
    basis = _numeric_varray(root, "basis", required=True)
    positions = _numeric_varray(root, "positions", required=True)
    assert basis is not None and positions is not None
    if basis.shape != (3, 3):
        raise ValueError("vasprun basis must contain exactly three vectors")
    species = tuple(document.species_labels)
    if species and len(species) != positions.shape[0]:
        raise ValueError(
            f"vasprun position count {positions.shape[0]} does not match atominfo {len(species)}"
        )
    if not species:
        species = tuple(UnknownSpecies(f"vasp:atom:{index + 1}") for index in range(len(positions)))
    forces = _numeric_varray(root, "forces", required=False)
    if forces is not None and forces.shape != positions.shape:
        raise ValueError("vasprun forces count does not match positions")
    velocities = _numeric_varray(root, "velocities", required=False)
    if velocities is not None and velocities.shape != positions.shape:
        raise ValueError("vasprun velocities count does not match positions")
    kind = str(reference.metadata.get("xml_kind", "calculation"))
    reported = frame.reported_index if frame.reported_index is not None else frame.index + 1
    return VaspSnapshot(
        name=f"vasprun.xml — {kind} {reported}",
        lattice=basis,
        species=species,
        fractional=positions,
        frame_index=frame.index,
        source=reference.source
        or SourceReference(document.source_name, "vasp-xml", f"{kind}:{reported}"),
        velocities=velocities,
        velocity_mode="cartesian" if velocities is not None else None,
        velocity_unit="angstrom/fs" if velocities is not None else None,
        forces=forces,
        force_unit="eV/angstrom" if forces is not None else None,
    )


__all__ = ["load_vasprun_snapshot", "parse_vasprun"]
