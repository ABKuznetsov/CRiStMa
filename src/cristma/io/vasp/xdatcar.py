"""Indexed lazy reader for VASP XDATCAR trajectories."""

from __future__ import annotations

import re

import numpy as np

from cristma.chemistry import UnknownSpecies
from cristma.diagnostics import Diagnostic, Severity
from cristma.io.result import ReadResult, SourceInfo
from cristma.structure import FrameReference, SourceReference, StructureSequence

from .document import VaspFrameSpan, VaspHeader, VaspScale, VaspSnapshot, XdatcarDocument
from .numeric import scaled_lattice


_CONFIGURATION = re.compile(r"^\s*Direct\s+configuration\s*=\s*(\d+)", re.IGNORECASE)


def _numbers(line: str, count: int, label: str) -> tuple[float, ...]:
    tokens = line.split()
    if len(tokens) < count:
        raise ValueError(f"{label} requires {count} numbers")
    try:
        values = tuple(float(value) for value in tokens[:count])
    except ValueError as error:
        raise ValueError(f"invalid number in {label}") from error
    if not np.isfinite(values).all():
        raise ValueError(f"non-finite number in {label}")
    return values


def _parse_header(lines: list[str]) -> tuple[VaspHeader, int, bool]:
    if len(lines) < 7:
        raise ValueError("XDATCAR header is incomplete")
    scale_tokens = lines[1].split()
    scale = VaspScale(tuple(float(value) for value in scale_tokens))
    lattice = np.asarray([_numbers(lines[index], 3, "lattice") for index in range(2, 5)])
    population_index = 5
    unresolved = False
    try:
        counts = tuple(int(value) for value in lines[population_index].split())
        species = None
        unresolved = True
    except ValueError:
        species = tuple(lines[population_index].split())
        population_index += 1
        counts = tuple(int(value) for value in lines[population_index].split())
    header = VaspHeader(
        title=lines[0].rstrip("\r\n"),
        scale=scale,
        raw_lattice=lattice,
        species_labels=species,
        counts=counts,
        coordinate_mode="direct",
    )
    return header, population_index + 1, unresolved


def parse_xdatcar(source: str, source_name: str | None = None) -> ReadResult:
    """Index complete configurations without materializing their structures."""

    lines = source.splitlines(keepends=True)
    offsets = [0]
    for line in lines:
        offsets.append(offsets[-1] + len(line))
    diagnostics: list[Diagnostic] = []
    frames: list[VaspFrameSpan] = []
    header: VaspHeader | None = None
    try:
        header, start_index, unresolved = _parse_header(lines)
        if unresolved:
            diagnostics.append(
                Diagnostic(
                    Severity.WARNING,
                    "vasp.map.species_unresolved",
                    "VASP 4 XDATCAR does not report chemical species",
                )
            )
        expected = sum(header.counts)
        index = start_index
        while index < len(lines):
            match = _CONFIGURATION.match(lines[index])
            if match is None:
                index += 1
                continue
            coordinate_start = index + 1
            valid = True
            for row_index in range(coordinate_start, coordinate_start + expected):
                if row_index >= len(lines) or _CONFIGURATION.match(lines[row_index]):
                    valid = False
                    break
                try:
                    _numbers(lines[row_index], 3, "XDATCAR coordinate")
                except ValueError:
                    valid = False
                    break
            if not valid:
                diagnostics.append(
                    Diagnostic(
                        Severity.WARNING,
                        "vasp.xdatcar.frame_incomplete",
                        f"Configuration {match.group(1)} is incomplete and was ignored",
                    )
                )
                index = coordinate_start
                continue
            frames.append(
                VaspFrameSpan(
                    index=len(frames),
                    start_offset=offsets[coordinate_start],
                    end_offset=offsets[coordinate_start + expected],
                    reported_index=int(match.group(1)),
                )
            )
            index = coordinate_start + expected
    except (IndexError, ValueError) as error:
        diagnostics.append(Diagnostic(Severity.ERROR, "vasp.xdatcar.parse_error", str(error)))

    document = XdatcarDocument(source, source_name, header, tuple(frames))
    references = []
    for position, frame in enumerate(frames):
        source_reference = SourceReference(
            source_name,
            "vasp-xdatcar",
            f"configuration:{frame.reported_index}",
            frame.start_offset,
            frame.end_offset,
        )
        references.append(
            FrameReference(
                frame.index,
                role="final" if position == len(frames) - 1 else "intermediate",
                source=source_reference,
                metadata={"configuration": frame.reported_index},
            )
        )

    def load(reference: FrameReference):
        from . import mapper

        return mapper.map_vasp_snapshot(load_xdatcar_snapshot(document, reference))

    newline = "\r\n" if "\r\n" in source else "\n"
    return ReadResult(
        document=document,
        structures=StructureSequence(tuple(references), load),
        diagnostics=tuple(diagnostics),
        source_info=SourceInfo(source_name, "vasp-xdatcar", "utf-8", newline),
    )


def load_xdatcar_snapshot(
    document: XdatcarDocument,
    reference: FrameReference,
) -> VaspSnapshot:
    """Load one indexed XDATCAR configuration into numerical arrays."""

    header = document.header
    if header is None:
        raise ValueError("XDATCAR header is unavailable")
    try:
        frame = document.frames[reference.index]
    except IndexError as error:
        raise IndexError("XDATCAR frame reference is outside the document") from error
    rows = document.raw_source[frame.start_offset : frame.end_offset].splitlines()
    coordinates = np.asarray([_numbers(line, 3, "XDATCAR coordinate") for line in rows])
    expected = sum(header.counts)
    if coordinates.shape != (expected, 3):
        raise ValueError("XDATCAR frame coordinate count changed after indexing")
    species = []
    if header.species_labels is None:
        for type_index, count in enumerate(header.counts, start=1):
            species.extend([UnknownSpecies(f"vasp:type:{type_index}")] * count)
    else:
        for label, count in zip(header.species_labels, header.counts, strict=True):
            species.extend([label] * count)
    reported = frame.reported_index if frame.reported_index is not None else frame.index + 1
    return VaspSnapshot(
        name=f"{header.title or 'XDATCAR'} — configuration {reported}",
        lattice=scaled_lattice(header.scale, header.raw_lattice),
        species=tuple(species),
        fractional=coordinates,
        frame_index=frame.index,
        source=reference.source
        or SourceReference(document.source_name, "vasp-xdatcar", f"configuration:{reported}"),
    )


__all__ = ["load_xdatcar_snapshot", "parse_xdatcar"]
