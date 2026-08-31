"""Indexed lazy reader for structural frames reported by VASP OUTCAR."""

from __future__ import annotations

import re

import numpy as np

from cristma.chemistry import UnknownSpecies
from cristma.diagnostics import Diagnostic, Severity
from cristma.io.result import ReadResult, SourceInfo
from cristma.structure import FrameReference, SourceReference, StructureSequence

from .document import OutcarDocument, VaspFrameSpan, VaspSnapshot
from .numeric import fractional_from_cartesian


_VRHFIN = re.compile(r"VRHFIN\s*=\s*([A-Z][a-z]?)\s*:")
_TITEL = re.compile(r"TITEL\s*=\s*\S+\s+([A-Z][a-z]?)(?:_|\s|$)")
_IONS = re.compile(r"ions\s+per\s+type\s*=\s*(.*)", re.IGNORECASE)


def _floats(line: str, count: int) -> tuple[float, ...]:
    tokens = line.split()
    if len(tokens) < count:
        raise ValueError("numeric OUTCAR row is incomplete")
    try:
        values = tuple(float(value) for value in tokens[:count])
    except ValueError as error:
        raise ValueError("invalid numeric OUTCAR row") from error
    if not np.isfinite(values).all():
        raise ValueError("non-finite numeric OUTCAR row")
    return values


def parse_outcar(source: str, source_name: str | None = None) -> ReadResult:
    """Index complete ionic POSITION blocks and their current lattices."""

    lines = source.splitlines(keepends=True)
    offsets = [0]
    for line in lines:
        offsets.append(offsets[-1] + len(line))
    diagnostics: list[Diagnostic] = []
    vrhfin_evidence: list[str] = []
    titel_evidence: list[str] = []
    counts: tuple[int, ...] = ()
    frames: list[VaspFrameSpan] = []
    lattices: list[np.ndarray] = []
    current_lattice: np.ndarray | None = None

    for line in lines:
        match = _VRHFIN.search(line)
        if match and match.group(1) not in vrhfin_evidence:
            vrhfin_evidence.append(match.group(1))
        titel = _TITEL.search(line)
        if titel and titel.group(1) not in titel_evidence:
            titel_evidence.append(titel.group(1))
        ions = _IONS.search(line)
        if ions:
            try:
                counts = tuple(int(value) for value in ions.group(1).split())
            except ValueError:
                diagnostics.append(
                    Diagnostic(Severity.ERROR, "vasp.outcar.ions_per_type_invalid", "Invalid ions per type")
                )

    evidence = vrhfin_evidence or titel_evidence
    inconsistent = bool(evidence and counts and len(evidence) != len(counts))
    if inconsistent:
        diagnostics.append(
            Diagnostic(
                Severity.ERROR,
                "vasp.outcar.species_count_inconsistent",
                "Species evidence does not match ions-per-type columns",
            )
        )
    elif counts and not evidence:
        diagnostics.append(
            Diagnostic(
                Severity.WARNING,
                "vasp.map.species_unresolved",
                "OUTCAR reports atom counts but no usable species labels",
            )
        )

    expected = sum(counts)
    index = 0
    while index < len(lines):
        folded = lines[index].casefold()
        if "direct lattice vectors" in folded:
            try:
                candidate = np.asarray([_floats(lines[index + row], 3) for row in range(1, 4)])
                if abs(float(np.linalg.det(candidate))) <= 1e-15:
                    raise ValueError("singular direct lattice vectors")
                current_lattice = candidate
            except (IndexError, ValueError) as error:
                diagnostics.append(
                    Diagnostic(Severity.ERROR, "vasp.outcar.lattice_invalid", str(error))
                )
            index += 4
            continue
        if "position" in folded and "total-force" in folded:
            row_start = index + 1
            while row_start < len(lines) and (
                not lines[row_start].strip() or set(lines[row_start].strip()) == {"-"}
            ):
                row_start += 1
            rows: list[tuple[float, ...]] = []
            for row_index in range(row_start, row_start + expected):
                try:
                    rows.append(_floats(lines[row_index], 6))
                except (IndexError, ValueError):
                    break
            if expected <= 0:
                diagnostics.append(
                    Diagnostic(
                        Severity.ERROR,
                        "vasp.outcar.atom_count_missing",
                        "No valid ions-per-type field precedes the structural frames",
                    )
                )
            elif len(rows) != expected:
                diagnostics.append(
                    Diagnostic(
                        Severity.WARNING,
                        "vasp.outcar.frame_incomplete",
                        f"POSITION block has {len(rows)} of {expected} declared atoms",
                    )
                )
            elif current_lattice is None:
                diagnostics.append(
                    Diagnostic(
                        Severity.ERROR,
                        "vasp.outcar.frame_without_lattice",
                        "POSITION block has no preceding direct lattice vectors",
                    )
                )
            elif not inconsistent:
                frames.append(
                    VaspFrameSpan(
                        len(frames),
                        offsets[row_start],
                        offsets[row_start + expected],
                        len(frames) + 1,
                    )
                )
                lattices.append(np.array(current_lattice, copy=True))
            index = row_start + max(len(rows), 1)
            continue
        index += 1

    labels = tuple(evidence) if evidence and not inconsistent else None
    document = OutcarDocument(
        source,
        source_name,
        tuple(frames),
        labels,
        counts,
        tuple(lattices),
    )
    references = tuple(
        FrameReference(
            frame.index,
            role="final" if position == len(frames) - 1 else "intermediate",
            source=SourceReference(
                source_name,
                "vasp-outcar",
                f"ionic-step:{frame.index + 1}",
                frame.start_offset,
                frame.end_offset,
            ),
            metadata={"ionic_step": frame.index + 1},
        )
        for position, frame in enumerate(frames)
    )

    def load(reference: FrameReference):
        from . import mapper

        return mapper.map_vasp_snapshot(load_outcar_snapshot(document, reference))

    newline = "\r\n" if "\r\n" in source else "\n"
    return ReadResult(
        document=document,
        structures=StructureSequence(references, load),
        diagnostics=tuple(diagnostics),
        source_info=SourceInfo(source_name, "vasp-outcar", "utf-8", newline),
    )


def load_outcar_snapshot(document: OutcarDocument, reference: FrameReference) -> VaspSnapshot:
    """Load one indexed OUTCAR POSITION/TOTAL-FORCE block."""

    frame = document.frames[reference.index]
    rows = document.raw_source[frame.start_offset : frame.end_offset].splitlines()
    values = np.asarray([_floats(line, 6) for line in rows])
    expected = sum(document.counts)
    if values.shape != (expected, 6):
        raise ValueError("OUTCAR frame row count changed after indexing")
    species = []
    if document.species_labels is None:
        for type_index, count in enumerate(document.counts, start=1):
            species.extend([UnknownSpecies(f"vasp:type:{type_index}")] * count)
    else:
        for label, count in zip(document.species_labels, document.counts, strict=True):
            species.extend([label] * count)
    lattice = document.lattices[reference.index]
    return VaspSnapshot(
        name=f"OUTCAR — ionic step {reference.index + 1}",
        lattice=lattice,
        species=tuple(species),
        fractional=fractional_from_cartesian(lattice, values[:, :3]),
        frame_index=reference.index,
        source=reference.source
        or SourceReference(document.source_name, "vasp-outcar", f"ionic-step:{reference.index + 1}"),
        forces=values[:, 3:6],
        force_unit="eV/angstrom",
    )


__all__ = ["load_outcar_snapshot", "parse_outcar"]
