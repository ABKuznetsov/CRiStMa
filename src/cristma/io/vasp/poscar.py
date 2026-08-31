"""Loss-preserving native POSCAR and CONTCAR reader."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from cristma.chemistry import UnknownSpecies
from cristma.diagnostics import Diagnostic, Severity
from cristma.io.result import ReadResult, SourceInfo
from cristma.structure import SourceReference

from .document import PoscarDocument, VaspAtomRow, VaspHeader, VaspScale, VaspSnapshot
from .numeric import fractional_from_cartesian, scaled_cartesian, scaled_lattice


def _diagnostic(severity: Severity, code: str, message: str) -> Diagnostic:
    return Diagnostic(severity=severity, code=code, message=message)


def _tokens(lines: Sequence[str], index: int) -> list[str]:
    return lines[index].strip().split() if index < len(lines) else []


def _floats(tokens: Sequence[str], count: int, label: str) -> tuple[float, ...]:
    if len(tokens) < count:
        raise ValueError(f"{label} requires {count} numbers")
    try:
        values = tuple(float(value) for value in tokens[:count])
    except ValueError as error:
        raise ValueError(f"invalid number in {label}") from error
    if not np.isfinite(values).all():
        raise ValueError(f"non-finite number in {label}")
    return values


def _mode(token: str) -> str:
    return "cartesian" if token[:1].lower() in {"c", "k"} else "direct"


def _is_int_row(tokens: Sequence[str]) -> bool:
    if not tokens:
        return False
    try:
        return all(int(token) >= 0 for token in tokens)
    except ValueError:
        return False


def parse_poscar(source: str, source_name: str | None = None) -> ReadResult:
    """Parse POSCAR text while retaining the exact source document."""

    lines = source.splitlines(keepends=True)
    diagnostics: list[Diagnostic] = []
    header: VaspHeader | None = None
    positions: list[VaspAtomRow] = []
    velocities: list[VaspAtomRow] = []
    velocity_mode: str | None = None
    trailing_start: int | None = None

    try:
        if len(lines) < 7:
            raise ValueError("POSCAR header is incomplete")
        title = lines[0].rstrip("\r\n")
        scale = VaspScale(_floats(_tokens(lines, 1), len(_tokens(lines, 1)), "scale"))
        raw_lattice = np.array(
            [_floats(_tokens(lines, index), 3, "lattice vector") for index in range(2, 5)],
            dtype=float,
        )
        index = 5
        first_population = _tokens(lines, index)
        if _is_int_row(first_population):
            species_labels = None
            counts = tuple(int(value) for value in first_population)
            diagnostics.append(
                _diagnostic(
                    Severity.WARNING,
                    "vasp.map.species_unresolved",
                    "VASP 4 population row does not report chemical species",
                )
            )
            index += 1
        else:
            species_labels = tuple(first_population)
            index += 1
            count_tokens = _tokens(lines, index)
            if not _is_int_row(count_tokens):
                raise ValueError("invalid POSCAR population counts")
            counts = tuple(int(value) for value in count_tokens)
            index += 1

        selective = False
        current = _tokens(lines, index)
        if current and current[0][:1].lower() == "s":
            selective = True
            index += 1
            current = _tokens(lines, index)
        if not current:
            raise ValueError("POSCAR coordinate mode is missing")
        coordinate_mode = _mode(current[0])
        index += 1
        header = VaspHeader(
            title=title,
            scale=scale,
            raw_lattice=raw_lattice,
            species_labels=species_labels,
            counts=counts,
            coordinate_mode=coordinate_mode,
            selective_dynamics=selective,
        )

        expected = sum(counts)
        for atom_index in range(expected):
            row_index = index + atom_index
            row_tokens = _tokens(lines, row_index)
            if len(row_tokens) < 3:
                diagnostics.append(
                    _diagnostic(
                        Severity.ERROR,
                        "vasp.poscar.positions_incomplete",
                        f"Expected {expected} position rows, found {len(positions)}",
                    )
                )
                break
            try:
                coordinates = _floats(row_tokens, 3, "atomic position")
            except ValueError:
                diagnostics.append(
                    _diagnostic(
                        Severity.ERROR,
                        "vasp.poscar.positions_incomplete",
                        f"Position row {atom_index + 1} is invalid",
                    )
                )
                break
            flags = None
            if selective:
                flag_tokens = row_tokens[3:6]
                if len(flag_tokens) != 3 or any(value.upper() not in {"T", "F"} for value in flag_tokens):
                    diagnostics.append(
                        _diagnostic(
                            Severity.ERROR,
                            "vasp.poscar.selective_flag_invalid",
                            f"Selective-dynamics flags are invalid on atom row {atom_index + 1}",
                        )
                    )
                else:
                    flags = tuple(value.upper() == "T" for value in flag_tokens)
            positions.append(VaspAtomRow(coordinates, flags))
        index += len(positions)

        # A velocity block is introduced by its own Direct/Cartesian mode line.
        while index < len(lines) and not lines[index].strip():
            index += 1
        next_tokens = _tokens(lines, index)
        if next_tokens and next_tokens[0][:1].lower() in {"d", "c", "k"}:
            velocity_mode = _mode(next_tokens[0])
            index += 1
            for velocity_index in range(expected):
                row_tokens = _tokens(lines, index + velocity_index)
                if len(row_tokens) < 3:
                    break
                try:
                    coordinates = _floats(row_tokens, 3, "velocity")
                except ValueError:
                    break
                velocities.append(VaspAtomRow(coordinates))
            index += len(velocities)
        if index < len(lines):
            trailing_start = sum(len(line) for line in lines[:index])
    except (IndexError, ValueError) as error:
        diagnostics.append(
            _diagnostic(Severity.ERROR, "vasp.poscar.parse_error", str(error))
        )

    document = PoscarDocument(
        raw_source=source,
        source_name=source_name,
        header=header,
        positions=tuple(positions),
        velocity_mode=velocity_mode,
        velocities=tuple(velocities),
        trailing_start=trailing_start,
    )
    newline = "\r\n" if "\r\n" in source else "\n"
    return ReadResult(
        document=document,
        diagnostics=tuple(diagnostics),
        source_info=SourceInfo(source_name, "vasp-poscar", "utf-8", newline),
    )


def poscar_snapshot(document: PoscarDocument) -> VaspSnapshot:
    """Convert a parsed POSCAR document to a numerical VASP snapshot."""

    header = document.header
    if header is None:
        raise ValueError("POSCAR header is unavailable")
    expected = sum(header.counts)
    if len(document.positions) != expected:
        raise ValueError(
            f"position count {len(document.positions)} does not match declared count {expected}"
        )
    lattice = scaled_lattice(header.scale, header.raw_lattice)
    reported = np.asarray([row.coordinates for row in document.positions], dtype=float)
    fractional = (
        reported
        if header.coordinate_mode == "direct"
        else fractional_from_cartesian(
            lattice,
            scaled_cartesian(header.scale, reported, header.raw_lattice),
        )
    )
    species = []
    if header.species_labels is None:
        for type_index, count in enumerate(header.counts, start=1):
            species.extend([UnknownSpecies(f"vasp:type:{type_index}")] * count)
    else:
        for label, count in zip(header.species_labels, header.counts, strict=True):
            species.extend([label] * count)
    selective = None
    if header.selective_dynamics:
        selective = np.asarray(
            [row.selective if row.selective is not None else (False, False, False) for row in document.positions],
            dtype=bool,
        )
    velocity_values = None
    velocity_unit = None
    if document.velocities:
        if len(document.velocities) != expected:
            raise ValueError("velocity count does not match declared atom count")
        velocity_values = np.asarray([row.coordinates for row in document.velocities], dtype=float)
        velocity_unit = (
            "angstrom/fs"
            if document.velocity_mode == "cartesian"
            else "direct_lattice_vector/timestep"
        )
    source_length = len(document.raw_source)
    return VaspSnapshot(
        name=header.title or document.source_name or "VASP structure",
        lattice=lattice,
        species=tuple(species),
        fractional=fractional,
        frame_index=0,
        source=SourceReference(document.source_name, "vasp-poscar", "frame:0", 0, source_length),
        selective_dynamics=selective,
        velocities=velocity_values,
        velocity_mode=document.velocity_mode if velocity_values is not None else None,
        velocity_unit=velocity_unit,
    )


__all__ = ["parse_poscar", "poscar_snapshot"]
