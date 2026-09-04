"""Fixed-column PDB parsing and validation."""

from __future__ import annotations

import re

from cristma.chemistry import normalize_element
from cristma.diagnostics import Diagnostic, Severity
from cristma.io.result import ReadResult, SourceInfo

from .document import PdbAtomRecord, PdbCryst1Record, PdbDocument


def _diagnostic(severity: Severity, code: str, message: str) -> Diagnostic:
    return Diagnostic(severity, code, message)


def _float(line: str, start: int, end: int, field: str) -> float:
    token = line[start:end].strip()
    try:
        return float(token)
    except ValueError as exc:
        raise ValueError(f"invalid {field} value {token!r}") from exc


def _optional_float(line: str, start: int, end: int, field: str) -> float | None:
    token = line[start:end].strip()
    return None if not token else _float(line, start, end, field)


def _element_from_name(raw_name: str) -> str:
    if raw_name[:1].isspace() or raw_name[:1].isdigit():
        token = raw_name.strip().lstrip("0123456789")[:1]
    else:
        token = raw_name.strip().lstrip("0123456789")
    match = re.match(r"[A-Za-z]{1,2}", token)
    if match is None:
        return token
    letters = match.group(0)
    for width in (2, 1):
        candidate = letters[:width].capitalize()
        try:
            return normalize_element(candidate)
        except ValueError:
            continue
    return letters


def _atom(line: str, line_number: int, model_number: int) -> PdbAtomRecord:
    padded = line.ljust(80)
    serial_token = padded[6:11].strip()
    try:
        serial = int(serial_token)
    except ValueError as exc:
        raise ValueError(f"invalid atom serial {serial_token!r}") from exc
    raw_name = padded[12:16]
    name = raw_name.strip()
    element = padded[76:78].strip() or _element_from_name(raw_name)
    occupancy = _optional_float(padded, 54, 60, "occupancy")
    return PdbAtomRecord(
        padded[:6].strip().upper(),
        serial,
        name or f"ATOM{serial}",
        padded[16:17].strip(),
        padded[17:20].strip(),
        padded[21:22].strip(),
        padded[22:26].strip(),
        (
            _float(padded, 30, 38, "x coordinate"),
            _float(padded, 38, 46, "y coordinate"),
            _float(padded, 46, 54, "z coordinate"),
        ),
        1.0 if occupancy is None else occupancy,
        _optional_float(padded, 60, 66, "B factor"),
        element,
        line_number,
        model_number,
    )


def _cryst1(line: str, line_number: int) -> PdbCryst1Record:
    padded = line.ljust(80)
    z_token = padded[66:70].strip()
    try:
        z = int(z_token) if z_token else None
    except ValueError as exc:
        raise ValueError(f"invalid CRYST1 Z value {z_token!r}") from exc
    return PdbCryst1Record(
        _float(padded, 6, 15, "cell a"),
        _float(padded, 15, 24, "cell b"),
        _float(padded, 24, 33, "cell c"),
        _float(padded, 33, 40, "cell alpha"),
        _float(padded, 40, 47, "cell beta"),
        _float(padded, 47, 54, "cell gamma"),
        padded[55:66].strip(),
        z,
        line_number,
    )


def parse_pdb(source: str, source_name: str | None = None) -> ReadResult:
    """Parse one PDB coordinate model and return canonical structure data."""

    lines = source.splitlines()
    diagnostics: list[Diagnostic] = []
    atoms: list[PdbAtomRecord] = []
    cryst1: PdbCryst1Record | None = None
    current_model = 1
    for line_number, line in enumerate(lines, start=1):
        record = line[:6].strip().upper()
        try:
            if record == "CRYST1":
                cryst1 = _cryst1(line, line_number)
            elif record == "MODEL":
                token = line[10:14].strip()
                try:
                    current_model = int(token)
                except ValueError as exc:
                    raise ValueError(f"invalid MODEL number {token!r}") from exc
            elif record in {"ATOM", "HETATM"}:
                atoms.append(_atom(line, line_number, current_model))
        except ValueError as exc:
            diagnostics.append(
                _diagnostic(
                    Severity.ERROR,
                    "pdb.parse.value_invalid",
                    f"Line {line_number}: {exc}",
                )
            )
    if not atoms:
        diagnostics.append(
            _diagnostic(
                Severity.ERROR,
                "pdb.map.atoms_missing",
                "PDB contains no usable atoms",
            )
        )
    document = PdbDocument(source, source_name, tuple(lines), cryst1, tuple(atoms))
    structures = ()
    if not any(item.severity is Severity.ERROR for item in diagnostics):
        from .mapper import map_pdb_document

        model_numbers = tuple(dict.fromkeys(atom.model_number for atom in atoms))
        structures = tuple(
            map_pdb_document(document, model_number=model_number)
            for model_number in model_numbers
        )
        if cryst1 is not None and any(
            getattr(structure.space_group, "provenance", None) == "identity_fallback"
            for structure in structures
        ):
            diagnostics.append(
                _diagnostic(
                    Severity.WARNING,
                    "pdb.map.space_group_unresolved",
                    f"PDB space group {cryst1.space_group!r} could not be resolved; "
                    "identity symmetry was retained explicitly",
                )
            )
    newline = "\r\n" if "\r\n" in source else "\r" if "\r" in source else "\n"
    return ReadResult(
        document,
        structures,
        tuple(diagnostics),
        SourceInfo(source_name, "pdb", "utf-8", newline),
    )


__all__ = ["parse_pdb"]
