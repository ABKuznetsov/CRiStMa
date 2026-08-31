"""Selected-frame parser and source validation for XYZ/extXYZ."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from cristma.chemistry import element_from_atomic_number, normalize_element
from cristma.diagnostics import Diagnostic, Severity
from cristma.structure import FrameReference, SourceReference

from .document import XyzDocument, XyzFrame, XyzFrameSpan, XyzPropertySpec
from .metadata import parse_xyz_metadata


_PLAIN_SCHEMA = (
    XyzPropertySpec("species", "S", 1),
    XyzPropertySpec("pos", "R", 3),
)


def _frame_lines(document: XyzDocument, span: XyzFrameSpan) -> tuple[str, list[str]]:
    comment = document.raw_source[span.comment_start_offset : span.comment_end_offset].rstrip("\r\n")
    rows = document.raw_source[span.atom_rows_start_offset : span.end_offset].splitlines()
    if len(rows) != span.atom_count:
        raise ValueError("XYZ frame atom-row count does not match its index")
    return comment, rows


def _convert(token: str, kind: str) -> object:
    try:
        if kind == "S":
            return token
        if kind == "I":
            return int(token)
        if kind == "R":
            value = float(token.replace("D", "E").replace("d", "e"))
            if not np.isfinite(value):
                raise ValueError
            return value
        if kind == "L":
            normalized = token.upper()
            if normalized in {"T", "TRUE"}:
                return True
            if normalized in {"F", "FALSE"}:
                return False
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid XYZ property value {token!r}") from exc
    raise ValueError(f"invalid XYZ property value {token!r}")


def _parse_rows(
    rows: Iterable[str],
    schema: tuple[XyzPropertySpec, ...],
    *,
    plain: bool,
) -> tuple[dict[str, list[object]], bool]:
    required_width = sum(item.width for item in schema)
    values = {item.name: [] for item in schema}
    ignored_plain = False
    for row_number, row in enumerate(rows, start=1):
        tokens = row.split()
        if (plain and len(tokens) < required_width) or (not plain and len(tokens) != required_width):
            raise ValueError(f"XYZ row {row_number} column count does not match declared schema")
        if plain and len(tokens) > required_width:
            ignored_plain = True
            tokens = tokens[:required_width]
        offset = 0
        for item in schema:
            converted = [_convert(token, item.kind) for token in tokens[offset : offset + item.width]]
            values[item.name].append(converted[0] if item.width == 1 else converted)
            offset += item.width
    return values, ignored_plain


def _arrays(values: dict[str, list[object]], schema: tuple[XyzPropertySpec, ...]) -> dict[str, np.ndarray]:
    dtypes = {"S": str, "I": np.int64, "R": np.float64, "L": np.bool_}
    return {item.name: np.asarray(values[item.name], dtype=dtypes[item.kind]) for item in schema}


def _source_for(document: XyzDocument, span: XyzFrameSpan) -> SourceReference:
    return SourceReference(document.source_name, "xyz", f"frame:{span.index}", span.start_offset, span.end_offset)


def _diagnostic(severity: Severity, code: str, message: str) -> Diagnostic:
    return Diagnostic(severity=severity, code=code, message=message)


def _semantic_diagnostics(
    metadata: object,
    schema: tuple[XyzPropertySpec, ...],
    values: dict[str, list[object]],
    *,
    plain: bool,
) -> list[Diagnostic]:
    diagnostics: list[Diagnostic] = []
    specs = {item.name: item for item in schema}
    species_spec = specs.get("species")
    number_spec = specs.get("Z")
    position_spec = specs.get("pos")
    if not plain and (
        (species_spec is None and number_spec is None)
        or species_spec is not None and (species_spec.kind, species_spec.width) != ("S", 1)
        or number_spec is not None and (number_spec.kind, number_spec.width) != ("I", 1)
        or position_spec is None
        or (position_spec.kind, position_spec.width) != ("R", 3)
    ):
        diagnostics.append(
            _diagnostic(
                Severity.ERROR,
                "xyz.map.schema_invalid",
                "extXYZ requires species:S:1 or Z:I:1 and exactly pos:R:3",
            )
        )

    normalized_species: list[str | None] = []
    for token in values.get("species", []):
        try:
            normalized_species.append(normalize_element(str(token)))
        except ValueError:
            normalized_species.append(None)
            diagnostics.append(
                _diagnostic(
                    Severity.WARNING,
                    "xyz.map.species_unresolved",
                    f"Unknown XYZ species {token!r} was retained explicitly",
                )
            )
    normalized_numbers: list[str | None] = []
    for token in values.get("Z", []):
        try:
            normalized_numbers.append(element_from_atomic_number(int(token)))
        except (TypeError, ValueError):
            normalized_numbers.append(None)
            diagnostics.append(
                _diagnostic(
                    Severity.WARNING,
                    "xyz.map.species_unresolved",
                    f"Unsupported XYZ atomic number {token!r} was retained explicitly",
                )
            )
    if normalized_species and normalized_numbers and any(
        symbol is not None and number is not None and symbol != number
        for symbol, number in zip(normalized_species, normalized_numbers, strict=True)
    ):
        diagnostics.append(
            _diagnostic(
                Severity.ERROR,
                "xyz.map.species_conflict",
                "XYZ species and atomic number columns disagree",
            )
        )

    lattice = getattr(metadata, "lattice")
    pbc = getattr(metadata, "pbc")
    if pbc is not None and any(pbc) and lattice is None:
        diagnostics.append(
            _diagnostic(
                Severity.ERROR,
                "xyz.map.lattice_required",
                "Periodic XYZ axes require an explicit Lattice",
            )
        )
    return diagnostics


def load_xyz_frame(document: XyzDocument, reference: FrameReference) -> XyzFrame:
    """Materialize one indexed XYZ frame and no other trajectory frame."""

    try:
        span = document.frames[reference.index]
    except IndexError as exc:
        raise IndexError("XYZ frame reference is outside the document") from exc
    comment, rows = _frame_lines(document, span)
    metadata = parse_xyz_metadata(comment)
    plain = not metadata.schema
    schema = _PLAIN_SCHEMA if plain else metadata.schema
    values, _ = _parse_rows(rows, schema, plain=plain)
    semantic = _semantic_diagnostics(metadata, schema, values, plain=plain)
    errors = [item for item in semantic if item.severity is Severity.ERROR]
    if errors:
        raise ValueError(errors[0].message)
    name_value = metadata.values.get("name") or metadata.values.get("label")
    name = str(name_value or (comment if plain and comment else f"frame {span.index + 1}"))
    return XyzFrame(
        name=name,
        atom_count=span.atom_count,
        comment=comment,
        metadata=metadata.values,
        schema=schema,
        columns=_arrays(values, schema),
        lattice=metadata.lattice,
        pbc=metadata.pbc,
        source=reference.source or _source_for(document, span),
    )


def validate_xyz_frame(document: XyzDocument, span: XyzFrameSpan) -> tuple[Diagnostic, ...]:
    """Validate one indexed frame without allocating NumPy property arrays."""

    diagnostics: list[Diagnostic] = []
    try:
        comment, rows = _frame_lines(document, span)
        metadata = parse_xyz_metadata(comment)
        plain = not metadata.schema
        schema = _PLAIN_SCHEMA if plain else metadata.schema
        values, ignored_plain = _parse_rows(rows, schema, plain=plain)
        diagnostics.extend(_semantic_diagnostics(metadata, schema, values, plain=plain))
        if ignored_plain:
            diagnostics.append(
                _diagnostic(
                    Severity.WARNING,
                    "xyz.map.uninterpreted_plain_columns",
                    "Plain XYZ trailing atom columns were preserved only in the source",
                )
            )
        if metadata.lattice is not None and metadata.pbc is None:
            diagnostics.append(
                _diagnostic(
                    Severity.WARNING,
                    "xyz.map.lattice_without_pbc",
                    "Lattice is reported without explicit pbc; the frame remains molecular",
                )
            )
    except ValueError as exc:
        message = str(exc)
        code = "xyz.row.column_count_mismatch" if "column count" in message else "xyz.map.value_invalid"
        diagnostics.append(_diagnostic(Severity.ERROR, code, message))
    return tuple(diagnostics)


__all__ = ["load_xyz_frame", "validate_xyz_frame"]
