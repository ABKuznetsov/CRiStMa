"""Map loss-preserving CIF documents to canonical CRiStMa crystals."""

from __future__ import annotations

from itertools import chain

from cristma.chemistry.elements import normalize_element
from cristma.core.cell import UnitCell
from cristma.core.structure import Crystal, IndependentSite, SiteComponent
from cristma.core.values import MeasuredValue, parse_measured_value
from cristma.io.diagnostics import Diagnostic, Severity
from cristma.symmetry.affine import parse_xyz_operation
from cristma.symmetry.orbit import SpaceGroupDefinition, expand_orbit

from . import names
from .document import CifBlock, CifDocument, CifLoop, CifScalar
from .tokens import CifToken


def _scalar(block: CifBlock, aliases: tuple[str, ...]) -> CifScalar | None:
    for alias in aliases:
        value = block.scalar(alias)
        if value is not None:
            return value
    return None


def _scalar_text(block: CifBlock, aliases: tuple[str, ...]) -> str | None:
    scalar = _scalar(block, aliases)
    if scalar is None or scalar.value in {"?", "."}:
        return None
    return scalar.value


def _find_loop(block: CifBlock, required: tuple[str, ...]) -> CifLoop | None:
    required_names = {item.casefold() for item in required}
    for loop in block.loops:
        present = {item.casefold() for item in loop.tags}
        if required_names <= present:
            return loop
    return None


def _column(loop: CifLoop, aliases: tuple[str, ...]) -> int | None:
    for alias in aliases:
        index = loop.column_index(alias)
        if index is not None:
            return index
    return None


def _token(
    row: tuple[CifToken, ...],
    loop: CifLoop,
    aliases: tuple[str, ...],
) -> CifToken | None:
    index = _column(loop, aliases)
    return row[index] if index is not None else None


def _looks_structural(block: CifBlock) -> bool:
    if any(_scalar(block, alias) is not None for alias in (names.CELL_A, names.ATOM_LABEL)):
        return True
    return any(
        loop.column_index(names.ATOM_LABEL[0]) is not None
        for loop in block.loops
    )


def _cell(
    block: CifBlock,
    diagnostics: list[Diagnostic],
) -> UnitCell | None:
    fields = (
        (names.CELL_A, "angstrom"),
        (names.CELL_B, "angstrom"),
        (names.CELL_C, "angstrom"),
        (names.CELL_ALPHA, "degree"),
        (names.CELL_BETA, "degree"),
        (names.CELL_GAMMA, "degree"),
    )
    scalars = [_scalar(block, aliases) for aliases, _unit in fields]
    if any(item is None for item in scalars):
        diagnostics.append(
            Diagnostic(
                Severity.ERROR,
                "cif.map.cell_missing",
                f"Block {block.name!r} does not report all six unit-cell parameters",
                block.data_token.span,
            )
        )
        return None
    try:
        values = tuple(
            parse_measured_value(scalar.raw_value, unit=unit)
            for scalar, (_aliases, unit) in zip(scalars, fields, strict=True)
        )
        if any(value.value is None for value in values):
            raise ValueError("missing cell value")
        return UnitCell(*values)
    except ValueError as exc:
        diagnostics.append(
            Diagnostic(
                Severity.ERROR,
                "cif.map.cell_invalid",
                f"Invalid unit cell in block {block.name!r}: {exc}",
                block.data_token.span,
            )
        )
        return None


def _symmetry(
    block: CifBlock,
    diagnostics: list[Diagnostic],
) -> SpaceGroupDefinition | None:
    operation_tokens: list[CifToken] = []
    for loop in block.loops:
        index = _column(loop, names.SYMMETRY_OPERATION)
        if index is not None:
            operation_tokens.extend(row[index] for row in loop.row_tokens)
    scalar_operation = _scalar(block, names.SYMMETRY_OPERATION)
    if scalar_operation is not None:
        operation_tokens.append(scalar_operation.value_token)

    provenance = "reported"
    if not operation_tokens:
        provenance = "identity_fallback"
        diagnostics.append(
            Diagnostic(
                Severity.WARNING,
                "cif.map.symmetry_operations_missing",
                "No symmetry operations reported; using identity only",
                block.data_token.span,
                recovery="x,y,z",
            )
        )
        operations = (parse_xyz_operation("x,y,z", operation_id="op:1"),)
    else:
        parsed = []
        for index, token in enumerate(operation_tokens, start=1):
            try:
                parsed.append(
                    parse_xyz_operation(token.value, operation_id=f"op:{index}")
                )
            except ValueError as exc:
                diagnostics.append(
                    Diagnostic(
                        Severity.ERROR,
                        "cif.map.symmetry_operation_invalid",
                        str(exc),
                        token.span,
                    )
                )
        if len(parsed) != len(operation_tokens):
            return None
        operations = tuple(parsed)

    number = None
    number_text = _scalar_text(block, names.IT_NUMBER)
    if number_text is not None:
        try:
            number = int(float(number_text))
        except ValueError:
            diagnostics.append(
                Diagnostic(
                    Severity.WARNING,
                    "cif.map.space_group_number_invalid",
                    f"Invalid reported space-group number: {number_text!r}",
                    _scalar(block, names.IT_NUMBER).value_token.span,
                )
            )

    return SpaceGroupDefinition(
        operations=operations,
        provenance=provenance,
        number=number,
        hm_symbol=_scalar_text(block, names.HM_SYMBOL),
        hall_symbol=_scalar_text(block, names.HALL_SYMBOL),
        setting=_scalar_text(block, names.SETTING),
        origin_choice=_scalar_text(block, names.ORIGIN_CHOICE),
    )


def _sites(
    block: CifBlock,
    diagnostics: list[Diagnostic],
) -> tuple[IndependentSite, ...] | None:
    atom_loop = _find_loop(
        block,
        (
            names.ATOM_LABEL[0],
            names.FRACT_X[0],
            names.FRACT_Y[0],
            names.FRACT_Z[0],
        ),
    )
    if atom_loop is None:
        diagnostics.append(
            Diagnostic(
                Severity.ERROR,
                "cif.map.atom_loop_missing",
                f"Block {block.name!r} has no complete fractional atom-site loop",
                block.data_token.span,
            )
        )
        return None

    sites: list[IndependentSite] = []
    block_failed = False
    for row_index, row in enumerate(atom_loop.row_tokens):
        label_token = _token(row, atom_loop, names.ATOM_LABEL)
        coordinate_tokens = tuple(
            _token(row, atom_loop, aliases)
            for aliases in (names.FRACT_X, names.FRACT_Y, names.FRACT_Z)
        )
        try:
            coordinates = tuple(
                parse_measured_value(token.raw)
                for token in coordinate_tokens
                if token is not None
            )
        except ValueError as exc:
            diagnostics.append(
                Diagnostic(
                    Severity.ERROR,
                    "cif.map.coordinate_invalid",
                    str(exc),
                    label_token.span if label_token is not None else block.data_token.span,
                )
            )
            block_failed = True
            continue
        if len(coordinates) != 3 or any(value.value is None for value in coordinates):
            diagnostics.append(
                Diagnostic(
                    Severity.ERROR,
                    "cif.map.coordinate_missing",
                    f"Atom row {row_index + 1} has incomplete fractional coordinates",
                    label_token.span if label_token is not None else block.data_token.span,
                )
            )
            block_failed = True
            continue

        label = label_token.value if label_token is not None else f"site{row_index + 1}"
        type_token = _token(row, atom_loop, names.ATOM_TYPE)
        try:
            element = normalize_element(type_token.value if type_token is not None else label)
        except ValueError as exc:
            diagnostics.append(
                Diagnostic(
                    Severity.ERROR,
                    "cif.map.element_invalid",
                    str(exc),
                    (type_token or label_token).span,
                )
            )
            block_failed = True
            continue

        occupancy_token = _token(row, atom_loop, names.OCCUPANCY)
        if occupancy_token is None:
            occupancy = MeasuredValue(1.0, None, None)
            diagnostics.append(
                Diagnostic(
                    Severity.WARNING,
                    "cif.map.occupancy_defaulted",
                    f"Occupancy is absent for {label}; using CIF default 1",
                    label_token.span,
                    recovery="1.0",
                )
            )
        else:
            try:
                occupancy = parse_measured_value(occupancy_token.raw)
            except ValueError as exc:
                diagnostics.append(
                    Diagnostic(
                        Severity.ERROR,
                        "cif.map.occupancy_invalid",
                        str(exc),
                        occupancy_token.span,
                    )
                )
                block_failed = True
                continue
            if occupancy.value is None:
                diagnostics.append(
                    Diagnostic(
                        Severity.ERROR,
                        "cif.map.occupancy_missing",
                        f"Occupancy is unknown for {label}",
                        occupancy_token.span,
                    )
                )
                block_failed = True
                continue

        try:
            sites.append(
                IndependentSite(
                    id=f"{block.name}:{label}:{row_index}",
                    label=label,
                    components=(
                        SiteComponent(
                            element,
                            occupancy,
                            metadata={
                                "reported_type_symbol": type_token.value
                                if type_token is not None
                                else None,
                            },
                        ),
                    ),
                    fractional=coordinates,
                    metadata={"source_row": row_index},
                )
            )
        except ValueError as exc:
            diagnostics.append(
                Diagnostic(
                    Severity.ERROR,
                    "cif.map.site_invalid",
                    f"Invalid atom site {label}: {exc}",
                    label_token.span,
                )
            )
            block_failed = True

    return None if block_failed else tuple(sites)


def _metadata(block: CifBlock) -> dict[str, object]:
    values: dict[str, object] = {"cif_block": block.name}
    for key, aliases in names.METADATA.items():
        value = _scalar_text(block, aliases)
        if value is not None:
            values[key] = value
    return values


def map_cif_structures(
    document: CifDocument,
) -> tuple[tuple[Crystal, ...], tuple[Diagnostic, ...]]:
    """Map every structural CIF block to a canonical asymmetric-unit crystal."""

    structures: list[Crystal] = []
    diagnostics: list[Diagnostic] = []
    for block in document.blocks:
        if not _looks_structural(block):
            continue
        cell = _cell(block, diagnostics)
        if cell is None:
            continue
        symmetry = _symmetry(block, diagnostics)
        if symmetry is None:
            continue
        sites = _sites(block, diagnostics)
        if sites is None:
            continue
        expanded = tuple(
            chain.from_iterable(
                expand_orbit(site, symmetry.operations)
                for site in sites
            )
        )
        structures.append(
            Crystal(
                name=block.name,
                cell=cell,
                sites=sites,
                space_group=symmetry,
                formula=_scalar_text(block, names.FORMULA),
                metadata=_metadata(block),
                expanded_sites=expanded,
            )
        )
    return tuple(structures), tuple(diagnostics)
