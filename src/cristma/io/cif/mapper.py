"""Map loss-preserving CIF documents to canonical CRiStMa crystals."""

from __future__ import annotations

from dataclasses import replace
from itertools import chain
import math

import numpy as np

from cristma.chemistry.elements import normalize_element
from cristma.core.cell import UnitCell
from cristma.structure import (
    CrystalStructure,
    DisplacementParameters,
    IndependentSite,
    SiteComponent,
)
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


def _optional_measured(
    token: CifToken | None,
    *,
    unit: str | None = None,
) -> MeasuredValue | None:
    if token is None or token.value in {"?", "."}:
        return None
    return parse_measured_value(token.raw, unit=unit)


def _optional_text(token: CifToken | None) -> str | None:
    if token is None or token.value in {"?", "."}:
        return None
    return token.value


def _optional_integer(token: CifToken | None) -> int | None:
    if token is None or token.value in {"?", "."}:
        return None
    value = float(token.value)
    if not value.is_integer():
        raise ValueError(f"Expected integer, got {token.value!r}")
    return int(value)


def _label_identity(site: IndependentSite) -> str:
    identity = site.label
    for component in site.components:
        if identity.casefold().startswith(component.element.casefold()):
            identity = identity[len(component.element) :]
            break
    return identity.casefold()


def _coincident(
    left: IndependentSite,
    right: IndependentSite,
    tolerance: float = 1e-8,
) -> bool:
    return all(
        abs((float(a.value) - float(b.value) + 0.5) % 1.0 - 0.5) <= tolerance
        for a, b in zip(left.fractional, right.fractional, strict=True)
    )


def _can_merge_mixed(group: list[IndependentSite], candidate: IndependentSite) -> bool:
    sites = [*group, candidate]
    components = [component for site in sites for component in site.components]
    if len({component.element for component in components}) != len(components):
        return False
    if any(
        component.occupancy.raw is None
        or component.occupancy.value is None
        or not 0 <= component.occupancy.value < 1
        for component in components
    ):
        return False
    if math.fsum(float(component.occupancy.value) for component in components) > 1.0 + 1e-6:
        return False
    if len({site.displacement for site in sites}) > 1:
        return False

    assemblies = {site.disorder_assembly for site in sites}
    groups = {site.disorder_group for site in sites}
    explicit_disorder = None not in assemblies and len(assemblies) == 1 and len(groups) == 1
    matching_identity = len({_label_identity(site) for site in sites}) == 1
    return explicit_disorder or matching_identity


def _same_explicit_disorder_model(
    group: list[IndependentSite],
    candidate: IndependentSite,
) -> bool:
    sites = [*group, candidate]
    assemblies = {site.disorder_assembly for site in sites}
    groups = {site.disorder_group for site in sites}
    return None not in assemblies and len(assemblies) == 1 and len(groups) == 1


def _combined_occupancy(group: list[IndependentSite], candidate: IndependentSite) -> float:
    return math.fsum(
        float(component.occupancy.value)
        for site in (*group, candidate)
        for component in site.components
    )


def _merge_coincident_sites(
    sites: tuple[IndependentSite, ...],
    diagnostics: list[Diagnostic],
) -> tuple[IndependentSite, ...] | None:
    consumed: set[int] = set()
    merged: list[IndependentSite] = []
    for index, site in enumerate(sites):
        if index in consumed:
            continue
        group = [site]
        for candidate_index in range(index + 1, len(sites)):
            if candidate_index in consumed:
                continue
            candidate = sites[candidate_index]
            if not _coincident(site, candidate):
                continue
            if (
                _same_explicit_disorder_model(group, candidate)
                and _combined_occupancy(group, candidate) > 1.0 + 1e-12
            ):
                diagnostics.append(
                    Diagnostic(
                        Severity.ERROR,
                        "cif.map.occupancy_total_exceeds_one",
                        "Explicit disorder components have total occupancy above one",
                    )
                )
                return None
            if _can_merge_mixed(group, candidate):
                group.append(candidate)
                consumed.add(candidate_index)
            else:
                diagnostics.append(
                    Diagnostic(
                        Severity.WARNING,
                        "cif.map.coincident_sites_unmerged",
                        f"Coincident sites {site.label} and {candidate.label} were not merged",
                    )
                )

        if len(group) == 1:
            merged.append(site)
            continue
        wyckoff_values = {item.wyckoff for item in group}
        multiplicities = {item.reported_multiplicity for item in group}
        source_rows = tuple(item.metadata["source_row"] for item in group)
        merged.append(
            replace(
                site,
                id=f"{site.id.rsplit(':', 1)[0]}:mixed:{','.join(map(str, source_rows))}",
                label="/".join(item.label for item in group),
                components=tuple(
                    component
                    for item in group
                    for component in item.components
                ),
                wyckoff=next(iter(wyckoff_values)) if len(wyckoff_values) == 1 else None,
                reported_multiplicity=(
                    next(iter(multiplicities)) if len(multiplicities) == 1 else None
                ),
                metadata={"source_rows": source_rows},
            )
        )
    return tuple(merged)


def _attach_anisotropic_displacements(
    block: CifBlock,
    sites: tuple[IndependentSite, ...],
    diagnostics: list[Diagnostic],
) -> tuple[IndependentSite, ...]:
    required = (
        names.ANISO_LABEL[0],
        names.ANISO_U11[0],
        names.ANISO_U22[0],
        names.ANISO_U33[0],
        names.ANISO_U12[0],
        names.ANISO_U13[0],
        names.ANISO_U23[0],
    )
    loop = _find_loop(block, required)
    if loop is None:
        return sites

    by_label: dict[str, tuple[tuple[MeasuredValue, ...], CifToken]] = {}
    for row in loop.row_tokens:
        label_token = _token(row, loop, names.ANISO_LABEL)
        if label_token is None:
            continue
        value_tokens = tuple(
            _token(row, loop, aliases)
            for aliases in (
                names.ANISO_U11,
                names.ANISO_U22,
                names.ANISO_U33,
                names.ANISO_U12,
                names.ANISO_U13,
                names.ANISO_U23,
            )
        )
        try:
            values = tuple(
                parse_measured_value(token.raw, unit="angstrom^2")
                for token in value_tokens
                if token is not None
            )
        except ValueError as exc:
            diagnostics.append(
                Diagnostic(
                    Severity.ERROR,
                    "cif.map.adp_invalid",
                    str(exc),
                    label_token.span,
                )
            )
            continue
        if len(values) != 6 or any(value.value is None for value in values):
            diagnostics.append(
                Diagnostic(
                    Severity.ERROR,
                    "cif.map.adp_incomplete",
                    f"Anisotropic tensor is incomplete for {label_token.value}",
                    label_token.span,
                )
            )
            continue
        by_label[label_token.value.casefold()] = (values, label_token)

    updated = []
    for site in sites:
        item = by_label.get(site.label.casefold())
        if item is None:
            updated.append(site)
            continue
        values, label_token = item
        u11, u22, u33, u12, u13, u23 = values
        tensor = (
            (u11, u12, u13),
            (u12, u22, u23),
            (u13, u23, u33),
        )
        numeric = np.array(
            [[float(value.value) for value in row] for row in tensor],
            dtype=float,
        )
        if float(np.linalg.eigvalsh(numeric).min()) < -1e-12:
            diagnostics.append(
                Diagnostic(
                    Severity.WARNING,
                    "cif.map.adp_not_positive_semidefinite",
                    f"Reported anisotropic tensor is not positive semidefinite for {site.label}",
                    label_token.span,
                )
            )
        updated.append(
            replace(
                site,
                displacement=DisplacementParameters(
                    kind="U_aniso",
                    tensor=tensor,
                    reported_kind="U",
                ),
            )
        )
    return tuple(updated)


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
            if occupancy.value > 1:
                diagnostics.append(
                    Diagnostic(
                        Severity.ERROR,
                        "cif.map.occupancy_total_exceeds_one",
                        f"Occupancy exceeds one for {label}: {occupancy.value}",
                        occupancy_token.span,
                    )
                )
                block_failed = True
                continue

        oxidation_token = _token(row, atom_loop, names.OXIDATION)
        multiplicity_token = _token(row, atom_loop, names.MULTIPLICITY)
        u_iso_token = _token(row, atom_loop, names.U_ISO)
        b_iso_token = _token(row, atom_loop, names.B_ISO)
        try:
            oxidation = _optional_measured(oxidation_token)
            reported_multiplicity = _optional_integer(multiplicity_token)
            u_iso = _optional_measured(u_iso_token, unit="angstrom^2")
            b_iso = _optional_measured(b_iso_token, unit="angstrom^2")
        except ValueError as exc:
            diagnostics.append(
                Diagnostic(
                    Severity.ERROR,
                    "cif.map.site_value_invalid",
                    f"Invalid reported value for {label}: {exc}",
                    label_token.span,
                )
            )
            block_failed = True
            continue

        displacement = None
        if u_iso is not None:
            displacement = DisplacementParameters(
                kind="U_iso",
                isotropic=u_iso,
                reported_kind="U",
            )
        elif b_iso is not None:
            displacement = DisplacementParameters(
                kind="B_iso",
                isotropic=b_iso,
                reported_kind="B",
            )

        try:
            sites.append(
                IndependentSite(
                    id=f"{block.name}:{label}:{row_index}",
                    label=label,
                    components=(
                        SiteComponent(
                            element,
                            occupancy,
                            oxidation_state=oxidation,
                            metadata={
                                "reported_type_symbol": type_token.value
                                if type_token is not None
                                else None,
                            },
                        ),
                    ),
                    fractional=coordinates,
                    wyckoff=_optional_text(_token(row, atom_loop, names.WYCKOFF)),
                    reported_multiplicity=reported_multiplicity,
                    disorder_assembly=_optional_text(
                        _token(row, atom_loop, names.DISORDER_ASSEMBLY)
                    ),
                    disorder_group=_optional_text(
                        _token(row, atom_loop, names.DISORDER_GROUP)
                    ),
                    displacement=displacement,
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

    if block_failed:
        return None
    with_anisotropic = _attach_anisotropic_displacements(
        block,
        tuple(sites),
        diagnostics,
    )
    return _merge_coincident_sites(with_anisotropic, diagnostics)


def _metadata(block: CifBlock) -> dict[str, object]:
    values: dict[str, object] = {"cif_block": block.name}
    for key, aliases in names.METADATA.items():
        value = _scalar_text(block, aliases)
        if value is not None:
            values[key] = value
    return values


def map_cif_structures(
    document: CifDocument,
) -> tuple[tuple[CrystalStructure, ...], tuple[Diagnostic, ...]]:
    """Map every structural CIF block to a canonical asymmetric-unit crystal."""

    structures: list[CrystalStructure] = []
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
        structure_id = f"cif:{block.name}"
        checked_sites: list[IndependentSite] = []
        expanded_by_site = []
        for site in sites:
            orbit = expand_orbit(
                site,
                symmetry.operations,
                structure_id=structure_id,
            )
            calculated_multiplicity = len(orbit)
            checked_site = replace(
                site,
                calculated_multiplicity=calculated_multiplicity,
            )
            if (
                site.reported_multiplicity is not None
                and site.reported_multiplicity != calculated_multiplicity
            ):
                diagnostics.append(
                    Diagnostic(
                        Severity.WARNING,
                        "cif.map.multiplicity_mismatch",
                        f"{site.label}: reported multiplicity {site.reported_multiplicity}, "
                        f"calculated {calculated_multiplicity}",
                    )
                )
            checked_sites.append(checked_site)
            expanded_by_site.append(orbit)
        sites = tuple(checked_sites)
        expanded = tuple(chain.from_iterable(expanded_by_site))
        structures.append(
            CrystalStructure(
                name=block.name,
                cell=cell,
                sites=sites,
                id=structure_id,
                space_group=symmetry,
                formula=_scalar_text(block, names.FORMULA),
                metadata=_metadata(block),
                expanded_sites=expanded,
            )
        )
    return tuple(structures), tuple(diagnostics)
