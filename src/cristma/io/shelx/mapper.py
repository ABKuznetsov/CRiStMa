"""Map native SHELX documents to canonical CRiStMa structures."""

from __future__ import annotations

import math
from pathlib import Path
import re

import numpy as np

from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue, parse_measured_value
from cristma.io.diagnostics import Diagnostic, Severity
from cristma.structure import (
    CrystalStructure,
    DisplacementParameters,
    IndependentSite,
    SiteComponent,
    SourceReference,
    StructureCollection,
    StructureProvenance,
)
from cristma.symmetry.displacement import SymmetryConsistencyError
from cristma.symmetry.orbit import SpaceGroupDefinition, expand_orbit

from .document import ShelxAtomRecord, ShelxDocument, ShelxInstructionRecord, ShelxRecord
from .occupancy import ShelxOccupancyExpression
from .records import (
    ShelxCellInstruction,
    ShelxFvarInstruction,
    ShelxLattInstruction,
    ShelxPartInstruction,
    ShelxResiInstruction,
    ShelxSymmInstruction,
    ShelxZerrInstruction,
)
from .sfac import ShelxScatteringEntry, extract_sfac_entries
from .symmetry import build_shelx_operations


def _active_records(document: ShelxDocument) -> tuple[ShelxRecord, ...]:
    return tuple(
        record
        for record in document.records
        if not record.after_hklf and not record.after_end
    )


def _first(records: tuple[ShelxRecord, ...], record_type: type):
    return next((record for record in records if isinstance(record, record_type)), None)


def _reported_title(
    records: tuple[ShelxRecord, ...],
    source_name: str | None,
) -> tuple[str, str | None, str | None]:
    title_record = next(
        (
            record
            for record in records
            if isinstance(record, ShelxInstructionRecord) and record.keyword == "TITL"
        ),
        None,
    )
    raw_title = " ".join(title_record.fields).strip() if title_record else ""
    match = re.fullmatch(r"(.*?)\s+in\s+(.+)", raw_title, flags=re.IGNORECASE)
    if match:
        return match.group(1).strip(), match.group(2).strip(), raw_title
    fallback = Path(source_name).stem if source_name else "SHELX structure"
    return raw_title or fallback, None, raw_title or None


def _cell_with_zerr(
    cell_record: ShelxCellInstruction,
    zerr: ShelxZerrInstruction | None,
) -> UnitCell:
    if zerr is None:
        return cell_record.cell
    reported = (
        cell_record.cell.a,
        cell_record.cell.b,
        cell_record.cell.c,
        cell_record.cell.alpha,
        cell_record.cell.beta,
        cell_record.cell.gamma,
    )
    values = tuple(
        MeasuredValue(value.value, uncertainty.value, value.raw, value.unit, value.missing)
        for value, uncertainty in zip(
            reported,
            zerr.cell_uncertainties,
            strict=True,
        )
    )
    return UnitCell(*values)


def _displacement(
    values: tuple[str, ...],
    record: ShelxAtomRecord,
    diagnostics: list[Diagnostic],
) -> DisplacementParameters | None:
    if not values:
        return None
    try:
        parsed = tuple(
            parse_measured_value(value, unit="angstrom^2") for value in values
        )
    except ValueError as error:
        raise ValueError(f"invalid displacement value: {error}") from error
    if any(value.value is None or not math.isfinite(value.value) for value in parsed):
        raise ValueError("displacement values must be finite")
    if len(parsed) == 1:
        if float(parsed[0].value) < 0:
            diagnostics.append(
                Diagnostic(
                    Severity.WARNING,
                    "shelx.map.displacement_dependency_unmapped",
                    f"Negative SHELX U code for {record.keyword} is retained but not "
                    "misreported as a physical displacement",
                    record.span,
                )
            )
            return None
        return DisplacementParameters(
            kind="U_iso",
            isotropic=parsed[0],
            reported_kind="SHELX_U",
        )
    if len(parsed) != 6:
        raise ValueError("atom displacement must contain one or six values")
    u11, u22, u33, u23, u13, u12 = parsed
    tensor = (
        (u11, u12, u13),
        (u12, u22, u23),
        (u13, u23, u33),
    )
    numeric = np.asarray(
        [[float(value.value) for value in row] for row in tensor],
        dtype=float,
    )
    if float(np.linalg.eigvalsh(numeric).min()) < -1e-12:
        diagnostics.append(
            Diagnostic(
                Severity.WARNING,
                "shelx.map.adp_not_positive_semidefinite",
                f"Reported anisotropic tensor is not positive semidefinite for {record.keyword}",
                record.span,
            )
        )
    return DisplacementParameters(
        kind="U_aniso",
        tensor=tensor,
        reported_kind="SHELX_U",
    )


def _map_atom(
    record: ShelxAtomRecord,
    *,
    structure_id: str,
    site_index: int,
    sfac: tuple[ShelxScatteringEntry, ...],
    free_variables: tuple[float, ...],
    part: int,
    residue: ShelxResiInstruction | None,
    diagnostics: list[Diagnostic],
) -> IndependentSite:
    if record.keyword is None or len(record.fields) < 6:
        raise ValueError("atom record requires SFAC, xyz, occupancy, and displacement")
    try:
        sfac_index = int(record.fields[0])
    except ValueError as error:
        raise ValueError(f"invalid SFAC index {record.fields[0]!r}") from error
    if not 1 <= sfac_index <= len(sfac):
        raise IndexError(f"SFAC index {sfac_index} is outside 1..{len(sfac)}")
    coordinates = tuple(parse_measured_value(value) for value in record.fields[1:4])
    if any(value.value is None for value in coordinates):
        raise ValueError("fractional coordinates must be present")
    expression = ShelxOccupancyExpression.parse(record.fields[4])
    occupancy = expression.evaluate(free_variables)
    displacement = _displacement(record.fields[5:], record, diagnostics)
    scattering = sfac[sfac_index - 1]
    residue_metadata = None
    if residue is not None:
        residue_metadata = {
            "number": residue.residue_number,
            "class": residue.residue_class,
        }
    return IndependentSite(
        id=f"{structure_id}:site:{site_index}:{record.keyword}",
        label=record.keyword,
        components=(
            SiteComponent(
                scattering.species,
                MeasuredValue(occupancy, None, record.fields[4]),
                metadata={
                    "shelx_occupancy": expression,
                    "shelx_sfac_label": scattering.source_label,
                    "shelx_sfac_index": sfac_index,
                },
            ),
        ),
        fractional=coordinates,
        disorder_group=f"part:{part}" if part else None,
        displacement=displacement,
        metadata={
            "shelx_part": part,
            "shelx_residue": residue_metadata,
            "source_line": record.span.start.line,
        },
    )


def map_shelx_structures(
    document: ShelxDocument,
) -> tuple[StructureCollection, tuple[Diagnostic, ...]]:
    """Map one SHELX instruction document to one canonical crystal."""

    diagnostics: list[Diagnostic] = []
    records = _active_records(document)
    cell_record = _first(records, ShelxCellInstruction)
    if cell_record is None:
        diagnostics.append(
            Diagnostic(
                Severity.ERROR,
                "shelx.map.cell_missing",
                "SHELX document has no valid CELL instruction",
            )
        )
        return StructureCollection(), tuple(diagnostics)

    try:
        sfac = extract_sfac_entries(records)
    except ValueError as error:
        diagnostics.append(
            Diagnostic(Severity.ERROR, "shelx.map.sfac_invalid", str(error))
        )
        return StructureCollection(), tuple(diagnostics)
    if not sfac:
        diagnostics.append(
            Diagnostic(Severity.ERROR, "shelx.map.sfac_missing", "SHELX document has no SFAC entries")
        )
        return StructureCollection(), tuple(diagnostics)

    fvar_records = tuple(
        record for record in records if isinstance(record, ShelxFvarInstruction)
    )
    free_variables = tuple(
        float(value.value)
        for record in fvar_records
        for value in record.values
        if value.value is not None
    )
    latt_record = _first(records, ShelxLattInstruction)
    latt = latt_record.code if latt_record is not None else 1
    if latt_record is None:
        diagnostics.append(
            Diagnostic(
                Severity.WARNING,
                "shelx.map.latt_defaulted",
                "LATT is absent; using the documented centrosymmetric primitive default",
                recovery="LATT 1",
            )
        )
    explicit = tuple(
        record.operation for record in records if isinstance(record, ShelxSymmInstruction)
    )
    try:
        operations = build_shelx_operations(latt, explicit)
    except ValueError as error:
        diagnostics.append(
            Diagnostic(
                Severity.ERROR,
                "shelx.map.symmetry_invalid",
                str(error),
                latt_record.span if latt_record is not None else None,
            )
        )
        return StructureCollection(), tuple(diagnostics)

    name, hm_symbol, raw_title = _reported_title(records, document.source_name)
    source_key = document.source_name or "memory"
    structure_id = f"shelx:{source_key}:0"
    part = 0
    residue: ShelxResiInstruction | None = None
    sites: list[IndependentSite] = []
    failed = False
    for record in records:
        if isinstance(record, ShelxPartInstruction):
            part = record.part
            continue
        if isinstance(record, ShelxResiInstruction):
            residue = (
                None
                if record.residue_number in {None, 0} and record.residue_class is None
                else record
            )
            continue
        if not isinstance(record, ShelxAtomRecord):
            continue
        try:
            sites.append(
                _map_atom(
                    record,
                    structure_id=structure_id,
                    site_index=len(sites),
                    sfac=sfac,
                    free_variables=free_variables,
                    part=part,
                    residue=residue,
                    diagnostics=diagnostics,
                )
            )
        except IndexError as error:
            diagnostics.append(
                Diagnostic(
                    Severity.ERROR,
                    "shelx.map.sfac_index_invalid",
                    str(error),
                    record.span,
                )
            )
            failed = True
        except ValueError as error:
            code = (
                "shelx.map.occupancy_invalid"
                if "occupancy" in str(error).casefold() or "fvar" in str(error).casefold()
                else "shelx.map.atom_invalid"
            )
            diagnostics.append(Diagnostic(Severity.ERROR, code, str(error), record.span))
            failed = True
    if not sites:
        diagnostics.append(
            Diagnostic(Severity.ERROR, "shelx.map.atoms_missing", "No valid atom records were found")
        )
        failed = True
    if failed:
        return StructureCollection(), tuple(diagnostics)

    cell = _cell_with_zerr(cell_record, _first(records, ShelxZerrInstruction))
    space_group = SpaceGroupDefinition(
        operations=operations,
        provenance="reported",
        hm_symbol=hm_symbol,
    )
    for site in sites:
        try:
            expand_orbit(site, operations, cell=cell, structure_id=structure_id)
        except (SymmetryConsistencyError, ValueError) as error:
            diagnostics.append(
                Diagnostic(
                    Severity.ERROR,
                    "shelx.map.symmetry_site_invalid",
                    str(error),
                )
            )
            return StructureCollection(), tuple(diagnostics)

    zerr = _first(records, ShelxZerrInstruction)
    crystal = CrystalStructure(
        name=name,
        cell=cell,
        sites=tuple(sites),
        id=structure_id,
        space_group=space_group,
        provenance=StructureProvenance(
            SourceReference(
                source_name=document.source_name,
                format="shelx",
                record_id="structure:0",
                start_offset=0,
                end_offset=len(document.raw_source),
            )
        ),
        metadata={
            "shelx_title": raw_title,
            "shelx_formula_units": (
                None if zerr is None else zerr.formula_units.value
            ),
        },
    )
    return StructureCollection.from_structures((crystal,)), tuple(diagnostics)


__all__ = ["map_shelx_structures"]
