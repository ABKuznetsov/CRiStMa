"""Map parsed PDB coordinate records into canonical structures."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from cristma.chemistry import ElementSpecies, UnknownSpecies, normalize_element
from cristma.core import MeasuredValue, UnitCell
from cristma.crystallography import SpaceGroupCatalog
from cristma.structure import (
    CrystalStructure,
    DisplacementParameters,
    IndependentSite,
    MolecularAtom,
    MolecularStructure,
    SiteComponent,
    SourceReference,
    StructureProvenance,
)
from cristma.symmetry import SpaceGroupDefinition, parse_xyz_operation

from .document import PdbAtomRecord, PdbDocument


def _reported(value: float, raw: str, unit: str) -> MeasuredValue:
    return MeasuredValue(float(value), None, raw.strip(), unit)


def _name(document: PdbDocument) -> str:
    if document.source_name:
        return Path(document.source_name).stem
    return "PDB structure"


def _source(document: PdbDocument, model_number: int) -> SourceReference:
    return SourceReference(
        document.source_name,
        "pdb",
        f"model:{model_number}",
        0,
        len(document.raw_source),
    )


def _species(record: PdbAtomRecord):
    try:
        return ElementSpecies(normalize_element(record.element))
    except ValueError:
        return UnknownSpecies(
            f"pdb:{record.element or record.name}",
            record.element or record.name,
        )


def _component(record: PdbAtomRecord) -> SiteComponent:
    return SiteComponent(
        _species(record),
        _reported(record.occupancy, f"{record.occupancy:.2f}", "fraction"),
    )


def _displacement(record: PdbAtomRecord) -> DisplacementParameters | None:
    if record.b_iso is None:
        return None
    return DisplacementParameters(
        "B_iso",
        isotropic=_reported(record.b_iso, f"{record.b_iso:.2f}", "angstrom^2"),
        reported_kind="B_iso",
    )


def _normalized_symbol(value: str) -> str:
    return "".join(value.casefold().split())


def _space_group(symbol: str) -> SpaceGroupDefinition:
    normalized = _normalized_symbol(symbol)
    matches = tuple(
        setting
        for setting in SpaceGroupCatalog.default().settings
        if normalized in {
            _normalized_symbol(setting.hm_short),
            _normalized_symbol(setting.hm_full),
        }
    )
    if len(matches) == 1:
        return matches[0].definition(provenance="reported")
    return SpaceGroupDefinition(
        operations=(parse_xyz_operation("x,y,z", operation_id="op:1"),),
        provenance="identity_fallback",
        hm_symbol=symbol or None,
    )


def map_pdb_document(
    document: PdbDocument,
    *,
    model_number: int = 1,
) -> CrystalStructure | MolecularStructure:
    """Map the parsed coordinate model without consulting application state."""

    selected_atoms = tuple(
        atom for atom in document.atoms if atom.model_number == model_number
    )
    provenance = StructureProvenance(_source(document, model_number))
    name = _name(document)
    if document.cryst1 is None:
        atoms = tuple(
            MolecularAtom(
                id=f"atom:{record.serial}",
                label=record.name,
                components=(_component(record),),
                cartesian=record.cartesian,
                metadata={
                    "pdb_record": record.record_name,
                    "residue_name": record.residue_name,
                    "chain_id": record.chain_id,
                    "residue_sequence": record.residue_sequence,
                    "alternate_location": record.alternate_location,
                },
            )
            for record in selected_atoms
        )
        return MolecularStructure(
            name,
            atoms,
            id=f"model:{model_number}",
            provenance=provenance,
        )

    cell_record = document.cryst1
    cell = UnitCell(
        _reported(cell_record.a, f"{cell_record.a:.3f}", "angstrom"),
        _reported(cell_record.b, f"{cell_record.b:.3f}", "angstrom"),
        _reported(cell_record.c, f"{cell_record.c:.3f}", "angstrom"),
        _reported(cell_record.alpha, f"{cell_record.alpha:.2f}", "degree"),
        _reported(cell_record.beta, f"{cell_record.beta:.2f}", "degree"),
        _reported(cell_record.gamma, f"{cell_record.gamma:.2f}", "degree"),
    )
    inverse = np.linalg.inv(cell.matrix)
    sites = tuple(
        IndependentSite(
            id=f"site:{record.serial}",
            label=record.name,
            components=(_component(record),),
            fractional=tuple(
                _reported(value, f"{value:.16g}", "fraction")
                for value in np.asarray(record.cartesian, dtype=float) @ inverse
            ),
            displacement=_displacement(record),
            metadata={
                "pdb_record": record.record_name,
                "residue_name": record.residue_name,
                "chain_id": record.chain_id,
                "residue_sequence": record.residue_sequence,
                "alternate_location": record.alternate_location,
            },
        )
        for record in selected_atoms
    )
    return CrystalStructure(
        name,
        cell,
        sites,
        id=f"model:{model_number}",
        space_group=_space_group(cell_record.space_group),
        provenance=provenance,
        metadata={
            "pdb_declared_space_group": cell_record.space_group,
            "pdb_z": cell_record.z,
        },
    )


__all__ = ["map_pdb_document"]
