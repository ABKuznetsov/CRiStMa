"""Map format-neutral XYZ frames into canonical CRiStMa structures."""

from __future__ import annotations

import math

import numpy as np

from cristma.chemistry import ElementSpecies, UnknownSpecies, element_from_atomic_number, normalize_element
from cristma.core import MeasuredValue, UnitCell
from cristma.structure import (
    AtomicProperty,
    AtomicPropertyTable,
    CrystalStructure,
    IndependentSite,
    MolecularAtom,
    MolecularStructure,
    PropertyProvenance,
    SiteComponent,
    StructureProvenance,
)

from .document import XyzFrame


_STRUCTURAL_PROPERTIES = frozenset({"species", "Z", "pos"})


def _reported(value: float, unit: str) -> MeasuredValue:
    return MeasuredValue(float(value), None, f"{float(value):.16g}", unit)


def _cell_from_lattice(lattice: np.ndarray) -> UnitCell:
    lengths = np.linalg.norm(lattice, axis=1)

    def angle(left: np.ndarray, right: np.ndarray) -> float:
        cosine = float(np.dot(left, right) / (np.linalg.norm(left) * np.linalg.norm(right)))
        return math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))

    return UnitCell(
        _reported(lengths[0], "angstrom"),
        _reported(lengths[1], "angstrom"),
        _reported(lengths[2], "angstrom"),
        _reported(angle(lattice[1], lattice[2]), "degree"),
        _reported(angle(lattice[0], lattice[2]), "degree"),
        _reported(angle(lattice[0], lattice[1]), "degree"),
    )


def _species(frame: XyzFrame, index: int):
    if "species" in frame.columns:
        reported = str(frame.columns["species"][index])
        try:
            return ElementSpecies(normalize_element(reported))
        except ValueError:
            if "Z" not in frame.columns:
                return UnknownSpecies(f"xyz:{reported}", reported)
    if "Z" in frame.columns:
        reported_number = int(frame.columns["Z"][index])
        try:
            return ElementSpecies(element_from_atomic_number(reported_number))
        except ValueError:
            return UnknownSpecies(f"xyz:Z:{reported_number}", str(reported_number))
    return UnknownSpecies(f"xyz:atom:{index + 1}", None)


def _properties(frame: XyzFrame) -> AtomicPropertyTable:
    source_name = frame.source.source_name
    properties = tuple(
        AtomicProperty(
            item.name,
            frame.columns[item.name],
            unit=None,
            source_name=source_name,
            provenance=PropertyProvenance(
                source_name,
                f"Properties:{item.name}",
                "reported",
            ),
        )
        for item in frame.schema
        if item.name not in _STRUCTURAL_PROPERTIES
    )
    return AtomicPropertyTable(frame.atom_count, properties)


def _metadata(frame: XyzFrame) -> dict[str, object]:
    metadata = dict(frame.metadata)
    metadata["xyz_comment"] = frame.comment
    if frame.lattice is not None:
        metadata["xyz_lattice"] = tuple(tuple(float(value) for value in row) for row in frame.lattice)
    if frame.pbc is not None:
        metadata["xyz_pbc"] = frame.pbc
    return metadata


def _component(frame: XyzFrame, index: int) -> SiteComponent:
    return SiteComponent(_species(frame, index), _reported(1.0, "fraction"))


def map_xyz_frame(frame: XyzFrame) -> MolecularStructure | CrystalStructure:
    """Map one parsed XYZ frame without consulting its original file format."""

    positions = np.asarray(frame.columns.get("pos"), dtype=float)
    if positions.shape != (frame.atom_count, 3):
        raise ValueError("XYZ frame requires exactly one pos:R:3 property")
    provenance = StructureProvenance(frame.source)
    properties = _properties(frame)
    metadata = _metadata(frame)

    if frame.pbc is None or not any(frame.pbc):
        atoms = tuple(
            MolecularAtom(
                id=f"atom:{index + 1}",
                label=f"{_species(frame, index).label}{index + 1}",
                components=(_component(frame, index),),
                cartesian=tuple(float(value) for value in positions[index]),
            )
            for index in range(frame.atom_count)
        )
        return MolecularStructure(
            frame.name,
            atoms,
            id=frame.source.record_id,
            provenance=provenance,
            properties=properties,
            metadata=metadata,
        )

    if frame.lattice is None:
        raise ValueError("Periodic XYZ frame requires an explicit Lattice")
    cell = _cell_from_lattice(frame.lattice)
    fractional = positions @ np.linalg.inv(frame.lattice)
    sites = tuple(
        IndependentSite(
            id=f"site:{index + 1}",
            label=f"{_species(frame, index).label}{index + 1}",
            components=(_component(frame, index),),
            fractional=tuple(_reported(value, "fraction") for value in fractional[index]),
        )
        for index in range(frame.atom_count)
    )
    return CrystalStructure.explicit(
        frame.name,
        cell,
        sites,
        id=frame.source.record_id,
        periodic=frame.pbc,
        provenance=provenance,
        properties=properties,
        metadata=metadata,
    )


__all__ = ["map_xyz_frame"]
