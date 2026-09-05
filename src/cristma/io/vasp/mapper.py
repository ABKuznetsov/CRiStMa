"""Map format-neutral VASP snapshots into canonical CrIStMa structures."""

from __future__ import annotations

import math

import numpy as np

from cristma.core import MeasuredValue, UnitCell
from cristma.structure import (
    AtomicProperty,
    AtomicPropertyTable,
    IndependentSite,
    PropertyProvenance,
    SiteComponent,
    StructureProvenance,
)

from .document import VaspSnapshot


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


def _canonical_rotation(source_lattice: np.ndarray, canonical_lattice: np.ndarray) -> np.ndarray:
    rotation = np.linalg.solve(source_lattice, canonical_lattice)
    if not np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-10, rtol=1e-10):
        raise ValueError("VASP lattice and canonical cell are not related by a rotation")
    return rotation


def _rotate_rows(rows: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    result = np.asarray(rows, dtype=float) @ rotation
    return np.where(np.abs(result) < 1e-12, 0.0, result)


def map_vasp_snapshot(snapshot: VaspSnapshot):
    """Build one identity-symmetry canonical crystal from a VASP frame."""

    cell = _cell_from_lattice(snapshot.lattice)
    rotation = _canonical_rotation(snapshot.lattice, cell.matrix)
    sites = tuple(
        IndependentSite(
            id=f"site:{index + 1}",
            label=f"{species.label}{index + 1}",
            components=(SiteComponent(species, _reported(1.0, "fraction")),),
            fractional=tuple(
                _reported(float(value), "fraction") for value in snapshot.fractional[index]
            ),
        )
        for index, species in enumerate(snapshot.species)
    )
    source_name = snapshot.source.source_name
    properties: list[AtomicProperty] = []
    if snapshot.selective_dynamics is not None:
        properties.append(
            AtomicProperty(
                "selective_dynamics",
                snapshot.selective_dynamics,
                source_name=source_name,
                provenance=PropertyProvenance(source_name, "selective_dynamics", "reported"),
            )
        )
    if snapshot.velocities is not None:
        values = snapshot.velocities
        if snapshot.velocity_mode == "cartesian":
            values = _rotate_rows(values, rotation)
        properties.append(
            AtomicProperty(
                "velocity",
                values,
                unit=snapshot.velocity_unit,
                source_name=source_name,
                provenance=PropertyProvenance(
                    source_name,
                    f"velocity:{snapshot.velocity_mode}",
                    "reported",
                ),
            )
        )
    if snapshot.forces is not None:
        properties.append(
            AtomicProperty(
                "force",
                _rotate_rows(snapshot.forces, rotation),
                unit=snapshot.force_unit,
                source_name=source_name,
                provenance=PropertyProvenance(source_name, "force", "reported"),
            )
        )

    from cristma.structure import CrystalStructure

    return CrystalStructure.explicit(
        snapshot.name,
        cell,
        sites,
        id=f"vasp:frame:{snapshot.frame_index}",
        provenance=StructureProvenance(snapshot.source),
        properties=AtomicPropertyTable(len(sites), tuple(properties)),
        metadata={"vasp_frame_index": snapshot.frame_index},
    )


__all__ = ["map_vasp_snapshot"]
