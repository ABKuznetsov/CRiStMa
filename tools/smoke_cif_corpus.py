"""Run the native crystal-chemistry hierarchy over a local structure corpus."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import cristma
from cristma.chemistry import ChemistryAnalyzer, Composition, InteractionLayer
from cristma.crystal_chemistry import (
    ContactAnalyzer,
    PeriodicConnectivityAnalyzer,
    PolyhedronOrbitBuilder,
    RingFinder,
    ShellRole,
    ShellResolutionPolicy,
    StructuralBlockFinder,
    StructuralGraphBuilder,
    StructuralRepresentationBuilder,
    StructuralSelectionPolicy,
    StructuralUnitBuilder,
)
from cristma.crystallography import SymmetryContext
from cristma.structure import CrystalStructure


_SUPPORTED_SUFFIXES = frozenset({".cif", ".ins", ".res"})
_SHELL_POLICY = ShellResolutionPolicy(1.60, 0.01, 0.08, 0.01, 2.0)
_SELECTION_POLICY = StructuralSelectionPolicy(
    included_layers=frozenset(
        {
            InteractionLayer.STRUCTURAL,
            InteractionLayer.COORDINATION,
            InteractionLayer.INTRA_SUBSYSTEM,
        }
    ),
    included_shell_roles=frozenset({ShellRole.PRIMARY}),
)


def discover_structure_files(directory: Path) -> tuple[Path, ...]:
    """Return supported files in deterministic, case-insensitive order."""

    if not directory.is_dir():
        raise NotADirectoryError(directory)
    return tuple(
        sorted(
            (
                item
                for item in directory.iterdir()
                if item.is_file() and item.suffix.casefold() in _SUPPORTED_SUFFIXES
            ),
            key=lambda item: (item.name.casefold(), item.name),
        )
    )


def analyze_file(path: Path, *, include_rings: bool = True) -> list[dict[str, object]]:
    """Calculate hierarchy records for every crystal structure in one file."""

    loaded = cristma.read(path)
    if not loaded.structures:
        raise ValueError("reader returned no structures")
    output: list[dict[str, object]] = []
    for index, structure in enumerate(loaded.structures):
        if not isinstance(structure, CrystalStructure):
            raise TypeError("corpus entry is not a CrystalStructure")
        if not structure.sites:
            output.append(
                {
                    "file": path.name,
                    "structure_index": index,
                    "status": "not_applicable",
                    "sites": 0,
                    "atoms": 0,
                    "contacts": 0,
                    "polyhedra": 0,
                    "units": 0,
                    "unit_orbits": 0,
                    "unit_geometries": 0,
                    "blocks": 0,
                    "rings": 0 if include_rings else None,
                    "ring_orbits": 0 if include_rings else None,
                    "diagnostics": sorted({item.code for item in loaded.diagnostics}),
                }
            )
            continue
        chemistry = ChemistryAnalyzer().analyze(Composition.from_structure(structure))
        context = SymmetryContext.from_definition(structure.space_group, structure.cell)
        resolution = ContactAnalyzer(_SHELL_POLICY).analyze(
            structure, context, chemistry.grammar,
        )
        polyhedra = PolyhedronOrbitBuilder().build(resolution)
        units = StructuralUnitBuilder().build(resolution, polyhedra)
        graph = StructuralGraphBuilder().build(resolution, polyhedra)
        representation = StructuralRepresentationBuilder(_SELECTION_POLICY).build(graph)
        connectivity = PeriodicConnectivityAnalyzer().analyze(representation)
        blocks = StructuralBlockFinder().find(representation, connectivity)
        rings = (
            RingFinder().find(representation, blocks)
            if include_rings
            else None
        )
        diagnostic_codes = sorted(
            {
                item.code
                for item in (
                    *loaded.diagnostics,
                    *resolution.diagnostics,
                    *units.diagnostics,
                    *blocks.diagnostics,
                    *((rings.diagnostics if rings is not None else ())),
                )
            }
        )
        output.append(
            {
                "file": path.name,
                "structure_index": index,
                "status": resolution.status.value,
                "sites": len(structure.sites),
                "atoms": len(structure.atomic_view().atoms),
                "contacts": len(resolution.contacts),
                "polyhedra": len(polyhedra.polyhedra),
                "units": len(units.unit_orbits),
                "unit_orbits": len(units.unit_orbits),
                "unit_geometries": sum(
                    unit.representative_geometry is not None
                    for unit in units.unit_orbits
                ),
                "blocks": len(blocks.blocks),
                "rings": None if rings is None else sum(
                    ring.multiplicity_in_reference_cell for ring in rings.ring_orbits
                ),
                "ring_orbits": None if rings is None else len(rings.ring_orbits),
                "diagnostics": diagnostic_codes,
            }
        )
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--skip-rings", action="store_true")
    arguments = parser.parse_args(argv)
    files = discover_structure_files(arguments.directory)
    failures = 0
    structures = 0
    for path in files:
        try:
            rows = analyze_file(path, include_rings=not arguments.skip_rings)
        except Exception as error:  # keep the remaining corpus observable
            failures += 1
            print(json.dumps({"file": path.name, "error": str(error)}, ensure_ascii=False))
            continue
        structures += len(rows)
        for row in rows:
            print(json.dumps(row, ensure_ascii=False, sort_keys=True))
    print(
        json.dumps(
            {
                "files": len(files),
                "structures": structures,
                "failures": failures,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
