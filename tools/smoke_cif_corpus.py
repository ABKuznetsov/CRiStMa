"""Run the native crystal-chemistry hierarchy over a local structure corpus."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import cristma
from cristma.chemistry import ChemistryAnalyzer, Composition, InteractionLayer
from cristma.crystal_chemistry import (
    ContactClassification,
    CoordinationShellResolver,
    PeriodicConnectivityAnalyzer,
    RingFinder,
    ShellResolutionPolicy,
    StructuralBlockFinder,
    StructuralGraphBuilder,
    StructuralRepresentationBuilder,
    StructuralSelectionPolicy,
    StructuralUnitBuilder,
)
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
    included_classifications=frozenset({ContactClassification.PRIMARY}),
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
        view = structure.atomic_view()
        chemistry = ChemistryAnalyzer().analyze(Composition.from_structure(structure))
        resolution = CoordinationShellResolver(_SHELL_POLICY).resolve(
            structure,
            chemistry.grammar,
        )
        units = StructuralUnitBuilder().build(
            resolution,
            resolution.polyhedra,
            structure=structure,
            atomic_view=view,
        )
        graph = StructuralGraphBuilder().build(units.units, resolution.contacts)
        representation = StructuralRepresentationBuilder(_SELECTION_POLICY).build(graph)
        connectivity = PeriodicConnectivityAnalyzer().analyze(representation)
        blocks = StructuralBlockFinder().find(
            representation,
            connectivity,
            structure=structure,
            atomic_view=view,
        )
        rings = (
            RingFinder().find(structure, view, representation, blocks)
            if include_rings
            else None
        )
        if len({unit.unit_id for orbit in units.unit_orbits for unit in orbit.units}) != len(
            units.units
        ):
            raise ValueError("structural-unit orbit coverage is inconsistent")
        if len({block.block_id for orbit in blocks.block_orbits for block in orbit.blocks}) != len(
            blocks.blocks
        ):
            raise ValueError("structural-block orbit coverage is inconsistent")
        if rings is not None and any(
            not ring.parent_block_orbit_id for ring in rings.rings
        ):
            raise ValueError("ring lacks its parent structural-block orbit")
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
                "atoms": len(view.atoms),
                "contacts": len(resolution.contacts),
                "polyhedra": len(resolution.polyhedra),
                "units": len(units.units),
                "unit_orbits": len(units.unit_orbits),
                "unit_geometries": sum(unit.geometry is not None for unit in units.units),
                "blocks": len(blocks.blocks),
                "block_orbits": len(blocks.block_orbits),
                "rings": None if rings is None else len(rings.rings),
                "ring_orbits": None if rings is None else len(rings.orbits),
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
