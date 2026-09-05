# CrIStMa

**Crystallographic Infrastructure for Structures and Materials.**

A compact, physics-first Python library for crystallography, crystal chemistry,
and periodic structure analysis.

> **Status:** public beta `0.1.0b1`<br>
> **Compatibility:** Python 3.11+<br>
> **Runtime dependency:** NumPy<br>
> **License:** BSD-3-Clause

## Overview

CrIStMa provides a common scientific foundation for programs that work with
crystal and molecular structures. It reads widely used structural formats,
maps them into one canonical model, and offers independent tools for symmetry,
geometry, crystal chemistry, and periodic topology.

The project exists because scientific logic is often coupled to a particular
file parser, graphical application, or large external framework. That makes
calculations difficult to reuse, compare, and audit. CrIStMa keeps these layers
separate: formats end at the I/O boundary, scientific operations receive
explicit inputs, and results retain diagnostics and provenance.

CrIStMa is not an end-user application and does not prescribe a workflow. It
is a reusable scientific library for scripts, notebooks, research software,
desktop applications, and automated data-processing systems.

## What it can do

- read CIF, SHELX RES/INS, VASP, PDB, XYZ, and extXYZ structures through one
  content-aware API;
- preserve source documents where supported or write a normalized structure;
- represent periodic crystals and non-periodic molecules as distinct physical
  models;
- expand crystallographic sites using exact symmetry operations;
- use a bundled catalog of all 530 Hall settings and their Wyckoff positions;
- build symmetry orbits and assign Wyckoff positions;
- calculate finite and periodic neighbour graphs and coordination
  environments;
- analyze composition, oxidation-state evidence, coordination shells, and
  coordination polyhedra;
- assemble structural units and classify periodic blocks as finite units,
  chains, layers, or frameworks;
- find translation-aware finite rings in periodic structural
  representations;
- report recoverable problems as structured diagnostics instead of hiding
  assumptions or silently changing the input.

## Scientific model

All supported formats converge on the same native structures:

```text
CIF / RES / INS / POSCAR / XDATCAR / OUTCAR / vasprun.xml / PDB / XYZ / extXYZ
                                |
                                v
             CrystalStructure | MolecularStructure
                                |
                                v
        symmetry / geometry / chemistry / periodic topology
```

The canonical structure is the source of truth for calculations. Parsed
documents remain available for source preservation and provenance, but
file-specific details do not control downstream scientific semantics.

Calculated objects are immutable results rather than hidden application state.
The caller decides calculation order, caching, storage, presentation, and user
interaction.

## Installation

CrIStMa requires Python 3.11 or newer.

After the public beta is published on PyPI:

```bash
python -m pip install --pre cristma
```

To install the current source checkout:

```bash
git clone https://github.com/ABKuznetsov/CrIStMa.git
cd CrIStMa
python -m pip install -e .
```

NumPy is the only runtime dependency. Optional development and reference-data
dependencies are kept outside the scientific runtime.

## Quick start

```python
import cristma
from cristma.geometry import CoordinationAnalyzer, NeighborFinder
from cristma.symmetry import expand_structure

result = cristma.read("sample.cif")

for diagnostic in result.diagnostics:
    print(diagnostic.severity.value, diagnostic.code, diagnostic.message)

if not result.ok or not result.structures:
    raise RuntimeError("The structure could not be read")

crystal = result.structures.primary or result.structures[0]
view = expand_structure(crystal)
neighbors = NeighborFinder(cutoff=3.0).find(view)
coordination = CoordinationAnalyzer().analyze(view, neighbors)

print(crystal.cell.volume)
print(len(view.atoms), len(neighbors.edges))
print(len(coordination.environments))
```

The same entry point reads other supported formats:

```python
structure = cristma.read("POSCAR").structures[0]
trajectory = cristma.read("XDATCAR").structures
model = cristma.read("molecule.pdb").structures[0]
```

## Native structure I/O

| Format | Reading | Writing | Notes |
| --- | --- | --- | --- |
| CIF 1.1 | Yes | Preserve and canonical | Source order, comments, unknown tags, and numeric text can be retained |
| SHELX RES/INS | Yes | Preserve and canonical | Canonical output requires an explicit wavelength |
| VASP POSCAR/CONTCAR | Yes | — | Selective Dynamics and reported velocities are retained |
| VASP XDATCAR | Yes | — | Frames are indexed and loaded lazily |
| VASP OUTCAR | Structural frames | — | Per-atom forces and units are retained |
| `vasprun.xml` | Structural frames | — | Trajectory-oriented structural parsing |
| PDB | Yes | — | Crystal and molecular coordinate models |
| XYZ/extXYZ | Yes | — | Typed properties and lazy trajectories |

CrIStMa implements these readers natively. Gemmi, pymatgen, PyXtal, CrysPy,
GSAS-II, SHELX, and graphical frameworks are not required at runtime.

## Design principles

- **Physics before interface.** Scientific meaning is not determined by a GUI
  or storage format.
- **One canonical model.** Every reader produces the same structure types for
  downstream calculations.
- **Explicit assumptions.** Policies, tolerances, limits, and incomplete
  searches are visible in inputs and results.
- **Traceable results.** Symmetry images, reference data, transformations, and
  diagnostics retain provenance.
- **Composable tools.** Calculators are independent and do not rely on a
  hidden current structure or global workflow.
- **Small runtime.** The core depends only on Python and NumPy.

## Beta status

`0.1.0b1` is the first public beta. The implemented scientific core is covered
by automated tests and is ready for evaluation and integration. Until the
first stable release, public APIs may still change when required to correct or
clarify scientific contracts.

The current beta covers structural I/O, canonical structure models, symmetry,
periodic geometry, crystal chemistry, structural representations, periodic
block classification, and ring analysis. It does not yet calculate diffraction
patterns or perform structure refinement.

## Roadmap

Planned scientific layers are developed as independent milestones:

1. reciprocal metrics, reflection generation, systematic absences, reciprocal
   symmetry orbits, multiplicity, and Friedel relations;
2. scattering contexts and structure-factor calculations;
3. radiation-aware powder lines and physical corrections;
4. calculated diffraction profiles on explicit grids;
5. additional structural transforms, hierarchy and topology tools, and
   refinement built over the same forward calculations.

The roadmap describes direction, not a compatibility or release-date promise.
CrIStMa will remain independent of any particular consuming application.

## Documentation

- [Inorganic crystal chemistry](https://github.com/ABKuznetsov/CrIStMa/blob/main/docs/inorganic-crystal-chemistry.md)
- [Beta release notes](https://github.com/ABKuznetsov/CrIStMa/blob/main/docs/releases/0.1.0b1.md)
- [Third-party notices](https://github.com/ABKuznetsov/CrIStMa/blob/main/THIRD_PARTY_NOTICES.md)

## License and reference data

Original CrIStMa code is distributed under the permissive
[BSD-3-Clause license](https://github.com/ABKuznetsov/CrIStMa/blob/main/LICENSE).
It may be used in open-source, commercial, and closed-source software subject
to the license notice requirements.

Bundled reference resources retain their own attribution and provenance:

- space-group and Wyckoff data normalized from pinned spglib 2.7.0 resources
  under BSD-3-Clause;
- Cordero covalent radii compiled from QCElemental resources under
  BSD-3-Clause;
- Shannon radii compiled from a pinned pymatgen artifact under MIT;
- selected Crystallography Open Database fixtures under CC0/public-domain
  terms;
- curated chemical-reference rules with their scientific literature recorded
  in the versioned resources.

Versions, commits, hashes, known provenance limitations, and redistribution
requirements are listed in
[THIRD_PARTY_NOTICES.md](https://github.com/ABKuznetsov/CrIStMa/blob/main/THIRD_PARTY_NOTICES.md).

## Author

Artem B. Kuznetsov<br>
[GitHub](https://github.com/ABKuznetsov)
