# CRiStMa

CRiStMa is a physics-first Python library for crystallography, crystal
chemistry, structure modelling, diffraction, and refinement.

The library has no GUI dependency. Its native scientific model is independent
of application frameworks and specialized crystallographic packages.

## Canonical scientific model

Every supported external structure format terminates at the CRiStMa I/O boundary:

```text
CIF / RES / INS / POSCAR / XDATCAR / OUTCAR / vasprun.xml / XYZ / extXYZ
        -> native reader and semantic mapper
        -> CrystalStructure | MolecularStructure
```

`CrystalStructure` and `MolecularStructure` are the canonical scientific inputs
for all structure-based calculations. File-specific documents remain available
for provenance, diagnostics, and preserve-mode writing, but never determine
symmetry, geometry, diffraction, topology, or refinement mathematics after
mapping.

Derived results such as `AtomicView`, `PeriodicNeighborGraph`, and future
polyhedron or hierarchy results remain calculations over the canonical
structure; they do not replace it as the source of truth. Refinement likewise
applies a parameterization to produce a new canonical structure and evaluates
that structure with the same calculators used outside refinement.

Tools are stateless with respect to scientific inputs and results:

```python
crystal = cristma.read("sample.cif").structures[0]
view = expand_structure(crystal)
neighbors = NeighborFinder(cutoff=3.0).find(view)
coordination = CoordinationAnalyzer().analyze(view, neighbors)
```

The caller owns `crystal`, derived results, calculation order, and caching.
CRiStMa has no hidden `current_structure` or `last_result` session state.

CRiStMa does not assign modules to applications. Finder, CRAFT, Rietveld, a
script, or any future consumer may import any public scientific tool it needs.
Applications choose their workflow; CRiStMa supplies the shared implementation.

## Inorganic crystal chemistry

The current crystal-chemistry slice resolves chemically requested contacts,
coordination shells, and three-dimensional coordination polyhedra without
compound-specific branches or expected coordination-number tables:

```python
from cristma.chemistry import ChemistryAnalyzer, Composition
from cristma.crystal_chemistry import (
    CoordinationShellResolver,
    PolyhedronBuilder,
    ShellResolutionPolicy,
)

structure = cristma.read("sample.cif").structures[0]
chemistry = ChemistryAnalyzer().analyze(Composition.from_structure(structure))
policy = ShellResolutionPolicy(1.60, 0.01, 0.08, 0.01, 2.00)
resolution = CoordinationShellResolver(policy).resolve(structure, chemistry.grammar)
view = structure.atomic_view()
polyhedra = tuple(
    PolyhedronBuilder().build(shell, view)
    for shell in resolution.coordination_shells
)
```

The policy is explicit configuration, not a hidden universal preset. CRAFT can
retain these immutable results and render them directly; all workflow and UI
state remain application-owned. Scientific semantics, statuses, provenance,
radii policy, and current limits are documented in
[`docs/inorganic-crystal-chemistry.md`](docs/inorganic-crystal-chemistry.md).

## Structural crystallography

The packaged crystallographic catalog contains all 530 Hall settings and their
Wyckoff positions. Calculated symmetry orbits are the source of multiplicity;
catalog records are used to identify and validate the corresponding Wyckoff
position.

```python
import cristma
from cristma.crystallography import (
    SpaceGroupCatalog,
    assign_wyckoff,
    build_orbit,
)

crystal = cristma.read("sample.cif").structures[0]
setting = SpaceGroupCatalog.default().by_hall("P -4 2ab")

orbit = build_orbit(crystal.sites[0], setting, cell=crystal.cell)
assignment = assign_wyckoff(orbit, setting)

print(orbit.multiplicity)
print(assignment.position.letter if assignment.position else assignment.status)
```

When a CIF omits symmetry operations but identifies one catalog setting
unambiguously, the CIF reader derives the operations and records that decision
as a diagnostic. Explicit source operations remain authoritative and a
catalog disagreement is reported rather than silently corrected.

The normalized catalog is compiled from pinned spglib 2.7.0 data during
development. spglib is not a runtime dependency; installed CRiStMa packages
use their own versioned JSON resources. Source hashes, attribution, and the
BSD-3-Clause notice are stored beside those resources in
`cristma/reference_data/resources/crystallography/`.

## Native structure I/O

The current native readers are:

- CIF 1.1 (`.cif`), with preserving and canonical writing;
- SHELX RES/INS (`.res`, `.ins`), with loss-preserving documents and canonical
  writing from `CrystalStructure`;
- VASP POSCAR/CONTCAR, XDATCAR, structural OUTCAR frames, and `vasprun.xml`.
- plain XYZ and extXYZ, including typed per-atom columns and lazy trajectories.

Other structural formats are added to the same registry as independent native
readers. Applications should never implement their own suffix switch or parser.

```python
from pathlib import Path

import cristma
from cristma.structure import CrystalStructure, StructureCollection, StructureSequence

result = cristma.read("structure.cif")
for diagnostic in result.diagnostics:
    print(diagnostic.severity.value, diagnostic.code, diagnostic.message)

structures: StructureCollection = result.structures
if result.ok and structures:
    # Readers may mark a primary model explicitly; otherwise retain source order.
    crystal: CrystalStructure = structures.primary or structures[0]
    # The asymmetric unit is primary; symmetry-expanded sites are derived.
    print(crystal.cell.volume, [site.label for site in crystal.sites])

    # Emit a normalized CIF from the canonical scientific model.
    cristma.write(crystal, "canonical.cif", mode="canonical")

# Parsed documents retain comments, unknown tags, ordering, and numeric text.
document = cristma.read(Path("structure.cif")).document
cristma.write(document, "preserved.cif", mode="preserve")

# Text sources use the same content-aware format registry.
memory_result = cristma.read_text("data_demo\n_cell_length_a 5\n", format="cif")
```

SHELX uses the same auto-detecting read call. Canonical output requires the
measurement wavelength explicitly because that information does not belong to
the structure alone:

```python
from cristma.io.shelx import ShelxWriteOptions

result = cristma.read("refinement.res")
crystal = result.structures[0]

# Untouched source, including comments, instructions, order, and line endings.
cristma.write(result.document, "preserved.res", mode="preserve")

# New normalized instruction file generated from the canonical structure.
cristma.write(
    crystal,
    "canonical.ins",
    format="shelx",
    mode="canonical",
    options=ShelxWriteOptions(wavelength=0.71073),
)
```

VASP structure sources use the same one-line API. Trajectories are indexed
without creating every `CrystalStructure`; a frame is parsed, mapped, and
cached only when requested:

```python
import cristma

result = cristma.read("XDATCAR")
trajectory = result.structures
final = trajectory.final
print(len(trajectory), final.cell.volume)
```

POSCAR/CONTCAR preserve Selective Dynamics and explicitly reported velocities.
OUTCAR and `vasprun.xml` retain complete structural frames and per-atom forces
with units. A VASP 4 source without species labels produces explicit unknown
species plus a diagnostic; CRiStMa does not guess elements from coordinates.

Native VASP reading requires only NumPy and the Python standard library. This
slice deliberately excludes electronic-result parsing and canonical VASP
writing.

XYZ structure sources use the same registry and lazy sequence contract:

```python
result = cristma.read("trajectory.extxyz")
trajectory = result.structures
final = trajectory.final
print(len(trajectory), final.name, tuple(final.atomic_view().properties))
```

Plain XYZ frames are molecular. An extXYZ frame becomes periodic only when it
reports both a valid `Lattice` and explicit `pbc`; a lattice without `pbc` is
retained as source metadata but is not treated as evidence of periodicity. All
declared `S`, `I`, `R`, and `L` atom columns are typed, and non-structural
columns retain their reported names and provenance without inferred units.
Trajectories may change schema, cell, or molecular/periodic type between
frames. Native XYZ reading uses only NumPy and the Python standard library;
writing and bond inference are outside this reader slice.

`StructureCollection` represents finite multi-model files without losing the
role of a primary or final model. Large trajectories use the same sequence
contract through lazy `StructureSequence`: indexing loads and caches only the
requested frame.

All applications consume the same canonical scientific objects from
`cristma.structure`.
Periodic `CrystalStructure` and non-periodic `MolecularStructure` remain
distinct physical models, while both expose an immutable numerical
`AtomicView` for diffraction, topology, visualization, and refinement code.

CIF, SHELX RES/INS, VASP, and XYZ/extXYZ structural reading are implemented
natively. Gemmi, pymatgen, PyXtal, CrysPy, GSAS-II, SHELX, and Qt are not
required dependencies.
