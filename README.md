# CRiStMa

CRiStMa is a physics-first Python library for crystallography, crystal
chemistry, structure modelling, diffraction, and refinement.

The library has no GUI dependency. Its native scientific model is independent
of application frameworks and specialized crystallographic packages.

## Canonical scientific model

Every supported external structure format terminates at the CRiStMa I/O boundary:

```text
CIF / RES / INS / future registered structure source
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

## Native structure I/O

The current native readers are:

- CIF 1.1 (`.cif`), with preserving and canonical writing;
- SHELX RES/INS (`.res`, `.ins`), with loss-preserving documents and canonical
  writing from `CrystalStructure`.

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

`StructureCollection` represents finite multi-model files without losing the
role of a primary or final model. Large trajectories use the same sequence
contract through lazy `StructureSequence`: indexing loads and caches only the
requested frame.

All applications consume the same canonical scientific objects from
`cristma.structure`.
Periodic `CrystalStructure` and non-periodic `MolecularStructure` remain
distinct physical models, while both expose an immutable numerical
`AtomicView` for diffraction, topology, visualization, and refinement code.

CIF and SHELX RES/INS reading and writing are implemented natively. Gemmi,
pymatgen, PyXtal, CrysPy, GSAS-II, SHELX, and Qt are not required dependencies.
