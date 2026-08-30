# CRiStMa

CRiStMa is a physics-first Python library for crystallography, crystal
chemistry, structure modelling, diffraction, and refinement.

The library has no GUI dependency. Its native scientific model is independent
of application frameworks and specialized crystallographic packages.

## Native CIF structure I/O

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

`StructureCollection` represents finite multi-model files without losing the
role of a primary or final model. Large trajectories use the same sequence
contract through lazy `StructureSequence`: indexing loads and caches only the
requested frame.

All applications consume the same canonical objects from `cristma.structure`.
Periodic `CrystalStructure` and non-periodic `MolecularStructure` remain
distinct physical models, while both expose an immutable numerical
`AtomicView` for diffraction, topology, visualization, and refinement code.

Ordinary CIF reading and writing are implemented natively. Gemmi, pymatgen,
PyXtal, CrysPy, GSAS-II, and Qt are not required dependencies.
