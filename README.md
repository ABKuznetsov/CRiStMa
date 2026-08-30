# CRiStMa

CRiStMa is a physics-first Python library for crystallography, crystal
chemistry, structure modelling, diffraction, and refinement.

The library has no GUI dependency. Its native scientific model is independent
of application frameworks and specialized crystallographic packages.

## Native CIF structure I/O

```python
from pathlib import Path

import cristma

result = cristma.read("structure.cif")
for diagnostic in result.diagnostics:
    print(diagnostic.severity.value, diagnostic.code, diagnostic.message)

if result.ok and result.structures:
    crystal = result.structures[0]
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

Ordinary CIF reading and writing are implemented natively. Gemmi, pymatgen,
PyXtal, CrysPy, GSAS-II, and Qt are not required dependencies.
