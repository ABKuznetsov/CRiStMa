# Native Structure Readers Roadmap

**Goal:** Complete the CRiStMa reader branch without coupling independent
format implementations or repeatedly running unrelated application tests.

**Spec:** `docs/superpowers/specs/2026-08-30-native-structure-readers-design.md`

## Execution rule

Each row below is an independently reviewable sub-project. Before its code is
changed, it receives a focused implementation plan in this directory. During a
sub-project only its tests and the previously established CRiStMa contracts are
run. The complete CRiStMa suite and wheel audit run only at the final gate.

| Order | Sub-project | Deliverable | Required focused verification |
| --- | --- | --- | --- |
| 1 | Structure model and reader infrastructure | `cristma.structure`, species, atomic properties, collections/sequences, source resolver, lazy registry | `tests/chemistry tests/structure tests/io/test_registry.py tests/io/test_result.py tests/io/test_source.py tests/io/cif` |
| 2 | SHELX | Native RES/INS document, exact symmetry, symbolic occupancy, canonical mapping | `tests/io/shelx tests/symmetry` |
| 3 | VASP | POSCAR/CONTCAR/XDATCAR plus structural output frames | `tests/io/vasp tests/structure` |
| 4 | Quantum ESPRESSO | `pw.x` input/output structural cards and frame sequence | `tests/io/qe tests/structure` |
| 5 | CASTEP | `.cell` and text-output structure frames | `tests/io/castep tests/structure` |
| 6 | CP2K | Nested input structure sections, safe includes, output frames | `tests/io/cp2k tests/io/test_source.py tests/structure` |
| 7 | ABINIT | Input datasets and text-output frames with exact unit handling | `tests/io/abinit tests/structure` |
| 8 | SIESTA | FDF/XV/STRUCT_OUT, safe includes, species mapping | `tests/io/siesta tests/io/test_source.py tests/structure` |
| 9 | XYZ/extXYZ | Basic and extended typed columns, periodic and molecular frames | `tests/io/xyz tests/structure` |
| 10 | PDB and PDBx/mmCIF | Models, hierarchy, altloc, cell/symmetry, bonds | `tests/io/pdb tests/io/cif/test_mmcif_mapper.py tests/structure` |
| 11 | MOL/SDF | V2000/V3000 documents, bonds, properties, multiple records | `tests/io/mol tests/structure` |
| 12 | Reader matrix | Real fixture provenance, detection matrix, built-wheel and dependency audit | `pytest -q` plus clean wheel checks |

## Dependency order

```text
structure model + source/registry
        |
        +-- SHELX
        +-- VASP -- QE -- CASTEP -- CP2K -- ABINIT -- SIESTA
        +-- XYZ
        +-- PDB + PDBx/mmCIF
        `-- MOL/SDF
                 |
                 `-- complete reader matrix
```

Format readers depend only on the public structure and I/O contracts from the
first sub-project. They do not import one another. Shared lexical utilities may
be promoted into `cristma.io.text` only when two implemented readers require
the same behavior and focused tests prove the shared contract.

The later topology subsystem defined in
`docs/superpowers/specs/2026-08-30-structural-hierarchy-design.md` depends on
stable site/atom identity and provenance produced by these readers. Reader
packages do not depend on topology and do not attempt structural grouping.

## Reader-branch completion gate

The branch is ready to finish only when:

- all twelve sub-projects have reviewable commits;
- every claimed format has analytic, malformed, multi-model where applicable,
  and provenance-recorded real fixtures;
- `pytest -q` passes only inside the CRiStMa repository;
- the built wheel imports from a clean environment;
- wheel metadata has no mandatory ASE, Gemmi, pymatgen, RDKit, Open Babel,
  CrysPy, GSAS-II, or Qt dependency;
- `git status --short` is empty.
