# CRiStMa Structural Crystallography Design

**Date:** 2026-09-01  
**Status:** approved  
**Scope:** space-group catalog, Wyckoff positions, crystallographic orbits,
site symmetry, multiplicity, and strictly geometric local structure

## 1. Scientific boundary

This milestone implements only the crystallography of an atomic structure:

```text
Crystallography
├── Symmetry
│   ├── SpaceGroupCatalog
│   ├── WyckoffPosition
│   ├── CrystallographicOrbit
│   ├── SiteSymmetry
│   └── CalculatedMultiplicity
│
└── Local Geometry
    ├── periodic neighbour enumeration
    ├── distances and angles
    └── geometric coordination shells
```

The milestone explicitly excludes:

- crystal morphology, faceting, BFDH and twins;
- confirmed chemical bonds and bond orders;
- bond-valence analysis, polyhedra and structural hierarchy;
- reflections, diffraction and refinement;
- application projects, UI state and workflows.

The boundary between scientific layers is strict:

```text
Chemistry
    proposes meaningful interaction kinds from composition

Crystallography
    determines symmetry and geometric neighbourhoods

Crystal Chemistry
    interprets those neighbourhoods as bonds, coordination
    polyhedra, BVS, blocks and topology
```

A geometric neighbour is not automatically a chemical bond.

## 2. Canonical data flow

```text
IndependentSite + SpaceGroupDefinition
                    |
                    v
            CrystallographicOrbit
            ├── expanded positions
            ├── calculated multiplicity
            ├── stabilizer
            ├── site symmetry
            ├── matched Wyckoff position
            └── diagnostics
```

The actual orbit is the source of truth for calculated multiplicity. Catalog
values and values reported by a source file are validation targets:

```text
actual unique orbit size
        vs
catalog Wyckoff multiplicity
        vs
reported CIF multiplicity
```

CRiStMa never mutates the immutable `IndependentSite` or
`CrystalStructure` while performing this analysis. Applications may construct
a new snapshot if they wish to persist derived values.

## 3. Reference data

Crystallographic tables belong to Reference Data; algorithms that use them
belong to Crystallography.

```text
cristma/reference_data/resources/crystallography/
├── space_groups.json
├── wyckoff_positions.json
└── SOURCE.md
```

Every packaged dataset contains or is accompanied by:

```text
dataset_id
version
source
license
compiled_date
checksum
```

The resources are inert, versioned data shipped inside the CRiStMa wheel.
There is no runtime network access and no runtime dependency on Gemmi,
pymatgen, spglib or CrysPy.

`SpaceGroupCatalog` loads, validates and indexes these resources. It does not
contain application policy and does not infer structure from experimental
data.

### 3.1 Provenance and external code

The primary source is the BSD-3-Clause licensed spglib `v2.7.0` release,
commit `12355c77fb7c505a55f52cae36341d73b781a065`:

```text
spglib/database/spg.csv
spglib/database/Wyckoff.csv
spglib database-generation scripts required to interpret them
```

A development-only compiler converts these sources into CRiStMa's normalized,
validated JSON schema. The generated resources retain the spglib copyright and
BSD-3-Clause notice. `SOURCE.md` records the upstream repository, exact commit,
input checksums, compiler command, compilation date, output checksums and
scientific conventions. The compiler is deterministic: identical pinned
inputs produce byte-identical JSON.

PyXtal is the secondary cross-check after the provenance of the relevant CSV
files is recorded. Bilbao Crystallographic Server and Gemmi are independent
scientific oracles for selected fixtures. Their databases are not repackaged.
International Tables are not copied or systematically reproduced.

Suitable pymatgen code may be adapted under its MIT license with the required
copyright and license notice. Gemmi is dual-licensed under MPL-2.0 or LGPLv3;
it is used as a comparator rather than copied into the CRiStMa implementation.

None of spglib, PyXtal, pymatgen or Gemmi becomes a runtime dependency. A
catalog is never copied from an undocumented installed package.

## 4. Space-group identity

A space-group number alone does not identify a usable catalog entry. The
canonical identity is setting-sensitive and origin-sensitive.

```text
SpaceGroupKey
├── hall_number
├── hall_symbol
└── choice
```

`choice` preserves spglib's setting, unique-axis, cell-choice and origin-choice
code without discarding distinctions. The Hall symbol is the primary
crystallographic lookup key; `hall_number` is the pinned dataset key. Number
and Hermann-Mauguin symbol are searchable aliases and may return several
choices. The caller must resolve ambiguity explicitly.

`SpaceGroupRecord` contains:

```text
key
number
hm_symbol
hall_symbol
setting
origin_choice
crystal_system
centering
operations
wyckoff_positions
reference provenance
```

The existing `SpaceGroupDefinition` remains the compact symmetry object stored
inside `CrystalStructure`. A catalog record produces a
`SpaceGroupDefinition`; it does not introduce a second competing structure
model.

```python
record = SpaceGroupCatalog.default().by_hall_symbol("P 1")
definition = record.definition(provenance="derived")
```

Lookup by number or Hermann-Mauguin symbol returns a tuple of records unless
the caller supplies enough setting/origin information to select exactly one.

## 5. Wyckoff positions

`WyckoffPosition` is an immutable typed record:

```text
WyckoffPosition
├── space_group_key
├── letter
├── multiplicity
├── site_symmetry_symbol
├── coordinate representatives
├── degrees of freedom
└── reference provenance
```

Coordinate representatives use exact affine expressions and CRiStMa's
existing `AffineOperation` machinery wherever possible. Fractions remain
exact until numerical comparison with measured coordinates is required.

Wyckoff matching is performed against the selected setting and origin. A
letter from another setting is not silently accepted as equivalent.

## 6. Crystallographic orbit

`CrystallographicOrbit` is a derived result, not an alternative structure:

```text
CrystallographicOrbit
├── source_site_id
├── atoms: tuple[ExpandedAtom, ...]
├── calculated_multiplicity
├── stabilizer: tuple[SymmetryImageProvenance, ...]
├── site_symmetry: SiteSymmetry
├── wyckoff_position: WyckoffPosition | None
└── diagnostics
```

`OrbitAnalyzer` is a configurable stateless tool:

```python
orbit = OrbitAnalyzer(tolerance=1e-6).analyze(
    site,
    space_group,
    catalog=SpaceGroupCatalog.default(),
    cell=crystal.cell,
    structure_id=crystal.id,
)
```

It reuses the established exact symmetry expansion. No second implementation
of position expansion is allowed.

### 6.1 Multiplicity

Calculated multiplicity is the number of unique positions in the reference
cell after applying all operations and merging periodically equivalent
positions within the configured tolerance.

### 6.2 Stabilizer

The stabilizer contains symmetry images that map the source position onto
itself modulo an integer lattice translation. Its provenance retains both the
operation ID and normalization translation. A periodic neighbour translation
is a different concept and is not used here.

### 6.3 Site symmetry

`SiteSymmetry` contains the stabilizer operations and the catalog symbol.
The first milestone validates a catalog symbol against the computed
stabilizer; it does not invent Hermann-Mauguin symbols from floating-point
matrices using an unverified heuristic.

### 6.4 Wyckoff matching

Matching proceeds by:

1. selecting the exact catalog record by Hall symbol/setting/origin;
2. filtering Wyckoff positions by calculated multiplicity;
3. comparing the orbit and stabilizer to catalog representatives under
   periodic equivalence and the explicit tolerance;
4. returning one match, unresolved, or ambiguous with diagnostics.

The analyzer does not choose an arbitrary candidate when several positions
remain indistinguishable.

## 7. Diagnostics

Scientific inconsistencies are explicit diagnostics, including:

```text
crystallography.space_group.catalog_entry_missing
crystallography.space_group.lookup_ambiguous
crystallography.orbit.reported_multiplicity_mismatch
crystallography.orbit.wyckoff_multiplicity_mismatch
crystallography.orbit.reported_wyckoff_mismatch
crystallography.orbit.wyckoff_unresolved
crystallography.orbit.wyckoff_ambiguous
crystallography.orbit.site_symmetry_mismatch
```

Invalid CRiStMa objects and invalid tool configuration raise exceptions.
Unexpected or incomplete source-file information encountered during mapping or
scientific comparison produces diagnostics.

## 8. Geometric local structure

The existing `NeighborFinder`, `PeriodicNeighborGraph` and coordination
objects remain the basis of local geometry. This milestone may add reusable
distance, angle and shell results, but only with geometric semantics:

```text
GeometricCoordinationShell
├── center
├── neighbours
├── distances
├── coordination_number
└── search provenance
```

No radius table, electronegativity rule or BVS cutoff is allowed to convert a
geometric shell into confirmed bonds inside Crystallography. Such
interpretation belongs to Crystal Chemistry.

## 9. Package layout

Only implemented modules are created:

```text
cristma/
├── crystallography/
│   ├── __init__.py
│   ├── space_group.py
│   ├── wyckoff.py
│   ├── catalog.py
│   ├── orbit.py
│   └── local_geometry.py       # only when additional shell API is needed
│
├── symmetry/
│   ├── affine.py               # exact operations remain here
│   ├── displacement.py
│   └── orbit.py                # low-level expansion remains here
│
└── reference_data/resources/crystallography/
```

`symmetry` remains low-level mathematical machinery.
`crystallography` exposes scientific catalog and analysis tools built on it.

## 10. Verification gate

Focused analytic fixtures precede broad real-file tests:

- P1 general position;
- P-1 general and inversion-center special positions;
- a group with nontrivial setting/origin alternatives;
- tetragonal P-42\u2081m positions used by the gehlenite samples;
- rounded special coordinates within tolerance and outside tolerance;
- reported multiplicity and Wyckoff mismatches;
- deterministic catalog lookup and checksum validation.

Selected results are compared independently with pymatgen and Gemmi where
they expose equivalent functionality. These packages are development-time
comparators only.

The milestone is complete when a real CIF can travel through:

```text
CRiStMa reader
    -> CrystalStructure
    -> exact catalog-resolved SpaceGroupDefinition
    -> CrystallographicOrbit for every IndependentSite
    -> calculated multiplicity / site symmetry / Wyckoff validation
    -> geometric AtomicView and PeriodicNeighborGraph
```

with no external crystallographic runtime dependency and with all scientific
ambiguity visible in diagnostics.
