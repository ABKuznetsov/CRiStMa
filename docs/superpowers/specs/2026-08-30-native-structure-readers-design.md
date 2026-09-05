# CrIStMa native structure readers: design specification

Date: 2026-08-30

Status: proposed for user review

## 1. Purpose

CrIStMa shall read common crystallographic, molecular, and electronic-structure
geometry files into its own physics-first scientific model. Applications work
with CrIStMa structures rather than Gemmi, ASE, pymatgen, RDKit, Open Babel, or
code-specific objects.

This reader design is the first input gateway for the shared domain contract
defined in `2026-08-30-cristma-domain-language-vision.md`. Readers do not merely
convert files: they establish stable scientific identity and provenance that
later structure, diffraction, profile, and refinement subsystems consume.

This branch is a reader slice. It includes source parsing, semantic mapping,
diagnostics, and canonical in-memory structures. Except for the already
implemented CIF writer, writing new formats is intentionally deferred.

The supported format families are:

- crystallographic CIF and PDBx/mmCIF;
- SHELX RES/INS;
- VASP POSCAR, CONTCAR, XDATCAR, and structural snapshots in text/XML output;
- Quantum ESPRESSO `pw.x` input and output;
- CASTEP `.cell` and structural blocks in text output;
- CP2K input and text output;
- ABINIT input and text output;
- SIESTA FDF, XV, and STRUCT_OUT;
- XYZ and extended XYZ;
- PDB;
- MOL and SDF, V2000 and V3000.

Only structure-bearing information is interpreted from electronic-structure
files. Energies, forces, stresses, k-point meshes, potentials, convergence
settings, and electronic results remain uninterpreted source records.

## 2. Design principles

1. Scientific meaning is not guessed silently.
2. Source syntax and the canonical internal scientific model are separate
   objects. File-specific representations terminate at the I/O boundary.
3. Periodic crystals and non-periodic molecules share an atomic interface but
   are not forced into one physically false representation.
4. The asymmetric unit remains primary when the source reports symmetry.
5. Explicit simulation atoms remain explicit when the source does not report
   crystallographic symmetry.
6. Multi-model documents preserve every complete model and its order.
7. Large trajectories are indexed and loaded lazily.
8. Common NumPy arrays are allowed; specialized structure libraries are not
   mandatory runtime dependencies.
9. Every derived or recovered value records provenance.
10. An unsupported or incomplete source produces diagnostics, not plausible
    fabricated data.

## 3. Lessons adopted from existing libraries

CrIStMa adopts proven boundaries without adopting another library's internal
model:

- ASE demonstrates a useful central atomic view with cell, per-axis periodic
  flags, typed per-atom arrays, format registration, and lazy frame iteration.
- pymatgen demonstrates the importance of separate periodic `Structure` and
  non-periodic `Molecule` concepts, species-aware sites, and lazy format
  registration.
- Gemmi demonstrates separation between a lossless-ish CIF document,
  small-structure model, and macromolecular hierarchy.
- OVITO demonstrates indexed frame sources, source-frame provenance, and
  explicit mapping of arbitrary file columns to typed particle properties.
- Open Babel demonstrates extensible format handlers and format-specific read
  options, but CrIStMa does not perform implicit chemistry transformations as
  part of reading.

The important CrIStMa addition is that every successful import retains both
the format document and the canonical scientific model.

`CrystalStructure` and `MolecularStructure` are the only canonical scientific
inputs produced by structure readers. Format documents support provenance,
diagnostics, and format-preserving output; scientific calculators must not
branch on them or consume their records as structural state.

## 4. Public structure API

The stable public namespace is:

```python
from cristma.structure import (
    AtomicView,
    CrystalStructure,
    MolecularStructure,
    Structure,
    StructureCollection,
    StructureSequence,
)
```

`Structure` is a protocol implemented by `CrystalStructure` and
`MolecularStructure`. It exposes identity, source provenance, atomic
properties, periodicity, and `atomic_view()`.

The early `cristma.core.structure.Crystal` name becomes a compatibility alias
for `CrystalStructure`. New application code uses `cristma.structure`.

### 4.1 CrystalStructure

`CrystalStructure` contains:

- `UnitCell` represented by a 3x3 lattice and conventional parameters;
- per-axis periodicity;
- `SpaceGroupDefinition` with provenance;
- ordered independent sites when an asymmetric unit is known;
- explicit sites for DFT structures with unreported symmetry;
- derived expanded sites with exact symmetry-operation provenance;
- structure-level and site-level source references;
- typed atomic properties.

Identity symmetry used for explicit DFT atoms has provenance
`unreported_identity`. It does not assert that the physical space group is P1.

### 4.2 MolecularStructure

`MolecularStructure` contains:

- ordered Cartesian atoms;
- explicit bonds, including order, aromaticity, and stereochemical fields when
  reported;
- components, chains, residues, and alternate-location groups when reported;
- conformer/model provenance;
- typed molecular and atomic properties.

No artificial periodic box is created for a molecule.

### 4.3 AtomicView

`AtomicView` is a derived, read-only numerical view used by geometry,
visualization, diffraction, and adapters. It provides:

- chemical species;
- Cartesian coordinates;
- fractional coordinates when a cell exists;
- occupancy;
- per-axis periodicity and optional cell;
- stable atom/site identifiers;
- source and symmetry provenance;
- access to typed per-atom properties.

Creating an `AtomicView` never changes the primary asymmetric-unit or molecular
representation.

### 4.4 Chemical species

Sites refer to a typed species hierarchy:

```text
ChemicalSpecies
|- ElementSpecies
|- IsotopeSpecies
|- ChargedSpecies
`- UnknownSpecies
```

`UnknownSpecies` is used when a format supplies a type index without enough
information to identify the element, as in a VASP 4 POSCAR without a companion
POTCAR. It has a stable source label and blocks element-dependent calculations
until resolved. Existing `.element` access remains available for known
elements.

### 4.5 Atomic properties

Frequently used properties have stable semantic keys:

- occupancy;
- isotropic and anisotropic displacement;
- formal and partial charge;
- scalar or vector magnetic moment;
- mass and isotope;
- selective-dynamics flags;
- velocity;
- residue, chain, component, and source labels.

Additional columns are stored in an `AtomicPropertyTable`. Every property
records dtype, shape, unit, missing-value mask, source name, and provenance.
Unknown extXYZ or simulation columns therefore survive without becoming
untyped dictionaries.

Canonical structures are immutable. Edits use explicit builders or
transformations so refinement history can retain prior states safely.

## 5. Multiple structures and lazy frames

`StructureCollection` is an immutable, tuple-like finite collection used for
CIF blocks, SDF records, and PDB models. It preserves source order and supports:

```python
collection[0]
collection.primary
collection.final
```

`StructureSequence` is a lazy, indexable frame source used for XDATCAR,
extXYZ trajectories, and multi-step DFT outputs. Initial parsing builds a
`FrameIndex` containing byte/character spans and frame metadata. Individual
structures are mapped only when requested.

Both implement a common sequence protocol and expose `primary`, `final`, frame
roles, source-frame numbers, and source spans. The last complete DFT geometry is
marked `final`; incomplete trailing geometries are retained in the format
document but are not returned as structures.

## 6. Read result and source document

The public flow remains:

```python
result = cristma.read(path, format=None)

result.document
result.structures
result.diagnostics
result.source_info
```

`result.structures` is a `StructureCollection` or `StructureSequence`. Both are
sequence-compatible so existing indexing and iteration remain natural.

Each format defines its own immutable `FormatDocument`. At minimum it retains:

- raw source text or binary source metadata;
- token/record order and source spans;
- comments and unknown records;
- original numeric tokens and units;
- format-specific blocks, sections, or frames;
- recoverable malformed records;
- a frame index where applicable.

The source document is not a generic dictionary. It preserves enough format
semantics for diagnostics and future preserve-mode writing.

## 7. Format registry

Each format package owns a lazily imported handler descriptor containing:

- primary registry name and aliases;
- suffix patterns and special basenames;
- content probe;
- text/binary and compression capabilities;
- single/multiple/lazy structure capabilities;
- typed read options;
- parser and mapper entry points.

Selection priority is:

1. explicit `format=`;
2. decisive content signature;
3. special basename such as POSCAR or CONTCAR;
4. suffix;
5. an ambiguity diagnostic if equally plausible handlers remain.

External packages may register adapters through a CrIStMa entry-point group,
but built-in readers do not require those packages. Standard-library gzip,
bzip2, and xz wrappers may be decoded transparently before content probing.

## 8. Format packages

Each format follows one pipeline:

```text
source -> lexer/records -> FormatDocument -> mapper
       -> StructureCollection or StructureSequence
```

### 8.1 CIF and PDBx/mmCIF

The existing CIF 1.1 lexer and document parser remain shared. Semantic mapping
is selected by categories:

- small/inorganic CIF maps cell, symmetry, asymmetric sites, disorder,
  occupancy, oxidation, and displacement parameters;
- PDBx/mmCIF maps atom-site identifiers, models, chains, residues, alternate
  locations, assemblies where structurally relevant, cell, and symmetry.

Multiple data blocks preserve their order. A block may yield zero, one, or
multiple structures and diagnostics.

### 8.2 SHELX RES/INS

One line-oriented parser supports both suffixes. It recognizes `TITL`, `CELL`,
`ZERR`, `LATT`, `SYMM`, `SFAC`, `UNIT`, `FVAR`, atoms, Q peaks, occupancy/free
variable expressions, isotropic and anisotropic displacement, `PART`, `RESI`,
`HKLF`, and `END`.

Constraints and refinement instructions including `AFIX`, `HFIX`, `EXYZ`,
`EADP`, `DFIX`, `DANG`, `SADI`, `SAME`, `FLAT`, `RIGU`, `SIMU`, and `DELU` are
preserved even before they receive canonical constraint objects. `LATT` and
`SYMM` become exact rational affine operations. Free-variable occupancy remains
symbolic as well as evaluated.

### 8.3 VASP

The VASP package supports POSCAR/CONTCAR versions 4 and 5, positive scale,
negative target volume, anisotropic scale where valid, Direct and Cartesian
coordinates, Selective dynamics, velocities when reported, and XDATCAR frames.

Species names are taken only from explicit source data or an explicitly
provided resolver. Missing VASP 4 names produce `UnknownSpecies`; atom counts
are never guessed into elements.

Text/XML output readers extract only complete reported geometries. XML parsing
uses the standard library and may stream large files.

### 8.4 Quantum ESPRESSO

The reader understands `pw.x` namelists and structure cards, including
`ibrav`, `celldm`, A/B/C and angle parameters, `CELL_PARAMETERS`,
`ATOMIC_SPECIES`, and `ATOMIC_POSITIONS` in `alat`, `bohr`, `angstrom`, and
`crystal` units. Output geometries form a lazy sequence.

The mapper records whether a cell is explicit or derived from Bravais-lattice
parameters. Contradictory definitions are diagnostics.

### 8.5 CASTEP

The reader supports `.cell` blocks for lattice vectors/parameters, fractional
or absolute positions, species, units, and common ionic constraints as atomic
properties. Text `.castep`, geometry, and molecular-dynamics structural blocks
are indexed as frames when present. Non-structural blocks remain source data.

### 8.6 CP2K

The reader parses nested sections sufficiently to locate `&CELL`, `&COORD`,
periodicity, units, scaled coordinates, and included structural fragments when
their included text is supplied through an explicit source resolver. Text
output geometries are indexed as frames. The reader does not execute CP2K
preprocessor commands or access arbitrary files implicitly.

### 8.7 ABINIT

The reader maps `acell`, `rprim`/`rprimd`, `angdeg`, `xred`, `xcart`, `xangst`,
`natom`, `ntypat`, `typat`, and `znucl`, including dataset suffixes and unit
conversion. Each complete dataset is a structure. Text output geometries are
indexed in source order.

NetCDF/HDF5 ABINIT files are not mandatory-reader inputs. Without an optional
adapter they produce an explicit unsupported-binary diagnostic.

### 8.8 SIESTA

The reader supports FDF scalar and block syntax, species labels, lattice
constant/vectors/parameters, atomic-coordinate formats, per-axis periodicity,
and `%include` only through an explicit source resolver. XV and STRUCT_OUT map
their reported structures directly. Recursive includes are detected and
reported.

### 8.9 XYZ and extended XYZ

The basic reader supports repeated XYZ frames. The extended reader parses
`Lattice`, `pbc`, `Properties`, quoted metadata, and typed arbitrary columns.
Legacy extra columns require explicit mapping when their meaning cannot be
determined.

A frame with a valid lattice and any periodic axis maps to
`CrystalStructure`; otherwise it maps to `MolecularStructure`.

### 8.10 PDB

The reader recognizes `MODEL/ENDMDL`, `CRYST1`, `SCALE`, `ATOM`, `HETATM`,
`ANISOU`, alternate locations, occupancy, B factors, chains, residues,
insertion codes, elements, charges, and `CONECT`. Unknown records are retained.

A physically valid periodic cell maps to `CrystalStructure` while retaining
macromolecular groups. Otherwise the model maps to `MolecularStructure`.
Alternate locations remain explicit grouped alternatives.

### 8.11 MOL and SDF

The molecular reader supports V2000 and V3000 atom/bond blocks, charges,
isotopes, radicals, stereochemical fields, coordinates, titles, and arbitrary
SDF property blocks. Every SDF record maps to one collection item. Missing bond
information is not perceived automatically during reading.

## 9. Read options and external context

Readers accept typed format-specific options rather than arbitrary keyword
bags. External context is opt-in through a `SourceResolver`, for example:

- associating a POTCAR with a VASP 4 POSCAR;
- resolving SIESTA or CP2K includes;
- selecting legacy XYZ column meanings.

Default reading never scans neighboring files or executes include directives
without permission. Every external value records its source path and span when
available.

ASE, Gemmi, pymatgen, RDKit, and Open Babel may be optional development
comparators or explicit interoperability adapters. No adapter object enters
the canonical model.

## 10. Diagnostics and recovery

Diagnostics retain stable code, severity, message, source span, and recovery.
The following are never silent:

- incomplete coordinate triples or frames;
- unknown species or invalid atomic numbers;
- atom-count mismatches;
- missing or nonphysical cells for periodic models;
- unknown coordinate modes or units;
- conflicting cell definitions;
- invalid symmetry operations;
- occupancy or displacement values outside supported physical domains;
- unresolved includes or external species mappings;
- unsupported binary containers;
- ambiguous format detection;
- loss of a reported atom property during mapping.

A parsed document may be returned with no canonical structure. A structure
containing `UnknownSpecies` may be visualized and edited, but element-dependent
operations such as scattering-factor calculation must reject it explicitly.

## 11. Testing and evidence

Tests are capability-focused and do not import application code. Every reader
has:

1. a minimal analytic fixture;
2. malformed and ambiguous fixtures;
3. unit and coordinate-mode cases;
4. multi-structure or multi-frame coverage where supported;
5. content-probe coverage independent of filename;
6. at least one representative real fixture with source path, copy date, and
   SHA-256 provenance;
7. mapping checks for every claimed semantic field;
8. checks that unknown records and custom properties remain in the document;
9. optional development-only comparisons with an independent reader.

External libraries are comparators, not the expected truth. Format
specifications and analytic fixtures decide disagreements.

Large-frame tests prove that initial read builds an index without materializing
all structures. Thread safety and deterministic repeated frame access are part
of `StructureSequence` acceptance.

## 12. Implementation slices

The reader branch is implemented and committed in independently verifiable
slices:

1. public structure model, species, properties, collections, and compatibility;
2. lazy registry and source/frame infrastructure;
3. SHELX;
4. VASP;
5. Quantum ESPRESSO;
6. CASTEP;
7. CP2K;
8. ABINIT;
9. SIESTA;
10. XYZ/extXYZ;
11. PDB and PDBx/mmCIF mapping;
12. MOL/SDF;
13. real-fixture matrix, clean installation, and dependency audit.

Each slice follows TDD and runs only its focused CrIStMa tests during
development. The complete CrIStMa suite runs at reader-branch completion.

## 13. Acceptance criteria

The reader branch is complete when:

- every format in section 1 can produce a source document and diagnostics;
- every valid supported structure maps to `CrystalStructure` or
  `MolecularStructure` without specialized mandatory dependencies;
- multi-model files preserve every complete model and mark primary/final roles;
- large trajectories are lazy and indexable;
- periodicity, cells, units, species, coordinates, occupancy, displacement,
  bonds, and claimed format-specific properties map correctly;
- unresolved scientific identity is explicit rather than guessed;
- source provenance reaches structures, atoms, and custom properties;
- readers preserve unknown source records for later writing;
- public imports use `cristma.structure` and old `Crystal` code remains
  temporarily compatible;
- the built wheel contains no Qt/application code and requires no Gemmi, ASE,
  pymatgen, RDKit, Open Babel, CrysPy, or GSAS-II package;
- focused and complete CrIStMa tests pass from the built wheel.

## 14. Deferred work

The following remain separate projects:

- writers for formats other than the existing CIF writer;
- calculation inputs beyond their structural geometry;
- energies, forces, stresses, trajectories as simulation-result objects;
- implicit bond perception, aromaticity perception, protonation, and 3D
  coordinate generation;
- mandatory parsing of binary NetCDF, HDF5, DCD, and proprietary containers;
- diffraction, refinement, structure building, and crystal-chemistry analyses
  that consume the new canonical structures.
