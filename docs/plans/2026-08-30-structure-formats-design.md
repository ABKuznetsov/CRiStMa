# CRiStMa structure formats: design specification

Date: 2026-08-30

Status: CIF foundation implemented; remaining reader scope superseded by
`docs/superpowers/specs/2026-08-30-native-structure-readers-design.md`

## 1. Purpose

The first CRiStMa vertical slice provides independent reading and writing of
standard crystal-structure files without requiring Gemmi, pymatgen, PyXtal,
CrysPy, GSAS-II, Qt, or application code.

The initial formats are:

- CIF 1.1 (`.cif`);
- SHELX instruction and result documents (`.ins`, `.res`).

The slice must preserve the original scientific document, build a canonical
CRiStMa crystal model, report every ambiguity, and support round-trip writing.

## 2. Architectural boundaries

CRiStMa is a scientific library. It contains no Qt types, dialogs, application
state, or GUI-specific behavior.

Common scientific packages such as NumPy and SciPy may be normal dependencies.
Specialized crystallographic packages are optional interoperability and test
adapters only. Ordinary CIF/RES/INS reading must work when none is installed.

The dependency direction is:

```text
applications -> CRiStMa
optional adapters -> CRiStMa
CRiStMa core -X-> applications / Gemmi / pymatgen / PyXtal
```

## 3. Two-level import model

Parsing and scientific interpretation are separate operations:

```text
source bytes
    -> format document
    -> semantic mapper
    -> Crystal / StructureCollection
```

The format document is loss-preserving. It retains records that are unknown or
not yet interpreted. The semantic model is normalized and suitable for
symmetry, geometry, diffraction, and refinement.

This separation prevents unsupported tags or instructions from disappearing
when a file is opened and saved.

## 4. Common public concepts

### 4.1 Source and diagnostics

Every read operation returns a `ReadResult`, not a bare structure. It contains:

- one parsed source document;
- zero or more semantic structures;
- warnings and errors with source locations;
- the detected format and format version;
- encoding and newline information needed for faithful writing.

Diagnostics have stable machine-readable codes, severity, message, source
span, and optional recovery information. Scientific recovery is never silent.

### 4.2 Values

Numeric values preserve:

- parsed numeric value;
- standard uncertainty, when present;
- original token;
- unit;
- missing state: absent, unknown (`?`), or inapplicable (`.`).

The original token and the normalized value must coexist. For example,
`7.6959(2)` remains available verbatim while exposing value `7.6959` and
standard uncertainty `0.0002`.

### 4.3 Canonical crystal model

The primary structural representation is an asymmetric unit with exact
symmetry provenance:

```text
Crystal
|- UnitCell
|- SpaceGroupDefinition
|- AsymmetricUnit
|  `- IndependentSite[]
`- derived representations
   |- SymmetryOrbit[]
   |- ExpandedStructure
   |- PeriodicNeighborGraph
   `- later scientific analyses
```

Expanded atoms are derived and cannot become independent editable parameters.
Each expanded atom retains its independent-site ID, symmetry-operation ID, and
integer cell translation.

## 5. Format registry

Readers and writers implement small format-neutral protocols and are registered
by extension, basename, and content probe. Format selection does not depend
only on a filename extension.

The initial registry contains CIF and SHELX. Later additions can include
POSCAR/CONTCAR, PDB, XYZ, MOL/SDF, and application-specific formats without
changing the canonical model or public read API.

## 6. CIF 1.1 document model

The native CIF layer implements the syntax itself. It consists of:

```text
formats/cif/
|- tokens
|- lexer
|- parser
|- document
|- structure mapper
|- writer
`- diagnostics
```

### 6.1 Syntax coverage

The parser supports at least:

- comments;
- `data_` blocks, including multiple blocks per file;
- tags and scalar values;
- `loop_` tables;
- unquoted, single-quoted, and double-quoted values;
- semicolon-delimited multiline text;
- `.` and `?` missing-value semantics;
- case-insensitive reserved words and data names where required;
- line and column tracking;
- preservation of unknown tags and loops.

Malformed input yields diagnostics. Recoverable errors may produce a partial
document; they must not be converted silently into plausible scientific data.

### 6.2 Structural mapping

The first mapper recognizes modern and legacy aliases for:

- unit-cell parameters and uncertainties;
- space-group number, Hermann-Mauguin and Hall symbols;
- explicit symmetry operations;
- atom labels and type symbols;
- fractional coordinates and uncertainties;
- Wyckoff symbol and site multiplicity;
- occupancy;
- isotropic U and B parameters;
- anisotropic displacement tensors;
- disorder assembly and group;
- oxidation states where reported;
- chemical formula and common publication metadata.

Reported symmetry operations are authoritative source data. If operations are
derived from a named group, the derivation and selected setting/origin are
recorded explicitly.

Coincident rows are not merged merely because their coordinates match. Mixed
sites are combined only when labels, disorder metadata, occupancies, and source
semantics justify one physical crystallographic position.

### 6.3 CIF writing

Two writing modes are required:

- preserve mode: edit known values while retaining unknown source content and
  original ordering as far as possible;
- canonical mode: emit a normalized CRiStMa CIF from the semantic model.

CIF 2.0 is a future compatible extension, not part of this first slice.

## 7. SHELX RES/INS document model

RES and INS share one line-oriented SHELX document parser. The parser preserves
instruction order, comments, continuation lines, original tokens, and unknown
instructions.

### 7.1 Initial semantic coverage

The first mapper understands:

- `TITL`, `CELL`, `ZERR`, `LATT`, and `SYMM`;
- `SFAC`, `UNIT`, and `FVAR`;
- atom records and Q peaks;
- occupancy expressions involving free variables;
- isotropic and anisotropic displacement parameters;
- `PART` and `RESI` context;
- `HKLF` and `END` boundaries.

The document layer also preserves, even before full semantic support:

- `AFIX`/`HFIX`;
- `EXYZ`/`EADP`;
- `DFIX`, `DANG`, `SADI`, `SAME`, and `FLAT`;
- `RIGU`, `SIMU`, and `DELU`;
- weighting and refinement instructions;
- `REM` comments and unknown commands.

These retained instructions can later map to CRiStMa constraints for molecular
fragments, rigid polyhedra, structural blocks, shared atoms, occupancies, and
ADPs.

### 7.2 Symmetry and occupancy

`LATT` centering and centrosymmetry, explicit `SYMM` operations, and cell
translations are converted into exact rational affine operations. Floating
point approximations are not the canonical symmetry representation.

SHELX occupancy/free-variable expressions are represented symbolically. Their
evaluated value is available, but the original dependency on `FVAR` must not be
lost.

### 7.3 SHELX writing

As with CIF, the writer provides preserve and canonical modes. Preserve mode is
the default for an imported refinement document; canonical mode is used when a
new CRiStMa structure is exported to SHELX.

## 8. Reuse from existing applications

Existing project code is an algorithm and fixture source, not the new public
API.

Useful behavior to retain from Crystal Blocks/Craft includes:

- CIF measured values and standard uncertainties;
- source scalar/loop preservation;
- modern and legacy symmetry tags;
- mixed-site and disorder handling;
- SHELX `CELL`, `LATT`, `SYMM`, `SFAC`, `FVAR`, occupancy, and Uiso handling;
- existing format probes and representative fixtures.

Useful Organic behavior includes MOL/SDF and SMILES workflows, but RDKit stays
optional and those formats are outside the first slice.

No Qt types, app services, pymatgen structures, Gemmi structures, or RDKit
molecules enter the CRiStMa canonical model.

## 9. Error handling and invariants

The following are errors or explicit warnings, never silent defaults:

- missing or invalid unit-cell parameters;
- incomplete atom coordinate triples;
- invalid loop row width;
- unknown or conflicting symmetry settings;
- disagreement between reported and generated multiplicity;
- occupancy outside an allowed physical range;
- mixed-site occupancies inconsistent with their source model;
- an invalid displacement tensor;
- a SHELX atom referencing an absent `SFAC` or `FVAR` entry.

A partially parsed source document may still be returned for inspection even
when no valid `Crystal` can be constructed.

## 10. Verification strategy

Tests are organized by capability rather than by application:

1. lexer and parser fixtures for valid and malformed syntax;
2. semantic mapping fixtures for cell, symmetry, sites, occupancy, disorder,
   and displacement parameters;
3. read-write-read equivalence tests;
4. preserve-mode tests proving unknown content survives;
5. analytic symmetry fixtures with known orbits and multiplicities;
6. representative real files collected from the existing applications;
7. optional development-only comparisons against Gemmi, pymatgen, and PyXtal.

External libraries are comparators, not the expected truth. Disagreements are
resolved against the format specification and analytic fixtures.

## 11. Acceptance criteria for the first slice

The slice is complete when:

- CIF 1.1 and RES/INS files can be read without specialized crystallographic
  dependencies;
- an imported file yields a source document, diagnostics, and canonical
  asymmetric-unit structure when scientifically possible;
- important known fields are mapped without discarding their raw form;
- unknown CIF tags and SHELX instructions survive preserve-mode round trips;
- symmetry expansion has exact provenance;
- mixed occupancy and displacement data survive read-write-read;
- failures identify the offending source location and do not fabricate data;
- the API has no Qt or application dependency.

## 12. Deferred work

The following are intentionally deferred while keeping extension points:

- CIF 2.0 and dictionary-driven validation;
- full refinement semantics for every SHELX instruction;
- POSCAR/CONTCAR, PDB, XYZ, MOL/SDF, and SMILES;
- periodic neighbor graphs and polyhedra;
- diffraction and refinement engines;
- optional interoperability adapters.
