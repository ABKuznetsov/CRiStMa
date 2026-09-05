# Native SHELX and Universal CRAFT Structure Input Design

Date: 2026-08-30

Status: proposed for implementation

## 1. Purpose

CrIStMa will provide a native, Qt-free SHELX RES/INS subsystem that reads,
preserves, interprets, and writes structural documents without Gemmi,
pymatgen, SHELX, or application code. CrIStMa is an independent toolbox of
reusable crystallographic data types, readers, writers, transformations, and
calculators. It does not own application workflows or interpretations.

CRAFT will not own a reader switch for CIF, SHELX, VASP, PDB, XYZ, or future
structure formats. Every structure source goes through one auto-detecting
`cristma.read(path)` boundary and becomes a canonical CrIStMa
`CrystalStructure` or `MolecularStructure`. CRAFT then projects that canonical
object into its temporary presentation model and owns display and
cross-structure comparison. XPFF remains application-owned because it is a
Finder project container rather than a single scientific structure format.

The subsystem is library functionality, not a CRAFT editing feature. Preserve
and canonical writing exist for future Builder, Organica, refinement, and other
applications.

## 2. Scope

This slice includes:

- one line-oriented parser for `.res` and `.ins`;
- a loss-preserving `ShelxDocument`;
- physical-line and logical-record provenance;
- semantic mapping to canonical `CrystalStructure`;
- exact `LATT`/`SYMM` symmetry;
- symbolic and evaluated FVAR occupancy;
- isotropic and anisotropic displacement parameters;
- PART/RESI context;
- preserved refinement, restraint, comment, and unknown instructions;
- preserve-mode and canonical-mode writing;
- format-registry integration and public API;
- the format-neutral CRAFT structure-input boundary;
- analytic, malformed, round-trip, and real-file verification.

The CRAFT switch is completed only after native CrIStMa readers cover every
structure format currently advertised by CRAFT: CIF/mmCIF, RES/INS,
POSCAR/CONTCAR/VASP, PDB, and XYZ. Their reader implementations remain separate
CrIStMa sub-projects from the structure-reader roadmap. This slice does not yet
implement reusable polyhedron, block, or topology tools. Cross-structure
comparison and comparison-table logic remain application responsibilities of
CRAFT.

## 3. Boundaries

```text
structure source
        |
        v
cristma.read(path)            registry auto-detects the format
        |
        +----> format document and diagnostics
        |
        v
CrystalStructure / MolecularStructure
        |
        +----> temporary CRAFT adapter -> viewer presentation model
        |
        `----> independent CrIStMa geometry/hierarchy tools
```

CrIStMa does not import CRAFT, Sci, Qt, Gemmi, pymatgen, or SHELX. CRAFT
depends directly on CrIStMa. Sci may install compatible packages but does not
participate in application logic.

The canonical internal-model invariant is:

```text
CIF / RES / INS / POSCAR / PDB / XYZ
        -> cristma.read(path)
        -> CrystalStructure | MolecularStructure
```

Adding another structural format means adding one independent reader,
descriptor, and mapper to CrIStMa. It must not require a new parser, suffix
branch, or scientific model in CRAFT or another consuming application.

The boundary is governed by two mandatory rules:

> **Applications own project formats; CrIStMa owns structural formats.**

> **No application-specific structural parser, mapper, writer, or format
> registry is permitted outside CrIStMa.**

An application may unpack its own project container, but every embedded
structural payload crosses the same boundary:

```text
XPFF / RMP / application package       application responsibility
        |
        v
embedded CIF / RES / other structure
        |
        v
cristma.read(...) or cristma.read_text(...)
        |
        v
CrystalStructure | MolecularStructure
```

These are canonical scientific models, not canonical file formats. Once a
source is mapped, its format may affect provenance and later export but must not
affect symmetry, geometry, crystal chemistry, diffraction, hierarchy, topology,
or refinement mathematics. Derived views/results never replace the canonical
structure as the source of truth.

CrIStMa tools receive canonical objects and explicit derived inputs as
arguments and return explicit results. They do not retain a current structure,
last result, application cache, or workflow session. Applications own result
reuse and invalidation; for example CRAFT may retain one calculated neighbor
graph and pass it to coordination and polyhedron tools.

The ownership contract is:

| CrIStMa provides | CRAFT owns |
| --- | --- |
| structure and derived scientific entities | windows, panels, tables, and 3D actors |
| parsing, normalization, and validation | file-selection and interaction workflow |
| neighbors, coordination, polyhedra, and blocks | visibility, colours, camera, and selection |
| per-structure topology and reusable geometry calculators | cross-structure matching and comparison |
| scientific diagnostics and provenance | interpretation and presentation of results |

This table separates ownership, not access. CRAFT, Finder, Rietveld Manager,
and future consumers may use any public CrIStMa geometry, chemistry, hierarchy,
topology, diffraction, refinement, transform, analysis, or I/O tool. CrIStMa
does not maintain per-application capability sets.

CrIStMa tools are independent and do not know which application calls them or
why. CRAFT composes the tools needed for its own workflow. It must not develop a
second implementation of a reusable structural operation already supplied by
CrIStMa, but it retains application-specific comparison, mechanics narrative,
UI state, and workflow logic.

## 4. Document model

`ShelxDocument` retains:

- the original decoded source;
- encoding and original newline style;
- every physical line in source order;
- blank lines and comments;
- continuation markers and continuation lines;
- logical records assembled from physical lines;
- original tokens and source spans;
- records before and after `HKLF` and `END`;
- unknown and currently uninterpreted instructions;
- non-destructive source edits for preserve writing.

The document distinguishes physical layout from logical meaning. For example,
an anisotropic atom split with `=` is one `ShelxAtomRecord` with two physical
source lines. Reformatting is never required merely to interpret it.

Initial logical record types include:

```text
ShelxInstructionRecord
ShelxAtomRecord
ShelxQPeakRecord
ShelxCommentRecord
ShelxBlankRecord
ShelxUnknownRecord
```

Instruction names are case-insensitive for interpretation. Original spelling
and spacing remain source data.

## 5. Scientific mapping

### 5.1 Cell and measurement information

`CELL` maps the six unit-cell parameters to `UnitCell`. Its first numeric value
is retained as a measured radiation wavelength on the `ShelxCellInstruction`:

```text
ShelxCellInstruction
|- wavelength: MeasuredValue
`- cell: UnitCell
```

Wavelength is not stored inside `CrystalStructure`, because it belongs to the
measurement rather than the material structure. The imported document retains
it for writing and future creation of a format-neutral diffraction `Radiation`
object.

`ZERR` is retained with its reported formula-unit count and cell uncertainties.
The first structure slice does not invent cell uncertainties when `ZERR` is
absent.

### 5.2 Symmetry

`LATT` and explicit `SYMM` records produce exact rational
`AffineOperation` values. The mapper implements SHELX centring codes and the
sign convention for centrosymmetric versus non-centrosymmetric lattices.

The identity operation is included exactly once. Centrosymmetric inversion and
centring translations are combined deterministically with explicit operations.
Equivalent operations are deduplicated by exact rational coefficients, not
floating-point coordinates.

### 5.3 Elements and atoms

`SFAC` supports both element lists and one-record-per-element scattering-factor
forms. Canonical sites require a resolvable element; missing or invalid SFAC
references are errors with source spans.

Atom records before the structural termination boundary map to
`IndependentSite`. Fractional coordinates are retained as measured values and
may lie outside `[0, 1)` in the reported document; canonical periodic geometry
normalizes only when generating expanded positions.

Q peaks remain typed `ShelxQPeakRecord` objects and never become chemical atoms
or independent crystallographic sites. Records following `END` remain in the
document but do not enter the canonical structure.

### 5.4 Occupancy and free variables

SHELX occupancy is represented symbolically and numerically:

```text
ShelxOccupancyExpression
|- raw: str
|- free_variable_index: int
|- multiplier: float
|- complement: bool
`- evaluated: float
```

Codes using fixed occupancy, an FVAR value, or its complement are decoded
without discarding their dependency. The evaluated physical occupancy enters
`SiteComponent`; the symbolic expression remains in component metadata and in
the source document for future refinement and canonical writing.

An absent referenced FVAR, a non-finite expression, or an evaluated value
outside the physical range is an error. Recovery is never silent.

### 5.5 Displacement and disorder

One displacement value maps to `U_iso`. Six SHELX anisotropic values are
interpreted in the reported order:

```text
U11 U22 U33 U23 U13 U12
```

and mapped to the symmetric CrIStMa tensor. Symmetry expansion continues to
transform anisotropic tensors through the existing Structure Core rules.

`PART` and `RESI` establish explicit context for following atom records. The
context is retained as typed record data and mapped to site disorder/group
metadata without pretending that every residue is a crystallographic rigid
body.

### 5.6 Preserved constraints and refinement instructions

The parser recognizes and preserves `AFIX`, `HFIX`, `EXYZ`, `EADP`, `DFIX`,
`DANG`, `SADI`, `SAME`, `FLAT`, `RIGU`, `SIMU`, `DELU`, weighting commands,
`REM`, and unknown commands. They are not converted to canonical CrIStMa
constraints until a separately specified constraint layer can represent their
semantics without loss.

## 6. Writing

### 6.1 Preserve mode

```python
result = cristma.read("model.res")
cristma.write(result.document, "copy.res", mode="preserve")
```

An unchanged document writes the original text exactly, including newline
style, comments, continuations, unknown instructions, and content after `END`.
Supported edits replace only their explicit source spans.

### 6.2 Canonical mode

```python
from cristma.io.shelx import ShelxWriteOptions

cristma.write(
    crystal,
    "model.ins",
    format="shelx",
    mode="canonical",
    options=ShelxWriteOptions(wavelength=0.71073),
)
```

Canonical output contains a minimal complete structural instruction document:

```text
TITL
CELL
LATT / SYMM
SFAC
UNIT
FVAR when symbolic dependencies are retained
atom records with Uiso or Uaniso
HKLF
END
```

`UNIT` is calculated from symmetry-expanded unit-cell contents. A bare
`CrystalStructure` does not contain a wavelength, so `ShelxWriteOptions`
requires it. When canonicalizing an imported `ShelxDocument`, its measured
wavelength can be reused explicitly by the caller or writer helper.

Unknown species, unavailable wavelength, unresolved FVAR expressions, invalid
ADP tensors, and structures that cannot be represented in the selected SHELX
form are explicit errors.

### 6.3 Public dispatch

Top-level `cristma.write` gains additive format and options keywords. Document
types select their preserve writer automatically. Canonical structure output
requires an explicit format when more than CIF is possible.

Dedicated functions and types remain available from `cristma.io.shelx` for
users who want format-specific control.

## 7. Diagnostics

Diagnostics use stable `shelx.lex.*`, `shelx.parse.*`, and `shelx.map.*` codes
with physical source spans. A parsed document may be returned when scientific
mapping fails. No canonical structure is returned for an invalid cell,
unresolved element, invalid FVAR dependency, invalid occupancy, or unusable
displacement model.

Recoverable conditions, such as a preserved but unsupported instruction, are
warnings or information only when they do not alter the scientific structure.

## 8. Universal CRAFT structure input

CRAFT uses one structure reader call:

```text
CIF / RES / INS / POSCAR / PDB / XYZ / future registered format
   -> cristma.read(path)
   -> canonical CrystalStructure / MolecularStructure
   -> temporary CrIStMa-to-CRAFT compatibility projection
   -> viewer presentation and currently unmigrated calculations
```

CRAFT does not pass `format="cif"` or `format="shelx"`; the CrIStMa registry
selects by content, suffix, or basename. CRAFT obtains supported extensions and
basenames from CrIStMa rather than maintaining a second format matrix. Its only
separate file route is XPFF project loading. If XPFF contains embedded CIF or
another structure source, CRAFT locates and extracts that payload but delegates
its parsing and mapping to `cristma.read_text`.

The private `_load_shelx`, `_load_vasp`, `_load_pdb`, `_load_xyz`, and
CIF-specific application routes become unreachable and are removed only after
equivalent real-file tests pass. No SHELX writer, editor, export command, or new
UI is added to CRAFT.

Source diagnostics appear through the existing CRAFT document-warning path.
The canonical CrIStMa structure and `ShelxDocument` provenance remain available
to future applications even though the current CRAFT compatibility model uses
only the fields required for display and analysis.

Existing CRAFT readers and reusable single-structure algorithms are migration
sources for independent CrIStMa tools. CRAFT keeps its application-specific
comparison and mechanics workflows, including the decisions about which
structures and derived entities to compare and how to explain the result.

## 9. Verification

Focused CrIStMa tests cover:

- physical lines, continuations, comments, blank lines, and unknown records;
- valid and malformed core instructions;
- exact P1, inversion, and all supported centring translations;
- SFAC variants and invalid indices;
- fixed, FVAR, and complement occupancy expressions;
- Uiso and correct Uaniso component order;
- PART/RESI and Q-peak exclusion;
- exact preserve round-trip;
- canonical write-read scientific equivalence;
- real `tests/fixtures/shelx/zdk288.res` provenance fixture;
- `RES -> CrystalStructure -> expanded atoms -> periodic neighbors -> coordination`.

Focused CRAFT tests prove:

- every advertised structure format uses the same auto-detecting CrIStMa call;
- CRAFT contains no format-specific structure parser or format dispatch;
- the real `zdk288.res` renders and enters progressive analysis;
- Q peaks are not displayed as atoms;
- molecular structures receive only an application display projection and are
  not converted into fake crystals inside CrIStMa;
- XPFF remains an explicit application-project route.

The full CrIStMa suite and built-wheel audit run at the library completion gate.
The full CRAFT suite runs once after its read-only integration.

## 10. Acceptance criteria

The slice is complete when:

- CrIStMa reads RES/INS without specialized crystallographic dependencies;
- the source document survives an unchanged preserve round-trip exactly;
- a scientifically valid source yields a canonical asymmetric-unit structure;
- symmetry, occupancy dependencies, disorder, and ADP retain provenance;
- canonical SHELX output reads back to an equivalent CrIStMa structure;
- the real SHELX fixture reaches coordination analysis;
- CRAFT opens all advertised structure formats through the same CrIStMa call
  and displays the resulting canonical structure;
- no Qt, Sci, or application dependency enters CrIStMa.

## 11. Subsequent CrIStMa roadmap

After SHELX and the complete reader branch establish stable canonical inputs,
reusable single-structure processing is added to CrIStMa in this order:

1. `PolyhedronBuilder` and typed polyhedron results;
2. structural entities, blocks, connectors, and hierarchy;
3. reusable per-structure topology and geometry descriptors;
4. additional independent transforms and calculators needed by multiple
   applications.

CRAFT owns entity matching between different structures, comparison metrics,
mechanical interpretation across a series, and the comparison table. It uses
CrIStMa structures and per-structure results as inputs. Another application may
compose the same tools differently without importing CRAFT or inheriting its
comparison semantics.
