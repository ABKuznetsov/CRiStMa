# Native SHELX and CRAFT Integration Design

Date: 2026-08-30

Status: proposed for implementation

## 1. Purpose

CRiStMa will provide a native, Qt-free SHELX RES/INS subsystem that reads,
preserves, interprets, and writes structural documents without Gemmi,
pymatgen, SHELX, or application code. CRAFT will consume only the read side:
RES/INS files become canonical CRiStMa `CrystalStructure` objects and then use
the same temporary CRAFT compatibility projection as CIF.

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
- read-only RES/INS integration in CRAFT;
- analytic, malformed, round-trip, and real-file verification.

This slice does not migrate polyhedron, block, comparison, or topology
algorithms. Those remain separate roadmap stages after structural file inputs
share one stable canonical model.

## 3. Boundaries

```text
SHELX source text
        |
        v
ShelxDocument                 loss-preserving format layer
        |
        v
ShelxMapper                   scientific interpretation
        |
        v
CrystalStructure             canonical CRiStMa structure
        |
        +----> CRAFT adapter -> current viewer/analysis model
        |
        `----> future CRiStMa geometry/hierarchy/comparison tools
```

CRiStMa does not import CRAFT, Sci, Qt, Gemmi, pymatgen, or SHELX. CRAFT
depends directly on CRiStMa. Sci may install compatible packages but does not
participate in application logic.

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

and mapped to the symmetric CRiStMa tensor. Symmetry expansion continues to
transform anisotropic tensors through the existing Structure Core rules.

`PART` and `RESI` establish explicit context for following atom records. The
context is retained as typed record data and mapped to site disorder/group
metadata without pretending that every residue is a crystallographic rigid
body.

### 5.6 Preserved constraints and refinement instructions

The parser recognizes and preserves `AFIX`, `HFIX`, `EXYZ`, `EADP`, `DFIX`,
`DANG`, `SADI`, `SAME`, `FLAT`, `RIGU`, `SIMU`, `DELU`, weighting commands,
`REM`, and unknown commands. They are not converted to canonical CRiStMa
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

## 8. CRAFT integration

CRAFT uses the SHELX subsystem only as a reader:

```text
.res/.ins
   -> cristma.read
   -> canonical CrystalStructure
   -> existing CRiStMa-to-CRAFT compatibility projection
   -> viewer and current analysis
```

The private `_load_shelx` route becomes unreachable and is removed after
equivalent real-file tests pass. No SHELX writer, editor, export command, or
new UI is added to CRAFT.

Source diagnostics appear through the existing CRAFT document-warning path.
The canonical CRiStMa structure and `ShelxDocument` provenance remain available
to future applications even though the current CRAFT compatibility model uses
only the fields required for display and analysis.

## 9. Verification

Focused CRiStMa tests cover:

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

- `.res` and `.ins` use CRiStMa exclusively;
- the real `zdk288.res` renders and enters progressive analysis;
- Q peaks are not displayed as atoms;
- no legacy SHELX parser is called;
- CIF and non-SHELX routes remain unchanged.

The full CRiStMa suite and built-wheel audit run at the library completion gate.
The full CRAFT suite runs once after its read-only integration.

## 10. Acceptance criteria

The slice is complete when:

- CRiStMa reads RES/INS without specialized crystallographic dependencies;
- the source document survives an unchanged preserve round-trip exactly;
- a scientifically valid source yields a canonical asymmetric-unit structure;
- symmetry, occupancy dependencies, disorder, and ADP retain provenance;
- canonical SHELX output reads back to an equivalent CRiStMa structure;
- the real SHELX fixture reaches coordination analysis;
- CRAFT opens RES/INS only through CRiStMa and displays the resulting structure;
- no Qt, Sci, or application dependency enters CRiStMa.

## 11. Subsequent CRiStMa roadmap

After SHELX and the reader branch establish stable canonical inputs, reusable
crystal-chemistry processing moves from CRAFT into CRiStMa in this order:

1. `PolyhedronBuilder` and typed polyhedron results;
2. structural entities, blocks, connectors, and hierarchy;
3. entity matching between related structures;
4. `StructureComparator` and typed `ComparisonResult`;
5. optional kinematic and series analysis over matched entities.

CRiStMa owns comparison calculations and scientific descriptors. CRAFT keeps
the comparison table UI, colours, filters, exports, selection state, and links
between table rows and 3D objects.

