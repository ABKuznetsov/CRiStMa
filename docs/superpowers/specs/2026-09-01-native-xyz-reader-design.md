# CrIStMa native XYZ and extXYZ reader design

Date: 2026-09-01

Status: approved in chat; written specification awaiting final user review

Primary format reference: https://github.com/libAtoms/extxyz

## 1. Purpose

CrIStMa shall read ordinary XYZ and schema-declared extended XYZ sources
without ASE, `extxyz`, Open Babel, or another chemistry package. All complete
frames map through the existing canonical structure model and become available
to every application through the unchanged call:

```python
result = cristma.read("trajectory.xyz")
structure = result.structures[0]
```

The source format ends at the I/O boundary. Geometry and later
crystal-chemical calculations consume `MolecularStructure` or
`CrystalStructure`, never XYZ parser records.

This slice includes:

- ordinary four-column XYZ;
- schema-declared extXYZ per-frame and per-atom data;
- finite, fully periodic, and partially periodic frames;
- lazy multi-frame access;
- exact source retention, diagnostics, typed properties, and provenance;
- registry detection and public `cristma.read` integration.

Canonical XYZ/extXYZ writing, inferred bonds, unit conversion, and
application-specific trajectory state are outside this slice.

## 2. Governing rules

1. Atom count and the following frame extent are explicit; atom-looking text
   outside an indexed frame is never interpreted as structure data.
2. Plain XYZ is molecular unless the source contains an extXYZ schema and an
   explicit periodicity declaration.
3. CrIStMa deliberately requires both `Lattice` and explicit `pbc` before it
   creates a periodic structure. The extXYZ reference implementation commonly
   defaults `pbc` to true when `Lattice` is present, but CrIStMa does not infer
   physical periodicity from a box alone.
4. Every extXYZ column declared by `Properties` is validated and preserved.
5. Unspecified units remain unspecified. Names such as `forces` or `charge`
   do not by themselves prove a unit.
6. Frames are independent scientific snapshots. One source may contain both
   molecular and periodic frames or changing cells and schemas.
7. Parsing and mapping do not mutate source documents or canonical structures.
8. Applications never branch on `.xyz`; the format registry owns detection and
   dispatch.

## 3. Public and internal data flow

```text
decoded XYZ/extXYZ text
        |
        v
XyzDocument + tuple[XyzFrameSpan]
        |
        v  selected lazily
XyzFrame
        |
        v
map_xyz_frame(...)
        |
        +--> MolecularStructure
        `--> CrystalStructure
```

`ReadResult.document` is an `XyzDocument`. A one-frame source may still return
an eager `StructureCollection`; a multi-frame source returns a lazy
`StructureSequence`. The implementation may use `StructureSequence` for all
frame counts if this keeps loading behavior uniform, but the public scientific
result must retain the existing collection/sequence contracts.

## 4. Source document and frame contracts

The XYZ package provides immutable records:

```text
XyzDocument
|- raw_source
|- source_name
`- frames: tuple[XyzFrameSpan, ...]

XyzFrameSpan
|- index
|- atom_count
|- start_offset
|- end_offset
|- comment_start_offset
|- comment_end_offset
`- atom_rows_start_offset

XyzPropertySpec
|- name
|- type: string | integer | real | logical
`- width: positive integer

XyzFrame
|- name
|- atom_count
|- comment
|- metadata
|- schema
|- columns
|- lattice
|- pbc
`- source
```

Offsets are half-open Python character offsets into `raw_source`, including
when Unicode occurs in a comment. The document can reproduce the original text
exactly. Frame indexing reads only the count line, comment extent, and declared
number of following atom rows. It does not allocate numeric arrays for every
frame.

## 5. Ordinary XYZ

An ordinary frame contains:

```text
N
arbitrary comment
species x y z
... N rows
```

The first atom token may be an element symbol or an integral atomic number.
Coordinates are Cartesian and retain no inferred unit. The whole comment line
is preserved as frame metadata under `comment`; it is not required to parse as
key/value syntax.

If `Properties` is absent, only the first four columns have defined ordinary
XYZ semantics. Additional tokens are retained by the source document but are
not converted into unnamed scientific properties. A diagnostic explains that
schema-free trailing columns were left uninterpreted.

Ordinary XYZ maps to `MolecularStructure` even when its comment happens to
contain a `Lattice`-looking token but lacks a valid extXYZ `Properties` schema.

## 6. extXYZ frame metadata

The second line is parsed as quoted or bare `key=value` entries following the
libAtoms extXYZ grammar. The parser retains reported spelling and values while
also exposing typed metadata values where the grammar is unambiguous:

- integer;
- real;
- logical;
- string;
- one- or two-dimensional arrays of those primitive types.

Unknown keys are retained. `Properties`, `Lattice`, and `pbc` have structural
semantics. Other keys such as `energy`, `time`, or `step` remain per-frame
metadata and do not become hidden application state.

Malformed quoting, duplicate special keys, or a value whose required shape is
wrong produces a diagnostic. The parser never silently falls back from a
broken extXYZ declaration to ordinary XYZ merely because that recovery would
produce four readable columns.

## 7. extXYZ per-atom schema

`Properties` is a colon-separated sequence of triplets:

```text
name : type : width
```

Supported reported types are:

- `S`: string;
- `I`: integer;
- `R`: real;
- `L`: logical.

Width must be a positive integer. Property names must be non-empty and unique.
The sum of widths is the exact required atom-row column count. Every row must
match it; otherwise that frame is invalid rather than partially shifted into a
different schema.

Structural columns are:

- `species:S:1`, or `Z:I:1`, for chemical identity;
- `pos:R:3` for Cartesian coordinates.

At least one chemical-identity column and exactly one `pos:R:3` property are
required. If both `species` and `Z` are reported, their normalized elements
must agree row by row. A contradiction invalidates the frame. An unknown
symbol or unsupported atomic number becomes `UnknownSpecies` and emits a
diagnostic; it is not guessed from coordinates or neighboring atoms.

All remaining schema properties become immutable `AtomicProperty` objects.
Widths greater than one produce `(N, width)` arrays; width one produces `(N,)`
arrays. Reported names remain property names. Known names may receive stable
semantic aliases only in a future separately specified compatibility layer;
this reader does not rename user columns or infer units.

Each property stores source name, source field, and method `reported` in
`PropertyProvenance`.

## 8. Periodicity and cell mapping

`Lattice` is interpreted as nine row-major real values representing three
Cartesian lattice vectors. `pbc` is exactly three logical values.

Mapping rules are:

```text
valid Lattice + explicit pbc with at least one true axis
    -> CrystalStructure(periodic=pbc)

valid Lattice + explicit pbc=(False, False, False)
    -> MolecularStructure, lattice retained in metadata

valid Lattice + no pbc
    -> MolecularStructure, lattice retained in metadata,
       diagnostic xyz.map.periodicity_unspecified

true pbc axis + absent or invalid Lattice
    -> invalid frame
```

For a periodic frame the mapper derives the six `UnitCell` parameters from the
reported lattice, then rotates reported Cartesian positions and polar-vector
properties into CrIStMa's canonical cell frame. A property is rotated only if
its scientific transformation semantics are explicitly declared by the
reader. This slice declares `pos` as a Cartesian position; arbitrary width-3
properties are not assumed to be vectors and remain numerically reported.

Fractional positions are computed from the reported lattice before mapping to
the canonical cell orientation. Identity-only symmetry is attached because
XYZ/extXYZ reports explicit atoms rather than an asymmetric unit and space
group.

## 9. Molecular atomic properties

`MolecularStructure` gains an optional immutable `AtomicPropertyTable`, with
the same atom-count validation already enforced by `CrystalStructure`.
`MolecularStructure.atomic_view()` passes that table into `AtomicView`.

This is a general canonical-model capability, not an XYZ-specific field. It is
needed for molecular forces, velocities, charges, labels, and future MOL/SDF
properties without creating parallel per-format containers.

## 10. Multi-frame behavior

The index contains only complete frames. A `FrameReference` records frame
index, source span, and basic count metadata. All complete frames except the
last have role `intermediate`; the last complete frame has role `final`.

Loading a frame parses only its indexed comment and atom rows, maps it to a
canonical structure, and caches the result through `StructureSequence`.
Repeated access returns the same object. Frames may change atom count, schema,
periodicity, or lattice independently.

An incomplete final frame is retained in `XyzDocument.raw_source`, excluded
from `StructureSequence`, and reported as `xyz.frame.incomplete`. A malformed
complete frame remains indexed so its source identity is visible, but loading
it raises a precise mapping error; structural corruption is never represented
as a partially valid canonical structure.

## 11. Diagnostics

Stable diagnostic families include:

```text
xyz.frame.count_invalid
xyz.frame.incomplete
xyz.frame.blank_between_frames
xyz.comment.syntax_invalid
xyz.schema.invalid
xyz.schema.duplicate_property
xyz.row.column_count_mismatch
xyz.row.value_invalid
xyz.map.position_missing
xyz.map.species_missing
xyz.map.species_conflict
xyz.map.species_unresolved
xyz.map.lattice_invalid
xyz.map.pbc_invalid
xyz.map.periodicity_unspecified
xyz.map.uninterpreted_plain_columns
```

Invalid source data encountered while producing objects yields diagnostics or
a frame-load error. Invalid canonical CrIStMa objects continue to raise their
own validation exceptions. Diagnostics contain source spans whenever the
offending token or row has an indexed location.

## 12. Registry and detection

The built-in descriptor is lazy and named `xyz`, with alias `extxyz` and
suffixes `.xyz` and `.extxyz`. The lightweight probe recognizes a coherent
integer count, comment line, and enough plausible atom rows. A valid extXYZ
`Properties` declaration receives higher confidence than ordinary XYZ.

Constructing `builtin_format_descriptors()` must not import the XYZ parser or
mapper. `XyzFormatHandler` owns parsing and mapping. No application contains a
suffix branch or XYZ-specific parser.

## 13. Verification

Focused analytic tests cover:

- ordinary element-symbol and atomic-number XYZ;
- exact document preservation and Unicode offsets;
- valid `S`, `I`, `R`, and `L` schema columns;
- arbitrary scalar and vector properties;
- `species`/`Z` agreement and conflict;
- molecular, fully periodic, and partially periodic mapping;
- the conservative `Lattice`-without-`pbc` rule;
- changing multi-frame schemas and cells;
- lazy access and caching;
- truncated frames and malformed rows;
- descriptor laziness and public `cristma.read`.

Provenance-recorded fixtures verify plain XYZ, molecular extXYZ, and periodic
extXYZ. Cross-format tests compare equivalent canonical geometry and
coordination with CIF/POSCAR where periodic, and later with MOL/SDF where
molecular.

The final gate runs the complete CrIStMa suite once, builds a wheel without
network dependencies, installs it into a clean temporary environment, and
reads all XYZ fixtures from outside the repository.

## 14. Dependencies and future work

The implementation uses only the Python standard library and NumPy already
required by CrIStMa. It does not import ASE, libAtoms `extxyz`, pymatgen,
Open Babel, RDKit, Qt, or application code.

MOL/SDF is the next independent structure-format project. Crystal-chemical
analysis follows the format projects and consumes only canonical structures,
their atomic views, and explicit graph/property results.
