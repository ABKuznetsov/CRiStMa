# Native VASP Structure Reader Design

**Status:** Approved design for implementation  
**Date:** 2026-08-31  
**Parent specification:** `2026-08-30-native-structure-readers-design.md`

## 1. Purpose

CrIStMa will read the complete family of text VASP structure sources without
ASE, pymatgen, Qt, or application code:

- POSCAR and CONTCAR;
- XDATCAR;
- OUTCAR structural frames;
- `vasprun.xml` structural frames.

All successfully mapped frames become canonical `CrystalStructure` snapshots.
POSCAR and CONTCAR yield a finite one-structure collection. Multi-step sources
yield a lazy `StructureSequence`; the final complete frame has role `final`.
Applications continue to call only:

```python
result = cristma.read(path)
```

File grammar, probing, frame discovery, species resolution, coordinate
conversion, units, diagnostics, and provenance belong only to CrIStMa.

In this milestone, laziness applies to parsing and mapping trajectory frames.
The existing text registry may still decode and retain the complete source in
memory. File-backed span reading is a later I/O infrastructure improvement and
does not delay the VASP reader.

## 2. Architectural boundary

The VASP package is one registered format family with independent parsers:

```text
cristma/io/vasp/
├── document.py
├── poscar.py
├── xdatcar.py
├── outcar.py
├── vasprun.py
├── mapper.py
├── probe.py
├── handler.py
└── __init__.py
```

The parsers share immutable format records and small numeric utilities, not
parser state. They communicate with the rest of CrIStMa only through existing
I/O, chemistry, structure, and provenance types.

No parser imports another scientific format. No application contains a VASP
suffix switch, parser, mapper, or VASP-specific canonical structure.

## 3. Format registration and detection

One lazy `FormatDescriptor` is registered as `vasp`, with aliases describing
the recognized source families. Detection combines content and source name:

- exact basenames `POSCAR`, `CONTCAR`, `XDATCAR`, and `OUTCAR`;
- basename `vasprun.xml` and characteristic VASP XML root/content;
- POSCAR-like lattice/count/coordinate grammar when the name is unavailable;
- XDATCAR configuration markers;
- OUTCAR version and structural block markers.

Content probing must recognize renamed files. A basename is evidence but must
not cause a high-confidence false positive for content belonging to another
registered format. Parser and mapper modules remain unloaded until the VASP
descriptor is selected.

The descriptor advertises text input, multiple structures, and lazy frames.

## 4. Common scientific rules

### 4.1 Canonical snapshots

Every mapped frame is an immutable `CrystalStructure`. It contains:

- a validated `UnitCell` in angstrom;
- explicit identity-only symmetry with provenance `unreported_identity`;
- independent sites with stable source-derived IDs;
- normalized element species when known;
- fractional coordinates as the canonical periodic coordinates;
- VASP-specific per-site values in typed `AtomicPropertyTable` columns or
  metadata only when no shared scientific field exists;
- source offsets, source frame index, and reported coordinate mode.

VASP files do not report a crystallographic space group. The reader must not
infer one during import.

### 4.2 Scale factor and cell

The POSCAR scale line is interpreted exactly:

- one positive scalar multiplies every lattice vector and Cartesian position;
- one negative scalar denotes the requested positive cell volume; the required
  isotropic multiplier is computed from the absolute reported volume and raw
  lattice determinant;
- three positive scalars independently scale the x, y, and z Cartesian
  components of every lattice vector and Cartesian position as defined by
  VASP;
- zero, non-finite values, an invalid negative-volume cell, or unsupported
  mixed-sign vector scaling produce errors rather than fabricated geometry.

Direct coordinates are fractional and are not multiplied by the scale factor.
Following VASP's grammar, coordinate-mode lines beginning with C, c, K, or k
select Cartesian mode; every other token selects Direct mode. An unrecognized
Direct-mode token is retained and diagnosed but is not reinterpreted.
Cartesian coordinates are converted to fractional coordinates after all cell
and Cartesian scaling has been applied.

### 4.3 Species and populations

VASP 5 element-name rows map directly to normalized elements. VASP 4 files
without species names may use an explicitly supplied format option or recover
species from a scientifically explicit source field in OUTCAR/XML. When no
identity is available, the reader creates `UnknownSpecies` values keyed by type
index and emits a diagnostic; it never guesses elements from the comment.

Each listed atom is a fully occupied independent geometric position. VASP
species counts are validated against coordinate-row counts. Zero-count species
remain visible in the document but create no sites.

### 4.4 Selective dynamics, velocities, and extra data

Selective-dynamics flags are stored as a typed per-atom boolean triplet. The
meaning is retained exactly: `True` means the coordinate is allowed to move and
`False` means constrained. Invalid flags produce diagnostics.

Velocity blocks, when present and complete, are stored in a typed per-atom
vector property together with the explicit reported coordinate convention.
Cartesian velocities are retained in angstrom/fs and are not multiplied by the
POSCAR scale factor. Direct velocities remain in direct-lattice-vector per
timestep units; they are not silently converted without a reported timestep.
`VaspSnapshot` therefore carries `velocity_mode` as well as `velocity_unit`.
Predictor-corrector or other unimplemented trailing sections remain in the source
document and are reported as preserved unsupported data; they do not alter the
canonical positions.

`UnitCell` reconstructs a canonical Cartesian orientation from the six metric
parameters, while a VASP lattice may be globally rotated. The mapper computes
the orthogonal change of frame from the reported lattice to the canonical cell.
Forces and Cartesian velocities are rotated into that canonical frame. Direct
velocities remain components along the lattice vectors and are not rotated.
This preserves vector directions relative to the canonical structure.

Unknown or custom trailing lines never become plausible atoms.

## 5. POSCAR and CONTCAR

`PoscarDocument` losslessly retains source text and typed spans for:

- comment/title;
- scale specification;
- three lattice vectors;
- optional element-name row;
- population counts;
- optional Selective Dynamics marker;
- Direct/Cartesian/KPOINTS-style coordinate marker where valid;
- exactly the declared coordinate rows;
- optional complete velocity rows;
- all remaining source sections.

Mapping returns one `CrystalStructure`. CONTCAR is scientifically the same
grammar; its basename affects provenance, not the structure model.

Malformed headers or incomplete declared coordinate blocks produce no partial
canonical structure. Source documents and diagnostics remain available.

## 6. XDATCAR

`XdatcarDocument` retains one POSCAR-like header and indexes every complete
`Direct configuration=` or supported Cartesian configuration block by source
span. Initial reading parses the header and builds `FrameReference` entries
without constructing site objects or coordinate arrays for every frame.

Loading `sequence[i]` parses only the indexed frame, maps it with the shared
header, and caches the resulting `CrystalStructure` through the existing
`StructureSequence` contract. Repeated deterministic access returns the cached
snapshot. Thread-safe access is inherited from `StructureSequence`.

All complete frames except the last are `intermediate`; the final complete
frame is `final`. An incomplete trailing configuration is retained in the
document and produces a diagnostic but is excluded from the sequence.

## 7. OUTCAR

`OutcarDocument` performs a source scan that records shared metadata and frame
spans without building structures. It recognizes only explicit VASP fields:

- VASP version markers;
- `ions per type` and `VRHFIN`/`TITEL` species evidence;
- direct lattice-vector blocks associated with ionic states;
- `POSITION ... TOTAL-FORCE` coordinate blocks;
- clearly reported selective constraints or velocities when available.

Each complete POSITION block must have a compatible atom count and a valid
cell. Cartesian angstrom coordinates are converted into fractional positions.
Forces are retained as a typed per-atom vector property with explicit units;
they do not modify the structure.

Electronic iterations without a new complete ionic structure do not create
frames. A truncated final block is diagnostic-only. The last complete ionic
frame has role `final`.

## 8. vasprun.xml

`VasprunDocument` uses Python's standard XML facilities. The initial pass
validates XML syntax, reads global atom/species metadata, and indexes complete
`<calculation>` structure spans. It does not create a `CrystalStructure` for
every calculation.

Frame loading parses the selected indexed calculation, obtaining its basis,
fractional positions, and available per-atom forces. Initial/final structures
outside calculation elements are indexed only when they represent distinct
complete snapshots. Duplicate final snapshots are deduplicated by source
identity, not approximate geometry.

Namespaces and harmless unknown XML elements are tolerated and preserved in
the document. Missing atom metadata, inconsistent arrays, malformed numbers,
or incomplete XML produce location-aware diagnostics and never guessed values.

For very large XML files, the document may retain decoded source text because
the current `FormatHandler` contract is text-based, but canonical structures
and numerical frame arrays remain lazy. Moving source storage to a mapped-byte
contract is a separate I/O optimization and is not required by this reader.

### 8.1 Atomic-property symmetry boundary

VASP snapshots use explicit identity-only symmetry, so their site-indexed
properties map one-to-one into `AtomicView`. This milestone does not define a
generic transformation law for arbitrary `AtomicProperty` values under
non-identity symmetry. Polar vectors, axial vectors, fractional vectors, and
tensors require different transformations. CrIStMa must reject such expansion
until a property declares its transformation semantics; it must never copy a
force or velocity unchanged through inversion or rotation.

## 9. Results, identity, and provenance

Single-structure input returns a `ReadResult` containing a
`StructureCollection`. Trajectory input returns a `ReadResult` containing a
`StructureSequence`.

Stable IDs include source identity, source frame index, type index, and atom
index. Equivalent frames from different files are scientifically comparable
but do not share identity automatically.

Every frame and site maps back to its source span. The source document remains
available even when semantic mapping fails. Format-specific fields terminate
at this boundary; geometry and later calculations consume only canonical
objects and typed properties.

## 10. Diagnostics and failure policy

Errors prevent only structures or frames that cannot be represented honestly.
Stable diagnostic codes cover at least:

- incomplete or invalid header;
- invalid scale or singular cell;
- missing or inconsistent species/counts;
- incomplete coordinate, velocity, force, or frame block;
- unrecognized coordinate-mode token interpreted by VASP's Direct fallback;
- malformed selective-dynamics flags;
- malformed XML;
- numeric values outside finite supported domains;
- preserved unsupported trailing sections;
- VASP 4 unresolved species identity.

An invalid canonical object raises inside the mapper and is converted into a
source diagnostic by the handler. Reader APIs do not return partially valid
`CrystalStructure` objects.

## 11. Testing strategy

Implementation follows TDD in independently reviewable slices:

1. POSCAR/CONTCAR document, scale rules, mapper, and handler;
2. XDATCAR frame index and lazy sequence;
3. OUTCAR structural frame index and typed forces;
4. `vasprun.xml` index and lazy sequence;
5. registry integration, real fixtures, cross-format equivalence, and wheel
   smoke test.

Focused tests cover:

- VASP 4 and 5 headers;
- positive, negative-volume, and three-component scaling;
- Direct and Cartesian coordinates;
- Selective Dynamics and velocity blocks;
- malformed and truncated files;
- zero-count species and unresolved identity;
- complete and incomplete trajectories;
- lazy indexing, cache behavior, deterministic access, and final roles;
- content detection independent of basename;
- one provenance-recorded real fixture for every source family;
- equivalent CIF/POSCAR structures producing equivalent canonical geometry,
  neighbors, and coordination within declared tolerances;
- lazy descriptor import and wheel installation without optional
  crystallographic packages.

During implementation, only VASP and directly shared contract tests run. The
complete CrIStMa suite and wheel audit run once at the final gate.

## 12. Out of scope

This milestone does not add:

- canonical POSCAR, XDATCAR, OUTCAR, or XML writers;
- electronic energies, DOS, bands, wavefunctions, charge density, POTCAR, or
  pseudopotential interpretation beyond explicit species evidence;
- automatic bond, oxidation-state, space-group, or topology inference;
- ASE/pymatgen compatibility as a mandatory runtime path;
- CRAFT UI or project-format changes.

## 13. Primary format references

The implemented grammar and scientific semantics are checked against the
official VASP documentation:

- <https://vasp.at/wiki/POSCAR>
- <https://vasp.at/wiki/XDATCAR>
- <https://vasp.at/wiki/Vasprun.xml>
- <https://vasp.at/wiki/Output_files>

After VASP, the independent reader roadmap continues with Quantum ESPRESSO.
The universal CRAFT loader cutover remains gated by native SHELX, VASP,
XYZ/extXYZ, and PDB/PDBx/mmCIF support.
