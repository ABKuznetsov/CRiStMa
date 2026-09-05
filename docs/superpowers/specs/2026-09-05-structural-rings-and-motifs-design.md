# CrIStMa Structural Rings and Motifs Design

**Date:** 2026-09-05
**Status:** approved in discussion; written specification awaiting final review
**Implementation scope:** periodic-unit references, locally shortest structural rings,
symmetry grouping, and the CRAFT presentation contract

## 1. Purpose

CrIStMa already converts a selected crystal-chemical representation into exact
periodic structural blocks:

```text
StructuralUnitGraph
        ↓
StructuralRepresentation
        ↓
PeriodicConnectivity
        ↓
StructuralBlock
    0D / 1D / 2D / 3D
```

A block describes a maximal connected subsystem, but it does not describe the
finite motifs nested inside that subsystem. A rank-three B–O or T–O framework may
contain finite rings; a zeolite framework may later contain rings, cages and
building units. Those objects do not replace the parent block.

This milestone adds the first nested-motif analyzer:

```text
StructuralBlock
        ↓
RingFinder
        ↓
StructuralRing
        ↓
StructuralRingOrbit
```

CrIStMa owns the scientific identification and symmetry grouping. CRAFT owns
only tree placement, visibility, selection and rendering.

## 2. Scientific boundary

`RingFinder` operates on one already selected `StructuralRepresentation` and
the `StructuralBlock` objects produced from it:

```text
StructuralUnitGraph
→ StructuralRepresentation
→ PeriodicConnectivity
→ StructuralBlock
→ RingFinder
```

It must not operate on the unfiltered full graph. Consequently,
interstitial/context contacts cannot create rings in a primary structural
framework.

`RingFinder` must not:

- repeat neighbour or distance calculations;
- infer new chemical bonds;
- add connections absent from the selected representation;
- combine rings from different representations or parent blocks into one
  symmetry orbit;
- use compound names or framework codes in production decisions;
- apply CRAFT-specific visibility or rendering rules.

### 2.1 No structure-specific hardcoding

Production behavior must be derived only from the supplied canonical scientific
objects, selected representation, graph topology, periodic translations,
symmetry and explicit tool policy. It must not contain branches keyed by:

- chemical formula or material name;
- element labels such as `B`, `Si` or `O`;
- site labels such as `B1` or `O3`;
- input filename or source format;
- known framework, mineral or structure-type identifiers;
- expected ring size, composition or multiplicity of a test structure.

Names such as `B₃O₇ ring` are calculated from the canonical atom union of the
identified motif. Labels such as `B–O framework` are presentation of calculated
composition and block dimensionality; they are not inputs to ring recognition.

Lithium triborate, zeolites and other named structures are acceptance fixtures
only. Passing them must demonstrate the general algorithm, never activate a
special-case path.

## 3. Data model

### 3.1 `PeriodicUnitRef`

A unit identity alone is insufficient in a periodic graph. A path may visit
different periodic images of the same quotient-graph unit. Ring membership is
therefore expressed by:

```text
PeriodicUnitRef
├── unit_id
└── cell_translation: tuple[int, int, int]
```

The type is immutable and validates a three-integer translation. It is the
unit-level counterpart of `PeriodicAtomRef`.

### 3.2 `StructuralRing`

Each ring instance is an immutable result tied to its scientific context:

```text
StructuralRing
├── ring_id
├── parent_block_id
├── representation_id
├── unit_refs: tuple[PeriodicUnitRef, ...]
├── connection_ids: tuple[str, ...]
├── connector_atom_refs: tuple[PeriodicAtomRef, ...]
├── size
├── composition
├── translation_sum
└── provenance
```

Invariants:

- `size == len(unit_refs) == len(connection_ids)`;
- `size >= 3`;
- all referenced units and connections belong to the supplied representation;
- all units belong to `parent_block_id` in the selected representation;
- `translation_sum == (0, 0, 0)`;
- `connector_atom_refs` are copied from the referenced
  `StructuralConnection.shared_atom_refs`, never rediscovered geometrically;
- duplicate atom references are removed canonically.

`composition` is calculated from the canonical union of atoms represented by
the member units. Shared and connector atoms are counted once, not once per
unit. Thus a ring assembled from connected borate polyhedra may be reported as
`B₃O₇`, rather than as the arithmetic sum of three independently counted
coordination polyhedra.

### 3.3 `StructuralRingOrbit`

An orbit groups symmetry-equivalent instances, not merely rings with equal
composition or topology:

```text
StructuralRingOrbit
├── orbit_id
├── parent_block_id
├── representation_id
├── representative_ring_id
├── ring_ids
├── multiplicity
├── composition
└── size
```

All members must have the same parent block and representation. Topologically
equal rings that are not related by crystallographic symmetry remain different
orbits. A broader `RingType` is deliberately outside this milestone.

### 3.4 `RingAnalysisResult`

```text
RingAnalysisResult
├── rings
├── orbits
├── status: COMPLETE | INCOMPLETE
├── diagnostics
└── provenance
```

An empty complete result means that the analyzed block contains no qualifying
rings. Any computational safety limit must return `INCOMPLETE` with a precise
diagnostic; it may not silently omit candidates.

## 4. Ring definition

A `StructuralRing` in this milestone is a finite zero-translation cycle of
structural units in one selected representation and one parent block. It must
be obtained as a shortest return cycle for at least one of its connections and
must additionally satisfy the chordless criterion below.

This definition deliberately describes **locally shortest rings**, not every
chordless cycle or the complete cycle space.

For every eligible connection:

```text
anchor its first endpoint at translation (0, 0, 0)
        ↓
lift its second endpoint using the connection translation
        ↓
remove that exact anchored edge instance in both directions
        ↓
find all equal-length shortest return paths
        ↓
combine each path with the removed connection
        ↓
require at least three units
        ↓
require net translation (0, 0, 0)
        ↓
require a chordless cycle
        ↓
canonicalize and deduplicate
```

### 4.1 Translation-aware shortest paths

Shortest-path states are `PeriodicUnitRef` values, never bare `unit_id` values.
For an oriented connection:

```text
A --t--> B
```

the anchored removed edge is:

```text
A(0, 0, 0) → B(t)
```

and every return path is searched between those exact lifted states. Traversing
another connection adds its directed lattice translation to the current state;
reverse traversal adds the negated translation. Interior states of a candidate
cycle must be distinct.

Consequently, zero net translation follows naturally when the return path
reaches the exact target state. `translation_sum == (0, 0, 0)` remains an
explicit result invariant and defensive validation.

The lifted periodic state graph can be unbounded. Search policy therefore has
explicit hop/state limits. Reaching a limit before the target is resolved must
produce an `INCOMPLETE` result; it must not be reported as proof that no ring
exists.

### 4.2 Equal shortest paths

The algorithm must retain every equal-length shortest return path. Selecting
one arbitrary predecessor would lose distinct minimal rings in symmetric
graphs.

### 4.3 Periodic closure

Returning to the same quotient node is not sufficient. Directed traversal of a
connection accumulates its lattice translation; reverse traversal accumulates
the negated translation. A ring is finite only when the total translation is
zero.

For example, a path that returns to unit type `A` in image `(1, 0, 0)` belongs
to a periodic chain and is not a finite ring.

### 4.4 Chordless cycles

A candidate is chordless relative to the selected representation: no confirmed
connection joins two nonconsecutive unit images in the candidate and creates a
shorter closure. A square with a diagonal therefore contributes its minimal
triangles, not an additional composite square.

### 4.5 Eligible connection kinds

All confirmed shared-unit relations are eligible:

```text
SHARED_VERTEX
SHARED_EDGE
SHARED_FACE
```

Their connector atoms are already explicit:

```text
shared vertex → one connector atom
shared edge   → two connector atoms
shared face   → all atoms of the common face
```

`DIRECT_CONTACT` is excluded from the first ring analyzer unless a later
scientific requirement explicitly defines a ring over such contacts.

## 5. Canonical identity and symmetry

Ring identity is invariant under:

- cyclic permutation of the path;
- reverse traversal;
- a common shift of every periodic unit image.

Thus `A–B–C–A`, `B–C–A–B` and `A–C–B–A` describe one canonical ring.

Canonicalization occurs before accumulation of results. Symmetry grouping is a
second, separate operation over canonical instances. The operation applies the
parent crystal symmetry to unit and atom provenance, then matches canonical
ring identities. Equal size and composition alone are insufficient.

The action of a space-group operation `g` on a periodic unit reference must be
translation-aware:

```text
g : (unit_id, cell_translation)
    → (mapped_unit_id,
       M_g · cell_translation + normalization_shift_g(unit_id))
```

Here `M_g` is the integer action of the symmetry rotation in the lattice basis,
and `normalization_shift_g(unit_id)` is the integer shift introduced when the
mapped representative unit is returned to the reference cell. The same action
is applied consistently to referenced connections and connector atoms.

After applying `g`, the complete ring is canonicalized again with an arbitrary
common periodic-image shift removed. Two instances belong to one
`StructuralRingOrbit` only when this full mapping succeeds.

## 6. Performance and completeness

The first implementation prioritizes correct lifted-state traversal over
quotient-graph pruning. Ordinary biconnected decomposition of the finite
quotient graph is not an allowed mandatory optimization: it may collapse
parallel connections with different translations or hide distinct periodic
images, so a quotient bridge is not automatically safe to discard.

Translation-aware pruning may be added later only with a correctness argument
and tests for periodic multigraphs. Until then, shortest-path predecessor DAGs
over `PeriodicUnitRef` states preserve all equal shortest paths. Candidates are
checked for closure and chords, then canonicalized immediately. The
implementation must avoid enumerating the complete cycle space.

Configurable safety limits may protect callers from pathological graphs, but
crossing a limit produces:

```text
status = INCOMPLETE
diagnostic = crystal_chemistry.rings.search_limit_reached
```

The result records the applied policy and tool version in provenance.

## 7. CRAFT presentation contract

CRAFT receives rings and orbits from CrIStMa. It does not search the graph,
reconstruct connector atoms or regroup rings chemically.

The hierarchy is contextual and nested:

```text
Structural Blocks
└─ B–O framework
   └─ Rings
      └─ B₃O₇ ring · ×4
         ├─ instance 1
         ├─ instance 2
         ├─ instance 3
         └─ instance 4
```

The orbit row toggles every symmetry-equivalent instance. Expanded instance
rows permit selecting or displaying one instance.

Tree nesting expresses scientific ownership but does not impose visibility.
The block and its ring orbits are independently displayable:

```text
block on,  rings off → complete block
block on,  rings on  → block plus connector-atom ring markers
block off, rings on  → only ring member polyhedra plus connector atoms
block off, rings off → neither representation
```

Ring member polyhedra retain the parent block color. A ring is not highlighted
by recoloring every polyhedron. CRAFT renders connector atoms in their normal
element colors. Central atoms inside member polyhedra are not rendered as part
of the ring representation.

When both block and ring are visible, existing block surfaces are not
duplicated; only the connector-atom overlay is added. When only the ring is
visible, CRAFT draws the required member-polyhedron surfaces and connector
atoms.

The publication legend reflects selected representations:

- block only: `B–O framework`;
- ring only: `B₃O₇ ring`;
- both: both entries.

## 8. Future motif and topology layers

`StructuralRing` is the first concrete nested motif, not the final structural
description. Future independent tools may add:

```text
StructuralBlock
├─ RingFinder         → rings
├─ CageFinder         → finite cages
└─ BuildingUnitFinder → SBU / CBU
```

Ring analysis always runs on the lowest selected structural-unit
representation that still preserves the connections defining the internal
motif. If a finite group such as `B₅O₁₀` is later collapsed into one higher-level
unit, its internal rings are inherited from the lower representation rather
than rediscovered inside the collapsed node.

For zeolites, the first `RingFinder` reports locally shortest chordless rings.
It does not claim exact equivalence with every IZA essential-ring, natural-
tiling or composite-building-unit convention. Curated framework databases may
serve as acceptance oracles without introducing framework-name branches or
runtime dependencies.

Void-space analysis is scientifically separate from structural motifs:

```text
StructuralBlock + geometry
→ ChannelFinder
→ PoreTopologyResult
   ├─ cavities
   ├─ channels
   └─ channel dimensionality
```

CRAFT may display pore topology beneath the same parent block, but CrIStMa must
not represent a channel as an atomic or polyhedral motif.

## 9. Acceptance tests

The implementation is complete when the following cases pass:

1. A finite triangle produces one three-unit locally shortest ring.
2. A square with a diagonal produces the minimal triangles and no composite
   square.
3. Distinct equal-length shortest return paths are all retained.
4. A quotient-graph loop with nonzero net translation is rejected.
5. Cycles containing fewer than three units are rejected.
6. Shared-vertex, shared-edge and shared-face rings preserve the exact connector
   atoms from their connections.
7. Cyclic permutations, reverse traversal and common periodic-image shifts
   deduplicate to one instance.
8. Symmetry-equivalent instances form one orbit with correct multiplicity;
   merely topology-equivalent instances do not.
9. Rings from different blocks or representations never share an orbit.
10. Lithium triborate yields one `B₃O₇` ring orbit with multiplicity four
    inside the B–O framework.
11. CRAFT displays the orbit below its parent block, supports independent block
    and ring visibility, and renders connector oxygen atoms without recoloring
    the member polyhedra.
12. An exceeded safety limit returns an explicit incomplete result and
    diagnostic.
13. Renaming sites and the input file without changing the canonical structure
    leaves the detected rings and orbits unchanged.
14. No production branch refers to a fixture formula, material name, element
    combination, expected ring size or expected multiplicity.
15. A path search with equal quotient `unit_id` values but different periodic
    translations closes only when the exact target `PeriodicUnitRef` is reached.
16. Parallel quotient edges carrying different translations remain distinct
    during search and are not removed by ordinary quotient-graph bridge
    pruning.
17. A space-group operation maps unit translations, normalization shifts,
    connections and connector atoms consistently before orbit matching.
18. Rings inside a later collapsed finite group remain attached to the lower
    representation that preserved their defining connections.

## 10. Out of scope for this milestone

- general topology-equivalent `RingType` grouping;
- zeolite essential-ring or natural-tiling equivalence;
- cages;
- SBU and CBU recognition;
- pore, cavity and channel analysis;
- motif-based mechanical rigidity;
- application project state and persistence;
- CRAFT-side scientific ring inference.
