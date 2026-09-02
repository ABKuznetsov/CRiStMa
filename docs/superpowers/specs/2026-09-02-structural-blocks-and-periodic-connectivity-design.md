# CRiStMa Structural Blocks and Periodic Connectivity Design

**Date:** 2026-09-02  
**Status:** approved in discussion; written specification awaiting final review  
**Scope:** inorganic structural units, interaction-layer graphs, periodic connectivity, structural blocks, and nested motifs

## 1. Purpose

This milestone turns already validated crystal-chemical contacts and
coordination polyhedra into higher structural objects. It does not repeat
distance-based neighbour detection and it does not infer mechanical rigidity.

The scientific pipeline is:

```text
ResolvedContact / CoordinationShell
                ↓
          StructuralUnit
                ↓
       selected representation
                ↓
        structural-unit graph
                ↓
       exact periodic connectivity
                ↓
 StructuralBlock / chain / layer / framework
                ↓
        nested structural motifs
```

The central rule is:

> A structural block is not a group of nearby atoms. It is a maximal
> chemically connected subsystem in an explicitly selected interaction
> representation.

## 2. Layer boundary

The input graph is a scientific result of the preceding CRiStMa layers:

```text
Crystallography
    geometric contacts

Chemistry
    CompositionGrammar and interaction meaning

Crystal Chemistry
    ResolvedContact, CoordinationShell, CoordinationPolyhedron
```

This milestone consumes those objects. It must not:

- rerun a neighbour search with a private cutoff;
- create a contact that is absent from `ResolvedContact`;
- turn a merely `AMBIGUOUS` or `INCOMPLETE` shell into a polyhedron;
- use compound names, filenames, expected coordination numbers, or
  element-specific production branches;
- declare a structural unit mechanically rigid.

## 3. Terminology

### 3.1 `StructuralUnit`

A `StructuralUnit` is a finite chemically meaningful graph node with a complete
mapping to canonical atoms:

```text
StructuralUnit
├── unit_id
├── unit_kind                 # atom, bonded group, polyhedron, ...
├── atom_refs
├── source_contact_ids
├── source_polyhedron_ids
├── composition
└── provenance
```

The initial builders produce:

- one polyhedral unit from each successfully built coordination polyhedron;
- one bonded-group unit when confirmed non-shell contacts define a finite
  chemical group;
- an atomic unit only where a representation deliberately operates at atomic
  resolution.

Units may share atoms. A bridging ligand is not copied and retains the same
canonical atom identity in every unit.

### 3.2 `StructuralConnection`

A `StructuralConnection` is a canonical relation between two units in the
finite quotient graph:

```text
StructuralConnection
├── first_unit_id
├── second_unit_id
├── lattice_translation
├── connection_kind           # shared_vertex, shared_edge, shared_face,
│                             # direct_bridge
├── shared_atom_refs
├── source_contact_ids
└── provenance
```

For polyhedral units, sharing is determined from canonical ligand identities:

```text
one shared ligand      → shared_vertex
two shared ligands     → shared_edge
three or more ligands  → shared_face
```

The relative lattice translation remains part of the relation. Forward and
reverse descriptions of the same periodic connection are canonicalized so the
physical relation occurs once.

`direct_bridge` is used only when a confirmed `ResolvedContact` joins members
of two units without shared membership. A geometric near contact cannot create
this relation.

### 3.3 `StructuralRepresentation`

A representation is an immutable, inspectable selection of units and
connections:

```text
StructuralRepresentation
├── representation_id
├── units
├── connections
├── selection_policy
├── excluded_connection_ids
├── diagnostics
└── provenance
```

It is not a workflow mode and CRiStMa does not choose a globally preferred
representation. The caller may construct several representations from the same
structure, for example:

```text
primary structural graph
all validated contacts graph
polyhedral graph
```

The initial reusable selectors use existing facts only:

- grammar interaction type;
- grammar priority;
- resolved primary/secondary classification;
- coordination-shell status;
- explicit caller selection.

Interstitial, secondary, and context interactions are not deleted from the
scientific result. They are excluded from a particular representation and
remain available in provenance.

### 3.4 `StructuralBlock`

A `StructuralBlock` is one maximal connected component of a
`StructuralRepresentation`:

```text
StructuralBlock
├── block_id
├── representation_id
├── unit_ids
├── atom_refs
├── connection_ids
├── periodic_rank             # 0, 1, 2, or 3
├── periodic_generators
├── classification
├── nested_motif_ids
├── completeness
├── diagnostics
└── provenance
```

Its boundary is therefore reproducible: change the selected interaction layer
and a different, equally explicit representation may yield different blocks.

The word `block` here is crystal-chemical. It does not imply a rigid body.

### 3.5 `StructuralMotif`

A motif is a finite or repeating subgraph inside a block:

```text
StructuralMotif
├── motif_id
├── parent_block_id
├── unit/image members
├── motif_kind                # initially shortest-path ring
├── atom_refs
├── completeness
└── provenance
```

A ring or visually compact fragment is not automatically promoted to a block.
It remains nested unless a separately selected interaction representation makes
it a maximal connected component.

## 4. Building the structural-unit graph

The graph builder receives explicit units and validated scientific relations.
It never receives raw coordinates as authority for connectivity.

```text
units
+ shared canonical membership
+ ResolvedContact relations
        ↓
StructuralGraphBuilder
        ↓
finite quotient multigraph with integer translation labels
```

Coordinates may be used only to calculate descriptive geometry or to validate
the orientation of an already established shared face. They cannot add an
edge.

Parallel relations are retained. Two units can have several crystallographically
distinct connections, including connections to translated images of the same
unit.

## 5. Exact periodic connectivity

Each connected component is analysed as a gain graph whose edge labels are
integer lattice translations.

For each component:

1. choose a spanning tree and assign an integer image offset to every node;
2. for every non-tree edge and periodic self-edge, calculate its closure
   translation;
3. discard zero closure translations for periodic-rank calculation;
4. calculate the exact rank of the integer closure vectors in three dimensions;
5. retain a deterministic independent set as periodic generators.

The rank is exact and belongs to the component:

```text
rank 0 → finite block
rank 1 → chain or ribbon
rank 2 → layer
rank 3 → framework
```

A cross-cell edge that belongs only to the spanning tree is an image choice,
not proof of infinite connectivity. Periodicity requires a non-zero closure.

The rank is invariant under node insertion order, graph traversal order, and
the choice of crystallographic reference cell. Generator coordinates may
change under a cell transformation, but their rank and represented periodic
subspace must not.

The base classification is determined from rank. Optional morphology refinement
may distinguish chain from ribbon or finite block from ring/cluster using graph
topology, but it cannot change the calculated rank.

## 6. Rank before motifs

Periodic-component recognition always precedes motif search:

```text
representation
↓
connected components
↓
periodic rank and StructuralBlock
↓
motifs inside each block
```

Consequently:

- a ring inside a layer remains a motif of the layer;
- a ring inside a framework remains a motif of the framework;
- a standalone zero-translation ring in a rank-0 component may describe the
  whole finite block;
- failure or truncation of motif search cannot split a chain, layer, or
  framework into artificial islands.

The initial ring detector operates in the lifted periodic graph:

1. enumerate bounded simple cycles;
2. accept only cycles with zero accumulated lattice translation;
3. reject a cycle when a shorter path exists between two of its members;
4. canonicalize rotation, reversal, starting image, and symmetry-equivalent
   duplicates where the required symmetry mapping is available;
5. retain exact unit images, shared atoms, and source connections.

Search limits are explicit configuration. Budget exhaustion produces an
incomplete motif result, never a claim that no rings exist.

## 7. Multiple representations and examples

### `CaMoO4`

```text
primary structural representation
    Mo–O structural contacts
    → isolated MoO4 blocks, rank 0

all validated representation
    Mo–O + Ca–O
    → broader coordination context
```

Ca coordination does not silently merge the primary Mo–O representation.

### `FeS2`

```text
Fe–S shells → FeS6 polyhedral units
S–S contacts → bonded sulfur units or direct structural relations
```

The representation records both kinds explicitly instead of forcing every
contact into a coordination polyhedron.

### `LiB3O5`

```text
BO3 / BO4 units
→ periodic B–O block with rank determined before ring search
→ B3O7-like rings retained as nested motifs when derived from the graph
```

No ring is selected because its formula or compound name is known in advance.

### Layered structure

```text
primary in-layer interactions
→ one rank-2 StructuralBlock

secondary/interstitial interlayer contacts
→ excluded relations recorded in representation provenance
```

The layer is a crystal-chemical block even though it may contain many future
mechanical blocks.

## 8. Error handling and uncertainty

- An unresolved input shell is retained in diagnostics and is not converted
  into a polyhedral unit.
- A missing unit mapping makes the affected connection incomplete rather than
  guessing atom identity.
- A non-integral periodic translation is invalid input for the quotient graph.
- Inconsistent forward/reverse connections produce a diagnostic and no silent
  duplicate.
- Incomplete contact resolution propagates into representation and block
  completeness when it could change connectivity.
- Motif-search truncation affects motif completeness only; it does not erase
  an already established block or rank.
- Alternative representations coexist. CRiStMa does not hide one as the
  universally correct interpretation.

## 9. Pre-release tool contracts

The initial internal API consists of independent immutable tools and results:

```python
unit_result = StructuralUnitBuilder().build(
    resolution=resolution,
    polyhedra=polyhedra,
)

connection_result = PolyhedronConnectionFinder().find(
    units=unit_result.units,
)

representation = StructuralGraphBuilder(
    selection=PrimaryStructuralSelection(),
).build(
    units=unit_result.units,
    connections=connection_result.connections,
    contacts=resolution.contacts,
)

connectivity = PeriodicConnectivityAnalyzer().analyze(representation)

blocks = StructuralBlockFinder().find(
    representation=representation,
    connectivity=connectivity,
)

motifs = PeriodicRingFinder(limits=RingSearchLimits()).find(
    representation=representation,
    blocks=blocks.blocks,
)
```

Tools store configuration only. They do not retain the current structure or
last result. Each result stores diagnostics and sufficient provenance to
reproduce its selection and calculation.

Convenience functions may be added over these contracts, but no mandatory
session or monolithic hierarchy manager is introduced.

These names remain pre-release while CRAFT and Finder integration exercise the
library. They do not constitute a frozen public PyPI API yet.

## 10. Mechanical interpretation is out of scope

The following concepts are deliberately excluded from this milestone:

- `RigidBlock` or a claim that an entire structural block is rigid;
- a scalar block-rigidity score;
- automatic merging solely because polyhedra share an edge or face;
- hinge degrees of freedom;
- transformations or refinement parameterization;
- evidence from a temperature, pressure, or composition series.

Later mechanical analysis may use this hierarchy as input:

```text
StructuralBlock and nested polyhedra
+ polyhedron distortion/BVS evidence
+ user constraints
+ series-preserved internal geometry
        ↓
MechanicalBlockCandidate / HingeCandidate
```

Corner sharing is useful evidence for a possible pivot, and edge/face sharing
is useful evidence for stronger coupling, but neither is a mechanical truth by
itself.

## 11. Acceptance criteria

### Graph construction

- Polyhedra sharing one, two, or at least three canonical ligands produce
  vertex-, edge-, and face-sharing relations.
- Periodic reverse descriptions collapse to one relation.
- A confirmed direct contact can connect units; a geometric near contact
  cannot.
- Excluded context interactions remain visible in representation provenance.

### Periodic connectivity

- A lone cross-cell tree edge remains rank 0.
- A periodic self-edge with non-zero translation produces rank 1.
- Synthetic quotient graphs produce ranks 0, 1, 2, and 3 exactly.
- Rank and component membership are invariant under insertion-order changes.
- Chain, layer, and framework classification occurs before motif search.

### Motifs

- A zero-translation shortest-path ring crossing the reference-cell boundary
  is retained.
- A winding cycle with non-zero accumulated translation is not a finite ring.
- A cycle containing a graph shortcut is not a shortest-path ring.
- A ring inside a periodic component is nested under that component.
- Search-budget exhaustion is reported as incomplete.

### Scientific fixtures

- `CaMoO4` keeps primary MoO4 blocks separate from Ca coordination context.
- `FeS2` retains FeS6 polyhedra and S–S structural contacts.
- `LiB3O5` determines the periodic rank of its B–O component before reporting
  nested rings.
- A layered fixture yields a rank-2 primary block while retained interlayer
  context contacts do not alter that representation.
- Reordered sites, renamed files, and equivalent periodic edge orientation do
  not change the scientific result.

## 12. Relationship to earlier CRiStMa design

This document narrows and supersedes the grouping/discovery portions of
`2026-08-30-structural-hierarchy-design.md` where they conflict. In particular:

- the first implementation is a small set of independent tools, not a general
  hierarchy manager;
- structural blocks are components of explicit interaction representations;
- motifs are nested descriptions, not automatic blocks;
- mechanical characterization remains a later independent layer.

The canonical atom identity, shared membership, finite quotient representation,
and separation of grouping from mechanics remain unchanged.
