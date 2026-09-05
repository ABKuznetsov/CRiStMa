# CrIStMa structural hierarchy: design specification

Date: 2026-08-30

Status: proposed for user review

## 1. Purpose

Given an ordinary atomistic or crystallographic structure, CrIStMa shall build
physically interpretable structural entities without assuming in advance that
any entity is rigid.

The analysis direction is:

```text
CrystalStructure
-> ExpandedStructure
-> PeriodicNeighborGraph
-> CoordinationEnvironments
-> StructuralEntity candidates
-> StructuralHierarchy
-> EntityTransform parameterizations
-> later mechanical characterization
```

The reverse direction resolves every entity and transformation back to the
same source atoms used by geometry, diffraction, visualization, and
refinement.

## 2. Governing principles

1. Connectivity is evidence-backed, not a single distance cutoff.
2. A polyhedron is a structural entity, not automatically a final block.
3. A block is a semantic entity resolving to atoms, not a copied atom list.
4. Structural hierarchy is a directed acyclic hypergraph, not a strict tree.
5. Shared atoms are stored once and may belong to multiple entities.
6. Periodic chains, layers, and frameworks use finite quotient
   representations rather than enumerating infinite atoms.
7. Grouping and mechanical characterization are separate stages.
8. A single structure may suggest candidate entities but cannot establish
   experimental rigidity by itself.
9. Entity transforms parameterize atomic coordinates; they never create a
   second authoritative structure.
10. Every upward grouping and downward resolution retains semantic identity
    and provenance.

## 3. Atom identity

The hierarchy operates on finite reference-cell expanded atoms:

```text
ExpandedAtom
|- expanded_atom_id
|- source_site_id
|- representative_image
|- equivalent_images
`- structure_id
```

An expanded atom is derived from one independent site and one or more exact
symmetry images normalized into the reference cell. Human-readable atom labels
are not used as unique identity. A lattice-translated image used by a periodic
contact is represented separately by `PeriodicAtomRef`.

private copy of an atom or its coordinate.
All structural entities resolve to sets of `ExpandedAtom`. No entity owns a
private copy of an atom or its coordinate.
private copy of an atom or its coordinate.

## 4. PeriodicNeighborGraph

The graph contains one representative periodic node for each relevant
expanded atom and evidence-bearing contacts:

```text
PeriodicContact
|- contact_id
|- atom_a_id
|- atom_b_id
|- image_translation
|- distance
|- expected_distance or interval
|- contact_kind
|- strength_estimate
|- confidence
|- Evidence[]
`- provenance
```

`image_translation` identifies the periodic image of atom B relative to atom
A. Reversing a contact negates this translation.

Independent contact providers may contribute:

- covalent and ionic radius intervals;
- oxidation-aware bond-valence expectations;
- coordination-number plausibility;
- reported bonds or restraints;
- chemical-type rules;
- geometric separation and competition;
- optional user declarations.

The graph stores the evidence and the selected contact policy. It does not
present heuristic contacts as experimentally observed bonds.

## 5. Coordination environments and polyhedra

A `CoordinationEnvironment` is centered on one atom and refers to selected
neighbor contacts. It records coordination number, neighbor species,
geometry classification, distances, angles, distortion descriptors, and
classification confidence.

A `PolyhedralEntity` wraps a coordination environment when a central atom and
ligand shell form a meaningful unit such as `SiO4`, `AlO6`, `BO3`, or `MoO4`.
It preserves the central atom separately from ligand members and retains every
supporting graph edge.

Polyhedra may overlap. A bridging ligand can belong to two or more polyhedra.
Corner-, edge-, and face-sharing relations are derived from the cardinality
and identity of shared ligand sets.

## 6. StructuralEntity

All higher structural concepts implement one immutable semantic contract:

```text
StructuralEntity
|- entity_id
|- entity_type
|- direct_atom_members
|- child_entity_ids
|- periodic_generators
|- dimensionality
|- symmetry_descriptor
|- Evidence[]
|- confidence
`- provenance
```

Entity types include polyhedron, ring, cluster, chain, layer, framework,
molecular fragment, and candidate block. The type vocabulary is extensible,
but every type must define its membership and periodic semantics.

`direct_atom_members` are atoms not supplied by children. Resolving an entity
computes the set union of direct members and all child members. Repeated
membership never duplicates an atom.

## 7. StructuralHierarchy as a DAG and hypergraph

The hierarchy contains entities and explicit relations:

```text
StructuralRelation
|- relation_id
|- relation_type
|- participant_entity_ids
|- shared_atom_ids
|- periodic_translation
|- Evidence[]
`- provenance
```

Relation types include:

- `contains`;
- `overlaps`;
- `shares_vertex`;
- `shares_edge`;
- `shares_face`;
- `connected_to`;
- `periodic_repeat_of`;
- `candidate_hinge`.

Containment edges form a DAG and are checked for cycles. Non-containment
relations are hyperedges and may connect any number of entities. An atom may
belong to multiple sibling entities, but its canonical coordinate remains
unique.

Example:

```text
Layer
|- Chain A
|  |- SiO4-1
|  `- SiO4-2
`- Chain B
   |- SiO4-3
   `- SiO4-4

SiO4-1 --shares_vertex(O17)-- SiO4-2
```

## 8. Periodic motifs

Infinite chains, layers, and frameworks are represented by a finite quotient
entity plus integer translation generators:

```text
PeriodicEntityDescriptor
|- representative_entity_ids
|- translation_generators
|- dimensionality: 0 | 1 | 2 | 3
|- quotient_relations
`- embedding_provenance
```

Dimensionality follows the rank of independent periodic generators, not a
visual bounding-box heuristic. A chain crossing the unit-cell boundary remains
one chain.

## 9. Candidate discovery and evidence

No single all-purpose classifier decides structural grouping. Independent
detectors propose `CandidateStructuralEntity` objects:

```text
CandidateStructuralEntity
|- proposed_entity
|- Evidence[]
|- confidence
|- alternatives
|- conflicts
`- detector_id and version
```

Evidence categories include:

- chemical: coordination, bond valence, contact strength;
- topological: connected components, rings, shared vertices/edges/faces;
- geometric: compactness, distortion, internal-coordinate regularity;
- symmetry: repeated equivalent motifs;
- prior knowledge: known strong coordination units;
- series: preservation of internal geometry across temperature, composition,
  pressure, or time.

Conflicting candidates coexist until a selection policy or user decision
chooses a hierarchy. Confidence is not converted into a physical stiffness.

## 10. Entity-local coordinates

Every finite entity can derive an `EntityReferenceFrame` from a reference
structure:

```text
EntityReferenceFrame
|- origin
|- orthonormal_axes
|- reference_atom_coordinates
|- frame_construction_method
`- degeneracy diagnostics
```

Principal axes may be used only when non-degenerate. Symmetric entities need a
chemically or symmetry anchored frame to avoid arbitrary axis flips.

An `EntityTransform` contains translation and orientation. Orientation uses a
normalized quaternion or rotation-vector representation internally; Euler
angles are presentation values only. Optional internal degrees of freedom are
explicit named coordinates.

## 11. Shared atoms and kinematic consistency

Independent transforms cannot be applied blindly to overlapping entities. If
two entities share an atom, separate transforms may predict different
coordinates for that atom.

`EntityKinematicModel` therefore owns transforms and relations together:

```text
EntityKinematicModel
|- reference_structure_id
|- active_entity_ids
|- EntityTransform[]
|- Joint[]
|- CoordinateDependency[]
|- consistency_constraints
`- diagnostics
```

A shared atom has one resolved coordinate. Compatible motion is expressed by a
joint, hinge, exact dependency, or explicit internal deformation. An
overconstrained or inconsistent transform set is rejected before diffraction
calculation.

## 12. Downward coordinate resolution

The reversible path is:

```text
StructuralHierarchy
-> selected entities and transforms
-> kinematic dependency resolution
-> unique ExpandedAtom coordinates
-> IndependentSite parameter mapping where valid
-> structure factors and calculated profile
```

Resolution returns both coordinates and a dependency trace. It never mutates
the reference structure. A new immutable structure/refinement state is
produced.

The inverse analysis path groups the same atom references upward:

```text
IndependentSite
-> ExpandedAtom
-> PolyhedralEntity
-> StructuralEntity
-> hierarchy relations
```

## 13. Mechanical characterization

Grouping does not assert rigidity. Mechanical evidence is attached after the
hierarchy exists.

For the initial CrIStMa model, formal rigidity scores are defined only for
polyhedra, consistent with the current scientific decision. They may use
internal distance/angle changes, distortion measures, bond-valence behavior,
and matched structures in a series.

Higher entities such as chains, layers, and candidate blocks report:

- internal RMSD after optimal rigid alignment;
- translation and rotation relative to a reference;
- hinge/joint changes;
- internal-coordinate changes;
- correspondence confidence across a series.

These are deformation and stability observations, not an assumed scalar block
rigidity. Experimental invariance is a result to report, not a theory inserted
into refinement.

## 14. Series correspondence

`EntityCorrespondence` matches entities across related structures using source
site identity where available, species, topology, symmetry, and geometry. It
records ambiguous alternatives and confidence.

Series evidence may show that a polyhedron keeps its internal geometry while a
hinge angle or entity orientation changes. Such evidence may justify a later
refinement parameterization, but it does not rewrite the observed structures.

## 15. Relationship to refinement

Refinement consumes selected entity parameterizations through the shared
CrIStMa domain contract:

```text
structural entities
-> EntityKinematicModel
-> semantic Parameter / Dependency / Constraint objects
-> unique atomic coordinates
-> forward diffraction model
```

Entity-level and atom-level refinement are alternative coordinate systems over
the same semantic atoms. They are not separate phase models.

The intended refinement order remains:

```text
structural blocks/motifs -> polyhedra -> independent atoms
```

but a higher-level stage is enabled only when its selected kinematic model is
well-defined and physically admissible.

## 16. Diagnostics and invariants

The subsystem rejects or reports:

- contact policies with invalid physical ranges;
- hierarchy containment cycles;
- references to absent atoms or entities;
- inconsistent periodic translations around a graph cycle;
- degenerate local frames without a fallback anchor;
- transforms assigning multiple coordinates to a shared atom;
- entity correspondence that changes chemical identity silently;
- mechanical claims unsupported by their recorded evidence.

Partial candidate sets and graphs remain inspectable when a full hierarchy
cannot be selected.

## 17. Testing strategy

Tests progress from analytic fixtures to real structures:

1. periodic graph translations and reverse-edge invariants;
2. known isolated polyhedra;
3. corner-, edge-, and face-sharing pairs;
4. rings and finite clusters;
5. chains crossing cell boundaries;
6. layers and 3D frameworks with known dimensionality;
7. overlapping entities and shared-atom resolution;
8. containment-cycle and inconsistent-transform rejection;
9. upward grouping followed by downward atom resolution;
10. matched temperature/composition series with controlled rotations and
    internal distortions;
11. representative real fixtures with provenance and independent manual
    expectations.

Every derived entity must resolve to exactly the expected set of stable atom
references.

## 18. Implementation order

This subsystem follows, rather than expands, the current native-reader branch:

1. `PeriodicNeighborGraph` and evidence-bearing contacts;
2. coordination environments and polyhedral entities;
3. `StructuralEntity`, relations, and DAG/hypergraph invariants;
4. periodic motif dimensionality;
5. candidate detectors and evidence aggregation;
6. entity reference frames and transforms;
7. shared-atom kinematic model and downward resolution;
8. series correspondence and polyhedral mechanical characterization;
9. refinement adapters over semantic parameters and dependencies.

Reader work must provide the stable atom/site identity and provenance required
by this design, but reader implementations do not import topology algorithms.

## 19. Acceptance criteria

The structural hierarchy foundation is complete when:

- graph edges retain periodic image translations and evidence;
- known polyhedra and their sharing modes are recovered analytically;
- hierarchy containment is acyclic while overlapping memberships are allowed;
- periodic chains/layers/frameworks use finite representations and correct
  dimensionality;
- every entity resolves back to unique source atom references;
- incompatible transforms of shared atoms are rejected;
- upward grouping followed by downward resolution preserves atomic identity;
- rigidity is not assumed during grouping;
- formal initial rigidity characterization is limited to polyhedra;
- all scientific recovery, selection, and confidence are inspectable.

## 20. Deferred work

- learned classifiers trained on external structure databases;
- DFT-derived elastic or phonon evidence;
- automatic adoption of candidate blocks into refinement without review;
- scalar rigidity claims for higher structural entities;
- UI for editing hierarchy candidates and kinematic joints.
