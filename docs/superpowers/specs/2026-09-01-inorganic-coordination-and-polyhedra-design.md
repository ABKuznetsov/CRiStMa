# CrIStMa Inorganic Coordination and Polyhedra Design

**Date:** 2026-09-01  
**Status:** approved design, pre-implementation  
**Scope:** inorganic contact resolution, coordination shells, and 3D coordination polyhedra

## 1. Purpose

This milestone connects the existing canonical structure, composition grammar,
periodic neighbour graph, and reference radii without turning CrIStMa into an
application workflow.

It answers three reusable scientific questions:

1. Which geometrically available inorganic contacts are chemically meaningful?
2. Which confirmed contacts form a coordination shell around a crystallographic
   centre?
3. Which resolved three-dimensional shells define coordination polyhedra?

The first milestone covers:

```text
oxides and oxysalts
halides
chalcogenides
nitrides
borides
carbides
pnictides
intermetallics
```

Organic and organometallic structures are deliberately excluded. Their primary
representation requires a molecular covalent graph and receives a separate
milestone.

## 2. Layer boundary

Crystallography reports geometry. Crystal Chemistry interprets it.

```text
Crystallography
AtomicView + NeighborGraph
        ↓
GeometricContact

Chemistry
CompositionGrammar

Reference Data
CovalentRadii
        ↓

Crystal Chemistry
CoordinationShellResolver
        ↓
CrystalChemistryResolution
├── ResolvedContact
└── CoordinationShell
        ↓
PolyhedronBuilder
        ↓
CoordinationPolyhedron
```

`NeighborFinder` remains a geometric tool. It does not read chemical grammar,
radii, oxidation states, BVS parameters, or expected coordination numbers.

`CoordinationShellResolver` belongs to `cristma.crystal_chemistry` because it
combines geometric contacts with chemical grammar and reference data.

## 3. Minimal scientific model

The milestone has three main contact/shell entities.

### 3.1 `GeometricContact`

`GeometricContact` is the canonical finite representation of one contact from a
finite or periodic neighbour graph:

```text
GeometricContact
├── contact_id
├── first_atom_id
├── second_atom_id
├── cell_translation          # None for finite geometry
├── distance
├── vector_cartesian
├── first_source_site_id
├── second_source_site_id
└── geometric_provenance
```

It contains no radii, normalized distances, bond labels, grammar priority, or
chemical classification.

For a periodic edge, these two descriptions identify the same contact:

```text
A --t--> B
B ---t-> A
```

For a periodic contact, `first_atom_id` is placed in the reference cell and
`cell_translation` locates the second endpoint. Canonicalization selects the
lexicographically smaller of the forward and reverse keys:

```text
(A, B, t)
(B, A, -t)
```

Consequently a physical periodic contact occurs once in the result. This does
not change the directed neighbour graph used to enumerate environments.

### 3.2 `ResolvedContact`

`ResolvedContact` adds explicit chemical interpretation to one
`GeometricContact`:

```text
ResolvedContact
├── geometric_contact
├── interaction_type
├── interaction_layer         # STRUCTURAL | INTERSTITIAL | COORDINATION | ...
├── grammar_priority          # PRIMARY | ALLOWED
├── contact_classification    # PRIMARY | SECONDARY
├── component_interpretations
├── normalized_distance_min
├── normalized_distance_max
├── occupancy_evidence
├── evidence
└── provenance
```

Grammar priority and calculated contact classification are independent facts.
A `PRIMARY` grammar request cannot make geometrically poor contact primary.
The interaction layer is independent from both and is preserved unchanged
from `CandidateInteraction`, so later graph builders need no second chemistry
or geometry pass.

One contact may have several component-pair interpretations. For example:

```text
M1 = Ca0.7 Sr0.3
X1 = O0.8 F0.2

M1—X1
├── Ca—O
├── Ca—F
├── Sr—O
└── Sr—F
```

Each interpretation records the two species, their reported occupancies, the
covalent-radius sum, normalized distance, applicable grammar evidence, and an
occupancy evidence weight. The product of reported occupancies is an
average-structure evidence weight, not a claim about the joint local-disorder
probability.

Mixed components never create additional geometric edges.

One geometric contact may produce several `ResolvedContact` records when its
component pair participates in distinct interaction contexts. Resolution scope
is `source_site_id + interaction + centre view`; matching never stops at the
first `CandidateInteraction`.

### 3.3 `CoordinationShell`

`CoordinationShell` groups resolved contacts around one expanded centre while
retaining the orbit-level decision:

```text
CoordinationShell
├── source_site_id
├── center_atom_id
├── contacts: tuple[ResolvedContact, ...]
├── geometric_CN
├── mean_occupied_neighbors
├── status
├── alternatives
├── evidence
├── diagnostics
└── provenance
```

The status values are:

```text
RESOLVED
AMBIGUOUS
INCOMPLETE
NOT_APPLICABLE
```

`geometric_CN` is the number of geometric neighbour positions in the shell.
It is never multiplied by the number of chemical components.

`mean_occupied_neighbors` is the sum of total occupancies of the neighbour
positions. Four ligand positions with occupancy `0.75` therefore give:

```text
geometric_CN = 4
mean_occupied_neighbors = 3.0
```

This does not assert that every local configuration has coordination number 3.

## 4. Result container

One resolver call returns a small immutable container:

```text
CrystalChemistryResolution
├── contacts: tuple[ResolvedContact, ...]
├── coordination_shells: tuple[CoordinationShell, ...]
├── diagnostics
└── provenance
```

No separate `ResolvedContactNetwork` exists. A periodic or finite graph is a
derived view over `resolution.contacts`.

Contacts produced for these operations may form coordination shells:

```text
CENTRE_LIGAND_SHELL
INTERSTITIAL_COORDINATION
MIXED_ANION_COORDINATION
```

Contacts produced for these operations remain confirmed contacts unless a
separate consumer groups them:

```text
COVALENT_NETWORK
INTRA_SUBSYSTEM_BONDS
METALLIC_COORDINATION
```

Examples:

```text
Fe—S → CoordinationShell around Fe → FeS6 polyhedron
S—S  → ResolvedContact

Si—C → ResolvedContact → periodic graph
Ni—Al / Ni—Ni / Al—Al → ResolvedContact → metallic graph
```

## 5. Normalized distance

For every grammar-allowed component pair:

```text
rho_ij = d_ij / (r_i + r_j)
```

The first milestone uses covalent radii as the universal normalization scale.
Ionic radii are not primary because they require oxidation state, coordination
number, and sometimes spin state; using them to discover coordination number
would create a circular dependency.

If an oxidation state is known independently, ionic radii may later add
validation evidence without replacing the primary scale.

CrIStMa packages the Shannon table as versioned `ReferenceData`, retaining
both ionic and crystal radii. Lookup is exact by element, oxidation state,
coordination label, and spin state; missing values are never inferred. Because
coordination belongs to the lookup key, Shannon radii are applied only after a
shell or other independent analysis has supplied the required coordination.

An optional `ShannonDistanceValidator` compares an observed distance with an
explicit caller-selected lower-bound ratio:

```text
minimum_distance = minimum_ratio * (r_first + r_second)
```

This is secondary evidence only. A short-distance contradiction is reported
but does not delete the geometric contact or alter the primary shell boundary.
The threshold is part of the validator configuration and result provenance;
CrIStMa does not hide a universal hard cutoff in the resolver.

Contacts with mixed ligands are compared in normalized space rather than by raw
angstrom distances.

## 6. Search completeness

`NeighborFinder` accepts one cutoff in angstroms. Before geometric search, the
resolver derives a guaranteed cutoff from every grammar-allowed component pair:

```text
search_cutoff =
    max(r_i + r_j for every allowed component pair)
    * candidate_rho_max
```

After enumeration, every interpretation is filtered by its own normalized
distance:

```text
rho_ij <= candidate_rho_max
```

A shell boundary is resolvable only when an outer distance group is observed.
If the candidate range ends at the selected shell, the result is:

```text
INCOMPLETE
diagnostic = crystal_chemistry.shell.search_boundary_not_observed
```

Increasing `candidate_rho_max` and rerunning is explicit application policy;
the resolver does not silently extend its search.

Missing radii, unknown species, invalid occupancy, unsupported disorder,
insufficient geometric candidates, or an incomplete orbit also produce
`INCOMPLETE`, not a guessed shell.

A missing radius remains a diagnostic in its interaction scope. A missing
`PRIMARY` component-pair radius makes that scope `INCOMPLETE`; the resolver may
not silently resolve from only the remaining components.

## 7. Orbit-level resolution

Shell boundaries are determined collectively for every crystallographic orbit:

```text
source_site_id
↓
all symmetry-equivalent ExpandedAtom centres
↓
common geometric contacts
↓
component-pair interpretations
↓
one normalized distance picture
↓
one boundary decision
↓
projection to every center_atom_id
```

Contact signatures never compare raw expanded-atom IDs or translation vectors.
The canonical orbit-aware signature uses:

```text
neighbour source_site_id
component-pair interpretation
interaction type
normalized-distance group
multiplicity of equivalent contacts
```

If these signatures disagree between symmetry-equivalent centres, every shell
projection for that orbit is `INCOMPLETE` with:

```text
crystal_chemistry.shell.symmetry_inconsistent
```

This is not `AMBIGUOUS`. It signals incomplete geometric search, an unsuitable
tolerance, incorrect symmetry expansion, special-position failure, or invalid
source data.

No shell is resolved from a subset of an orbit and copied to the missing
members.

## 8. Resolution algorithm

The resolver does not calculate one opaque bond-probability score. It applies a
traceable sequence:

```text
hard filters
↓
distance grouping
↓
candidate boundaries
↓
lexicographic comparison
↓
secondary validation
↓
RESOLVED / AMBIGUOUS / INCOMPLETE / NOT_APPLICABLE
```

### 8.1 Hard filters

A component-pair interpretation enters comparison only when:

- `CompositionGrammar` permits the pair and interaction type;
- both covalent radii are available and valid;
- its normalized distance is within the candidate range;
- all symmetry-equivalent centres have complete geometric input;
- source occupancy and species data are usable.

No expected coordination number is used. Rules such as “Si must be CN4” or
“Fe must be CN6” are forbidden inputs to the resolver.

### 8.2 Distance groups

Contacts are ordered by normalized distance. Symmetry-equivalent contacts and
numerically close normalized distances are placed in common groups using the
dimensionless `distance_group_tolerance`.

Grouping is an internal calculation recorded as evidence. `DistanceGroup` is
not a top-level scientific result.

### 8.3 Candidate boundaries

Every boundary followed by an observed outer group is considered. For boundary
group `k`, the relative gap is:

```text
relative_gap_k = (rho_first_outside - rho_last_inside) / rho_last_inside
```

`minimum_shell_gap` is a dimensionless minimum relative gap in normalized
distance space. It is not an angstrom cutoff.

Internal spread is dimensionless and defined as:

```text
internal_spread =
    (rho_max_inside - rho_min_inside) / rho_median_inside
```

It is zero for a single-valued shell. All input values are retained in method
provenance.

A gap is physically significant when it is at least `minimum_shell_gap`.
Boundaries below that value may be retained as alternatives, but cannot alone
produce `RESOLVED`. `strong_contacts_outside` is true when the immediately
following group contains a grammar-permitted interpretation within this
minimum relative-gap interval.

### 8.4 Lexicographic comparison

Candidates are compared without a weighted sum:

1. all hard filters must pass;
2. prefer the larger physically significant relative gap;
3. when gaps are equivalent, prefer smaller internal spread;
4. then prefer no strong contact immediately outside the boundary;
5. only then use grammar priority.

`ambiguity_tolerance` applies separately to each compared dimensionless
criterion. For example, two gaps whose difference lies within this tolerance
are considered equivalent and comparison proceeds to internal spread.

If the criteria do not establish dominance, the shell is `AMBIGUOUS` and all
non-dominated alternatives are retained in its evidence.

### 8.5 Secondary validation

Bond-valence and coordination-geometry analyses return independent evidence:

```text
BVS evidence
SUPPORTIVE | NEUTRAL | CONTRADICTORY | NOT_AVAILABLE

Geometry evidence
SUPPORTIVE | NEUTRAL | CONTRADICTORY | NOT_APPLICABLE
```

They do not create a contact absent from geometric candidates and do not enter
the primary comparison through hidden weights. Contradictory secondary evidence
may prevent an otherwise close decision from being declared resolved, but it
does not silently replace the recorded lexicographic outcome.

Shannon-radius distance checks follow the same separation: they report
`SUPPORTIVE` or `CONTRADICTORY` evidence for an explicitly supplied threshold,
but never select or exclude a contact themselves.

## 9. Resolution policy and provenance

```text
ShellResolutionPolicy
├── candidate_rho_max
├── distance_group_tolerance
├── minimum_shell_gap
└── ambiguity_tolerance
```

Every field is finite, positive, and dimensionless. The resolver receives a
policy explicitly; there is no hidden mutable global policy. A named preset may
be added only after calibration against the acceptance corpus and must carry a
method version.

Every `CrystalChemistryResolution` records:

```text
policy values
derived search cutoff in angstroms
maximum observed normalized distance
whether every accepted boundary had an outer group
CompositionGrammar method/version
ReferenceData version/checksum
resolver method/version
input structure identity
```

This provenance is sufficient to reproduce why a structure was resolved,
ambiguous, incomplete, or not applicable.

## 10. Polyhedron construction

`PolyhedronBuilder` consumes only a `RESOLVED CoordinationShell`.

```text
RESOLVED CoordinationShell
↓
one local Cartesian vertex per geometric neighbour position
↓
affine-rank check
↓
rank 3 → CoordinationPolyhedron
rank < 3 → NOT_APPLICABLE
```

```text
CoordinationPolyhedron
├── polyhedron_id
├── source_site_id
├── center_atom_id
├── shell_provenance
├── vertex_contacts
├── local_vertices
├── faces                         # maximal coplanar polygons
├── volume
├── geometric_centroid
├── center_offset
└── diagnostics
```

`PolyhedronBuildResult` holds either one polyhedron or a non-success status and
diagnostics. Linear, planar, insufficient, ambiguous, and incomplete shells are
not converted into artificial volume-bearing polyhedra.

Occupancy does not remove a mean-structure vertex. Mixed chemical components
remain accessible through the vertex's `ResolvedContact`.

Faces are the convex hull of ligand positions around the centre. Equivalent
centres must have the same face topology under an orbit-aware signature.
Public faces are ordered maximal coplanar polygons (`tuple[int, ...]`), not an
internal triangular tessellation. A cube has six quadrilateral faces.
Distortion, BVS, and future rigidity analysis describe a constructed
polyhedron; they do not determine whether its shell contacts exist.

Failure status is preserved: linear and planar resolved shells return
`NOT_APPLICABLE`; ambiguous and incomplete shells return `AMBIGUOUS` and
`INCOMPLETE` respectively.

## 11. Package organization

Only reusable scientific operations enter CrIStMa:

```text
src/cristma/
├── crystallography/
│   └── local_geometry.py
│       GeometricContact and graph-to-contact conversion
│
└── crystal_chemistry/
    ├── __init__.py
    ├── contacts.py
    ├── policy.py
    ├── resolver.py
    └── polyhedra.py
```

Internal distance-group and candidate-boundary helpers remain private to the
resolver until another independent scientific consumer requires them.

CRAFT, Finder, and other applications own result caching, selection, display,
and workflow. CrIStMa owns the scientific objects and calculations only.

## 12. Diagnostics and failure semantics

Invalid canonical objects and invalid tool configuration raise exceptions.
Scientific uncertainty or incomplete source/reference information returns a
status and machine-readable diagnostics.

Initial stable diagnostic codes include:

```text
crystal_chemistry.contact.radius_missing
crystal_chemistry.contact.species_unknown
crystal_chemistry.contact.occupancy_invalid
crystal_chemistry.shell.search_boundary_not_observed
crystal_chemistry.shell.candidates_insufficient
crystal_chemistry.shell.mixed_occupancy_disagreement
crystal_chemistry.shell.symmetry_inconsistent
crystal_chemistry.shell.boundary_ambiguous
crystal_chemistry.polyhedron.shell_not_resolved
crystal_chemistry.polyhedron.not_three_dimensional
crystal_chemistry.polyhedron.symmetry_inconsistent
```

Mixed occupancy alone is not an error and does not lower status. It yields
`AMBIGUOUS` only when valid component interpretations imply different shell
boundaries.

## 13. Acceptance corpus

The implementation is accepted against committed, attributed structure
fixtures covering the shared algorithm rather than material-specific code:

```text
NaF       ionic coordination shell
SiC       covalent resolved contacts and periodic graph view
Si3N4     resolved Si-centred SiN4 shells
FeS2      FeS6 shell plus independent S–S contacts
CaN2      Ca–N shell plus N–N contacts when present in the structure
Na3P      resolved Na–P ionic shell
Bi2Te3    Bi–Te shell with retained secondary candidates
CaMoO4    MoO4 and Ca coordination shells
LiB3O5    BO3 and BO4 environments
anorthite resolved Al/SiO4 polyhedra, or explicit AMBIGUOUS/INCOMPLETE shells
          with no forced polyhedron when the generic resolver cannot prove them
```

Acceptance tests assert scientifically established outcomes, but those values
are never supplied to the resolver as expected coordination numbers.

Additional analytic tests cover:

- canonical periodic contact identity and reversal;
- mixed occupancy on one and both ends;
- vacancy-weighted mean neighbour counts;
- equivalent-orbit projection and deliberate symmetry inconsistency;
- missing radii and unknown species;
- absent outer distance groups;
- competing shell boundaries and retained alternatives;
- planar and linear non-polyhedral shells;
- polyhedron face topology under symmetry-equivalent centres;
- absence of family- or compound-specific Python branches.

## 14. Out of scope

This milestone does not add:

- molecular graph construction for organic structures;
- organometallic coordination;
- intermolecular contacts;
- expected-CN lookup tables;
- Voronoi as a source of truth;
- ionic radii as the primary shell scale;
- BVS-created contacts;
- structural blocks, chains, layers, rings, or topology;
- polyhedral distortion and rigidity models;
- application caches, project state, UI, or CRAFT rendering code.

It also adds no runtime dependency beyond the package's existing NumPy
requirement. Convex-hull construction for the small coordination environments
is implemented within CrIStMa rather than adding SciPy solely for polyhedra.

The next milestone may use `CoordinationPolyhedron` and resolved contact graphs
to build structural units and hierarchy without changing the contracts above.
