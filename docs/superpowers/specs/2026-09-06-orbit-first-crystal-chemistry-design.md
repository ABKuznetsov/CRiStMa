# Orbit-first crystal chemistry design

## Purpose

CrIStMa shall calculate direct-space crystal chemistry from independent sites
and validated crystallographic symmetry. Symmetry orbits are the primary
scientific objects. Expanded contacts are created only at the outward
materialization boundary for scene construction, export, or other consumers.

This replaces the expanded-first route that expands the reference cell, finds
thousands of individual contacts, interprets each contact, and subsequently
tries to recover symmetry equivalence. The replacement follows the same
general principle as crystallographic asymmetric-unit pair tables: generate
independent periodic pair relations, quotient them by endpoint stabilizers,
and resolve their chemistry once per orbit.

The model is frozen by this design. New public entities are added during
implementation only when a concrete mathematical invariant cannot be
represented by the types defined here.

## Scope

This milestone includes:

- validated direct-space `SymmetryContext` construction;
- exact affine-periodic relation algebra;
- asymmetric-unit site mappings and stabilizers;
- complete cutoff-bounded asymmetric-unit pair search;
- geometric contact orbits;
- orbit-first chemical interpretations and oriented incidences;
- orbit-first coordination-shell resolution;
- local polyhedron realization from shell orbits;
- materialization of contacts and shells for consumers;
- migration of structural units, blocks, rings, and later motif analysis to
  orbit-first inputs;
- removal of the expanded-first crystal-chemistry route after acceptance.

It does not include presentation, colours, visibility, scene meshes,
experimental comparison, refinement, or Finder/CRAFT policy.

Diffraction remains outside this input contract. Diffraction v1 continues to
require a catalogued `SpaceGroupSetting` because reciprocal-basis identity,
reflection IDs, and extinction provenance require an unambiguous setting.

## Architecture

The scientific data flow is:

```text
IndependentSite + SymmetryContext
                ↓
AsymmetricUnitMapping
                ↓
SymmetryPairFinder
                ↓
SymmetryPairTable
                ↓
SymmetryContactOrbit
                ↓
ResolvedContactOrbit
                ↓
ContactIncidenceOrbit
                ↓
CoordinationShellOrbit
                ↓
PolyhedronOrbit
                ↓
Structural units / blocks / rings / FBB / motifs

                ↓ outward only
Materializers
                ↓
ResolvedContact / materialized shell / consumers
```

No scientific stage after contact resolution may consume materialized
`ResolvedContact` objects.

## Symmetry context

`SymmetryContext` is the only public symmetry input accepted by orbit-first
direct-space analyses:

```text
SymmetryContext
├─ operations: tuple[AffineOperation, ...]
├─ operation_keys
├─ basis_convention
├─ cell_fingerprint
├─ symmetry_action_fingerprint
├─ setting_id: str | None
├─ source_kind
├─ status
├─ diagnostics
└─ provenance
```

The existing exact `AffineOperation`, whose matrix and translation entries are
`Fraction`, is used. No duplicate exact-operation class is introduced.

Contexts are created through explicit constructors:

```text
SymmetryContext.from_setting(setting, cell)
SymmetryContext.from_definition(definition, cell)
SymmetryContext.from_operations(operations, cell, provenance=...)
```

Contact analysis never accepts `SpaceGroupDefinition` or a raw operation list
directly. Each constructor performs the same validation and canonicalization
before returning a usable context.

The validation gate requires:

- one identity operation;
- unique normalized operations;
- closure modulo lattice translations;
- an inverse for every operation;
- integer rotational entries;
- crystallographically valid rotational determinants;
- exact rational translations normalized modulo the lattice;
- compatibility of each rotational part with the numerical cell metric.

Group algebra is exact. Metric compatibility is numerical and uses one
documented metric tolerance. An invalid operation group does not produce a
working context and raises a symmetry invariant error. Identity fallback is
possible only through an explicitly constructed fallback context carrying its
own diagnostic; analysis never silently replaces invalid symmetry.

Operations are normalized and sorted by their exact `(R, t mod 1)` descriptor.
`operation_key` is derived from that descriptor and never from CIF row order,
source text, or a transient operation number.

`source_kind` distinguishes `CATALOG_SETTING`, `VALID_EXPLICIT_OPERATIONS`, and
an explicitly requested identity fallback. A complete valid explicit group is
not scientifically incomplete merely because `setting_id` is `None`. It may
emit `symmetry.setting_unresolved` as an identification diagnostic.

The full context fingerprint includes the cell fingerprint and all provenance
required to validate downstream inputs. The separate
`symmetry_action_fingerprint` contains the fractional-basis convention and
canonical operation set but excludes metric cell parameters. Scientific orbit
IDs use the action fingerprint, so a metric-only cell change preserves IDs
provided independent-site IDs and the symmetry action remain unchanged.

## Periodic symmetry relations

`PeriodicSymmetryRelation` is an element of the complete affine-periodic
action, not an unrelated operation/vector pair:

```text
PeriodicSymmetryRelation
├─ operation_key
└─ lattice_translation: tuple[int, int, int]
```

If normalized operation `g` has `(R_g, t_g)`, relation `(g, n)` acts as:

```text
x ↦ R_g x + t_g + n
```

The type supplies exact `compose`, `inverse`, and `normalize` operations.
Composition rotates the second lattice translation by the first rotational
part and retains the exact integer carry created when the composed fractional
translation is normalized. Inversion likewise retains its normalization
carry. All group-law invariants are checked with exact arithmetic.

## Asymmetric-unit mapping

`AsymmetricUnitMapping` contains one mapping per independent site:

```text
SiteOrbitMapping
├─ independent_site_id
├─ stabilizer_relations
└─ reference_cell_images
   └─ SiteImage
      ├─ image_id
      ├─ representative_relation
      ├─ equivalent_relations
      ├─ fractional_position
      └─ normalization_translation
```

A stabilizer contains complete `PeriodicSymmetryRelation` objects because an
operation may fix a special position only modulo an integer lattice
translation. An operation key alone is insufficient.

Operations producing the same special-position image are represented by one
`SiteImage`. Its equivalent relations retain the complete coset evidence.
Image identity is deterministic and independent of operation ordering.

The mapping has its own fingerprint derived from independent-site identities,
coordinates, occupations relevant to site identity, the symmetry action, and
the mapping convention. Numerical fractional-coordinate tolerance is explicit
and recorded in provenance.

## Complete asymmetric-unit pair search

`SymmetryPairFinder` generates the complete finite set of symmetry and lattice
images needed for a requested cutoff. It must not use a fixed `-1..1`
supercell. Translation bounds are derived from the numerical unit-cell metric
and cutoff and are valid for skewed and triclinic cells. Every accepted pair
passes a final numerical metric-distance check.

The search pipeline is:

```text
independent sites
→ complete symmetry/lattice cutoff buffer
→ Cartesian spatial bins
→ candidate pairs
→ numerical distance check with documented tolerance
→ periodic-relation canonicalization
→ endpoint-stabilizer quotient
→ SymmetryContactOrbit
```

Symmetry relations and lattice arithmetic are exact. Distances and cutoff
membership are numerical. Distance tolerance applies only to metric boundary
decisions and never to group identity.

For endpoints `A` and `B`, raw pair relations are quotiented by the endpoint
stabilizers as a double coset in the full affine-periodic group:

```text
H_A \ relation / H_B
```

For undirected contacts, endpoint exchange adds:

```text
(A, B, relation) ~ (B, A, inverse(relation))
```

The minimum exact descriptor after both equivalences is the canonical pair
relation.

## Pair table and geometric orbit

`SymmetryPairTable` is a reproducible geometric result:

```text
SymmetryPairTable
├─ contact_orbits: tuple[SymmetryContactOrbit, ...]
├─ symmetry_context_fingerprint
├─ asymmetric_unit_mapping_fingerprint
├─ cutoff
├─ distance_tolerance
├─ status
├─ diagnostics
└─ provenance
```

`SymmetryContactOrbit` is a pure geometric fact:

```text
SymmetryContactOrbit
├─ geometry_orbit_id
├─ first_independent_site_id
├─ second_independent_site_id
├─ canonical_relation
├─ equivalent_relations
├─ endpoint_stabilizers
├─ representative_distance
├─ representative_vector_cartesian
├─ multiplicity_in_reference_cell
├─ status
├─ diagnostics
└─ provenance
```

`geometry_orbit_id` is derived from independent-site IDs, canonical periodic
relation, fractional-basis convention, and symmetry-action identity. It does
not contain cutoff, distance, metric cell parameters, or enumeration order.
Thus a metric-only cell change preserves an orbit ID when independent-site IDs
and symmetry action are preserved, while the calculated distance changes.

`multiplicity_in_reference_cell` is the number of unique undirected instances
owned by one reference cell under the single ownership rule defined below. It
is not a coordination number.

## Chemical contact resolution

One geometric pair orbit may support multiple chemical interpretations. They
are preserved instead of being destructively coalesced:

```text
ResolvedContactOrbit
├─ resolved_contact_orbit_id
├─ geometry_orbit_id
├─ interpretations
│  └─ ContactInterpretation
│     ├─ interpretation_id
│     ├─ interaction_type
│     ├─ interaction_layer
│     ├─ grammar_priority
│     ├─ orientation_mode
│     ├─ endpoint_roles
│     ├─ component_pair_interpretations
│     ├─ normalized_distance_range
│     ├─ status
│     └─ evidence
├─ status
├─ diagnostics
└─ provenance
```

Pair-level chemistry states that a geometric pair admits an interpretation.
`orientation_mode` and `endpoint_roles` state whether the interpretation is
undirected or has scientifically meaningful centre/ligand roles. Mixed-site
component alternatives remain explicit.

## Oriented contact incidences

`ContactIncidenceOrbit` connects undirected pair topology to local
coordination around one independent centre:

```text
ContactIncidenceOrbit
├─ incidence_orbit_id
├─ resolved_contact_orbit_id
├─ interpretation_id
├─ center_independent_site_id
├─ ligand_independent_site_id
├─ oriented_periodic_relation
├─ incidence_multiplicity_per_center
├─ effective_neighbor_occupancy
├─ status
└─ evidence
```

Its identity exists before shell resolution. It therefore contains no
`PRIMARY` or `SECONDARY` classification and remains stable when shell policy
changes.

Incidence multiplicity is calculated through the centre stabilizer. A global
pair orbit of multiplicity one may contribute two local incidences to each
centre, as in a simple periodic `A—A—A` chain. More complicated special
positions may yield multiple incidence orbits of multiplicity one instead.

`effective_neighbor_occupancy` includes only components participating in the
selected interpretation. For a mixed `O0.6 F0.4` ligand and an `M—O`
interpretation, it is `0.6`, not `1.0`. Centre occupancy is not multiplied into
coordination counts because a shell is defined per occupied centre.

## Coordination-shell orbits

Shell selection belongs to alternatives, not incidence identity:

```text
CoordinationShellAlternative
├─ alternative_id
├─ primary_incidences
├─ secondary_incidences
├─ geometric_CN
├─ mean_occupied_neighbors
├─ boundary_evidence
└─ status

CoordinationShellOrbit
├─ shell_orbit_id
├─ center_independent_site_id
├─ selected_alternative: alternative_id | None
├─ alternatives
├─ status
├─ diagnostics
└─ provenance
```

For one alternative:

```text
geometric_CN = Σ incidence_multiplicity_per_center over primary incidences

mean_occupied_neighbors =
    Σ incidence_multiplicity_per_center × effective_neighbor_occupancy
      over primary incidences
```

Distance groups contain incidence orbits rather than duplicated contact
instances. Group populations, geometric CN, occupied-neighbour counts, and
internal statistics retain incidence-multiplicity weights. A resolved shell
selects one alternative. An ambiguous shell preserves all genuine alternatives
and has no selected alternative.

## Polyhedra and structural hierarchy

A representative coordination polyhedron is constructed directly from a
`CoordinationShellOrbit`. The builder may realize the finite local neighbour
positions around one representative centre through oriented relations and
incidence multiplicities. This is local geometric realization, not global
expanded-contact analysis.

Polyhedron orbits, structural units, blocks, rings, FBBs, and motif topology
consume orbit-first results. They never receive the compatibility
`result.contacts` view. The silicate and borate hierarchy is subsequently
expressed over these orbit objects:

```text
full structural unit
→ framework / layer / chain / finite component
→ repeat unit or FBB
→ ring
→ constituent polyhedra
```

A three-dimensional framework may contain scientifically defined child layers,
chains, rings, or FBBs. Arbitrary subgraphs are not reported as motifs. Pure
coordination packing such as NaCl must not create structural rings, while an
isolated borate ring such as the tested K7 structure remains detectable.

## Materialization boundary

Materialization belongs to CrIStMa but lies outside the scientific-analysis
pipeline:

```text
ContactMaterializer.materialize(
    result,
    region,
    contact_orbit_ids=None,
    interpretation_ids=None,
    shell_alternative_ids=None,
) -> tuple[ResolvedContact, ...]
```

Regions are scientific periodic ranges:

```text
ReferenceCell
CellRange(a_min, a_max, b_min, b_max, c_min, c_max)
```

`ReferenceCell` means all contacts whose canonical owner belongs to the
reference cell. It does not require both endpoints to lie in `[0,1)`. A bond
from an atom in the reference cell to an atom in cell `(+1,0,0)` is retained
once when its owner is in the reference cell. `CellRange` applies the same rule
to a larger set of owner cells.

The canonical ownership rule is shared by pair multiplicity and
materialization. An undirected instance is normalized under global lattice
translation and endpoint exchange; the exact minimum endpoint/relation
descriptor owns the instance. Translating that canonical owner through the
requested owner cells generates the region. Consequently:

```text
len(materialize(ReferenceCell, one_orbit))
    == one_orbit.multiplicity_in_reference_cell
```

`ResolvedContact` is a materialized record:

```text
ResolvedContact
├─ contact_id
├─ resolved_contact_orbit_id
├─ first_atom_ref
├─ second_atom_ref
├─ distance
├─ vector_cartesian
├─ interpretations
├─ shell_memberships
└─ provenance
```

Shell membership, rather than one global contact classification, preserves
ambiguous shell semantics:

```text
ShellMembership
├─ shell_orbit_id
├─ alternative_id
└─ role: PRIMARY | SECONDARY
```

An instance ID is derived from its resolved orbit ID and canonical periodic
endpoint references. Output order is deterministic.

`ContactAnalysisResult.contacts` remains temporarily available as a read-only
compatibility property equivalent to materialization in `ReferenceCell`. It is
not stored as independent state and cannot be an input to another CrIStMa
scientific calculator. Explicit materializers are used for larger regions and
filtered projections.

## Result and status propagation

The primary result is:

```text
ContactAnalysisResult
├─ pair_table
├─ contact_orbits
├─ contact_incidence_orbits
├─ coordination_shell_orbits
├─ status
├─ diagnostics
├─ configuration
└─ provenance
```

Symmetry source status and scientific result status are distinct.
`CATALOG_SETTING` and `VALID_EXPLICIT_OPERATIONS` describe provenance and do
not imply completeness differences.

Pair-table status is `COMPLETE` unless an explicit, controlled search/resource
limit truncates the mathematically required region. Such truncation produces
`INCOMPLETE`; invalid symmetry produces an exception instead of an incomplete
table.

Status propagates downward:

```text
PairTable INCOMPLETE
→ ContactAnalysisResult INCOMPLETE
→ dependent shells cannot be RESOLVED
```

With complete geometry:

- all applicable shells resolved produces `RESOLVED`;
- at least one genuine unresolved alternative produces `AMBIGUOUS`;
- missing required search space or evidence produces `INCOMPLETE`;
- `NOT_APPLICABLE` means no applicable chemical analysis for that centre or
  request and never means geometric failure.

Diagnostics remain data and are never parsed to infer status.

## Migration and deletion

The new route is developed behind the existing production entry point until
its scientific acceptance criteria pass. Development proceeds in dependency
order:

```text
SymmetryContext
→ PeriodicSymmetryRelation
→ AsymmetricUnitMapping
→ SymmetryPairFinder / SymmetryPairTable
→ contact interpretations and incidences
→ shell orbits
→ local polyhedron realization
→ structural hierarchy consumers
→ materializers
→ public cutover
→ legacy deletion
```

The following expanded-first objects and stages are removed at cutover rather
than retained as a second API:

- `CoordinationShellResolver`;
- `CrystalChemistryResolution`;
- `build_contact_orbits`;
- expanded-contact shell analysis;
- downstream scientific dependencies on materialized contacts.

The new `ContactAnalysisResult.contacts` property is the only compatibility
view. It is derived from orbit state by the CrIStMa materializer.

No intermediate orbit-first API is released as a stable beta contract before
the cutover is complete.

## Scientific and performance acceptance

Migration acceptance is not blind equality with the legacy route. Results are
checked against:

- hand-derived periodic chains and special-position fixtures;
- analytical multiplicity and incidence expectations;
- exact group, stabilizer, double-coset, ownership, and materialization
  invariants;
- curated real structures;
- legacy values only where the legacy calculation is independently known to
  be correct.

When orbit-first and legacy results disagree, the discrepancy is diagnosed.
Correction of a demonstrated legacy error is acceptance, not regression.

Required real-structure coverage includes at least:

- K7CaY2(B5O10)3, retaining its isolated ring orbit;
- NaCl, producing no structural rings;
- natrolite, retaining framework ring analysis;
- La2Zr2O7 and La0.5Zr0.5O1.75;
- LiB3O5, Na6Mo11O36, NaLiMoO, grossular, and gehlenite;
- organic, mixed-occupancy, split-position, and special-position structures;
- the maintained Desktop CIF corpus and Finder CIF corpus.

The La2Zr2O7 benchmark that currently performs approximately 1.57 million
atom-image lookups and 89.6 million periodic-coordinate comparisons must show
an order-of-magnitude reduction in contact-analysis time on the same machine.
Wall-clock thresholds are benchmark gates, not scientific model fields.

Key invariants include:

- exact relation composition and inversion satisfy the group laws;
- operation and orbit identities are independent of input operation order;
- every site image retains all equivalent coset relations;
- every pair orbit is unique under both endpoint stabilizers and endpoint
  exchange;
- `geometry_orbit_id` is independent of cutoff and distance;
- shell CN uses local incidence multiplicity, never pair multiplicity;
- ambiguous shells retain alternatives and select none;
- mixed occupancy uses interpretation-specific participating occupancy;
- materialized reference-cell count equals pair-orbit reference-cell
  multiplicity;
- materialization never feeds back into scientific analysis;
- controlled incompleteness propagates to every dependent result.
