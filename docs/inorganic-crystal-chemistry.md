# Inorganic crystal chemistry

CrIStMa separates geometric enumeration, chemical interpretation, shell
resolution, and polyhedron construction. Applications own the order and cache
the immutable results they need.

```python
import cristma
from cristma.chemistry import ChemistryAnalyzer, Composition
from cristma.crystal_chemistry import (
    CoordinationShellResolver,
    PeriodicConnectivityAnalyzer,
    PolyhedronBuilder,
    ResolutionStatus,
    ShellResolutionPolicy,
    StructuralBlockFinder,
    StructuralGraphBuilder,
    StructuralRepresentationBuilder,
    StructuralSelectionPolicy,
    StructuralUnitBuilder,
)
from cristma.chemistry import InteractionLayer
from cristma.crystal_chemistry import ContactClassification

structure = cristma.read("sample.cif").structures[0]
chemistry = ChemistryAnalyzer().analyze(Composition.from_structure(structure))

# Explicit example calibration, not a universal preset.
policy = ShellResolutionPolicy(
    candidate_rho_max=1.60,
    distance_group_tolerance=0.01,
    minimum_shell_gap=0.08,
    ambiguity_tolerance=0.01,
    search_rho_max=2.00,
)
resolution = CoordinationShellResolver(policy).resolve(
    structure,
    chemistry.grammar,
)

view = structure.atomic_view()
polyhedron_results = tuple(
    PolyhedronBuilder().build(shell, view)
    for shell in resolution.coordination_shells
)
polyhedra = tuple(
    item.polyhedron
    for item in polyhedron_results
    if item.status is ResolutionStatus.RESOLVED and item.polyhedron is not None
)

unit_result = StructuralUnitBuilder().build(resolution, polyhedra)
unit_graph = StructuralGraphBuilder().build(
    unit_result.units,
    resolution.contacts,
)

selection = StructuralSelectionPolicy(
    included_layers=frozenset({InteractionLayer.STRUCTURAL}),
    included_classifications=frozenset({ContactClassification.PRIMARY}),
)
representation = StructuralRepresentationBuilder(selection).build(unit_graph)
connectivity = PeriodicConnectivityAnalyzer().analyze(representation)
blocks = StructuralBlockFinder().find(representation, connectivity)
```

Callers retain `structure`, `resolution`, and `polyhedra`, then decide how to
display, select, compare, or cache them. CrIStMa keeps the reusable symmetry
expansion, contact resolution, and convex-hull mathematics independent of that
workflow.

## Structural units and their periodic graph

The next derived layer consumes only results already calculated by crystal
chemistry:

```text
CrystalChemistryResolution + CoordinationPolyhedron
    -> StructuralUnitBuilder
    -> StructuralUnit
    -> StructuralGraphBuilder
    -> StructuralUnitGraph
```

A polyhedron becomes one unit whose centre is in the reference cell and whose
ligands retain their exact integer lattice translations. A resolved-contact
endpoint that belongs to no polyhedron remains available as an atomic unit.
The graph then records two independent kinds of evidence:

- common unit membership gives `shared_vertex`, `shared_edge`, or
  `shared_face`;
- a supplied `ResolvedContact` gives `direct_contact`.

Every connection is canonical under reversal, so `(A, B, t)` and
`(B, A, -t)` describe one physical relation. Source contact IDs,
`InteractionLayer`, and `ContactClassification` remain attached to the graph
connection. The graph never accepts raw distances and never reruns Chemistry,
neighbour search, shell resolution, or polyhedron construction.

`StructuralUnitGraph` itself remains an unclassified finite periodic quotient
graph. Classification is an explicit later operation:

```text
StructuralUnitGraph
    -> StructuralRepresentationBuilder
    -> StructuralRepresentation
    -> PeriodicConnectivityAnalyzer
    -> PeriodicConnectivityResult
    -> StructuralBlockFinder
    -> StructuralBlockResult
```

The selection policy reuses the `InteractionLayer` and
`ContactClassification` already assigned by Chemistry and crystal chemistry;
it does not classify the composition again. Excluded unit and connection IDs
remain in the representation, so applications can explain why an interstitial
or secondary interaction did not affect a particular block model.

Periodic rank is calculated from exact integer cycle closures in the selected
gain graph. A lattice-crossing tree edge only chooses an image and remains
rank 0. Independent non-zero closures produce:

```text
rank 0  finite_block
rank 1  one_periodic
rank 2  layer
rank 3  framework
```

Rank 1 is intentionally called `one_periodic`: distinguishing a chain from a
ribbon requires a later morphology analysis. Rings and other motifs are also
not searched in this slice. A crystal-chemical block makes no claim of
mechanical rigidity.

Current scientific fixtures give four isolated rank-0 MoO4 blocks for
CaMoO4, one rank-3 B-O framework for LiB3O5, and one rank-3 framework for FeS2
when its Fe-S coordination and S-S subsystem contacts are selected together.

## Scientific stages

1. `ChemistryAnalyzer` classifies the composition and compiles interaction
   requests from Reference Data `grammar_templates`, without reading
   coordinates or using material-family Python branches.
2. `NeighborFinder` enumerates periodic geometry using a cutoff derived from
   the interaction grammar and the complete Cordero covalent-radius catalog.
3. `CoordinationShellResolver` groups normalized distances and chooses a shell
   boundary by explicit hard rules and lexicographic criteria.
4. The decision is made collectively for every crystallographic orbit and
   projected onto its expanded centres.
5. `PolyhedronBuilder` builds a convex coordination polyhedron only for a
   resolved three-dimensional shell.

The search horizon is deliberately separate from the maximum plausible
first-shell distance. An outer group must be observed before a boundary can be
called resolved. Contacts after the selected boundary are retained as
`SECONDARY`, allowing viewers to expose long contacts without adding them to
the geometric coordination number.

Three independent meanings are retained on the way to the result:

```text
GrammarOperation       chemical mode of the requested interaction
InteractionLayer       structural / interstitial / coordination / ... role
ContactClassification  primary or secondary geometric-shell membership
```

`ResolvedContact` keeps all three. A later structural graph can therefore be
built directly from resolved contacts without rerunning Chemistry or the
neighbour search. Current acceptance examples are:

```text
CaMoO4  Mo-O structural;    Ca-O interstitial
LiB3O5  B-O structural;     Li-O interstitial
FeS2    Fe-S coordination;  S-S intra-subsystem
```

## Result statuses

- `RESOLVED`: one boundary is supported by the configured rules;
- `AMBIGUOUS`: more than one boundary remains comparably plausible;
- `INCOMPLETE`: the available neighborhood or reference data cannot establish
  a boundary, or symmetry-equivalent centres disagree;
- `NOT_APPLICABLE`: a valid shell exists but does not define a volumetric
  polyhedron, for example a planar BO3 group.

Mixed occupancy remains one geometric position. `geometric_CN` counts
positions, while `mean_occupied_neighbors` records their statistical
occupation.

## Radii and current limits

Cordero covalent radii are used for oxidation- and CN-independent candidate
enumeration. CrIStMa packages 101 source records for 96 elements and preserves
the published carbon and spin variants; no missing element receives a guessed
fallback.

Shannon radii are indexed exactly by element, oxidation state, coordination,
and spin. They are secondary ionic-distance evidence after a candidate shell
exists, because using a CN-dependent radius to discover CN would be circular.
The current `ShannonDistanceValidator` is an independent tool; automatic
oxidation-state assignment and its integration into shell evidence are a
later crystal-chemistry slice.

BVS and coordination-shape evidence are currently reported as unavailable or
not applicable rather than simulated. They will be added as independent
analyzers. No expected-CN table, material-name branch, Voronoi dependency, or
weighted hidden bond score is used.

Every result records the resolution policy, search cutoff, observed normalized
distance range, structure identity, grammar/reference version, and method
version in provenance.
