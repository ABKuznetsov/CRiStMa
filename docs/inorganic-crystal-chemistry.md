# Inorganic crystal chemistry

CRiStMa separates geometric enumeration, chemical interpretation, shell
resolution, and polyhedron construction. Applications own the order and cache
the immutable results they need.

```python
import cristma
from cristma.chemistry import ChemistryAnalyzer, Composition
from cristma.crystal_chemistry import (
    CoordinationShellResolver,
    PolyhedronBuilder,
    ResolutionStatus,
    ShellResolutionPolicy,
)

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
```

This is also the intended CRAFT boundary: CRAFT retains `structure`,
`resolution`, and `polyhedra`, then decides how to display, select, compare, or
cache them. CRAFT does not reproduce symmetry expansion, contact resolution,
or convex-hull mathematics.

## Scientific stages

1. `ChemistryAnalyzer` classifies the composition and emits interaction
   requests without reading coordinates.
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
enumeration. CRiStMa packages 101 source records for 96 elements and preserves
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
