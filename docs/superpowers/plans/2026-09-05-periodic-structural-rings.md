# Periodic Structural Rings Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a stateless, translation-aware `RingFinder` that identifies locally shortest finite rings inside existing CrIStMa structural blocks and groups crystallographically equivalent instances.

**Architecture:** Public immutable result types live in `rings.py`; lifted periodic shortest-path traversal lives in `_ring_search.py`; symmetry mapping lives in `_ring_symmetry.py`; `RingFinder` composes those pieces without repeating chemistry, neighbour finding, block detection, or CRAFT presentation. Search operates on `PeriodicUnitRef` states and consumes only an existing `StructuralRepresentation`, matching `StructuralBlockResult`, canonical `CrystalStructure`, and its expanded `AtomicView`.

**Tech Stack:** Python 3.11+, standard library, NumPy, pytest; no NetworkX or new runtime dependency.

**Spec:** `docs/superpowers/specs/2026-09-05-structural-rings-and-motifs-design.md`

## Global Constraints

- Production behavior must not branch on formula, material name, element combination, site label, filename, source format, expected ring size, or expected multiplicity.
- `StructuralRing` means a finite zero-translation chordless cycle obtained as a shortest return cycle for at least one connection; it does not mean the complete cycle space.
- Search states are `(unit_id, cell_translation)` and never bare quotient `unit_id` values.
- Ordinary biconnected pruning of the quotient graph is forbidden in this milestone.
- `DIRECT_CONTACT` connections are excluded; `SHARED_VERTEX`, `SHARED_EDGE`, and `SHARED_FACE` are eligible.
- Any exceeded search limit yields `RingAnalysisStatus.INCOMPLETE` plus a diagnostic; truncation is never silent.
- CrIStMa remains Qt-free and stateless, with NumPy as its only runtime dependency.
- Run only the task-specific test file during each red/green cycle; run the crystal-chemistry and inorganic integration slices once at the final gate.

## File Map

- Create `src/cristma/crystal_chemistry/rings.py`: immutable public ring records, status, and search policy.
- Create `src/cristma/crystal_chemistry/_ring_search.py`: directed lifted adjacency and all-shortest-return-path traversal.
- Create `src/cristma/crystal_chemistry/_ring_symmetry.py`: translation-aware atom/unit mapping and orbit grouping.
- Create `src/cristma/crystal_chemistry/ring_finder.py`: validation, candidate filtering, canonicalization, composition, connectors, and public tool orchestration.
- Modify `src/cristma/crystal_chemistry/__init__.py`: export only the public contracts and `RingFinder`.
- Create `tests/crystal_chemistry/test_ring_types.py`: immutable-contract validation.
- Create `tests/crystal_chemistry/test_ring_search.py`: lifted-state traversal and limits.
- Create `tests/crystal_chemistry/test_ring_finder.py`: ring definition, composition, connectors, and no-hardcode invariants.
- Create `tests/crystal_chemistry/test_ring_symmetry.py`: space-group action and orbit grouping.
- Modify `tests/integration/test_inorganic_crystal_chemistry.py`: lithium-triborate end-to-end acceptance.

---

### Task 1: Public ring contracts and deterministic policy

**Files:**
- Create: `src/cristma/crystal_chemistry/rings.py`
- Test: `tests/crystal_chemistry/test_ring_types.py`

**Interfaces:**
- Consumes: `cristma.chemistry.Composition`, `cristma.diagnostics.Diagnostic`, `cristma.structure.PeriodicAtomRef`.
- Produces: `PeriodicUnitRef`, `RingAnalysisStatus`, `RingSearchPolicy`, `StructuralRing`, `StructuralRingOrbit`, `RingAnalysisResult`.

- [ ] **Step 1: Write failing value-object tests**

```python
from cristma.chemistry import Composition
from cristma.crystal_chemistry.rings import (
    PeriodicUnitRef,
    RingAnalysisResult,
    RingAnalysisStatus,
    RingSearchPolicy,
    StructuralRing,
)
from cristma.structure import PeriodicAtomRef


def test_structural_ring_requires_a_finite_closed_cycle() -> None:
    ring = StructuralRing(
        ring_id="ring:abc",
        parent_block_id="block:framework",
        representation_id="representation:structural",
        unit_refs=(
            PeriodicUnitRef("unit:A", (0, 0, 0)),
            PeriodicUnitRef("unit:B", (0, 0, 0)),
            PeriodicUnitRef("unit:C", (0, 0, 0)),
        ),
        connection_ids=("edge:AB", "edge:BC", "edge:CA"),
        connector_atom_refs=(PeriodicAtomRef("O1", (0, 0, 0)),),
        composition=Composition.from_mapping({"A": 1, "B": 1, "C": 1}),
        translation_sum=(0, 0, 0),
        provenance=(("method", "test"),),
    )

    assert ring.size == 3


def test_ring_policy_is_cloneable_and_rejects_invalid_limits() -> None:
    policy = RingSearchPolicy(maximum_ring_size=12, maximum_states_per_connection=50000)
    assert policy.clone(maximum_ring_size=8).maximum_ring_size == 8
    assert policy.get_config()["maximum_ring_size"] == 12


def test_incomplete_result_requires_a_diagnostic() -> None:
    with pytest.raises(ValueError, match="diagnostic"):
        RingAnalysisResult((), (), RingAnalysisStatus.INCOMPLETE)
```

- [ ] **Step 2: Run the contract test and confirm the missing module failure**

Run: `pytest tests/crystal_chemistry/test_ring_types.py -q`

Expected: FAIL during import because `cristma.crystal_chemistry.rings` does not exist.

- [ ] **Step 3: Implement immutable public records**

```python
class RingAnalysisStatus(StrEnum):
    COMPLETE = "complete"
    INCOMPLETE = "incomplete"


@dataclass(frozen=True, slots=True)
class PeriodicUnitRef:
    unit_id: str
    cell_translation: tuple[int, int, int]


@dataclass(frozen=True, slots=True)
class RingSearchPolicy:
    maximum_ring_size: int = 12
    maximum_states_per_connection: int = 50_000
    maximum_paths_per_connection: int = 4_096

    def get_config(self) -> dict[str, int]:
        return {
            "maximum_ring_size": self.maximum_ring_size,
            "maximum_states_per_connection": self.maximum_states_per_connection,
            "maximum_paths_per_connection": self.maximum_paths_per_connection,
        }

    def clone(self, **changes: object) -> "RingSearchPolicy":
        return replace(self, **changes)


@dataclass(frozen=True, slots=True)
class StructuralRing:
    ring_id: str
    parent_block_id: str
    representation_id: str
    unit_refs: tuple[PeriodicUnitRef, ...]
    connection_ids: tuple[str, ...]
    connector_atom_refs: tuple[PeriodicAtomRef, ...]
    composition: Composition
    translation_sum: tuple[int, int, int]
    provenance: tuple[tuple[str, object], ...] = ()

    @property
    def size(self) -> int:
        return len(self.unit_refs)


@dataclass(frozen=True, slots=True)
class StructuralRingOrbit:
    orbit_id: str
    parent_block_id: str
    representation_id: str
    representative_ring_id: str
    ring_ids: tuple[str, ...]
    multiplicity: int
    composition: Composition
    size: int


@dataclass(frozen=True, slots=True)
class RingAnalysisResult:
    rings: tuple[StructuralRing, ...]
    orbits: tuple[StructuralRingOrbit, ...]
    status: RingAnalysisStatus
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()
```

Implement strict validation: non-empty IDs, exactly three integer translation entries, three or more distinct interior unit states, equal unit/connection counts, zero translation sum, unique connector refs and ring IDs, orbit members present in `rings`, one block/representation per orbit, and at least one diagnostic for `INCOMPLETE`.

- [ ] **Step 4: Run the contract tests**

Run: `pytest tests/crystal_chemistry/test_ring_types.py -q`

Expected: PASS.

- [ ] **Step 5: Commit the contracts**

```bash
git add src/cristma/crystal_chemistry/rings.py tests/crystal_chemistry/test_ring_types.py
git commit -m "Add structural ring contracts"
```

---

### Task 2: Translation-aware all-shortest-return-path search

**Files:**
- Create: `src/cristma/crystal_chemistry/_ring_search.py`
- Test: `tests/crystal_chemistry/test_ring_search.py`

**Interfaces:**
- Consumes: `StructuralRepresentation`, `StructuralBlock`, `StructuralConnection`, `PeriodicUnitRef`, `RingSearchPolicy`.
- Produces internally: `_LiftedStep`, `_LiftedPath`, `_ReturnPathResult`, and `find_shortest_return_paths(...)` for `RingFinder`.

- [ ] **Step 1: Write failing lifted-state tests**

Build tiny `StructuralRepresentation` fixtures directly from atomic `StructuralUnit` objects and shared-vertex `StructuralConnection` objects. Assert these cases:

```python
def test_return_path_targets_the_exact_periodic_image() -> None:
    representation, block, removed = periodic_image_fixture()
    result = find_shortest_return_paths(representation, block, removed, RingSearchPolicy())
    assert result.paths
    assert all(path.states[0] == PeriodicUnitRef("unit:A", (0, 0, 0)) for path in result.paths)
    assert all(path.states[-1] == PeriodicUnitRef("unit:B", (1, 0, 0)) for path in result.paths)


def test_equal_shortest_return_paths_are_all_retained() -> None:
    result = find_shortest_return_paths(*diamond_fixture(), RingSearchPolicy())
    assert len(result.paths) == 2


def test_parallel_edges_with_different_translations_remain_distinct() -> None:
    result = find_shortest_return_paths(*periodic_parallel_edge_fixture(), RingSearchPolicy())
    assert {path.states[-1].cell_translation for path in result.paths} == {(1, 0, 0)}


def test_state_budget_returns_explicit_incomplete_search() -> None:
    policy = RingSearchPolicy(maximum_states_per_connection=2)
    result = find_shortest_return_paths(*wide_fixture(), policy)
    assert not result.complete
    assert result.limit_name == "maximum_states_per_connection"
```

- [ ] **Step 2: Run only the lifted search tests**

Run: `pytest tests/crystal_chemistry/test_ring_search.py -q`

Expected: FAIL because `_ring_search` is missing.

- [ ] **Step 3: Implement directed lifted adjacency**

```python
@dataclass(frozen=True, slots=True)
class _LiftedStep:
    connection_id: str
    source: PeriodicUnitRef
    target: PeriodicUnitRef
    first_unit_image: tuple[int, int, int]


def _steps_from(
    state: PeriodicUnitRef,
    connections: tuple[StructuralConnection, ...],
) -> tuple[_LiftedStep, ...]:
    """Translate each quotient connection into exact neighboring lifted states."""
```

Forward traversal adds `connection.lattice_translation`; reverse traversal adds its negation. `first_unit_image` always records the image of `connection.first_unit_id`, which later anchors `shared_atom_refs` correctly for either traversal direction.

- [ ] **Step 4: Implement bounded BFS and predecessor-DAG backtracking**

```python
def find_shortest_return_paths(
    representation: StructuralRepresentation,
    block: StructuralBlock,
    removed: StructuralConnection,
    policy: RingSearchPolicy,
) -> _ReturnPathResult:
    source = PeriodicUnitRef(removed.first_unit_id, (0, 0, 0))
    target = PeriodicUnitRef(removed.second_unit_id, removed.lattice_translation)
    adjacency = _block_adjacency(representation, block)
    predecessors, target_depth, limit_name = _breadth_first_predecessors(
        adjacency,
        source,
        target,
        removed.connection_id,
        policy,
    )
    if limit_name is not None:
        return _ReturnPathResult((), False, limit_name)
    paths, path_limit = _backtrack_all_shortest_paths(
        predecessors,
        source,
        target,
        target_depth,
        policy.maximum_paths_per_connection,
    )
    return _ReturnPathResult(
        paths,
        path_limit is None,
        path_limit,
    )
```

Restrict quotient connections to `block.connection_ids`; do not collapse parallel IDs or translations. A path may not repeat an interior `PeriodicUnitRef`. If state, path, or ring-size limits prevent a conclusion, return `complete=False` and the exact limit name.

- [ ] **Step 5: Run the lifted search tests**

Run: `pytest tests/crystal_chemistry/test_ring_search.py -q`

Expected: PASS.

- [ ] **Step 6: Commit lifted traversal**

```bash
git add src/cristma/crystal_chemistry/_ring_search.py tests/crystal_chemistry/test_ring_search.py
git commit -m "Add periodic-state ring traversal"
```

---

### Task 3: Candidate validation, canonical identity, composition, and connectors

**Files:**
- Create: `src/cristma/crystal_chemistry/ring_finder.py`
- Test: `tests/crystal_chemistry/test_ring_finder.py`

**Interfaces:**
- Consumes: Task 1 records, Task 2 return paths, `CrystalStructure`, `AtomicView[ExpandedAtom]`, `StructuralRepresentation`, `StructuralBlockResult`.
- Produces initially: complete ring instances with empty orbit tuple through `RingFinder.find_instances(...)`; Task 4 adds orbit grouping to `RingFinder.find(...)`.

- [ ] **Step 1: Write failing scientific ring tests**

```python
def test_square_with_diagonal_keeps_triangles_not_composite_square() -> None:
    result = find_instances(square_with_diagonal_fixture())
    assert {ring.size for ring in result.rings} == {3}
    assert len(result.rings) == 2


def test_nonzero_translation_walk_is_not_a_ring() -> None:
    result = find_instances(winding_chain_fixture())
    assert result.rings == ()


def test_ring_connector_atoms_come_from_structural_connections() -> None:
    result = find_instances(shared_edge_triangle_fixture())
    assert result.rings[0].connector_atom_refs == expected_shared_refs()


def test_ring_composition_counts_shared_atoms_once() -> None:
    result = find_instances(three_borate_units_fixture())
    assert result.rings[0].composition.normalized_formula == "B3O7"


def test_direct_contacts_do_not_create_structural_rings() -> None:
    result = find_instances(direct_contact_triangle_fixture())
    assert result.rings == ()
```

- [ ] **Step 2: Run the ring-finder tests and confirm failure**

Run: `pytest tests/crystal_chemistry/test_ring_finder.py -q`

Expected: FAIL because `RingFinder` does not exist.

- [ ] **Step 3: Implement canonical cycle identity**

```python
def _canonical_cycle_key(path: _LiftedPath) -> tuple[object, ...]:
    """Choose the minimum forward/reverse cyclic variant after common shift."""


def _is_chordless(
    unit_refs: tuple[PeriodicUnitRef, ...],
    representation: StructuralRepresentation,
    block: StructuralBlock,
) -> bool:
    """Reject an exact lifted edge between nonconsecutive cycle states."""
```

The canonical key alternates normalized `(unit_id, relative_translation)` tokens with connection IDs so parallel edges stay distinct. Generate every cyclic rotation in forward and reverse order, subtract the first unit translation from each variant, and take the lexicographic minimum.

- [ ] **Step 4: Materialize composition and connector atoms without geometry inference**

For each ring unit state, translate every `StructuralUnit.atom_refs` into the ring frame and take the canonical set union. Look up each unique atom in the supplied `AtomicView`, sum its `SiteComponent` occupancies into `Composition.from_mapping`, and never sum already-aggregated unit formulae.

For each traversed connection, translate `shared_atom_refs` by the image of `connection.first_unit_id` recorded in `_LiftedStep`; deduplicate the resulting `PeriodicAtomRef` values. Do not search coordinates for a bridge atom.

- [ ] **Step 5: Implement instance discovery and incomplete diagnostics**

```python
@dataclass(frozen=True, slots=True)
class RingFinder:
    policy: RingSearchPolicy = field(default_factory=RingSearchPolicy)

    def get_config(self) -> dict[str, int]:
        return self.policy.get_config()

    def clone(self, **changes: object) -> "RingFinder":
        return replace(self, **changes)

    def find_instances(
        self,
        structure: CrystalStructure,
        atomic_view: AtomicView[ExpandedAtom],
        representation: StructuralRepresentation,
        blocks: StructuralBlockResult,
    ) -> RingAnalysisResult:
        _validate_inputs(structure, atomic_view, representation, blocks)
        found: dict[tuple[object, ...], StructuralRing] = {}
        diagnostics: list[Diagnostic] = []
        complete = True
        connection_by_id = {
            item.connection_id: item for item in representation.connections
        }
        for block in blocks.blocks:
            for connection_id in block.connection_ids:
                connection = connection_by_id[connection_id]
                if connection.connection_kind is StructuralConnectionKind.DIRECT_CONTACT:
                    continue
                search = find_shortest_return_paths(
                    representation, block, connection, self.policy
                )
                if not search.complete:
                    complete = False
                    diagnostics.append(_limit_diagnostic(connection, search.limit_name))
                for path in search.paths:
                    if _is_chordless(path.states, representation, block):
                        key = _canonical_cycle_key(path)
                        found.setdefault(
                            key,
                            _materialize_ring(
                                key, path, block, representation, atomic_view
                            ),
                        )
        status = (
            RingAnalysisStatus.COMPLETE
            if complete
            else RingAnalysisStatus.INCOMPLETE
        )
        return RingAnalysisResult(
            tuple(found[key] for key in sorted(found)),
            (),
            status,
            tuple(diagnostics),
            (("method", "cristma.ring_finder:1"),),
        )
```

Validate matching representation IDs and known units/connections before search. Search only eligible connections. Stable `ring_id` values come from SHA-256 of the canonical key, never Python's randomized `hash()`. Aggregate limit diagnostics with code `crystal_chemistry.rings.search_limit_reached` and return `INCOMPLETE` if any connection search was truncated.

- [ ] **Step 6: Run ring-finder tests**

Run: `pytest tests/crystal_chemistry/test_ring_finder.py -q`

Expected: PASS.

- [ ] **Step 7: Commit instance discovery**

```bash
git add src/cristma/crystal_chemistry/ring_finder.py tests/crystal_chemistry/test_ring_finder.py
git commit -m "Find locally shortest structural rings"
```

---

### Task 4: Translation-aware crystallographic ring orbits

**Files:**
- Create: `src/cristma/crystal_chemistry/_ring_symmetry.py`
- Modify: `src/cristma/crystal_chemistry/ring_finder.py`
- Test: `tests/crystal_chemistry/test_ring_symmetry.py`

**Interfaces:**
- Consumes: `AffineOperation`, `ExpandedAtom` symmetry provenance, `AtomicView`, representation units/connections, canonical ring instances.
- Produces internally: `map_periodic_atom_ref(...)`, `map_periodic_unit_ref(...)`, `build_ring_orbits(...)`; updates `RingFinder.find(...)` to return rings plus orbits.

- [ ] **Step 1: Write failing translation-aware symmetry tests**

```python
def test_symmetry_maps_periodic_atom_image_with_rotation_and_wrap_shift() -> None:
    mapped = map_periodic_atom_ref(operation, PeriodicAtomRef("atom:A", (1, 0, 0)), view)
    assert mapped == PeriodicAtomRef("atom:A-prime", expected_image_translation)


def test_unit_mapping_uses_transformed_atom_membership_not_unit_name() -> None:
    mapped = map_periodic_unit_ref(operation, PeriodicUnitRef("arbitrary-name", (0, 0, 0)), view, representation)
    assert mapped == PeriodicUnitRef("another-arbitrary-name", expected_shift)


def test_symmetry_equivalent_rings_form_one_orbit() -> None:
    orbits = build_ring_orbits(structure, view, representation, symmetry_ring_fixture())
    assert len(orbits) == 1
    assert orbits[0].multiplicity == 2


def test_topology_equal_but_not_symmetry_related_rings_stay_separate() -> None:
    orbits = build_ring_orbits(identity_structure, view, representation, two_equal_rings())
    assert len(orbits) == 2
```

- [ ] **Step 2: Run the symmetry tests and confirm failure**

Run: `pytest tests/crystal_chemistry/test_ring_symmetry.py -q`

Expected: FAIL because `_ring_symmetry` is missing.

- [ ] **Step 3: Implement periodic atom mapping with the correct sign convention**

For `x' = R x + τ`, transform the full periodic fractional coordinate `fractional + cell_translation`. Match the wrapped coordinate to an `ExpandedAtom` with the same `source_site_id`. If `_wrap_with_translation` returns additive normalization `n` such that `wrapped = raw + n`, the mapped periodic image is `-n`; do not confuse that value with `SymmetryImageProvenance.normalization_translation`.

```python
def map_periodic_atom_ref(
    operation: AffineOperation,
    atom_ref: PeriodicAtomRef,
    view: AtomicView[ExpandedAtom],
    tolerance: float = DEFAULT_FRACTIONAL_TOLERANCE,
) -> PeriodicAtomRef:
    atom_by_id = {atom.id: atom for atom in view.atoms}
    source = atom_by_id[atom_ref.atom_id]
    full = np.asarray(source.fractional) + np.asarray(atom_ref.cell_translation)
    rotation = np.asarray(operation.rotation, dtype=float)
    offset = np.asarray(operation.translation, dtype=float)
    raw = rotation @ full + offset
    additive_normalization = -np.floor(raw + tolerance).astype(int)
    wrapped = raw + additive_normalization
    target = _match_expanded_atom(
        wrapped,
        source.source_site_id,
        view,
        tolerance,
    )
    return PeriodicAtomRef(
        target.id,
        tuple(int(-value) for value in additive_normalization),
    )
```

- [ ] **Step 4: Implement unit mapping through canonical transformed atom membership**

Transform every atom ref of the source unit. Normalize the transformed membership by one common lattice shift and match `(unit kind, normalized atom membership, interaction layers, contact classifications)` against units in the same representation. Return the matched `unit_id` plus common shift. Missing or non-unique mappings produce an explicit symmetry diagnostic rather than name parsing.

- [ ] **Step 5: Group ring instances by actual space-group action**

Apply every `structure.space_group.operations` item to each ring's unit refs, connector refs, and connection endpoints; map back to representation IDs; canonicalize the transformed ring again; and union only identities that match an already discovered ring in the same block and representation.

```python
def build_ring_orbits(
    structure: CrystalStructure,
    view: AtomicView[ExpandedAtom],
    representation: StructuralRepresentation,
    rings: tuple[StructuralRing, ...],
) -> tuple[tuple[StructuralRingOrbit, ...], tuple[Diagnostic, ...]]:
    links, diagnostics = _symmetry_links(
        structure,
        view,
        representation,
        rings,
    )
    groups = _connected_ring_components(
        tuple(ring.ring_id for ring in rings),
        links,
    )
    ring_by_id = {ring.ring_id: ring for ring in rings}
    orbits = tuple(
        _orbit_from_group(group, ring_by_id)
        for group in groups
    )
    return orbits, diagnostics
```

Orbit IDs use a stable digest of sorted member ring IDs. `multiplicity == len(ring_ids)`. Equal formula and size never create an orbit without a successful symmetry mapping.

- [ ] **Step 6: Compose orbit grouping in the public tool**

```python
def find(
    self,
    structure: CrystalStructure,
    atomic_view: AtomicView[ExpandedAtom],
    representation: StructuralRepresentation,
    blocks: StructuralBlockResult,
) -> RingAnalysisResult:
    instances = self.find_instances(structure, atomic_view, representation, blocks)
    orbits, diagnostics = build_ring_orbits(
        structure, atomic_view, representation, instances.rings
    )
    return replace(
        instances,
        orbits=orbits,
        diagnostics=instances.diagnostics + diagnostics,
    )
```

- [ ] **Step 7: Run symmetry and instance tests**

Run: `pytest tests/crystal_chemistry/test_ring_symmetry.py tests/crystal_chemistry/test_ring_finder.py -q`

Expected: PASS.

- [ ] **Step 8: Commit symmetry grouping**

```bash
git add src/cristma/crystal_chemistry/_ring_symmetry.py src/cristma/crystal_chemistry/ring_finder.py tests/crystal_chemistry/test_ring_symmetry.py
git commit -m "Group structural rings by crystal symmetry"
```

---

### Task 5: Public API and lithium-triborate end-to-end acceptance

**Files:**
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Modify: `tests/integration/test_inorganic_crystal_chemistry.py`
- Test: `tests/crystal_chemistry/test_ring_types.py`
- Test: `tests/crystal_chemistry/test_ring_search.py`
- Test: `tests/crystal_chemistry/test_ring_finder.py`
- Test: `tests/crystal_chemistry/test_ring_symmetry.py`

**Interfaces:**
- Consumes: all earlier tasks and the existing inorganic calculation helpers in `tests/integration/test_inorganic_crystal_chemistry.py`.
- Produces: stable imports from `cristma.crystal_chemistry` and an end-to-end CrIStMa result ready for CRAFT.

- [ ] **Step 1: Write the failing end-to-end acceptance**

Extend the integration helper to return `structure`, `view`, `representation`, and `blocks`, then add:

```python
def test_lithium_triborate_has_b3o7_ring_orbit_inside_bo_framework() -> None:
    calculation = build_structural_calculation(
        "LiB3O5_3000122.cif",
        {InteractionLayer.STRUCTURAL},
    )
    result = RingFinder().find(
        calculation.structure,
        calculation.view,
        calculation.representation,
        calculation.blocks,
    )

    assert result.status is RingAnalysisStatus.COMPLETE
    assert len(result.orbits) == 1
    orbit = result.orbits[0]
    assert orbit.composition.normalized_formula == "B3O7"
    assert orbit.size == 3
    assert orbit.multiplicity == 4
    assert {
        block.periodic_rank
        for block in calculation.blocks.blocks
        if block.block_id == orbit.parent_block_id
    } == {3}
```

Add a second generic regression that rebuilds an isomorphic synthetic graph with unrelated unit/site labels and different chemical components; assert unchanged ring sizes/counts while only calculated composition changes. Add a static guard:

Run: `! rg -n "LiB3O5|B3O7|lithium_triborate" src/cristma`

Expected: exit status 0 from shell negation because production sources contain none of those fixture names.

- [ ] **Step 2: Export the public API**

Add these names to `cristma.crystal_chemistry.__all__` and imports:

```python
PeriodicUnitRef
RingAnalysisResult
RingAnalysisStatus
RingFinder
RingSearchPolicy
StructuralRing
StructuralRingOrbit
```

- [ ] **Step 3: Run the lithium-triborate acceptance test**

Run: `pytest tests/integration/test_inorganic_crystal_chemistry.py::test_lithium_triborate_has_b3o7_ring_orbit_inside_bo_framework -q`

Expected: PASS with one orbit, size 3, multiplicity 4, and formula `B3O7`.

- [ ] **Step 4: Run the focused ring suite**

Run: `pytest tests/crystal_chemistry/test_ring_types.py tests/crystal_chemistry/test_ring_search.py tests/crystal_chemistry/test_ring_finder.py tests/crystal_chemistry/test_ring_symmetry.py -q`

Expected: PASS.

- [ ] **Step 5: Run the affected scientific regression slice**

Run: `pytest tests/crystal_chemistry tests/integration/test_inorganic_crystal_chemistry.py -q`

Expected: PASS. Do not run unrelated I/O, diffraction, or packaging tests for this change.

- [ ] **Step 6: Verify import and source independence**

Run: `python -c "from cristma.crystal_chemistry import RingFinder, RingSearchPolicy, StructuralRing; print(RingFinder(RingSearchPolicy()).get_config())"`

Expected: a dictionary containing the three deterministic search limits.

Run: `! rg -n "LiB3O5|B3O7|lithium_triborate|networkx|PySide6" src/cristma/crystal_chemistry`

Expected: success with no matches.

- [ ] **Step 7: Commit the completed CrIStMa slice**

```bash
git add src/cristma/crystal_chemistry/__init__.py tests/integration/test_inorganic_crystal_chemistry.py
git commit -m "Expose structural ring analysis"
```

Record the CrIStMa commit hash for the dependent CRAFT integration plan.
