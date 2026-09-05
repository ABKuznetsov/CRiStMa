# Periodic Connectivity and Structural Blocks Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Select an explicit interaction representation, calculate exact periodic rank for every graph component, and expose finite blocks, one-periodic systems, layers, and frameworks without motif or rigidity inference.

**Architecture:** `StructuralRepresentationBuilder` filters the existing unit graph using semantic interaction roles already present in CrIStMa results. `PeriodicConnectivityAnalyzer` treats the selected graph as an integer gain graph and derives exact cycle-closure rank. `StructuralBlockFinder` converts each analysed component into an immutable crystal-chemical block while preserving canonical unit, atom, connection, and translation mappings.

**Tech Stack:** Python 3.11, frozen dataclasses, integer/Fraction linear algebra, existing CrIStMa structural unit graph, pytest.

**Spec:** `docs/superpowers/specs/2026-09-02-structural-blocks-and-periodic-connectivity-design.md`

## Global Constraints

- CrIStMa remains Qt-free, application-independent, and free of runtime Gemmi, pymatgen, or spglib dependencies.
- Do not rerun Chemistry, neighbour search, contact resolution, or polyhedron construction.
- Connectivity authority is limited to `StructuralUnit`, `StructuralConnection`, and their preserved semantic provenance.
- Periodic rank uses exact integer arithmetic; coordinates and floating tolerances cannot alter it.
- A cross-cell spanning-tree edge alone is not evidence of periodicity; only a non-zero cycle closure contributes to rank.
- The caller chooses the interaction representation. CrIStMa does not select one globally preferred structural interpretation.
- Do not implement rings, motifs, morphology refinement, mechanical rigidity, hinges, or refinement parameterization.

---

### Task 1: Explicit structural representations

**Files:**
- Modify: `src/cristma/crystal_chemistry/structural_units.py`
- Create: `src/cristma/crystal_chemistry/representation.py`
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Modify: `tests/crystal_chemistry/test_structural_units.py`
- Create: `tests/crystal_chemistry/test_structural_representation.py`

**Interfaces:**
- Extends: `StructuralUnit.interaction_layers` and `StructuralUnit.contact_classifications`.
- Produces: `StructuralSelectionPolicy`, `StructuralRepresentation`, and `StructuralRepresentationBuilder.build(graph)`.

- [ ] **Step 1: Write failing tests for unit semantics and representation selection**

```python
def test_unit_builder_preserves_semantics_from_source_contacts():
    unit = StructuralUnitBuilder().build(resolution, (polyhedron,)).units[0]
    assert unit.interaction_layers == (InteractionLayer.STRUCTURAL,)
    assert unit.contact_classifications == (ContactClassification.PRIMARY,)


def test_structural_representation_excludes_interstitial_units_and_edges():
    policy = StructuralSelectionPolicy(
        included_layers=frozenset({InteractionLayer.STRUCTURAL}),
        included_classifications=frozenset({ContactClassification.PRIMARY}),
    )
    representation = StructuralRepresentationBuilder(policy).build(graph)
    assert tuple(unit.unit_id for unit in representation.units) == ("unit:MoO4",)
    assert representation.excluded_unit_ids == ("unit:CaO8",)
    assert representation.excluded_connection_ids == ("connection:Ca-Mo",)
```

- [ ] **Step 2: Run focused tests and verify RED**

Run:

```bash
python3 -m pytest -p no:cacheprovider \
  tests/crystal_chemistry/test_structural_units.py \
  tests/crystal_chemistry/test_structural_representation.py -q
```

Expected: new semantic fields and representation contracts are absent.

- [ ] **Step 3: Preserve unit semantics and implement immutable selection**

```python
@dataclass(frozen=True, slots=True)
class StructuralSelectionPolicy:
    included_layers: frozenset[InteractionLayer]
    included_classifications: frozenset[ContactClassification]


@dataclass(frozen=True, slots=True)
class StructuralRepresentation:
    representation_id: str
    units: tuple[StructuralUnit, ...]
    connections: tuple[StructuralConnection, ...]
    selection_policy: StructuralSelectionPolicy
    excluded_unit_ids: tuple[str, ...]
    excluded_connection_ids: tuple[str, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()
```

Populate each unit's interaction layers and classifications from its source `ResolvedContact` records. A unit is selected when both semantic sets intersect the policy. A connection is selected only when both endpoints are selected and its own layer and classification intersect the same policy. Preserve deterministic ordering and every excluded ID.

- [ ] **Step 4: Run focused tests and verify GREEN**

Run:

```bash
python3 -m pytest -p no:cacheprovider \
  tests/crystal_chemistry/test_structural_units.py \
  tests/crystal_chemistry/test_structural_representation.py -q
```

Expected: all focused tests pass.

- [ ] **Step 5: Commit the representation slice**

```bash
git add src/cristma/crystal_chemistry/structural_units.py \
  src/cristma/crystal_chemistry/representation.py \
  src/cristma/crystal_chemistry/__init__.py \
  tests/crystal_chemistry/test_structural_units.py \
  tests/crystal_chemistry/test_structural_representation.py
git commit -m "Add explicit structural representations"
```

### Task 2: Exact periodic connectivity

**Files:**
- Create: `src/cristma/crystal_chemistry/periodic_connectivity.py`
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Create: `tests/crystal_chemistry/test_periodic_connectivity.py`

**Interfaces:**
- Consumes: `PeriodicConnectivityAnalyzer.analyze(representation: StructuralRepresentation)`.
- Produces: `PeriodicComponent` and `PeriodicConnectivityResult` with exact ranks and deterministic generators.

- [ ] **Step 1: Write failing exact-rank tests**

```python
def test_cross_cell_tree_edge_is_finite():
    result = analyzer_for_edges((edge("A", "B", (1, 0, 0)),))
    assert result.components[0].periodic_rank == 0
    assert result.components[0].periodic_generators == ()


def test_periodic_self_edge_has_rank_one():
    result = analyzer_for_edges((edge("A", "A", (2, 0, 0)),))
    assert result.components[0].periodic_rank == 1
    assert result.components[0].periodic_generators == ((2, 0, 0),)


@pytest.mark.parametrize(
    ("translations", "expected_rank"),
    (
        (((1, 0, 0),), 1),
        (((1, 0, 0), (0, 1, 0)), 2),
        (((1, 0, 0), (0, 1, 0), (0, 0, 1)), 3),
    ),
)
def test_periodic_rank_is_exact(translations, expected_rank):
    result = analyzer_for_periodic_self_edges(translations)
    assert result.components[0].periodic_rank == expected_rank
```

Also build equivalent graphs with reversed unit and connection order and assert identical component membership, rank, and generators.

- [ ] **Step 2: Run the connectivity tests and verify RED**

Run: `python3 -m pytest -p no:cacheprovider tests/crystal_chemistry/test_periodic_connectivity.py -q`

Expected: import failure because the connectivity contracts do not exist.

- [ ] **Step 3: Implement gain-graph traversal and exact rank**

```python
@dataclass(frozen=True, slots=True)
class PeriodicComponent:
    component_id: str
    unit_ids: tuple[str, ...]
    connection_ids: tuple[str, ...]
    image_offsets: tuple[tuple[str, tuple[int, int, int]], ...]
    closure_translations: tuple[tuple[int, int, int], ...]
    periodic_rank: int
    periodic_generators: tuple[tuple[int, int, int], ...]


@dataclass(frozen=True, slots=True)
class PeriodicConnectivityResult:
    representation_id: str
    components: tuple[PeriodicComponent, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()
```

For each undirected component, assign integer image offsets through a deterministic spanning tree. For every edge, calculate `offset[first] + edge_translation - offset[second]`; retain non-zero closures. Calculate rank by exact Gaussian elimination over `fractions.Fraction`. Sort closure vectors deterministically and retain the first vectors that increase rank as generators without dividing their magnitude.

- [ ] **Step 4: Run connectivity tests and verify GREEN**

Run: `python3 -m pytest -p no:cacheprovider tests/crystal_chemistry/test_periodic_connectivity.py -q`

Expected: all rank, tree-edge, self-edge, and order-invariance tests pass.

- [ ] **Step 5: Commit exact connectivity**

```bash
git add src/cristma/crystal_chemistry/periodic_connectivity.py \
  src/cristma/crystal_chemistry/__init__.py \
  tests/crystal_chemistry/test_periodic_connectivity.py
git commit -m "Calculate exact periodic connectivity"
```

### Task 3: Crystal-chemical blocks and scientific acceptance

**Files:**
- Create: `src/cristma/crystal_chemistry/structural_blocks.py`
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Create: `tests/crystal_chemistry/test_structural_blocks.py`
- Modify: `tests/integration/test_inorganic_crystal_chemistry.py`
- Modify: `docs/inorganic-crystal-chemistry.md`

**Interfaces:**
- Consumes: `StructuralBlockFinder.find(representation, connectivity)`.
- Produces: `StructuralBlockClassification`, `StructuralBlock`, and `StructuralBlockResult`.

- [ ] **Step 1: Write failing block-classification tests**

```python
@pytest.mark.parametrize(
    ("rank", "classification"),
    (
        (0, StructuralBlockClassification.FINITE_BLOCK),
        (1, StructuralBlockClassification.ONE_PERIODIC),
        (2, StructuralBlockClassification.LAYER),
        (3, StructuralBlockClassification.FRAMEWORK),
    ),
)
def test_block_classification_follows_exact_rank(rank, classification):
    result = StructuralBlockFinder().find(representation, connectivity_for_rank(rank))
    assert result.blocks[0].classification is classification
```

Add integration assertions:

```python
def test_camo4_primary_structural_representation_keeps_moo4_finite():
    blocks = build_blocks("CaMoO4_9009632.cif", {InteractionLayer.STRUCTURAL})
    assert blocks.blocks
    assert {block.periodic_rank for block in blocks.blocks} == {0}


def test_lib3o5_primary_boron_oxygen_component_is_periodic():
    blocks = build_blocks("LiB3O5_3000122.cif", {InteractionLayer.STRUCTURAL})
    assert max(block.periodic_rank for block in blocks.blocks) >= 1
```

The LiB3O5 assertion first verifies genuine periodicity without forcing rank 3 if the current finite graph lacks a required connection; any lower-than-known rank must remain an explicit scientific follow-up, not be repaired with a fixture-specific branch.

- [ ] **Step 2: Run focused block and integration tests and verify RED**

Run:

```bash
python3 -m pytest -p no:cacheprovider \
  tests/crystal_chemistry/test_structural_blocks.py \
  tests/integration/test_inorganic_crystal_chemistry.py -q
```

Expected: block contracts are absent.

- [ ] **Step 3: Implement block records as a projection of connectivity**

```python
class StructuralBlockClassification(StrEnum):
    FINITE_BLOCK = "finite_block"
    ONE_PERIODIC = "one_periodic"
    LAYER = "layer"
    FRAMEWORK = "framework"


@dataclass(frozen=True, slots=True)
class StructuralBlock:
    block_id: str
    representation_id: str
    unit_ids: tuple[str, ...]
    atom_refs: tuple[PeriodicAtomRef, ...]
    connection_ids: tuple[str, ...]
    periodic_rank: int
    periodic_generators: tuple[tuple[int, int, int], ...]
    classification: StructuralBlockClassification
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()


@dataclass(frozen=True, slots=True)
class StructuralBlockResult:
    representation_id: str
    blocks: tuple[StructuralBlock, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()
```

Use the component image offsets to express every member atom in one deterministic lifted frame. Deduplicate canonical `PeriodicAtomRef` values. Reject connectivity results belonging to another representation rather than guessing correspondence.

- [ ] **Step 4: Document the boundary and run complete verification**

Document:

```text
StructuralUnitGraph
    -> explicit StructuralRepresentation
    -> exact PeriodicConnectivityResult
    -> StructuralBlockResult
```

State that rank-1 remains `one_periodic` until an independent morphology tool distinguishes chain from ribbon, and that motifs and mechanical rigidity remain absent.

Run:

```bash
python3 -m pytest -p no:cacheprovider -q
python3 -m pip wheel . --no-deps --no-build-isolation \
  --wheel-dir /tmp/cristma-periodic-block-dist
unzip -l /tmp/cristma-periodic-block-dist/cristma-0.1.0.dev0-py3-none-any.whl \
  | rg 'representation|periodic_connectivity|structural_blocks'
git diff --check
```

Expected: zero test failures, wheel build exit code 0, all three modules packaged, and no whitespace errors.

- [ ] **Step 5: Commit the completed block slice**

```bash
git add src/cristma/crystal_chemistry/structural_blocks.py \
  src/cristma/crystal_chemistry/__init__.py \
  tests/crystal_chemistry/test_structural_blocks.py \
  tests/integration/test_inorganic_crystal_chemistry.py \
  docs/inorganic-crystal-chemistry.md
git commit -m "Classify periodic structural blocks"
```
