# Structural Units and Graph Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert existing polyhedra and semantically resolved contacts into immutable structural units and a finite periodic quotient graph, without repeating chemistry or neighbour search.

**Architecture:** `StructuralUnitBuilder` maps polyhedra and unrepresented contact endpoints to canonical unit membership. `StructuralGraphBuilder` then derives shared-atom and direct-contact connections exclusively from `ResolvedContact`, preserving lattice translations and interaction layers. Periodic rank, blocks, and motifs are intentionally excluded.

**Tech Stack:** Python 3.11, frozen dataclasses, existing CRiStMa identity/contact/polyhedron types, pytest.

**Spec:** `docs/superpowers/specs/2026-09-02-structural-blocks-and-periodic-connectivity-design.md`

## Global Constraints

- CRiStMa remains Qt-free and application-independent.
- Do not run Chemistry, `NeighborFinder`, or contact resolution inside this layer.
- Use only canonical atom IDs, integer lattice translations, `ResolvedContact`, and `CoordinationPolyhedron` as connectivity evidence.
- Preserve `GrammarOperation`, `InteractionLayer`, and `ContactClassification` through source contacts.
- Inputs remain immutable; tools store configuration only and return explicit immutable results.
- Do not implement periodic rank, blocks, chains, layers, frameworks, rings, motifs, or mechanical rigidity in this plan.

---

### Task 1: Canonical structural-unit records

**Files:**
- Create: `src/cristma/crystal_chemistry/structural_units.py`
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Test: `tests/crystal_chemistry/test_structural_units.py`

**Interfaces:**
- Consumes: `PeriodicAtomRef`, `ResolvedContact`, `CoordinationPolyhedron`.
- Produces: `StructuralUnitKind`, `StructuralUnit`, `StructuralUnitBuildResult`, `StructuralUnitBuilder.build(resolution, polyhedra)`.

- [ ] **Step 1: Write failing model and builder tests**

```python
def test_polyhedron_unit_retains_periodic_atom_membership_and_sources():
    result = StructuralUnitBuilder().build(resolution, (polyhedron,))
    unit = result.units[0]
    assert unit.kind is StructuralUnitKind.POLYHEDRON
    assert unit.atom_refs == (
        PeriodicAtomRef("Mo", (0, 0, 0)),
        PeriodicAtomRef("O", (1, 0, 0)),
    )
    assert unit.source_contact_ids == ("contact:Mo|O|1,0,0",)
    assert unit.source_polyhedron_id == polyhedron.polyhedron_id


def test_unrepresented_contact_endpoint_gets_one_atomic_unit():
    result = StructuralUnitBuilder().build(resolution, ())
    assert {(unit.kind, unit.atom_refs[0].atom_id) for unit in result.units} == {
        (StructuralUnitKind.ATOM, "S1"),
        (StructuralUnitKind.ATOM, "S2"),
    }
```

- [ ] **Step 2: Run the focused tests and verify RED**

Run: `python3 -m pytest -p no:cacheprovider tests/crystal_chemistry/test_structural_units.py -q`

Expected: import failure because the structural-unit contracts do not exist.

- [ ] **Step 3: Implement immutable records and periodic membership mapping**

```python
class StructuralUnitKind(StrEnum):
    POLYHEDRON = "polyhedron"
    ATOM = "atom"


@dataclass(frozen=True, slots=True)
class StructuralUnit:
    unit_id: str
    kind: StructuralUnitKind
    atom_refs: tuple[PeriodicAtomRef, ...]
    source_contact_ids: tuple[str, ...]
    source_polyhedron_id: str | None
    provenance: tuple[tuple[str, object], ...] = ()


@dataclass(frozen=True, slots=True)
class StructuralUnitBuildResult:
    units: tuple[StructuralUnit, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()
```

For each polyhedron, place its centre at translation `(0, 0, 0)`. Derive each ligand image from the orientation of its `GeometricContact`: use the stored translation when the centre is the first endpoint and its negation when the centre is the second endpoint. Add one atomic unit per contact endpoint not represented by any polyhedral unit. Sort and deduplicate references and IDs deterministically.

- [ ] **Step 4: Run the focused tests and verify GREEN**

Run: `python3 -m pytest -p no:cacheprovider tests/crystal_chemistry/test_structural_units.py -q`

Expected: all tests pass.

- [ ] **Step 5: Commit the unit contracts**

```bash
git add src/cristma/crystal_chemistry/structural_units.py src/cristma/crystal_chemistry/__init__.py tests/crystal_chemistry/test_structural_units.py
git commit -m "Add canonical structural units"
```

### Task 2: Structural connections from membership and contacts

**Files:**
- Create: `src/cristma/crystal_chemistry/structural_graph.py`
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Test: `tests/crystal_chemistry/test_structural_graph.py`

**Interfaces:**
- Consumes: `StructuralUnitBuildResult.units`, `CrystalChemistryResolution.contacts`.
- Produces: `StructuralConnectionKind`, `StructuralConnection`, `StructuralUnitGraph`, `StructuralGraphBuilder.build(units, contacts)`.

- [ ] **Step 1: Write failing tests for sharing and direct contacts**

```python
@pytest.mark.parametrize(
    ("shared_count", "kind"),
    ((1, "shared_vertex"), (2, "shared_edge"), (3, "shared_face")),
)
def test_shared_membership_classifies_polyhedron_connection(shared_count, kind):
    graph = StructuralGraphBuilder().build(units_with_shared_atoms(shared_count), ())
    assert graph.connections[0].connection_kind.value == kind


def test_periodic_reverse_direct_contacts_collapse_to_one_connection():
    graph = StructuralGraphBuilder().build(atomic_units, (forward_contact, reverse_contact))
    assert len(graph.connections) == 1
    assert graph.connections[0].lattice_translation == (1, 0, 0)
```

- [ ] **Step 2: Run the graph tests and verify RED**

Run: `python3 -m pytest -p no:cacheprovider tests/crystal_chemistry/test_structural_graph.py -q`

Expected: import failure because graph contracts do not exist.

- [ ] **Step 3: Implement the finite quotient graph**

```python
class StructuralConnectionKind(StrEnum):
    SHARED_VERTEX = "shared_vertex"
    SHARED_EDGE = "shared_edge"
    SHARED_FACE = "shared_face"
    DIRECT_CONTACT = "direct_contact"


@dataclass(frozen=True, slots=True)
class StructuralConnection:
    connection_id: str
    first_unit_id: str
    second_unit_id: str
    lattice_translation: tuple[int, int, int]
    connection_kind: StructuralConnectionKind
    shared_atom_refs: tuple[PeriodicAtomRef, ...]
    source_contact_ids: tuple[str, ...]
    interaction_layers: tuple[InteractionLayer, ...]
    contact_classifications: tuple[ContactClassification, ...]
    provenance: tuple[tuple[str, object], ...] = ()


@dataclass(frozen=True, slots=True)
class StructuralUnitGraph:
    units: tuple[StructuralUnit, ...]
    connections: tuple[StructuralConnection, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()
```

Canonicalize `(A, B, t)` against `(B, A, -t)` and retain one physical relation. Determine shared vertex/edge/face only from canonical membership. Add `DIRECT_CONTACT` only from a supplied `ResolvedContact`; raw distances are not accepted by this API. Preserve every source contact ID, interaction layer, and geometric classification needed by later selection.

- [ ] **Step 4: Run graph tests and verify GREEN**

Run: `python3 -m pytest -p no:cacheprovider tests/crystal_chemistry/test_structural_graph.py -q`

Expected: all tests pass.

- [ ] **Step 5: Commit graph construction**

```bash
git add src/cristma/crystal_chemistry/structural_graph.py src/cristma/crystal_chemistry/__init__.py tests/crystal_chemistry/test_structural_graph.py
git commit -m "Build structural unit graph"
```

### Task 3: Scientific-fixture acceptance and documentation

**Files:**
- Modify: `tests/integration/test_inorganic_crystal_chemistry.py`
- Modify: `docs/inorganic-crystal-chemistry.md`

**Interfaces:**
- Consumes: existing CIF → Chemistry → resolver → polyhedron pipeline plus Tasks 1–2.
- Produces: verified end-to-end unit graphs for CaMoO4, LiB3O5, and FeS2.

- [ ] **Step 1: Write failing end-to-end assertions**

```python
def test_camo4_structural_graph_keeps_interstitial_contacts_semantic():
    graph = build_unit_graph("CaMoO4_9009632.cif")
    assert any(unit.kind is StructuralUnitKind.POLYHEDRON for unit in graph.units)
    assert {layer for edge in graph.connections for layer in edge.interaction_layers} >= {
        InteractionLayer.STRUCTURAL,
        InteractionLayer.INTERSTITIAL,
    }


def test_fes2_graph_retains_coordination_and_sulfur_subsystem_contacts():
    graph = build_unit_graph("FeS2_9000594.cif")
    assert {layer for edge in graph.connections for layer in edge.interaction_layers} >= {
        InteractionLayer.COORDINATION,
        InteractionLayer.INTRA_SUBSYSTEM,
    }
```

- [ ] **Step 2: Run acceptance tests and verify RED or expose missing mapping**

Run: `python3 -m pytest -p no:cacheprovider tests/integration/test_inorganic_crystal_chemistry.py -q`

Expected: the new graph assertions fail until all source-contact mappings are retained.

- [ ] **Step 3: Complete only missing mappings and document the API boundary**

Document:

```text
CrystalChemistryResolution + CoordinationPolyhedron
    -> StructuralUnitBuilder
    -> StructuralGraphBuilder
    -> StructuralUnitGraph
```

State explicitly that the graph is still unclassified: no rank, block, chain, layer, framework, or motif claim is made in this slice.

- [ ] **Step 4: Run full verification and build the wheel**

Run:

```bash
python3 -m pytest -p no:cacheprovider -q
python3 -m pip wheel . --no-deps --no-build-isolation --wheel-dir /tmp/cristma-structural-graph-dist
unzip -l /tmp/cristma-structural-graph-dist/cristma-0.1.0.dev0-py3-none-any.whl | rg structural_
git diff --check
```

Expected: zero test failures, wheel build exit code 0, both new modules packaged, and no whitespace errors.

- [ ] **Step 5: Commit the completed slice**

```bash
git add tests/integration/test_inorganic_crystal_chemistry.py docs/inorganic-crystal-chemistry.md
git commit -m "Verify structural unit graph pipeline"
```
