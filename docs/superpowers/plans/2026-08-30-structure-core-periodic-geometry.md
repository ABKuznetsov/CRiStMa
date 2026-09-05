# Structure Core Periodic Geometry Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement CrIStMa's immutable symmetry-expanded atomic identity, generic atomic views, finite and periodic neighbor graphs, and basic coordination analysis.

**Architecture:** Canonical structures retain independent sites only. Symmetry expansion creates a finite `AtomicView[ExpandedAtom]` with complete symmetry-image provenance; periodicity is introduced only through `PeriodicAtomRef` and periodic graph edges. Molecular structures use the same generic `AtomicView` contract while returning a distinct finite `NeighborGraph`.

**Tech Stack:** Python 3.11+, frozen slotted dataclasses, `typing.Protocol`/`Generic`, NumPy, pytest, setuptools/build.

**Spec:** `docs/superpowers/specs/2026-08-30-structure-core-periodic-identity-design.md`

## Global Constraints

- Inputs and result objects are immutable; NumPy arrays exposed by public data objects are read-only.
- `ExpandedAtom` is one unique symmetry-expanded position in the reference cell, never a periodic copy.
- Mixed occupancy remains one geometric position with multiple `SiteComponent` values.
- `SymmetryImageProvenance.normalization_translation` and `PeriodicAtomRef.cell_translation` are distinct concepts with distinct types.
- Anisotropic displacement tensors are transformed under symmetry and checked for special-position consistency.
- Invalid direct object construction raises; invalid source data becomes a namespaced reader diagnostic and does not create a partial canonical structure.
- Periodic cutoff search is complete for arbitrary valid triclinic cells and all supported partial-periodicity masks.
- This plan does not implement the SHELX parser. SHELX is phase 2 of `docs/superpowers/plans/2026-08-30-native-structure-readers-roadmap.md` and must later traverse the same public Structure Core path.
- No PyPI upload or public release is performed. The package version remains a development version until Viewer and Finder integrations are complete.
- Run targeted tests after every change; run the full suite once at the final integration task.

## File map

**New scientific modules**

- `src/cristma/diagnostics.py`: shared diagnostic severity, source location, and diagnostic value types.
- `src/cristma/structure/occupation.py`: shared chemical occupation components.
- `src/cristma/structure/position.py`: minimal `AtomicPosition` capability protocol.
- `src/cristma/symmetry/displacement.py`: symmetry transformation and consistency comparison for ADPs.
- `src/cristma/geometry/neighbors.py`: finite/periodic neighbor references, edges, graphs, and shared protocols.
- `src/cristma/geometry/finder.py`: exact cutoff neighbor enumeration.
- `src/cristma/geometry/coordination.py`: coordination result types and analyzer.
- `src/cristma/geometry/__init__.py`: stable geometry exports.

**Existing modules changed**

- `src/cristma/structure/crystal.py`: occupancy invariants, `ExpandedAtom` integration, and crystal atomic views.
- `src/cristma/structure/identity.py`: symmetry-image provenance, expanded identity, and periodic references.
- `src/cristma/structure/molecular.py`: component-based molecular positions and generic views.
- `src/cristma/structure/view.py`: generic `AtomicView[TAtom]` using `UnitCell`.
- `src/cristma/structure/__init__.py`: public structure exports.
- `src/cristma/symmetry/orbit.py`: deterministic finite expansion and provenance merging.
- `src/cristma/symmetry/__init__.py`: public symmetry exports.
- `src/cristma/core/cell.py`: read-only matrix contract.
- `src/cristma/io/cif/mapper.py`: source diagnostics and derived-expansion handling.
- `src/cristma/__init__.py`, `pyproject.toml`: development version and public imports.

**New tests**

- `tests/structure/test_position.py`
- `tests/symmetry/test_displacement.py`
- `tests/geometry/test_finite_neighbors.py`
- `tests/geometry/test_periodic_neighbors.py`
- `tests/geometry/test_coordination.py`
- `tests/integration/test_structure_core.py`

---

### Task 1: Enforce occupation invariants and development-version status

**Files:**
- Create: `src/cristma/structure/occupation.py`
- Modify: `src/cristma/structure/crystal.py`
- Modify: `src/cristma/structure/molecular.py`
- Modify: `src/cristma/structure/__init__.py`
- Modify: `src/cristma/io/cif/mapper.py`
- Modify: `tests/structure/test_crystal.py`
- Modify: `tests/io/cif/test_mapper_advanced.py`
- Modify: `src/cristma/__init__.py`
- Modify: `pyproject.toml`
- Modify: `tests/test_public_api.py`

**Interfaces:**
- Consumes: existing `MeasuredValue`, `SiteComponent`, `IndependentSite`, `Diagnostic`, and `ReadResult` contracts.
- Produces: shared `SiteComponent`, strict component/total occupancy validation, `IndependentSite.total_occupancy`, `IndependentSite.vacancy_fraction`, diagnostic codes `cif.map.occupancy_out_of_range` and `cif.map.occupancy_total_exceeds_one`, and version `0.1.0.dev0`.

- [ ] **Step 1: Write failing model tests**

Add to `tests/structure/test_crystal.py`:

```python
def test_component_occupancy_above_one_is_invalid() -> None:
    with pytest.raises(ValueError, match="occupancy must lie between zero and one"):
        SiteComponent("Ca", number(1.01))


def test_site_reports_total_occupancy_and_vacancy_fraction() -> None:
    site = IndependentSite(
        id="site:M1",
        label="M1",
        components=(
            SiteComponent("Ca", number(0.6)),
            SiteComponent("Sr", number(0.2)),
        ),
        fractional=(number(0), number(0), number(0)),
    )

    assert site.total_occupancy == pytest.approx(0.8)
    assert site.vacancy_fraction == pytest.approx(0.2)
```

Add an import-contract test proving `SiteComponent` is defined by
`cristma.structure.occupation` and re-exported unchanged by
`cristma.structure`.

- [ ] **Step 2: Run the model tests and verify failure**

Run:

```bash
pytest tests/structure/test_crystal.py -q
```

Expected: the occupancy-above-one test does not raise and the two properties are absent.

- [ ] **Step 3: Implement strict occupation validation and derived properties**

In `src/cristma/structure/crystal.py`, validate every component in `[0, 1]`, retain `math.fsum` for site totals, reject totals above `1.0 + 1e-12`, and add:

```python
@property
def total_occupancy(self) -> float:
    return math.fsum(float(item.occupancy.value) for item in self.components)

@property
def vacancy_fraction(self) -> float:
    return max(0.0, 1.0 - self.total_occupancy)
```

- [ ] **Step 4: Add a failing CIF source-diagnostic test**

Add to `tests/io/cif/test_mapper_advanced.py` a complete one-site CIF with `_atom_site_occupancy 1.2`, then assert:

```python
structures, diagnostics = map_cif_structures(parse_cif(source).document)

assert not structures
assert "cif.map.occupancy_out_of_range" in {
    item.code for item in diagnostics
}
```

Add a source with occupancy `-0.2` and assert the same
`cif.map.occupancy_out_of_range` code. Add a third source with two coincident rows in the same explicit disorder
assembly/group and occupancies `0.7` and `0.6`. Assert that the block also
produces no structure and emits `cif.map.occupancy_total_exceeds_one`.

- [ ] **Step 5: Run the reader test and verify failure**

Run:

```bash
pytest tests/io/cif/test_mapper_advanced.py -q
```

Expected: FAIL because the mapper currently emits only `cif.map.site_invalid`.

- [ ] **Step 6: Map invalid source occupancy to the specific diagnostic**

Move `SiteComponent` unchanged from `crystal.py` to
`structure/occupation.py`; import it from both crystal and molecular modules and
re-export it from `cristma.structure`.

In `_sites()` in `src/cristma/io/cif/mapper.py`, detect a parsed occupancy outside `[0, 1]`, append an error `Diagnostic` with code `cif.map.occupancy_out_of_range`, mark the block failed, and skip canonical site construction. Preserve the token span and reported value in the message.

Change `_merge_coincident_sites()` to return `None` when rows explicitly marked
as one disorder assembly/group have a combined occupancy above one. Emit
`cif.map.occupancy_total_exceeds_one` and make `_sites()` fail the complete CIF
block rather than returning separate coincident canonical sites. Preserve the
existing warning behavior for coincident full sites that are not declared as
one disorder model.

- [ ] **Step 7: Mark the package as an internal development build**

Set both version declarations to the same PEP 440 value:

```toml
# pyproject.toml
version = "0.1.0.dev0"
```

```python
# src/cristma/__init__.py
__version__ = "0.1.0.dev0"
```

Update `tests/test_public_api.py` to assert that exact value.

- [ ] **Step 8: Run targeted tests**

Run:

```bash
pytest tests/structure/test_crystal.py tests/io/cif/test_mapper_advanced.py tests/test_public_api.py -q
```

Expected: PASS.

- [ ] **Step 9: Commit**

```bash
git add pyproject.toml src/cristma/__init__.py src/cristma/structure/occupation.py src/cristma/structure/crystal.py src/cristma/structure/molecular.py src/cristma/structure/__init__.py src/cristma/io/cif/mapper.py tests/structure/test_crystal.py tests/io/cif/test_mapper_advanced.py tests/test_public_api.py
git commit -m "feat: enforce structure occupation invariants"
```

---

### Task 2: Replace ExpandedAtomRef with complete symmetry-expanded identity

**Files:**
- Modify: `src/cristma/structure/identity.py`
- Modify: `src/cristma/symmetry/orbit.py`
- Create: `src/cristma/symmetry/displacement.py`
- Modify: `src/cristma/structure/crystal.py`
- Modify: `src/cristma/io/cif/mapper.py`
- Modify: `src/cristma/structure/__init__.py`
- Modify: `src/cristma/symmetry/__init__.py`
- Modify: `tests/symmetry/test_orbit.py`
- Create: `tests/symmetry/test_displacement.py`
- Modify: `tests/io/cif/test_mapper_advanced.py`

**Interfaces:**
- Consumes: `IndependentSite`, `SiteComponent`, `DisplacementParameters`, `UnitCell`, `AffineOperation`, and `MeasuredValue`.
- Produces: `SymmetryImageProvenance`, `ExpandedAtom`, `SymmetryConsistencyError`, `transform_displacement()`, and `expand_orbit(site, operations, *, cell, structure_id=None, tolerance=1e-8)`.

- [ ] **Step 1: Replace orbit tests with the approved provenance contract**

In `tests/symmetry/test_orbit.py`, pass `cell=UnitCell.cubic(number(4.0))` to every expansion and assert the normalization sign convention:

```python
def test_orbit_records_operation_and_normalization_translation() -> None:
    site = site_at(1.12, -0.04, 0.35)
    operation = parse_xyz_operation("x,y,z", operation_id="op:identity")

    atom = expand_orbit(
        site,
        (operation,),
        cell=UnitCell.cubic(number(4.0)),
        structure_id="structure:test",
    )[0]

    assert atom.fractional == pytest.approx((0.12, 0.96, 0.35))
    assert atom.representative_image == SymmetryImageProvenance(
        "op:identity", (-1, 1, 0)
    )
    assert atom.representative_image in atom.equivalent_images
```

Update the special-position test to assert two `SymmetryImageProvenance` values rather than operation-ID strings. Add a mixed-site test asserting one `ExpandedAtom` contains both original components.

Add a split-position test with two `IndependentSite` objects at distinct
fractional coordinates and complementary occupancies; assert expansion returns
two `ExpandedAtom` objects even when their disorder assembly/group labels
match.

- [ ] **Step 2: Run orbit tests and verify failure**

Run:

```bash
pytest tests/symmetry/test_orbit.py -q
```

Expected: FAIL because `SymmetryImageProvenance`, `ExpandedAtom`, and the required cell argument do not exist.

- [ ] **Step 3: Implement symmetry-image and expanded-atom records**

Replace `ExpandedAtomRef` in `src/cristma/structure/identity.py` with frozen slotted records equivalent to:

```python
@dataclass(frozen=True, slots=True)
class SymmetryImageProvenance:
    operation_id: str
    normalization_translation: tuple[int, int, int]


@dataclass(frozen=True, slots=True)
class ExpandedAtom:
    id: str
    structure_id: str | None
    source_site_id: str
    fractional: tuple[float, float, float]
    cartesian: tuple[float, float, float]
    components: tuple[SiteComponent, ...]
    displacement: DisplacementParameters | None
    representative_image: SymmetryImageProvenance
    equivalent_images: tuple[SymmetryImageProvenance, ...]
```

Use `TYPE_CHECKING` imports or move only the annotations needed to avoid a runtime cycle between `identity.py` and `crystal.py`. Do not retain an `ExpandedAtomRef` compatibility alias because no public release exists.

- [ ] **Step 4: Implement deterministic orbit identity and normalization provenance**

In `src/cristma/symmetry/orbit.py`:

- define normalization so `canonical = raw + normalization_translation`;
- merge periodically equal positions using the declared tolerance;
- construct Cartesian coordinates as `np.asarray(fractional) @ cell.matrix`;
- retain the complete ordered tuple of generating images;
- derive IDs from `structure_id`, `source_site_id`, and the tolerance-quantized canonical position, never from a representative operation ID.

Use the quantized key:

```python
position_key = tuple(round(value / tolerance) for value in fractional)
atom_id = (
    f"expanded:{structure_id or 'unassigned'}:{site.id}:"
    + ",".join(str(value) for value in position_key)
)
```

- [ ] **Step 5: Write failing anisotropic ADP tests**

Create `tests/symmetry/test_displacement.py` with:

```python
def test_anisotropic_u_tensor_rotates_with_symmetry_operation() -> None:
    displacement = anisotropic_u(((1.0, 0.2, 0.0), (0.2, 2.0, 0.0), (0.0, 0.0, 3.0)))
    operation = parse_xyz_operation("-y,x,z", operation_id="op:rotate")

    transformed = transform_displacement(displacement, operation)

    rotation = np.asarray(operation.rotation, dtype=float)
    expected = rotation @ numeric_tensor(displacement) @ rotation.T
    assert np.allclose(numeric_tensor(transformed), expected)
    assert transformed.reported_kind == displacement.reported_kind


def test_special_position_rejects_inconsistent_equivalent_adp_images() -> None:
    site = site_with_anisotropic_u_at_origin(diagonal=(1.0, 2.0, 3.0))
    operations = (
        parse_xyz_operation("x,y,z", operation_id="op:1"),
        parse_xyz_operation("y,x,z", operation_id="op:2"),
    )

    with pytest.raises(SymmetryConsistencyError, match="anisotropic displacement"):
        expand_orbit(site, operations, cell=cubic_cell(4.0))
```

The local fixture helpers construct `MeasuredValue` tensors and return numeric arrays explicitly; do not compare raw CIF strings.

- [ ] **Step 6: Run displacement tests and verify failure**

Run:

```bash
pytest tests/symmetry/test_displacement.py -q
```

Expected: FAIL because the transformation module and consistency error are absent.

- [ ] **Step 7: Implement ADP transformation and consistency checking**

In `src/cristma/symmetry/displacement.py`:

- return `None` unchanged;
- preserve `U_iso` and `B_iso` values;
- for `U_aniso`, calculate `R @ U @ R.T` using the exact operation rotation converted to float;
- create derived `MeasuredValue` entries with `raw=None`;
- treat `(U11, U22, U33, U12, U13, U23)` as six independent reported values;
- build the exact `6 x 6` linear transform by applying `R @ basis @ R.T` to the six symmetric tensor basis matrices;
- calculate transformed values with that matrix and propagate independent standard uncertainties as `sqrt(sum((coefficient * sigma) ** 2))`; use `None` for an output uncertainty when any contributing nonzero coefficient refers to a missing input uncertainty;
- preserve `kind` and `reported_kind`;
- compare transformed tensors with `np.allclose(..., rtol=0.0, atol=tolerance)`;
- raise `SymmetryConsistencyError` if equivalent images of a merged position disagree.

In `expand_orbit()`, transform ADP for every generating image and attach the agreed result to `ExpandedAtom`.

- [ ] **Step 8: Convert expansion failures from CIF sources into diagnostics**

In `src/cristma/io/cif/mapper.py`, call `expand_orbit(..., cell=cell)` and catch `SymmetryConsistencyError`. Append an error diagnostic with code `cif.map.adp_symmetry_inconsistent`, omit the affected CIF block from `structures`, and retain the parsed source document in `ReadResult`.

Add a CIF mapper test with a special position and inconsistent anisotropic tensor, asserting no canonical structure and the exact diagnostic code.

- [ ] **Step 9: Run targeted symmetry and mapper tests**

Run:

```bash
pytest tests/symmetry/test_orbit.py tests/symmetry/test_displacement.py tests/io/cif/test_mapper_advanced.py -q
```

Expected: PASS.

- [ ] **Step 10: Commit**

```bash
git add src/cristma/structure/identity.py src/cristma/structure/crystal.py src/cristma/structure/__init__.py src/cristma/symmetry/orbit.py src/cristma/symmetry/displacement.py src/cristma/symmetry/__init__.py src/cristma/io/cif/mapper.py tests/symmetry/test_orbit.py tests/symmetry/test_displacement.py tests/io/cif/test_mapper_advanced.py
git commit -m "feat: model symmetry-expanded atomic identity"
```

---

### Task 3: Make AtomicView generic across crystals and molecules

**Files:**
- Create: `src/cristma/structure/position.py`
- Modify: `src/cristma/structure/view.py`
- Modify: `src/cristma/structure/crystal.py`
- Modify: `src/cristma/structure/molecular.py`
- Modify: `src/cristma/core/cell.py`
- Modify: `src/cristma/io/cif/mapper.py`
- Modify: `src/cristma/structure/__init__.py`
- Create: `tests/structure/test_position.py`
- Modify: `tests/structure/test_atomic_view.py`
- Modify: `tests/structure/test_molecular.py`
- Modify: `tests/io/cif/test_end_to_end.py`

**Interfaces:**
- Consumes: `ExpandedAtom`, `SiteComponent`, `MolecularAtom`, `UnitCell`, `AtomicPropertyTable`, and `expand_orbit()` from Task 2.
- Produces: `AtomicPosition` protocol, generic `AtomicView[TAtom]`, `expand_structure(crystal, tolerance=1e-8)`, and component-based `MolecularAtom`.

- [ ] **Step 1: Write failing protocol and generic-view tests**

Create `tests/structure/test_position.py`:

```python
def test_crystal_view_retains_expanded_atom_objects() -> None:
    crystal = one_site_p1_crystal()

    view = expand_structure(crystal)

    assert isinstance(view.atoms[0], ExpandedAtom)
    assert view.atoms[0].source_site_id == crystal.sites[0].id
    assert view.cell is crystal.cell
    assert view.cell_matrix.flags.writeable is False


def test_molecular_view_retains_molecular_atom_objects() -> None:
    atom = MolecularAtom(
        id="atom:C1",
        label="C1",
        components=(SiteComponent("C", number(1.0)),),
        cartesian=(1.0, 2.0, 3.0),
    )

    view = MolecularStructure("methane-fragment", atoms=(atom,)).atomic_view()

    assert view.atoms == (atom,)
    assert view.fractional is None
    assert view.cell is None
    assert view.periodic == (False, False, False)
```

Add a validation test asserting that any periodic axis with `cell=None` or `fractional=None` raises `ValueError`.

- [ ] **Step 2: Run view tests and verify failure**

Run:

```bash
pytest tests/structure/test_position.py tests/structure/test_atomic_view.py tests/structure/test_molecular.py -q
```

Expected: FAIL because views expose parallel IDs/species arrays instead of concrete atomic positions.

- [ ] **Step 3: Implement AtomicPosition and generic AtomicView**

Create `src/cristma/structure/position.py`:

```python
class AtomicPosition(Protocol):
    id: str
    cartesian: tuple[float, float, float]
    components: tuple[SiteComponent, ...]
```

In `src/cristma/structure/view.py`, define `TAtom = TypeVar("TAtom", bound=AtomicPosition)` and make `AtomicView(Generic[TAtom])` hold:

```python
atoms: tuple[TAtom, ...]
cell: UnitCell | None
periodic: tuple[bool, bool, bool]
properties: AtomicPropertyTable
```

In `__post_init__`, derive the read-only Cartesian array from
`atom.cartesian`. Derive the fractional array only when every atom provides a
fractional position; otherwise set it to `None`. Validate uniqueness of atom
IDs, finite coordinates, and the cell/fractional requirement for periodic
axes. Remove the parallel constructor inputs `ids`, `species`,
`source_site_ids`, `cartesian`, and `fractional`. Provide read-only properties
only when they are mathematically unambiguous:

```python
@property
def ids(self) -> tuple[str, ...]:
    return tuple(atom.id for atom in self.atoms)

@property
def cell_matrix(self) -> np.ndarray | None:
    return None if self.cell is None else self.cell.matrix
```

- [ ] **Step 4: Make UnitCell matrices read-only**

In `UnitCell.matrix`, set `matrix.flags.writeable = False` before returning it. Add a test that mutation raises `ValueError` and that `metric` and `volume` remain numerically unchanged.

- [ ] **Step 5: Convert MolecularAtom to the shared component contract**

Replace `species` and scalar `occupancy` storage in `MolecularAtom` with:

```python
components: tuple[SiteComponent, ...]
```

Require at least one component and a total occupancy no greater than one. Update molecular fixtures to construct a single full-occupancy component. `MolecularStructure.atomic_view()` returns `AtomicView[MolecularAtom]` with its original atoms and Cartesian array.

- [ ] **Step 6: Build crystal views from finite expansion results**

Add `expand_structure()` to `src/cristma/symmetry/orbit.py`. It expands every site with the crystal's cell and operations, concatenates atoms, rejects duplicate expanded IDs, and returns `AtomicView[ExpandedAtom]`.

Change `CrystalStructure.atomic_view(expanded=True)` to delegate to `expand_structure(self)`. Remove `CrystalStructure.expanded_sites`; it is derived state and must not live in the canonical snapshot. In the CIF mapper, retain calculated multiplicity but stop storing an expanded-atom cache on the structure.

- [ ] **Step 7: Update existing view and CIF end-to-end assertions**

Replace direct `structure.expanded_sites` assertions with:

```python
view = expand_structure(structure)
assert len(view.atoms) == expected_multiplicity
```

For mixed occupancy, assert `view.atoms[0].components` contains both species and no `UnknownSpecies` placeholder is created.

For a site with total occupancy `0.8`, assert the crystal view still contains
exactly one geometric atom and `vacancy_fraction == pytest.approx(0.2)` on its
source occupation model.

- [ ] **Step 8: Run targeted structure and CIF tests**

Run:

```bash
pytest tests/structure tests/io/cif/test_end_to_end.py tests/io/cif/test_mapper_basic.py tests/io/cif/test_mapper_advanced.py -q
```

Expected: PASS.

- [ ] **Step 9: Commit**

```bash
git add src/cristma/core/cell.py src/cristma/structure src/cristma/symmetry/orbit.py src/cristma/io/cif/mapper.py tests/structure tests/io/cif/test_end_to_end.py tests/io/cif/test_mapper_basic.py tests/io/cif/test_mapper_advanced.py
git commit -m "feat: add generic atomic views"
```

---

### Task 4: Add finite molecular neighbor graphs

**Files:**
- Create: `src/cristma/diagnostics.py`
- Modify: `src/cristma/io/diagnostics.py`
- Modify: `tests/io/test_diagnostics.py`
- Create: `src/cristma/geometry/__init__.py`
- Create: `src/cristma/geometry/neighbors.py`
- Create: `src/cristma/geometry/finder.py`
- Create: `tests/geometry/test_finite_neighbors.py`

**Interfaces:**
- Consumes: `AtomicPosition`, `AtomicView[TAtom]`, and read-only Cartesian arrays.
- Produces: `Neighbor`, `NeighborGraph[TAtom]`, `NeighborLike`, `NeighborGraphLike[TAtom, TNeighbor]`, and configurable `NeighborFinder(cutoff, tolerance=1e-12)`.

- [ ] **Step 1: Write failing finite-neighbor tests**

Create `tests/geometry/test_finite_neighbors.py`:

```python
def test_finite_cutoff_graph_contains_two_directed_edges() -> None:
    view = molecular_view(
        ("A", (0.0, 0.0, 0.0)),
        ("B", (1.0, 0.0, 0.0)),
        ("C", (3.0, 0.0, 0.0)),
    )

    graph = NeighborFinder(cutoff=1.1).find(view)

    assert isinstance(graph, NeighborGraph)
    assert [(edge.target_atom_id, edge.distance) for edge in graph.neighbors("A")] == [
        ("B", pytest.approx(1.0))
    ]
    assert graph.neighbors("B")[0].vector_cartesian == pytest.approx((-1.0, 0.0, 0.0))
    assert graph.neighbors("C") == ()


def test_finite_graph_rejects_nonpositive_cutoff() -> None:
    with pytest.raises(ValueError, match="cutoff must be positive"):
        NeighborFinder(cutoff=0.0)


def test_finder_configuration_is_inspectable_and_clone_is_immutable() -> None:
    finder = NeighborFinder(cutoff=3.0, tolerance=1e-12)

    clone = finder.clone(cutoff=2.5)

    assert finder.get_config() == {"cutoff": 3.0, "tolerance": 1e-12}
    assert clone.get_config() == {"cutoff": 2.5, "tolerance": 1e-12}
    assert finder.cutoff == 3.0
```

Also test deterministic neighbor ordering by `(distance, target_atom_id)` and
no self-edge at zero distance. Add two distinct molecular positions at the same
coordinates and assert that no zero-length edge is created and the graph
contains diagnostic code `geometry.coincident_positions`.

- [ ] **Step 2: Run the finite tests and verify failure**

Run:

```bash
pytest tests/geometry/test_finite_neighbors.py -q
```

Expected: FAIL because `cristma.geometry` does not exist.

- [ ] **Step 3: Promote diagnostics to a toolbox-neutral module**

Move the existing `Severity`, `SourcePosition`, `SourceSpan`, and `Diagnostic`
definitions to `src/cristma/diagnostics.py`. Make
`src/cristma/io/diagnostics.py` import and re-export those same class objects so
existing imports and `isinstance` behavior remain valid. Extend
`tests/io/test_diagnostics.py` with:

```python
from cristma.diagnostics import Diagnostic as SharedDiagnostic
from cristma.io import Diagnostic as IoDiagnostic


def test_io_reexports_shared_diagnostic_type() -> None:
    assert IoDiagnostic is SharedDiagnostic
```

- [ ] **Step 4: Implement immutable finite graph records**

In `src/cristma/geometry/neighbors.py`, define frozen slotted generic graph records. `NeighborGraph` validates unique node IDs, known edge endpoints, finite positive distances, and row ownership. Store adjacency behind immutable tuples and expose:

```python
def neighbors(self, atom_id: str) -> tuple[Neighbor, ...]: ...
```

Define `NeighborLike` with `source_atom_id`, `target_atom_id`, `distance`, and `vector_cartesian`. Define a two-parameter `NeighborGraphLike[TAtom, TNeighbor]` protocol with `atoms` and `neighbors()`. Graph diagnostics use `cristma.diagnostics.Diagnostic`, never an I/O-owned type.

- [ ] **Step 5: Implement exact finite cutoff enumeration**

In `src/cristma/geometry/finder.py`, validate positive finite `cutoff` and
`tolerance`. Implement `get_config()` and `clone(**changes)` directly on the
frozen tool-class without introducing a superclass. Branch on
`any(view.periodic)` and initially implement the non-periodic branch with all
unordered atom pairs. Emit two directed `Neighbor` edges per accepted pair,
with opposite vectors. Sort each adjacency tuple by `(distance,
target_atom_id)`. For distinct IDs separated by at most `self.tolerance`, emit
one `geometry.coincident_positions` warning and no edge.

- [ ] **Step 6: Run finite geometry tests**

Run:

```bash
pytest tests/io/test_diagnostics.py tests/geometry/test_finite_neighbors.py -q
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add src/cristma/diagnostics.py src/cristma/io/diagnostics.py src/cristma/geometry tests/io/test_diagnostics.py tests/geometry/test_finite_neighbors.py
git commit -m "feat: add finite neighbor graphs"
```

---

### Task 5: Add complete periodic cutoff enumeration

**Files:**
- Modify: `src/cristma/structure/identity.py`
- Modify: `src/cristma/structure/__init__.py`
- Modify: `src/cristma/geometry/neighbors.py`
- Modify: `src/cristma/geometry/finder.py`
- Modify: `src/cristma/geometry/__init__.py`
- Create: `tests/geometry/test_periodic_neighbors.py`

**Interfaces:**
- Consumes: `AtomicView[ExpandedAtom]`, `UnitCell.matrix`, and finite graph interfaces from Task 4.
- Produces: `PeriodicAtomRef`, `PeriodicNeighbor`, `PeriodicNeighborGraph[TAtom]`, and the periodic overload of `NeighborFinder.find()`.

Add concrete overloads for
`AtomicView[MolecularAtom] -> NeighborGraph[MolecularAtom]` and
`AtomicView[ExpandedAtom] -> PeriodicNeighborGraph[ExpandedAtom]`, followed by
the union implementation signature. Runtime dispatch remains based on
`view.periodic`.

- [ ] **Step 1: Write failing periodic identity and boundary tests**

Create `tests/geometry/test_periodic_neighbors.py` with a P1 cubic structure containing atoms at fractional `(0.05, 0, 0)` and `(0.95, 0, 0)` in a 10 Å cell:

```python
def test_periodic_graph_finds_boundary_crossing_image_once() -> None:
    view = expand_structure(boundary_crossing_crystal())

    graph = NeighborFinder(cutoff=1.1).find(view)

    assert isinstance(graph, PeriodicNeighborGraph)
    edge = graph.neighbors(view.atoms[0].id)[0]
    assert edge.target.atom_id == view.atoms[1].id
    assert edge.target.cell_translation == (-1, 0, 0)
    assert edge.distance == pytest.approx(1.0)
    assert edge.vector_cartesian == pytest.approx((-1.0, 0.0, 0.0))
```

Add tests for a periodic self-neighbor with nonzero translation, exact reverse-edge translation/vector, and no zero-distance `(same atom, (0,0,0))` edge.

- [ ] **Step 2: Write the skewed triclinic completeness test**

Construct a valid triclinic `UnitCell` with `a=4`, `b=4`, `c=4`, `alpha=90`, `beta=90`, `gamma=20`, place a one-site P1 atom, and use a `3.0` Å cutoff. The image at translation `(2, -2, 0)` is about `2.78` Å away even though a naive `ceil(cutoff / edge) == 1` bound would miss it. Compute the complete expected image set by exhaustive enumeration over `range(-4, 5)` in the test only, then assert the production graph returns that exact set.

The oracle code is explicit:

```python
expected = {
    translation
    for translation in itertools.product(range(-4, 5), repeat=3)
    if translation != (0, 0, 0)
    and np.linalg.norm(np.asarray(translation) @ cell.matrix) <= cutoff + 1e-12
}
actual = {
    edge.target.cell_translation
    for edge in graph.neighbors(atom.id)
}
assert actual == expected
```

- [ ] **Step 3: Run periodic tests and verify failure**

Run:

```bash
pytest tests/geometry/test_periodic_neighbors.py -q
```

Expected: FAIL because periodic reference and graph types are absent.

- [ ] **Step 4: Implement periodic reference and edge records**

Add frozen slotted `PeriodicAtomRef(atom_id, cell_translation)` to `src/cristma/structure/identity.py`. Add `PeriodicNeighbor` and generic `PeriodicNeighborGraph` to `src/cristma/geometry/neighbors.py`. `PeriodicNeighbor.target_atom_id` is a property delegating to `target.atom_id`, allowing it to satisfy `NeighborLike` without erasing periodic data.

- [ ] **Step 5: Implement reciprocal/metric translation bounds**

For each source-target fractional difference `delta`, use row-basis matrix `A = cell.matrix` and `inverse = np.linalg.inv(A)`. For Cartesian vectors with norm at most `cutoff`, each fractional component obeys:

```python
component_bound = cutoff * np.linalg.norm(inverse[:, axis])
lower = math.ceil(-delta[axis] - component_bound - tolerance)
upper = math.floor(-delta[axis] + component_bound + tolerance)
```

Enumerate `range(lower, upper + 1)` on periodic axes and `(0,)` on non-periodic axes. For every candidate, calculate:

```python
vector_fractional = delta + np.asarray(translation, dtype=float)
vector_cartesian = vector_fractional @ cell.matrix
distance = np.linalg.norm(vector_cartesian)
```

Accept `tolerance < distance <= cutoff + tolerance`, key edges by `(source_id, target_id, translation)`, and sort deterministically. Exclude the trivial `(same atom, (0, 0, 0))` candidate without warning. For any other candidate at or below tolerance, emit one `geometry.coincident_positions` warning and no edge. This derivation is the completeness guarantee; do not replace it with bounds based only on `a`, `b`, and `c`.

- [ ] **Step 6: Add partial-periodicity coverage**

Add a test with `periodic=(True, True, False)` asserting every returned translation has `translation[2] == 0` while boundary-crossing neighbors along `a` and `b` remain present. Implement the per-axis mask in the enumerator.

- [ ] **Step 7: Run all neighbor tests**

Run:

```bash
pytest tests/geometry/test_finite_neighbors.py tests/geometry/test_periodic_neighbors.py -q
```

Expected: PASS.

- [ ] **Step 8: Commit**

```bash
git add src/cristma/structure/identity.py src/cristma/structure/__init__.py src/cristma/geometry tests/geometry/test_periodic_neighbors.py
git commit -m "feat: add periodic neighbor geometry"
```

---

### Task 6: Add graph-neutral coordination analysis

**Files:**
- Create: `src/cristma/geometry/coordination.py`
- Modify: `src/cristma/geometry/__init__.py`
- Create: `tests/geometry/test_coordination.py`

**Interfaces:**
- Consumes: `AtomicView[TAtom]`, `NeighborGraphLike`, `Neighbor`, and `PeriodicNeighbor`.
- Produces: `CoordinationEnvironment`, `CoordinationResult`, and stateless configurable `CoordinationAnalyzer`.

- [ ] **Step 1: Write failing coordination tests for both graph kinds**

Create `tests/geometry/test_coordination.py`:

```python
@pytest.mark.parametrize("case", ["molecule", "crystal"])
def test_coordination_number_comes_from_geometric_neighbors(case: str) -> None:
    view, graph, center_id, expected = coordination_case(case)

    result = CoordinationAnalyzer().analyze(view, graph)
    environment = result.by_atom(center_id)

    assert environment.center_atom_id == center_id
    assert environment.coordination_number == expected
    assert len(environment.neighbors) == expected
```

Add a mixed Ca/Sr center test asserting one environment is returned, its coordination number is not multiplied by two, and `environment.center_components` retains both components.

- [ ] **Step 2: Run coordination tests and verify failure**

Run:

```bash
pytest tests/geometry/test_coordination.py -q
```

Expected: FAIL because coordination result types do not exist.

- [ ] **Step 3: Implement explicit coordination results**

In `src/cristma/geometry/coordination.py`, add frozen slotted records:

```python
@dataclass(frozen=True, slots=True)
class CoordinationEnvironment:
    center_atom_id: str
    center_components: tuple[SiteComponent, ...]
    neighbors: tuple[Neighbor | PeriodicNeighbor, ...]

    @property
    def coordination_number(self) -> int:
        return len(self.neighbors)


@dataclass(frozen=True, slots=True)
class CoordinationResult:
    environments: tuple[CoordinationEnvironment, ...]
    diagnostics: tuple[Diagnostic, ...] = ()

    def by_atom(self, atom_id: str) -> CoordinationEnvironment: ...
```

`by_atom()` raises `KeyError(atom_id)` for an unknown atom.

- [ ] **Step 4: Implement graph-neutral analysis**

`CoordinationAnalyzer.analyze(view, graph)` validates that graph and view atom IDs match exactly, then creates one environment per view atom in view order. It takes graph edges as-is and does not expand `SiteComponent` values into geometric neighbors.

- [ ] **Step 5: Run geometry tests**

Run:

```bash
pytest tests/geometry -q
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/cristma/geometry tests/geometry/test_coordination.py
git commit -m "feat: add coordination analysis"
```

---

### Task 7: Verify the complete CIF-to-coordination slice and package API

**Files:**
- Create: `tests/integration/test_structure_core.py`
- Modify: `tests/test_public_api.py`
- Modify: `src/cristma/__init__.py`
- Modify: `src/cristma/structure/__init__.py`
- Modify: `src/cristma/symmetry/__init__.py`
- Modify: `src/cristma/geometry/__init__.py`
- Modify: `docs/superpowers/specs/2026-08-30-cristma-domain-language-vision.md`
- Modify: `docs/superpowers/specs/2026-08-30-structural-hierarchy-design.md`

**Interfaces:**
- Consumes: all Structure Core interfaces from Tasks 1–6 and the existing CIF registry.
- Produces: stable internal imports for the completed slice and evidence that a built wheel works in a clean environment.

- [ ] **Step 1: Write the real-fixture integration test**

Create `tests/integration/test_structure_core.py`:

```python
def test_real_cif_traverses_structure_core() -> None:
    result = cristma.read("tests/fixtures/cif/cod_3000098_barium_borate.cif")

    assert result.ok
    crystal = result.structures[0]
    view = expand_structure(crystal)
    graph = NeighborFinder(cutoff=2.6).find(view)
    coordination = CoordinationAnalyzer().analyze(view, graph)

    assert view.atoms
    assert len(graph.atoms) == len(view.atoms)
    assert len(coordination.environments) == len(view.atoms)
    assert all(atom.source_site_id for atom in view.atoms)
```

Add a second integration test for `tests/fixtures/cif/mixed_disorder.cif` proving one mixed expanded position and no zero-length neighbor edge.

- [ ] **Step 2: Run integration tests and verify any missing exports**

Run:

```bash
pytest tests/integration/test_structure_core.py -q
```

Expected before export cleanup: collection succeeds scientifically, while one or more intended package-level imports may be absent.

- [ ] **Step 3: Finalize namespace exports**

Expose:

```python
from cristma.structure import (
    AtomicPosition,
    AtomicView,
    ExpandedAtom,
    PeriodicAtomRef,
)
from cristma.symmetry import (
    SymmetryImageProvenance,
    expand_structure,
)
from cristma.geometry import (
    CoordinationAnalyzer,
    CoordinationEnvironment,
    CoordinationResult,
    Neighbor,
    NeighborFinder,
    NeighborGraph,
    PeriodicNeighbor,
    PeriodicNeighborGraph,
)
```

Keep top-level `cristma` limited to broadly used entry points and version metadata; do not re-export every geometry record there.

- [ ] **Step 4: Update stale architectural names**

Replace conceptual `ExpandedAtomRef` references in the architecture and hierarchy specs with `ExpandedAtom` or `PeriodicAtomRef` according to meaning. Do not change historical implementation plans, which document the API that existed when they were executed.

- [ ] **Step 5: Run targeted integration and the full suite once**

Run:

```bash
pytest tests/integration/test_structure_core.py -q
pytest -q
```

Expected: both commands PASS; the full suite has no stale constructor or import failures.

- [ ] **Step 6: Build the wheel**

Run:

```bash
python -m build
```

Expected: one `cristma-0.1.0.dev0-*.whl` and one source distribution are created in `dist/`.

- [ ] **Step 7: Verify installation in a clean temporary environment**

Run from the repository root:

```bash
CRISTMA_VERIFY_DIR="$(mktemp -d)"
python -m venv "$CRISTMA_VERIFY_DIR/venv"
"$CRISTMA_VERIFY_DIR/venv/bin/python" -m pip install --quiet dist/cristma-0.1.0.dev0-*.whl
"$CRISTMA_VERIFY_DIR/venv/bin/python" -c "from cristma.structure import AtomicView, ExpandedAtom, PeriodicAtomRef; from cristma.symmetry import expand_structure; from cristma.geometry import NeighborFinder, PeriodicNeighborGraph, CoordinationAnalyzer; print('installed Structure Core imports OK')"
```

Expected output: `installed Structure Core imports OK`.

- [ ] **Step 8: Record the outstanding SHELX integration gate**

In `docs/superpowers/plans/2026-08-30-native-structure-readers-roadmap.md`, add the exact post-SHELX verification command:

```bash
pytest tests/io/shelx tests/integration/test_structure_core_shelx.py -q
```

The future `test_structure_core_shelx.py` must read `tests/fixtures/shelx/zdk288.res`, obtain a canonical structure, call `expand_structure()`, construct its neighbor graph, and run `CoordinationAnalyzer`. Do not create a skipped placeholder test in this plan.

- [ ] **Step 9: Commit**

```bash
git add src/cristma tests docs/superpowers/specs/2026-08-30-cristma-domain-language-vision.md docs/superpowers/specs/2026-08-30-structural-hierarchy-design.md docs/superpowers/plans/2026-08-30-native-structure-readers-roadmap.md
git commit -m "feat: complete structure core geometry slice"
```

## Completion gate

This implementation plan is complete when Tasks 1–7 pass, the full existing
suite passes once, and the development wheel imports in a clean environment.
That completes the CIF-backed Structure Core geometry slice. It does not make
CrIStMa public and does not close the separate SHELX reader phase.
