# Native VASP Structure Reader Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add dependency-free native reading of POSCAR, CONTCAR, XDATCAR, OUTCAR, and `vasprun.xml`, mapping every complete structure into canonical CRiStMa snapshots and exposing trajectories lazily.

**Architecture:** One lazily registered `vasp` format family owns four independent source parsers. Each parser produces a loss-preserving document and either one normalized `VaspSnapshot` or indexed lazy frame loaders; one mapper converts snapshots into `CrystalStructure`. A small general extension carries typed independent-site properties through symmetry expansion so Selective Dynamics, velocities, and forces remain scientific arrays rather than metadata.

**Tech Stack:** Python 3.11+, standard library (`dataclasses`, `xml.parsers.expat`, `xml.etree.ElementTree`), NumPy, existing CRiStMa structure/I/O contracts, pytest.

**Spec:** `docs/superpowers/specs/2026-08-31-native-vasp-reader-design.md`

## Global Constraints

- No ASE, pymatgen, Gemmi, Qt, VASP installation, or application dependency.
- `CrystalStructure` is the only periodic scientific snapshot returned by VASP readers.
- POSCAR/CONTCAR return `StructureCollection`; trajectories return lazy `StructureSequence`.
- VASP input never implies a reported space group; use identity symmetry with `unreported_identity` provenance.
- Invalid or incomplete source data produces diagnostics, never plausible partial structures.
- Unknown source content remains in the document; semantic calculations never consume parser records.
- During development run only the task's focused tests and established shared contracts; run all CRiStMa tests once at the final gate.

---

### Task 1: Carry typed site properties into expanded atomic views

**Files:**
- Modify: `src/cristma/structure/crystal.py`
- Modify: `src/cristma/symmetry/orbit.py`
- Test: `tests/structure/test_properties.py`
- Test: `tests/symmetry/test_orbit.py`

**Interfaces:**
- Consumes: `AtomicProperty`, `AtomicPropertyTable`, `CrystalStructure.sites`.
- Produces: `CrystalStructure.properties: AtomicPropertyTable`; `expand_structure(crystal)` copies each independent-site property row to every symmetry-expanded atom derived from that site.

- [ ] **Step 1: Write failing property ownership and expansion tests**

```python
def test_crystal_property_rows_follow_independent_sites() -> None:
    table = AtomicPropertyTable(
        1,
        (AtomicProperty("selective_dynamics", np.array([[True, False, True]])),),
    )
    crystal = CrystalStructure.explicit("demo", cell(), (site(),), properties=table)
    assert crystal.properties["selective_dynamics"].values.tolist() == [[True, False, True]]


def test_site_properties_expand_with_symmetry_images() -> None:
    crystal = inversion_crystal(
        properties=AtomicPropertyTable(
            1,
            (AtomicProperty("force", np.array([[1.0, 2.0, 3.0]]), unit="eV/angstrom"),),
        )
    )
    view = crystal.atomic_view()
    assert view.properties["force"].values.tolist() == [[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]
```

- [ ] **Step 2: Run the focused tests and confirm the constructor/API failure**

Run: `python3 -m pytest -p no:cacheprovider tests/structure/test_properties.py tests/symmetry/test_orbit.py -q`

Expected: FAIL because `CrystalStructure` does not accept `properties` and expanded views are empty.

- [ ] **Step 3: Add validated properties to the canonical crystal**

```python
@dataclass(frozen=True, slots=True)
class CrystalStructure:
    # existing fields ...
    properties: AtomicPropertyTable | None = field(default=None, compare=False)

    def __post_init__(self) -> None:
        # existing validation ...
        properties = self.properties or AtomicPropertyTable(len(self.sites))
        if properties.atom_count != len(self.sites):
            raise ValueError("property table atom count does not match independent sites")
        object.__setattr__(self, "properties", properties)
```

In `expand_structure`, build `site_index_by_id`, derive one expanded row index for every `ExpandedAtom.source_site_id`, and construct new `AtomicProperty` arrays with `np.take(prop.values, row_indices, axis=0)`. Copy `unit`, `missing`, `source_name`, and `provenance` unchanged.

- [ ] **Step 4: Run focused tests**

Run: `python3 -m pytest -p no:cacheprovider tests/structure/test_properties.py tests/symmetry/test_orbit.py tests/integration/test_structure_core.py -q`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/cristma/structure/crystal.py src/cristma/symmetry/orbit.py tests/structure/test_properties.py tests/symmetry/test_orbit.py
git commit -m "feat(structure): propagate typed site properties"
```

### Task 2: Define loss-preserving VASP documents and exact scale mathematics

**Files:**
- Create: `src/cristma/io/vasp/__init__.py`
- Create: `src/cristma/io/vasp/document.py`
- Create: `src/cristma/io/vasp/numeric.py`
- Test: `tests/io/vasp/test_document.py`
- Test: `tests/io/vasp/test_numeric.py`

**Interfaces:**
- Produces: `VaspScale`, `VaspHeader`, `VaspAtomRow`, `VaspFrameSpan`, `VaspSnapshot`, `PoscarDocument`, `XdatcarDocument`, `OutcarDocument`, `VasprunDocument`.
- Produces: `scaled_lattice(scale, raw_lattice) -> ndarray`, `scaled_cartesian(scale, rows) -> ndarray`, `fractional_from_cartesian(cell, rows) -> ndarray`.

- [ ] **Step 1: Write failing exact-scale tests**

```python
def test_negative_scalar_reconstructs_requested_volume() -> None:
    lattice = scaled_lattice(VaspScale((-64.0,)), np.eye(3))
    assert np.linalg.det(lattice) == pytest.approx(64.0)


def test_three_scalars_scale_cartesian_components() -> None:
    raw = np.array([[1.0, 1.0, 1.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])
    assert scaled_lattice(VaspScale((2.0, 3.0, 4.0)), raw).tolist() == [
        [2.0, 3.0, 4.0], [0.0, 6.0, 0.0], [0.0, 0.0, 12.0]
    ]


@pytest.mark.parametrize("values", [(0.0,), (-1.0, 2.0, 3.0), (1.0, 2.0)])
def test_invalid_scale_shape_or_domain_is_rejected(values) -> None:
    with pytest.raises(ValueError):
        VaspScale(values)
```

- [ ] **Step 2: Run tests and confirm missing-module failure**

Run: `python3 -m pytest -p no:cacheprovider tests/io/vasp/test_document.py tests/io/vasp/test_numeric.py -q`

Expected: FAIL because `cristma.io.vasp` does not exist.

- [ ] **Step 3: Implement immutable document contracts and numeric helpers**

```python
@dataclass(frozen=True, slots=True)
class VaspScale:
    values: tuple[float, ...]

    def __post_init__(self) -> None:
        if len(self.values) not in {1, 3}:
            raise ValueError("VASP scale requires one or three numbers")
        if not all(math.isfinite(value) for value in self.values):
            raise ValueError("VASP scale values must be finite")
        if len(self.values) == 3 and not all(value > 0 for value in self.values):
            raise ValueError("three VASP scale values must be positive")
        if len(self.values) == 1 and self.values[0] == 0:
            raise ValueError("VASP scale must not be zero")
```

`VaspSnapshot` contains `name`, a numeric `(3, 3)` lattice in angstrom, ordered species, fractional positions, optional selective flags, velocities/unit, forces/unit, `frame_index`, and `SourceReference`. Every array is copied and marked read-only in `__post_init__`.

```python
@dataclass(frozen=True, slots=True)
class VaspSnapshot:
    name: str
    lattice: np.ndarray
    species: tuple[ChemicalSpecies, ...]
    fractional: np.ndarray
    frame_index: int
    source: SourceReference
    selective_dynamics: np.ndarray | None = None
    velocities: np.ndarray | None = None
    velocity_unit: str | None = None
    forces: np.ndarray | None = None
    force_unit: str | None = None
```

- [ ] **Step 4: Run focused tests**

Run: `python3 -m pytest -p no:cacheprovider tests/io/vasp/test_document.py tests/io/vasp/test_numeric.py -q`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/cristma/io/vasp tests/io/vasp/test_document.py tests/io/vasp/test_numeric.py
git commit -m "feat(io): define native VASP source contracts"
```

### Task 3: Parse and map POSCAR/CONTCAR completely

**Files:**
- Create: `src/cristma/io/vasp/poscar.py`
- Create: `src/cristma/io/vasp/mapper.py`
- Test: `tests/io/vasp/test_poscar.py`
- Test: `tests/io/vasp/test_mapper.py`

**Interfaces:**
- Produces: `parse_poscar(source: str, source_name: str | None = None) -> ReadResult` with `PoscarDocument` and diagnostics.
- Produces: `poscar_snapshot(document: PoscarDocument) -> VaspSnapshot`.
- Produces: `map_vasp_snapshot(snapshot: VaspSnapshot) -> CrystalStructure`.

- [ ] **Step 1: Write failing parser tests for VASP 5, VASP 4, coordinates, constraints, and velocities**

```python
def test_vasp5_cartesian_selective_and_velocity_sections() -> None:
    result = parse_poscar(POSCAR_CARTESIAN_SELECTIVE_VELOCITIES, "CONTCAR")
    document = result.document
    assert document.header.species == ("Si", "O")
    assert document.header.counts == (1, 2)
    assert document.positions[0].selective == (True, False, True)
    assert len(document.velocities) == 3
    assert document.render_preserved() == POSCAR_CARTESIAN_SELECTIVE_VELOCITIES


def test_vasp4_keeps_unknown_species_explicit() -> None:
    result = parse_poscar(VASP4_POSCAR, "POSCAR")
    snapshot = poscar_snapshot(result.document)
    crystal = map_vasp_snapshot(snapshot)
    assert crystal.sites[0].components[0].element is None
    assert any(item.code == "vasp.map.species_unresolved" for item in result.diagnostics)
```

Also assert negative-volume and three-factor Cartesian coordinates against hand-calculated cells/fractional positions, zero-count species, incomplete coordinates, invalid flags, lattice-velocity preservation, and unknown trailing text.

- [ ] **Step 2: Run tests and confirm missing parser/mapper failure**

Run: `python3 -m pytest -p no:cacheprovider tests/io/vasp/test_poscar.py tests/io/vasp/test_mapper.py -q`

Expected: FAIL because the parser and mapper do not exist.

- [ ] **Step 3: Implement the POSCAR state machine**

Parse by declared counts, never by guessing atom-looking lines. Coordinate mode follows VASP exactly:

```python
def coordinate_mode(token: str) -> Literal["cartesian", "direct"]:
    return "cartesian" if token[:1] in {"C", "c", "K", "k"} else "direct"
```

Determine whether the row before counts is species or counts by requiring every count token to be a non-negative integer. Retain physical spans and all trailing sections in `PoscarDocument`. Return a document even when semantic diagnostics prevent mapping.

- [ ] **Step 4: Implement canonical snapshot mapping**

Create one `IndependentSite` per declared coordinate row with occupancy 1.0 and stable ID `vasp:{source}:frame:{frame}:site:{index}`. Use `CrystalStructure.explicit`, attach `AtomicProperty("selective_dynamics", dtype=bool)`, velocities with reported unit/convention, and source provenance. Convert all final coordinates to fractional values before constructing sites.

- [ ] **Step 5: Run focused tests**

Run: `python3 -m pytest -p no:cacheprovider tests/io/vasp/test_document.py tests/io/vasp/test_numeric.py tests/io/vasp/test_poscar.py tests/io/vasp/test_mapper.py tests/structure/test_properties.py -q`

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/cristma/io/vasp/poscar.py src/cristma/io/vasp/mapper.py tests/io/vasp/test_poscar.py tests/io/vasp/test_mapper.py
git commit -m "feat(io): read POSCAR and CONTCAR structures"
```

### Task 4: Register the public VASP handler

**Files:**
- Create: `src/cristma/io/vasp/probe.py`
- Create: `src/cristma/io/vasp/handler.py`
- Modify: `src/cristma/io/vasp/__init__.py`
- Modify: `src/cristma/io/formats.py`
- Test: `tests/io/vasp/test_handler.py`
- Modify: `tests/io/test_builtin_formats.py`
- Modify: `tests/test_public_api.py`

**Interfaces:**
- Produces: `probe_vasp(source: str) -> float`.
- Produces: `VaspFormatHandler.read_text(...) -> ReadResult`.
- Extends: `cristma.read(path)` and `cristma.read_text(source, format="vasp")` without changing their signatures.

- [ ] **Step 1: Write failing registry and public API tests**

```python
def test_poscar_basename_and_content_select_lazy_vasp_descriptor() -> None:
    registry = FormatRegistry(builtin_format_descriptors())
    descriptor = registry.select(POSCAR_TEXT, basename="renamed.data")
    assert descriptor.name == "vasp"
    assert descriptor.capabilities.multiple
    assert descriptor.capabilities.lazy_frames


def test_public_read_maps_poscar(tmp_path: Path) -> None:
    path = tmp_path / "POSCAR"
    path.write_text(POSCAR_TEXT)
    result = cristma.read(path)
    assert result.ok
    assert result.structures[0].name == "Silicon"
```

Add a subprocess assertion that calling `builtin_format_descriptors()` does not import `cristma.io.vasp.poscar`, `outcar`, or `vasprun`.

- [ ] **Step 2: Run focused tests and verify selection failure**

Run: `python3 -m pytest -p no:cacheprovider tests/io/vasp/test_handler.py tests/io/test_builtin_formats.py tests/test_public_api.py -q`

Expected: FAIL because no VASP descriptor is registered.

- [ ] **Step 3: Implement probe, handler, lazy exports, and descriptor**

The descriptor uses basenames `POSCAR`, `CONTCAR`, `XDATCAR`, `OUTCAR`, `vasprun.xml`, aliases `("poscar", "contcar", "xdatcar", "outcar", "vasprun")`, and no misleading generic suffix except `.xml` at content-qualified confidence. The handler dispatches internally by content and basename to the four native parser functions; applications never perform this dispatch.

- [ ] **Step 4: Run focused tests**

Run: `python3 -m pytest -p no:cacheprovider tests/io/vasp tests/io/test_builtin_formats.py tests/io/test_registry.py tests/test_public_api.py -q`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/cristma/io/vasp src/cristma/io/formats.py tests/io/vasp/test_handler.py tests/io/test_builtin_formats.py tests/test_public_api.py
git commit -m "feat(io): register native VASP structure reader"
```

### Task 5: Index and lazily load XDATCAR trajectories

**Files:**
- Create: `src/cristma/io/vasp/xdatcar.py`
- Test: `tests/io/vasp/test_xdatcar.py`
- Modify: `src/cristma/io/vasp/handler.py`

**Interfaces:**
- Produces: `parse_xdatcar(source, source_name=None) -> ReadResult` containing `XdatcarDocument` and `StructureSequence`.
- Produces: `load_xdatcar_snapshot(document, reference) -> VaspSnapshot`.

- [ ] **Step 1: Write failing lazy-frame tests**

```python
def test_xdatcar_indexes_complete_frames_without_mapping_them(monkeypatch) -> None:
    calls = []
    monkeypatch.setattr(mapper, "map_vasp_snapshot", lambda frame: calls.append(frame) or sentinel())
    result = parse_xdatcar(THREE_FRAME_XDATCAR, "XDATCAR")
    assert len(result.structures) == 3
    assert calls == []
    assert result.structures[-1].name.endswith("configuration 3")
    assert len(calls) == 1
    assert result.structures.final is result.structures[-1]


def test_incomplete_trailing_configuration_is_diagnostic_only() -> None:
    result = parse_xdatcar(TWO_COMPLETE_ONE_TRUNCATED, "XDATCAR")
    assert len(result.structures) == 2
    assert result.structures.references[-1].role == "final"
    assert any(item.code == "vasp.xdatcar.frame_incomplete" for item in result.diagnostics)
```

- [ ] **Step 2: Run focused tests and verify missing implementation**

Run: `python3 -m pytest -p no:cacheprovider tests/io/vasp/test_xdatcar.py -q`

Expected: FAIL because `parse_xdatcar` does not exist.

- [ ] **Step 3: Implement header reuse, span indexing, and lazy loader**

Scan only configuration markers and declared row counts. Store `FrameReference(index, role, source, metadata={"configuration": reported_number})`; the closure passed to `StructureSequence` calls `load_xdatcar_snapshot` and `map_vasp_snapshot` only on access. A complete final frame receives role `final`.

- [ ] **Step 4: Run focused tests**

Run: `python3 -m pytest -p no:cacheprovider tests/io/vasp/test_xdatcar.py tests/structure/test_sequence.py tests/io/vasp/test_mapper.py -q`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/cristma/io/vasp/xdatcar.py src/cristma/io/vasp/handler.py tests/io/vasp/test_xdatcar.py
git commit -m "feat(io): load XDATCAR frames lazily"
```

### Task 6: Index OUTCAR ionic structures and typed forces

**Files:**
- Create: `src/cristma/io/vasp/outcar.py`
- Test: `tests/io/vasp/test_outcar.py`
- Modify: `src/cristma/io/vasp/handler.py`

**Interfaces:**
- Produces: `parse_outcar(source, source_name=None) -> ReadResult` with `OutcarDocument` and lazy `StructureSequence`.
- Produces: `load_outcar_snapshot(document, reference) -> VaspSnapshot` with forces in `eV/angstrom`.

- [ ] **Step 1: Write failing frame/species/force tests**

```python
def test_outcar_uses_explicit_species_cell_positions_and_forces() -> None:
    result = parse_outcar(TWO_IONIC_STEP_OUTCAR, "OUTCAR")
    assert len(result.structures) == 2
    final = result.structures.final
    assert [site.components[0].element for site in final.sites] == ["Na", "Cl"]
    assert final.properties["force"].unit == "eV/angstrom"
    assert final.atomic_view().properties["force"].values.shape == (2, 3)


def test_electronic_iterations_do_not_create_structure_frames() -> None:
    result = parse_outcar(ONE_IONIC_MANY_ELECTRONIC_OUTCAR, "OUTCAR")
    assert len(result.structures) == 1
```

Add cases for changing lattice vectors, inconsistent `ions per type`, missing species, and a truncated final POSITION block.

- [ ] **Step 2: Run focused tests and confirm missing implementation**

Run: `python3 -m pytest -p no:cacheprovider tests/io/vasp/test_outcar.py -q`

Expected: FAIL because the OUTCAR parser does not exist.

- [ ] **Step 3: Implement an explicit-block source index**

Record `VRHFIN`/`TITEL`, `ions per type`, each complete `direct lattice vectors` block, and each complete `POSITION ... TOTAL-FORCE` block. Associate each POSITION frame with the most recent valid lattice. Do not parse atom-looking rows elsewhere. Frame loaders convert Cartesian positions to fractional and add forces as an immutable property.

- [ ] **Step 4: Run focused tests**

Run: `python3 -m pytest -p no:cacheprovider tests/io/vasp/test_outcar.py tests/io/vasp/test_mapper.py tests/structure/test_properties.py -q`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/cristma/io/vasp/outcar.py src/cristma/io/vasp/handler.py tests/io/vasp/test_outcar.py
git commit -m "feat(io): read OUTCAR ionic structures lazily"
```

### Task 7: Index vasprun.xml structures with standard-library XML

**Files:**
- Create: `src/cristma/io/vasp/vasprun.py`
- Test: `tests/io/vasp/test_vasprun.py`
- Modify: `src/cristma/io/vasp/handler.py`

**Interfaces:**
- Produces: `parse_vasprun(source, source_name=None) -> ReadResult` with `VasprunDocument` and lazy `StructureSequence`.
- Produces: `load_vasprun_snapshot(document, reference) -> VaspSnapshot`.

- [ ] **Step 1: Write failing XML correctness and laziness tests**

```python
def test_vasprun_indexes_ionic_steps_and_final_structure_once() -> None:
    result = parse_vasprun(THREE_STEP_XML, "vasprun.xml")
    assert len(result.structures) == 3
    assert result.structures.references[-1].role == "final"
    final = result.structures.final
    assert final.properties["force"].values.shape == (2, 3)


def test_truncated_xml_returns_document_and_diagnostic_without_fake_tail() -> None:
    result = parse_vasprun(TRUNCATED_AFTER_ONE_COMPLETE_STEP, "vasprun.xml")
    assert len(result.structures) == 1
    assert any(item.code == "vasp.vasprun.xml_incomplete" for item in result.diagnostics)
```

Also test namespaces, initial-only single-point data, atominfo/count mismatches, malformed numeric arrays, and repeated access caching.

- [ ] **Step 2: Run focused tests and verify missing implementation**

Run: `python3 -m pytest -p no:cacheprovider tests/io/vasp/test_vasprun.py -q`

Expected: FAIL because `parse_vasprun` does not exist.

- [ ] **Step 3: Implement XML indexing without materializing structures**

Use `xml.parsers.expat.ParserCreate()` callbacks and `CurrentByteIndex` to record complete `<calculation>` and named `<structure>` byte spans. Precompute UTF-8 byte-to-character boundaries once so `SourceReference` offsets address `raw_source`. Parse global `<atominfo>` once. The lazy loader wraps and parses only the selected fragment with `xml.etree.ElementTree`, extracting `basis`, fractional `positions`, optional `velocities`, and `forces`.

Deduplicate a separately reported `finalpos` only when its source identity explicitly refers to the already indexed final calculation; do not compare rounded geometry.

- [ ] **Step 4: Run focused tests**

Run: `python3 -m pytest -p no:cacheprovider tests/io/vasp/test_vasprun.py tests/structure/test_sequence.py tests/io/vasp/test_mapper.py -q`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/cristma/io/vasp/vasprun.py src/cristma/io/vasp/handler.py tests/io/vasp/test_vasprun.py
git commit -m "feat(io): read vasprun XML structures lazily"
```

### Task 8: Verify real fixtures and cross-format scientific equivalence

**Files:**
- Create: `tests/fixtures/vasp/POSCAR`
- Create: `tests/fixtures/vasp/XDATCAR`
- Create: `tests/fixtures/vasp/OUTCAR`
- Create: `tests/fixtures/vasp/vasprun.xml`
- Create: `tests/fixtures/vasp/PROVENANCE.md`
- Create: `tests/io/vasp/test_real_fixtures.py`
- Create: `tests/integration/test_structure_core_vasp.py`

**Interfaces:**
- Verifies: public `cristma.read`, all four VASP source families, canonical geometry, typed properties, lazy roles, neighbors, and coordination.

- [ ] **Step 1: Add provenance-recorded fixtures and integration acceptance tests**

`PROVENANCE.md` records for each fixture: original source path or official example URL, copy date `2026-08-31`, SHA-256, VASP family, expected atom/frame count, and any intentional truncation.

```python
@pytest.mark.parametrize(
    ("name", "expected_frames"),
    [("POSCAR", 1), ("XDATCAR", 3), ("OUTCAR", 2), ("vasprun.xml", 2)],
)
def test_real_vasp_fixture_reads(name, expected_frames) -> None:
    result = cristma.read(FIXTURES / name)
    assert result.ok
    assert len(result.structures) == expected_frames


def test_equivalent_cif_and_poscar_have_equal_geometry_and_coordination() -> None:
    cif = cristma.read(CIF_EQUIVALENT).structures[0]
    poscar = cristma.read(FIXTURES / "POSCAR").structures[0]
    assert np.allclose(cif.cell.matrix, poscar.cell.matrix)
    assert np.allclose(cif.atomic_view().fractional, poscar.atomic_view().fractional)
    assert coordination_signature(cif) == coordination_signature(poscar)
```

- [ ] **Step 2: Run the integration tests and capture any real-format gaps**

Run: `python3 -m pytest -p no:cacheprovider tests/io/vasp/test_real_fixtures.py tests/integration/test_structure_core_vasp.py -q`

Expected: PASS when the analytic implementation covers the real sources. Any
failure is evidence of one specific unsupported construct and must be reduced
before changing production code.

- [ ] **Step 3: Reduce and correct a fixture gap only when Step 2 fails**

Use systematic debugging to identify the exact source record. Add a minimal
regression such as the following to the responsible parser test, run it red,
implement only that grammar rule, and run it green:

```python
def test_outcar_accepts_version_specific_position_separator() -> None:
    result = parse_outcar(MINIMAL_SOURCE_WITH_REPORTED_SEPARATOR, "OUTCAR")
    assert result.ok
    assert len(result.structures) == 1
```

If Step 2 passes, perform no production edit in this step. Do not add
electronic-property parsing or dependencies outside the specification.

- [ ] **Step 4: Run all VASP and shared contract tests**

Run: `python3 -m pytest -p no:cacheprovider tests/io/vasp tests/integration/test_structure_core_vasp.py tests/structure/test_properties.py tests/structure/test_sequence.py tests/io/test_builtin_formats.py tests/test_public_api.py -q`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add tests/fixtures/vasp tests/io/vasp/test_real_fixtures.py tests/integration/test_structure_core_vasp.py src/cristma/io/vasp
git commit -m "test(io): verify real VASP structure workflows"
```

### Task 9: Document and verify the installable VASP slice

**Files:**
- Modify: `README.md`
- Modify: `docs/superpowers/plans/2026-08-30-native-structure-readers-roadmap.md`
- Verify: `pyproject.toml`

**Interfaces:**
- Documents: one-line `cristma.read(...)`, supported VASP sources, lazy trajectory access, and unresolved VASP 4 species behavior.

- [ ] **Step 1: Update user-facing documentation and roadmap status**

Add this runnable example:

```python
result = cristma.read("XDATCAR")
trajectory = result.structures
final = trajectory.final
print(len(trajectory), final.cell.volume)
```

State explicitly that VASP reading needs only NumPy and the standard library; electronic results and canonical VASP writing are outside this milestone.

- [ ] **Step 2: Scan the implementation for placeholders and dependency leaks**

Run: `rg -n "TODO|FIXME|NotImplementedError|pymatgen|ase|PySide|Qt" src/cristma/io/vasp README.md`

Expected: no implementation placeholders and no runtime imports of forbidden dependencies; documentation may mention them only to state they are unnecessary.

- [ ] **Step 3: Run the complete CRiStMa test suite once**

Run: `python3 -m pytest -p no:cacheprovider -q`

Expected: all collected tests PASS.

- [ ] **Step 4: Build a wheel without network isolation and install it into a clean temporary environment**

```bash
python3 -m pip wheel . --no-deps --no-build-isolation --wheel-dir /private/tmp/cristma-vasp-wheel
python3 -m venv --system-site-packages /private/tmp/cristma-vasp-smoke/venv
/private/tmp/cristma-vasp-smoke/venv/bin/python -m pip install --no-deps /private/tmp/cristma-vasp-wheel/cristma-0.1.0.dev0-py3-none-any.whl
```

From outside the repository, use the installed wheel to read all four real fixtures, access the final lazy frames, expand atoms, and assert the imported `cristma.__file__` resides inside the temporary environment.

- [ ] **Step 5: Verify the final tree and commit documentation**

Run: `git diff --check`

Then run: `git status --short`

Expected: only the intended README/roadmap changes before committing, then an empty status.

```bash
git add README.md docs/superpowers/plans/2026-08-30-native-structure-readers-roadmap.md
git commit -m "docs(io): document native VASP structure input"
```
