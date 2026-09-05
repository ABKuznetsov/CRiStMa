# Native XYZ and extXYZ Reader Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add dependency-free native XYZ/extXYZ reading that maps complete frames into canonical CrIStMa structures and exposes trajectories lazily through `cristma.read(...)`.

**Architecture:** A lightweight probe and frame index retain the exact source without allocating all trajectory arrays. The selected frame is parsed into immutable XYZ records, mapped independently to `MolecularStructure` or identity-symmetry `CrystalStructure`, and cached by `StructureSequence`. Molecular and periodic structures share the existing immutable atomic-property contract.

**Tech Stack:** Python 3.11+, standard library, NumPy, pytest, setuptools.

**Spec:** `docs/superpowers/specs/2026-09-01-native-xyz-reader-design.md`

## Global Constraints

- No ASE, libAtoms `extxyz`, pymatgen, Open Babel, RDKit, Qt, or application dependency.
- `Lattice` without explicit `pbc` never enables periodicity.
- Units are never inferred from property names.
- All declared extXYZ columns are validated and preserved with source provenance.
- Plain XYZ is molecular and interprets only chemical identity plus Cartesian position.
- Multi-frame numerical data remains lazy and each successfully loaded frame is cached once.
- Applications use only the unchanged public `cristma.read(...)` boundary.
- Use focused tests during Tasks 1–6; run the complete suite only in Task 7.

---

### Task 1: Add canonical prerequisites for molecular properties and atomic numbers

**Files:**
- Modify: `src/cristma/chemistry/elements.py`
- Modify: `src/cristma/chemistry/__init__.py`
- Modify: `src/cristma/structure/molecular.py`
- Test: `tests/chemistry/test_elements.py`
- Modify: `tests/structure/test_properties.py`

**Interfaces:**
- Produces: `element_from_atomic_number(value: int) -> str`.
- Extends: `MolecularStructure(..., properties: AtomicPropertyTable | None = None)`.
- Guarantees: `MolecularStructure.atomic_view().properties is structure.properties`.

- [ ] **Step 1: Write failing atomic-number and molecular-property tests**

```python
def test_atomic_number_maps_to_iupac_symbol() -> None:
    assert element_from_atomic_number(1) == "H"
    assert element_from_atomic_number(14) == "Si"
    assert element_from_atomic_number(118) == "Og"


@pytest.mark.parametrize("value", [0, 119, -1, True])
def test_invalid_atomic_number_is_rejected(value) -> None:
    with pytest.raises((TypeError, ValueError)):
        element_from_atomic_number(value)


def test_molecular_properties_reach_atomic_view() -> None:
    molecule = MolecularStructure(
        "water",
        atoms,
        properties=AtomicPropertyTable(
            3,
            (AtomicProperty("charge", np.array([-0.8, 0.4, 0.4])),),
        ),
    )
    assert molecule.atomic_view().properties is molecule.properties
```

- [ ] **Step 2: Run the focused tests and verify RED**

Run: `python3 -m pytest -p no:cacheprovider tests/chemistry/test_elements.py tests/structure/test_properties.py -q`

Expected: FAIL because `element_from_atomic_number` and molecular properties do not exist.

- [ ] **Step 3: Implement the canonical capabilities**

Keep an ordered periodic-table tuple alongside `ELEMENT_SYMBOLS`, validate that atomic numbers are non-boolean integers in `[1, 118]`, and return the one-based element. Add `properties` to `MolecularStructure`, default it to an empty table of the exact atom count, reject mismatched leading dimensions, and pass the same table to `AtomicView`.

```python
def element_from_atomic_number(value: int) -> str:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError("atomic number must be an integer")
    if not 1 <= value <= len(ELEMENT_SYMBOLS_BY_ATOMIC_NUMBER):
        raise ValueError("atomic number must lie between 1 and 118")
    return ELEMENT_SYMBOLS_BY_ATOMIC_NUMBER[value - 1]
```

- [ ] **Step 4: Run the focused tests and verify GREEN**

Run: `python3 -m pytest -p no:cacheprovider tests/chemistry/test_elements.py tests/structure/test_properties.py tests/structure/test_atomic_view.py -q`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/cristma/chemistry src/cristma/structure/molecular.py tests/chemistry/test_elements.py tests/structure/test_properties.py
git commit -m "feat(structure): support molecular atomic properties"
```

### Task 2: Define immutable XYZ source contracts and index complete frames

**Files:**
- Create: `src/cristma/io/xyz/__init__.py`
- Create: `src/cristma/io/xyz/document.py`
- Create: `src/cristma/io/xyz/index.py`
- Test: `tests/io/xyz/test_document.py`
- Test: `tests/io/xyz/test_index.py`

**Interfaces:**
- Produces: `XyzPropertySpec`, `XyzFrameSpan`, `XyzDocument`, `XyzFrame`.
- Produces: `index_xyz(source: str, source_name: str | None = None) -> tuple[XyzDocument, tuple[Diagnostic, ...]]`.
- Consumes later: `XyzFrameSpan` character offsets and declared atom count.

- [ ] **Step 1: Write failing document validation tests**

```python
def test_property_spec_requires_known_type_and_positive_width() -> None:
    assert XyzPropertySpec("forces", "R", 3).width == 3
    with pytest.raises(ValueError):
        XyzPropertySpec("", "R", 3)
    with pytest.raises(ValueError):
        XyzPropertySpec("forces", "Q", 3)
    with pytest.raises(ValueError):
        XyzPropertySpec("forces", "R", 0)


def test_document_preserves_exact_source() -> None:
    document, diagnostics = index_xyz(TWO_FRAME_SOURCE, "trajectory.xyz")
    assert diagnostics == ()
    assert document.render_preserved() == TWO_FRAME_SOURCE
    assert len(document.frames) == 2
```

- [ ] **Step 2: Write failing frame-index tests**

Cover zero-atom frames, Unicode comments, CRLF input, changing atom counts,
trailing blank lines, a negative/non-integral count, a blank line between
frames, and a final frame declaring more rows than remain.

```python
def test_truncated_tail_is_not_indexed() -> None:
    document, diagnostics = index_xyz(ONE_COMPLETE_ONE_TRUNCATED, "run.xyz")
    assert len(document.frames) == 1
    assert any(item.code == "xyz.frame.incomplete" for item in diagnostics)


def test_unicode_offsets_slice_exact_frame() -> None:
    document, _ = index_xyz(UNICODE_SOURCE, "модель.xyz")
    span = document.frames[0]
    assert document.raw_source[span.start_offset:span.end_offset] == UNICODE_SOURCE
```

- [ ] **Step 3: Run the focused tests and verify RED**

Run: `python3 -m pytest -p no:cacheprovider tests/io/xyz/test_document.py tests/io/xyz/test_index.py -q`

Expected: FAIL because `cristma.io.xyz` does not exist.

- [ ] **Step 4: Implement immutable records and the count-driven index**

`XyzFrameSpan` stores `index`, `atom_count`, frame/comment/row character offsets.
`XyzDocument` stores exact decoded text, optional source name, and frame spans.
`XyzFrame` stores parsed frame values but no application state. Arrays are
copied and marked read-only in `__post_init__`.

The index advances strictly by `2 + atom_count` physical lines per complete
frame. It accepts trailing blank lines but diagnoses blank lines between
frames. A malformed count stops further indexing because subsequent frame
boundaries are unknowable. A truncated tail remains in `raw_source` and is not
added to `frames`.

- [ ] **Step 5: Run the focused tests and verify GREEN**

Run: `python3 -m pytest -p no:cacheprovider tests/io/xyz/test_document.py tests/io/xyz/test_index.py -q`

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/cristma/io/xyz tests/io/xyz/test_document.py tests/io/xyz/test_index.py
git commit -m "feat(io): index XYZ frames without materializing arrays"
```

### Task 3: Parse extXYZ metadata, schemas, and selected atom rows

**Files:**
- Create: `src/cristma/io/xyz/metadata.py`
- Create: `src/cristma/io/xyz/parser.py`
- Test: `tests/io/xyz/test_metadata.py`
- Test: `tests/io/xyz/test_parser.py`

**Interfaces:**
- Produces: `parse_xyz_metadata(comment: str) -> XyzMetadata`.
- Produces: `parse_property_schema(value: str) -> tuple[XyzPropertySpec, ...]`.
- Produces: `load_xyz_frame(document: XyzDocument, reference: FrameReference) -> XyzFrame`.
- Produces: `validate_xyz_frame(document: XyzDocument, span: XyzFrameSpan) -> tuple[Diagnostic, ...]`.
- Consumes: one indexed frame only; does not map canonical structures.

- [ ] **Step 1: Write failing metadata and schema tests**

```python
def test_extxyz_metadata_parses_special_and_unknown_values() -> None:
    metadata = parse_xyz_metadata(
        'Lattice="2 0 0 0 2 0 0 0 2" '
        'Properties=species:S:1:pos:R:3:forces:R:3 '
        'pbc="T T F" energy=-1.25 label="relaxed cell"'
    )
    assert metadata.lattice.tolist() == [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    assert metadata.pbc == (True, True, False)
    assert metadata.values["energy"] == -1.25
    assert metadata.values["label"] == "relaxed cell"


def test_schema_rejects_duplicate_names_or_wrong_triplets() -> None:
    with pytest.raises(ValueError, match="duplicate"):
        parse_property_schema("species:S:1:pos:R:3:pos:R:3")
    with pytest.raises(ValueError, match="triplets"):
        parse_property_schema("species:S:1:broken")
```

Include quoted strings, escaped quotes, integers, reals, logical spellings,
one-dimensional arrays, malformed `Lattice`, malformed `pbc`, and duplicate
special keys.

- [ ] **Step 2: Write failing selected-frame row tests**

```python
def test_all_declared_property_types_are_typed() -> None:
    frame = load_xyz_frame(document, document_reference)
    assert frame.columns["species"].dtype.kind in {"U", "O"}
    assert frame.columns["Z"].dtype.kind == "i"
    assert frame.columns["forces"].shape == (2, 3)
    assert frame.columns["fixed"].dtype.kind == "b"


def test_row_width_must_equal_schema_width() -> None:
    with pytest.raises(ValueError, match="column count"):
        load_xyz_frame(broken_document, broken_reference)
```

Also verify ordinary XYZ implicit columns, ignored schema-free trailing tokens,
invalid numeric/logical cells, `species` plus `Z`, and source provenance.

`validate_xyz_frame` walks the selected text without constructing NumPy arrays.
It reports malformed metadata/schema/row values, unknown species,
`species`/`Z` conflicts, uninterpreted plain columns, and the conservative
`Lattice`-without-`pbc` warning. This keeps all source diagnostics available in
the immutable `ReadResult` while canonical structures remain lazy.

- [ ] **Step 3: Run focused tests and verify RED**

Run: `python3 -m pytest -p no:cacheprovider tests/io/xyz/test_metadata.py tests/io/xyz/test_parser.py -q`

Expected: FAIL because metadata and frame parsers do not exist.

- [ ] **Step 4: Implement a dependency-free metadata lexer**

Use a small character scanner rather than `str.split()` so quoted values,
escaped quotes, spaces around `=`, and arrays remain intact. Preserve raw
entries and return typed special values. Parse `Properties` as exact triplets;
convert row slices according to `S/I/R/L`; use only the selected source span.

Plain XYZ creates implicit columns `species`, `pos`; its arbitrary comment is
stored under `comment`. Extra row tokens are ignored scientifically and emit
`xyz.map.uninterpreted_plain_columns` during indexing or frame loading.

- [ ] **Step 5: Run focused tests and verify GREEN**

Run: `python3 -m pytest -p no:cacheprovider tests/io/xyz/test_metadata.py tests/io/xyz/test_parser.py tests/io/xyz/test_index.py -q`

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/cristma/io/xyz tests/io/xyz/test_metadata.py tests/io/xyz/test_parser.py
git commit -m "feat(io): parse plain and extended XYZ frames"
```

### Task 4: Map XYZ frames into canonical molecular and periodic structures

**Files:**
- Create: `src/cristma/io/xyz/mapper.py`
- Test: `tests/io/xyz/test_mapper.py`

**Interfaces:**
- Consumes: `map_xyz_frame(frame: XyzFrame) -> Structure`.
- Produces: `MolecularStructure` or identity-only `CrystalStructure`.
- Produces: immutable `AtomicPropertyTable` excluding structural columns.

- [ ] **Step 1: Write failing molecular mapping tests**

```python
def test_plain_xyz_maps_to_molecule() -> None:
    structure = map_xyz_frame(plain_frame)
    assert isinstance(structure, MolecularStructure)
    assert structure.periodic == (False, False, False)
    assert structure.atoms[0].components[0].element == "O"


def test_lattice_without_pbc_remains_molecular() -> None:
    structure = map_xyz_frame(lattice_only_frame)
    assert isinstance(structure, MolecularStructure)
    assert structure.metadata["xyz_lattice"] == lattice_tuple
```

Cover element symbols, atomic numbers, unknown species, `species`/`Z`
agreement and conflict, frame naming, metadata, stable atom IDs, and property
provenance with `unit is None`.

- [ ] **Step 2: Write failing periodic mapping tests**

```python
def test_explicit_periodic_extxyz_maps_to_crystal() -> None:
    structure = map_xyz_frame(periodic_frame)
    assert isinstance(structure, CrystalStructure)
    assert structure.periodic == (True, True, False)
    assert structure.space_group.provenance == "unreported_identity"
    assert np.allclose(structure.atomic_view().cartesian, canonical_positions)


def test_true_pbc_without_lattice_is_rejected() -> None:
    with pytest.raises(ValueError, match="lattice"):
        map_xyz_frame(pbc_without_lattice)
```

Use an analytically rotated lattice to prove that fractional positions are
computed in the reported frame and then reconstructed in CrIStMa's canonical
cell orientation. Prove arbitrary width-three properties are not rotated.

- [ ] **Step 3: Run focused tests and verify RED**

Run: `python3 -m pytest -p no:cacheprovider tests/io/xyz/test_mapper.py -q`

Expected: FAIL because the mapper does not exist.

- [ ] **Step 4: Implement canonical mapping**

Resolve each atom to `ElementSpecies` or `UnknownSpecies`; occupancy is one.
For molecular frames create `MolecularAtom` rows in reported Cartesian space.
For periodic frames derive `UnitCell`, compute fractional coordinates as
`cartesian @ inv(reported_lattice)`, create independent sites, and use
`CrystalStructure.explicit(..., periodic=frame.pbc)`.

Build `AtomicPropertyTable` from every non-structural schema property. Preserve
reported arrays exactly and attach `PropertyProvenance(source_name,
f"Properties:{name}", "reported")`. Do not infer units or vector transformation
semantics.

- [ ] **Step 5: Run focused tests and verify GREEN**

Run: `python3 -m pytest -p no:cacheprovider tests/io/xyz/test_mapper.py tests/structure/test_properties.py tests/symmetry/test_orbit.py -q`

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/cristma/io/xyz/mapper.py tests/io/xyz/test_mapper.py
git commit -m "feat(io): map XYZ frames to canonical structures"
```

### Task 5: Expose lazy XYZ trajectories through the format registry

**Files:**
- Create: `src/cristma/io/xyz/probe.py`
- Create: `src/cristma/io/xyz/handler.py`
- Modify: `src/cristma/io/xyz/parser.py`
- Modify: `src/cristma/io/xyz/__init__.py`
- Modify: `src/cristma/io/formats.py`
- Create: `tests/io/xyz/test_handler.py`
- Modify: `tests/io/test_builtin_formats.py`
- Modify: `tests/test_public_api.py`

**Interfaces:**
- Produces: `probe_xyz(source: str) -> float`.
- Produces: `parse_xyz(source: str, source_name: str | None = None) -> ReadResult`.
- Produces: `XyzFormatHandler.read_text(...) -> ReadResult`.
- Extends: `cristma.read(path)` and `cristma.read_text(source, format="xyz")`.

- [ ] **Step 1: Write failing lazy-sequence tests**

```python
def test_multiframe_xyz_maps_only_requested_frame(monkeypatch) -> None:
    calls = []
    monkeypatch.setattr(mapper, "map_xyz_frame", lambda frame: calls.append(frame) or sentinel(frame))
    result = parse_xyz(THREE_FRAME_XYZ, "trajectory.xyz")
    assert len(result.structures) == 3
    assert calls == []
    assert result.structures.final.name == "frame 3"
    assert len(calls) == 1
    assert result.structures.final is result.structures[-1]
```

Verify each `FrameReference` source span and role, changing schemas, caching,
and retention of complete frames before an incomplete tail.

- [ ] **Step 2: Write failing registry/public API tests**

```python
def test_xyz_descriptor_is_lazy_and_content_aware() -> None:
    descriptor = FormatRegistry(builtin_format_descriptors()).select(XYZ_TEXT)
    assert descriptor.name == "xyz"
    assert descriptor.aliases == ("extxyz",)
    assert descriptor.capabilities.multiple
    assert descriptor.capabilities.lazy_frames


def test_public_read_maps_xyz(tmp_path: Path) -> None:
    path = tmp_path / "molecule.xyz"
    path.write_text(WATER_XYZ)
    result = cristma.read(path)
    assert result.ok
    assert isinstance(result.structures[0], MolecularStructure)
```

Add a subprocess assertion that `builtin_format_descriptors()` does not import
`cristma.io.xyz.parser` or `cristma.io.xyz.mapper`.

- [ ] **Step 3: Run focused tests and verify RED**

Run: `python3 -m pytest -p no:cacheprovider tests/io/xyz/test_handler.py tests/io/test_builtin_formats.py tests/test_public_api.py -q`

Expected: FAIL because XYZ is not registered.

- [ ] **Step 4: Implement lazy parsing, probe, handler, and descriptor**

`parse_xyz` indexes source, calls `validate_xyz_frame` for each span to collect
source diagnostics without allocating trajectory arrays, and creates
`FrameReference` values whose loader calls `load_xyz_frame` and
`map_xyz_frame` on demand. Complete frames except the last are `intermediate`;
the last is `final`. Frames with validation errors remain source-addressable in
the sequence but raise the same precise error if loaded. The descriptor is:

```python
FormatDescriptor(
    name="xyz",
    aliases=("extxyz",),
    suffixes=(".xyz", ".extxyz"),
    basenames=(),
    probe=probe_xyz,
    factory=_xyz_handler,
    capabilities=FormatCapabilities(text=True, multiple=True, lazy_frames=True),
)
```

The probe validates a non-negative count, a comment line, and the first few
declared atom rows. Valid `Properties=` receives confidence `0.98`; coherent
plain XYZ receives `0.80`. Package exports use lazy `__getattr__` so descriptor
construction does not import parser/mapper modules.

- [ ] **Step 5: Run focused tests and verify GREEN**

Run: `python3 -m pytest -p no:cacheprovider tests/io/xyz tests/io/test_builtin_formats.py tests/io/test_registry.py tests/test_public_api.py tests/structure/test_sequence.py -q`

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/cristma/io/xyz src/cristma/io/formats.py tests/io/xyz tests/io/test_builtin_formats.py tests/test_public_api.py
git commit -m "feat(io): register native XYZ structure reader"
```

### Task 6: Verify fixtures and cross-format scientific equivalence

**Files:**
- Create: `tests/fixtures/xyz/water.xyz`
- Create: `tests/fixtures/xyz/water.extxyz`
- Create: `tests/fixtures/xyz/silicon.extxyz`
- Create: `tests/fixtures/xyz/trajectory.xyz`
- Create: `tests/fixtures/xyz/PROVENANCE.md`
- Create: `tests/io/xyz/test_reference_fixtures.py`
- Create: `tests/integration/test_structure_core_xyz.py`

**Interfaces:**
- Verifies: public reading, lazy frames, typed atomic properties, molecular geometry, periodic geometry, neighbors, and coordination.

- [ ] **Step 1: Add provenance-recorded reference fixtures**

Use hand-reduced examples whose grammar follows the libAtoms extXYZ
specification. Record the official URL, reduction date, SHA-256, atom/frame
counts, periodicity, and the fact that they are reference fixtures rather than
claimed physical-calculation outputs.

- [ ] **Step 2: Write fixture and integration tests**

```python
@pytest.mark.parametrize(
    ("name", "frames"),
    [("water.xyz", 1), ("water.extxyz", 1), ("silicon.extxyz", 1), ("trajectory.xyz", 3)],
)
def test_xyz_fixture_reads_through_public_api(name, frames) -> None:
    result = cristma.read(FIXTURES / name)
    assert result.ok
    assert len(result.structures) == frames


def test_plain_and_extended_water_have_equal_geometry() -> None:
    plain = cristma.read(FIXTURES / "water.xyz").structures[0]
    extended = cristma.read(FIXTURES / "water.extxyz").structures[0]
    assert np.allclose(plain.atomic_view().cartesian, extended.atomic_view().cartesian)


def test_extxyz_and_poscar_silicon_have_equal_periodic_coordination() -> None:
    xyz = cristma.read(FIXTURES / "silicon.extxyz").structures[0]
    poscar = cristma.read(VASP_FIXTURES / "POSCAR").structures[0]
    assert np.allclose(xyz.cell.matrix, poscar.cell.matrix)
    assert coordination_signature(xyz) == coordination_signature(poscar) == (6,)
```

Also prove extXYZ forces and user columns reach `AtomicView` with `unit is None`.

- [ ] **Step 3: Run the integration slice**

Run: `python3 -m pytest -p no:cacheprovider tests/io/xyz/test_reference_fixtures.py tests/integration/test_structure_core_xyz.py -q`

Expected: PASS. If a fixture exposes a format gap, reduce it into the responsible parser test before changing production code.

- [ ] **Step 4: Run all XYZ and shared contract tests**

Run: `python3 -m pytest -p no:cacheprovider tests/io/xyz tests/integration/test_structure_core_xyz.py tests/chemistry/test_elements.py tests/structure/test_properties.py tests/structure/test_sequence.py tests/io/test_builtin_formats.py tests/test_public_api.py -q`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add tests/fixtures/xyz tests/io/xyz/test_reference_fixtures.py tests/integration/test_structure_core_xyz.py
git commit -m "test(io): verify XYZ structure workflows"
```

### Task 7: Document and verify the installable XYZ slice

**Files:**
- Modify: `README.md`
- Modify: `docs/superpowers/plans/2026-08-30-native-structure-readers-roadmap.md`
- Verify: `pyproject.toml`

**Interfaces:**
- Documents: ordinary/extXYZ support, explicit-periodicity rule, arbitrary typed columns, and lazy trajectories.
- Verifies: complete suite and installed wheel outside the repository.

- [ ] **Step 1: Update user-facing documentation and roadmap**

Add a runnable example:

```python
result = cristma.read("trajectory.extxyz")
trajectory = result.structures
final = trajectory.final
print(len(trajectory), final.name, tuple(final.atomic_view().properties))
```

State that XYZ reading uses only NumPy and the standard library; `Lattice`
without explicit `pbc` remains molecular; writer and bond inference are outside
this slice. Mark roadmap row 9 implemented on the feature branch.

- [ ] **Step 2: Scan for placeholders and dependency leaks**

Run: `rg -n "TODO|FIXME|NotImplementedError|import ase|import extxyz|pymatgen|openbabel|rdkit|PySide|Qt" src/cristma/io/xyz README.md`

Expected: no implementation placeholders or forbidden runtime imports.

- [ ] **Step 3: Run the complete CrIStMa suite exactly once**

Run: `python3 -m pytest -p no:cacheprovider -q`

Expected: all collected tests PASS.

- [ ] **Step 4: Build and install the wheel without network dependencies**

```bash
python3 -m pip wheel . --no-deps --no-build-isolation --wheel-dir /private/tmp/cristma-xyz-wheel
python3 -m venv --system-site-packages /private/tmp/cristma-xyz-smoke/venv
/private/tmp/cristma-xyz-smoke/venv/bin/python -m pip install --no-deps /private/tmp/cristma-xyz-wheel/cristma-0.1.0.dev0-py3-none-any.whl
```

From `/private/tmp`, use the installed wheel to read every XYZ fixture, access
the final lazy frame, inspect typed properties, and assert `cristma.__file__`
resides inside the temporary environment.

- [ ] **Step 5: Verify the final tree and commit documentation**

Run: `git diff --check`

Run: `git status --short`

Expected: only intended README/roadmap edits before commit, then a clean tree.

```bash
git add README.md docs/superpowers/plans/2026-08-30-native-structure-readers-roadmap.md
git commit -m "docs(io): document native XYZ structure input"
```
