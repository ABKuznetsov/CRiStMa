# Crystallographic Catalog and Orbits Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship a self-contained, setting-aware space-group and Wyckoff catalog and use it to calculate and validate crystallographic orbits from CRiStMa structures.

**Architecture:** A development-only compiler converts the pinned spglib 2.7.0 database into deterministic CRiStMa JSON resources. Small immutable catalog records load those resources; `OrbitAnalyzer` reuses the existing exact `expand_orbit()` implementation to calculate multiplicity and stabilizers, then matches catalog Wyckoff positions and emits explicit diagnostics. CIF mapping may resolve missing operations from the catalog, but file-specific representations terminate at the I/O boundary.

**Tech Stack:** Python 3.11+, standard library, NumPy 1.26+, pytest 8+; development-only `spglib==2.7.0`; no new runtime dependency.

**Spec:** `docs/superpowers/specs/2026-09-01-cristma-structural-crystallography-design.md`

## Global Constraints

- CRiStMa remains independent of Qt and all application packages.
- Runtime dependencies remain `numpy>=1.26` only.
- The upstream source is spglib tag `v2.7.0`, commit `12355c77fb7c505a55f52cae36341d73b781a065`, under BSD-3-Clause.
- Input SHA-256 values are `4457df1042b14a65ea62af0bad7b5b609a4fc33592df245802cfe005b221f95e` for `database/spg.csv` and `d3d786a1f0187e5c6d69a3ade35648ffab34fd1b977d61ad84d8b0434b8b7ca0` for `database/Wyckoff.csv`.
- Generated resources retain the spglib copyright and BSD-3-Clause notice.
- Gemmi and Bilbao are comparators/oracles only; their databases are not packaged.
- International Tables content is not copied.
- `CrystalStructure`, `IndependentSite` and source documents remain immutable.
- Existing `cristma.symmetry.expand_orbit()` is the only position-expansion implementation.
- Calculated orbit multiplicity is authoritative; catalog and reported values are validation inputs.
- Focused tests run after every task; the full suite and wheel installation run only at the final gate.
- This plan does not add bonds, BVS, polyhedra, hierarchy, morphology, faceting, twins, diffraction or refinement.
- Existing periodic geometry is verified for compatibility but not redesigned here; additional geometric-shell APIs receive a separate plan only when an application needs them.

---

## File map

```text
tools/
├── __init__.py
└── compile_spglib_crystallography.py
    Development-only deterministic compiler. Reads pinned spg.csv and
    Wyckoff.csv, asks spglib 2.7.0 for operations, and writes normalized JSON.

src/cristma/reference_data/resources/crystallography/
├── space_groups.json
├── wyckoff_positions.json
├── SOURCE.md
└── SPGLIB_LICENSE.txt
    Generated runtime resources and mandatory provenance.

src/cristma/crystallography/
├── __init__.py
├── space_group.py
│   Immutable keys and records; conversion to SpaceGroupDefinition.
├── wyckoff.py
│   Exact affine Wyckoff representatives and typed positions.
├── catalog.py
│   Resource loading, validation and lookup indexes.
└── orbit.py
    OrbitAnalyzer, SiteSymmetry, CrystallographicOrbit and diagnostics.

tests/tools/
└── test_compile_spglib_crystallography.py

tests/crystallography/
├── test_models.py
├── test_catalog.py
├── test_orbit_analyzer.py
└── test_cif_integration.py

tests/fixtures/spglib/
├── spg_minimal.csv
└── Wyckoff_minimal.csv
```

Existing files modified:

```text
pyproject.toml
src/cristma/reference_data/facade.py
src/cristma/reference_data/__init__.py
src/cristma/io/cif/mapper.py
src/cristma/io/cif/handler.py
tests/io/cif/test_mapper_advanced.py
```

---

### Task 1: Deterministic spglib database compiler

**Files:**
- Create: `tools/__init__.py`
- Create: `tools/compile_spglib_crystallography.py`
- Create: `tests/tools/test_compile_spglib_crystallography.py`
- Create: `tests/fixtures/spglib/spg_minimal.csv`
- Create: `tests/fixtures/spglib/Wyckoff_minimal.csv`
- Modify: `pyproject.toml`

**Interfaces:**
- Consumes: spglib `spg.csv`, `Wyckoff.csv`, and `spglib.get_symmetry_from_database(hall_number)` from exactly version 2.7.0.
- Produces: `compile_catalog(spg_path: Path, wyckoff_path: Path, output_dir: Path, *, upstream_commit: str, compiled_date: str) -> tuple[Path, Path]` and two canonical JSON resources.

- [ ] **Step 1: Add minimal upstream-format fixtures**

Use real spglib rows for Hall entries 1, 2 and 390. `spg_minimal.csv` contains the unmodified corresponding rows from `spg.csv`. `Wyckoff_minimal.csv` contains the complete blocks beginning with `1:`, `2:` and `390:` and terminates with `end of data`.

The P-42₁m block must include these records:

```text
390:P -4 2sub1 m:::::::
::8:f:1:(x,y,z):(-x,-y,z):(y,-x,-z):(-y,x,-z)
:::::(-x+1/2,y+1/2,-z):(x+1/2,-y+1/2,-z):(-y+1/2,-x+1/2,z):(y+1/2,x+1/2,z)
::4:e:..m:(x,x+1/2,z):(-x,-x+1/2,z):(x+1/2,-x,-z):(-x+1/2,x,-z)
::4:d:2..:(0,0,z):(0,0,-z):(1/2,1/2,-z):(1/2,1/2,z)
::2:c:2.mm:(0,1/2,z):(1/2,0,-z)::
::2:b:-4..:(0,0,1/2):(1/2,1/2,1/2)::
::2:a:-4..:(0,0,0):(1/2,1/2,0)::
```

- [ ] **Step 2: Write the failing compiler tests**

```python
def test_compiler_emits_canonical_catalogs(tmp_path: Path) -> None:
    space_groups, wyckoffs = compile_catalog(
        FIXTURES / "spg_minimal.csv",
        FIXTURES / "Wyckoff_minimal.csv",
        tmp_path,
        upstream_commit=SPGLIB_COMMIT,
        compiled_date="2026-09-01",
    )

    groups = json.loads(space_groups.read_text(encoding="utf-8"))
    positions = json.loads(wyckoffs.read_text(encoding="utf-8"))

    assert [record["hall_number"] for record in groups["records"]] == [1, 2, 390]
    assert groups["records"][2]["hall_symbol"] == "P -4 2ab"
    assert len(groups["records"][2]["operations"]) == 8
    assert positions["records"]["390"][0]["letter"] == "f"
    assert positions["records"]["390"][-1]["letter"] == "a"
    assert positions["records"]["390"][-1]["multiplicity"] == 2


def test_compiler_output_is_byte_deterministic(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    kwargs = {"upstream_commit": SPGLIB_COMMIT, "compiled_date": "2026-09-01"}
    compile_catalog(SPG, WYCKOFF, first, **kwargs)
    compile_catalog(SPG, WYCKOFF, second, **kwargs)
    assert (first / "space_groups.json").read_bytes() == (second / "space_groups.json").read_bytes()
    assert (first / "wyckoff_positions.json").read_bytes() == (second / "wyckoff_positions.json").read_bytes()


def test_compiler_rejects_wrong_spglib_version(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(spglib, "__version__", "2.6.0")
    with pytest.raises(RuntimeError, match="spglib 2.7.0 is required"):
        compile_catalog(SPG, WYCKOFF, Path("unused"), upstream_commit=SPGLIB_COMMIT, compiled_date="2026-09-01")
```

- [ ] **Step 3: Run the focused tests and verify failure**

Run:

```bash
pytest tests/tools/test_compile_spglib_crystallography.py -q
```

Expected: collection fails because `tools.compile_spglib_crystallography` does not exist.

- [ ] **Step 4: Implement the safe source parsers and canonical writer**

The compiler must:

1. use `csv.reader` for `spg.csv` and the colon-delimited continuation logic from spglib's BSD-licensed `make_Wyckoff_db.py` for `Wyckoff.csv`;
2. validate Hall numbers are exactly the expected fixture subset or `1..530` for a full build;
3. obtain operations from `spglib.get_symmetry_from_database(hall_number)`;
4. convert every translation to a reduced rational pair using `Fraction(float(value)).limit_denominator(24)` and reject an approximation error above `1e-12`;
5. parse each Wyckoff coordinate into an exact affine parameter map, supporting `x`, `y`, `z`, signs, integer coefficients and rational constants without `eval`;
6. expand every representative by the conventional centering translations used by spglib, exactly as rational fractions:

```python
CENTERING = {
    "P": ((0, 0, 0),),
    "A": ((0, 0, 0), (0, 12, 12)),
    "B": ((0, 0, 0), (12, 0, 12)),
    "C": ((0, 0, 0), (12, 12, 0)),
    "I": ((0, 0, 0), (12, 12, 12)),
    "F": ((0, 0, 0), (0, 12, 12), (12, 0, 12), (12, 12, 0)),
    "R": ((0, 0, 0),),
    "H": ((0, 0, 0), (16, 8, 8), (8, 16, 16)),
}
# integer values above are divided by 24 before addition
```

7. verify expanded representative count equals the reported multiplicity;
8. serialize fractions as `[numerator, denominator]`;
9. serialize with `json.dumps(..., sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n"`;
10. sort space-group records by Hall number while preserving the upstream general-to-special Wyckoff order;
11. write schema metadata into both files:

```python
metadata = {
    "dataset_id": "cristma.crystallography.spglib",
    "schema_version": "1.0.0",
    "upstream": "spglib",
    "upstream_version": "2.7.0",
    "upstream_commit": upstream_commit,
    "license": "BSD-3-Clause",
    "compiled_date": compiled_date,
}
```

The normalized space-group record is:

```python
{
    "hall_number": 390,
    "number": 113,
    "choice": "",
    "hall_symbol": "P -4 2ab",
    "hm_short": "P -4 2_1 m",
    "hm_full": "P -4 2_1 m",
    "point_group": "-42m",
    "centering": "P",
    "crystal_system": "tetragonal",
    "operations": [
        {
            "rotation": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            "translation": [[0, 1], [0, 1], [0, 1]],
        }
    ],
}
```

The normalized Wyckoff record is:

```python
{
    "letter": "a",
    "multiplicity": 2,
    "site_symmetry": "-4..",
    "representatives": [
        {
            "parameter_matrix": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            "translation": [[0, 1], [0, 1], [0, 1]],
            "source": "(0,0,0)",
        },
        {
            "parameter_matrix": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            "translation": [[1, 2], [1, 2], [0, 1]],
            "source": "(1/2,1/2,0)",
        },
    ],
}
```

- [ ] **Step 5: Add the development-only dependency group**

Add to `pyproject.toml`:

```toml
[project.optional-dependencies]
test = [
    "build>=1.2",
    "pytest>=8",
]
reference-build = [
    "spglib==2.7.0",
]
```

Do not add spglib to `[project].dependencies`.

- [ ] **Step 6: Run the compiler tests**

Run:

```bash
pytest tests/tools/test_compile_spglib_crystallography.py -q
```

Expected: all tests pass.

- [ ] **Step 7: Commit**

```bash
git add pyproject.toml tools/__init__.py tools/compile_spglib_crystallography.py tests/tools/test_compile_spglib_crystallography.py tests/fixtures/spglib
git commit -m "feat: add deterministic spglib catalog compiler"
```

---

### Task 2: Generate and document the complete packaged catalog

**Files:**
- Create: `src/cristma/reference_data/resources/crystallography/space_groups.json`
- Create: `src/cristma/reference_data/resources/crystallography/wyckoff_positions.json`
- Create: `src/cristma/reference_data/resources/crystallography/SOURCE.md`
- Create: `src/cristma/reference_data/resources/crystallography/SPGLIB_LICENSE.txt`
- Modify: `pyproject.toml`
- Test: `tests/tools/test_compile_spglib_crystallography.py`

**Interfaces:**
- Consumes: Task 1 compiler, pinned `/private/tmp` inputs or newly downloaded inputs whose SHA-256 values match Global Constraints.
- Produces: complete immutable runtime resources containing 530 Hall settings and all corresponding Wyckoff records.

- [ ] **Step 1: Add a failing full-resource metadata test**

```python
def test_packaged_catalog_has_complete_pinned_source() -> None:
    root = files("cristma.reference_data").joinpath("resources/crystallography")
    groups = json.loads(root.joinpath("space_groups.json").read_text(encoding="utf-8"))
    wyckoffs = json.loads(root.joinpath("wyckoff_positions.json").read_text(encoding="utf-8"))

    assert groups["metadata"]["upstream_version"] == "2.7.0"
    assert groups["metadata"]["upstream_commit"] == SPGLIB_COMMIT
    assert len(groups["records"]) == 530
    assert {record["number"] for record in groups["records"]} == set(range(1, 231))
    assert set(wyckoffs["records"]) == {str(number) for number in range(1, 531)}
    assert all(wyckoffs["records"][str(number)] for number in range(1, 531))
    assert root.joinpath("SPGLIB_LICENSE.txt").is_file()
```

- [ ] **Step 2: Run the test and verify failure**

Run:

```bash
pytest tests/tools/test_compile_spglib_crystallography.py::test_packaged_catalog_has_complete_pinned_source -q
```

Expected: FAIL because resources do not exist.

- [ ] **Step 3: Generate the complete resources**

Verify inputs first:

```bash
shasum -a 256 /private/tmp/cristma-spg-v270.csv /private/tmp/cristma-Wyckoff-v270.csv
```

Expected hashes are the exact values in Global Constraints. Then run:

```bash
python tools/compile_spglib_crystallography.py \
  --spg /private/tmp/cristma-spg-v270.csv \
  --wyckoff /private/tmp/cristma-Wyckoff-v270.csv \
  --output src/cristma/reference_data/resources/crystallography \
  --upstream-commit 12355c77fb7c505a55f52cae36341d73b781a065 \
  --compiled-date 2026-09-01
```

Copy the complete unmodified spglib `COPYING` text to
`SPGLIB_LICENSE.txt`. `SOURCE.md` must list the two upstream raw URLs pinned to
the commit, both input hashes, both generated output hashes, the exact command
above, schema version, compilation date, and this statement:

```text
The normalized JSON is generated from spglib 2.7.0 database files under the
BSD-3-Clause license. No International Tables pages or Bilbao database records
are copied into this package.
```

- [ ] **Step 4: Include nested resources in wheels**

Replace the package-data rule with:

```toml
[tool.setuptools.package-data]
cristma = [
    "reference_data/resources/*.json",
    "reference_data/resources/crystallography/*.json",
    "reference_data/resources/crystallography/*.md",
    "reference_data/resources/crystallography/*.txt",
]
```

- [ ] **Step 5: Run the full-resource compiler tests**

Run:

```bash
pytest tests/tools/test_compile_spglib_crystallography.py -q
```

Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml src/cristma/reference_data/resources/crystallography tests/tools/test_compile_spglib_crystallography.py
git commit -m "data: add pinned spglib crystallography catalog"
```

---

### Task 3: Typed space-group and Wyckoff records

**Files:**
- Create: `src/cristma/crystallography/__init__.py`
- Create: `src/cristma/crystallography/space_group.py`
- Create: `src/cristma/crystallography/wyckoff.py`
- Create: `tests/crystallography/test_models.py`

**Interfaces:**
- Consumes: existing `AffineOperation`, `SpaceGroupDefinition` and exact `Fraction` values.
- Produces: `SpaceGroupKey`, `SpaceGroupRecord`, `AffineCoordinateMap`, `WyckoffPosition`.

- [ ] **Step 1: Write failing model tests**

```python
def test_space_group_record_builds_existing_definition() -> None:
    identity = parse_xyz_operation("x,y,z", operation_id="hall:1:op:1")
    record = SpaceGroupRecord(
        key=SpaceGroupKey(1, "P 1", ""),
        number=1,
        hm_short="P 1",
        hm_full="P 1",
        point_group="1",
        centering="P",
        crystal_system="triclinic",
        operations=(identity,),
        wyckoff_positions=(),
    )
    definition = record.definition(provenance="derived")
    assert isinstance(definition, SpaceGroupDefinition)
    assert definition.number == 1
    assert definition.hall_symbol == "P 1"
    assert definition.operations == (identity,)


def test_affine_coordinate_map_evaluates_exact_parameterization() -> None:
    representative = AffineCoordinateMap.from_xyz("x,x+1/2,-z")
    assert representative.evaluate((Fraction(1, 4), Fraction(0), Fraction(1, 3))) == (
        Fraction(1, 4), Fraction(3, 4), Fraction(-1, 3)
    )


def test_wyckoff_position_rejects_wrong_representative_count() -> None:
    with pytest.raises(ValueError, match="multiplicity"):
        WyckoffPosition(
            space_group_key=SpaceGroupKey(390, "P -4 2ab", ""),
            letter="a",
            multiplicity=2,
            site_symmetry_symbol="-4..",
            representatives=(AffineCoordinateMap.from_xyz("0,0,0"),),
        )
```

- [ ] **Step 2: Run tests and verify failure**

Run:

```bash
pytest tests/crystallography/test_models.py -q
```

Expected: import failure because `cristma.crystallography` does not exist.

- [ ] **Step 3: Implement immutable typed records**

Use these contracts:

```python
@dataclass(frozen=True, slots=True)
class SpaceGroupKey:
    hall_number: int
    hall_symbol: str
    choice: str


@dataclass(frozen=True, slots=True)
class AffineCoordinateMap:
    parameter_matrix: tuple[tuple[Fraction, Fraction, Fraction], ...]
    translation: tuple[Fraction, Fraction, Fraction]
    source: str | None = None

    @classmethod
    def from_xyz(cls, text: str) -> "AffineCoordinateMap": ...

    def evaluate(
        self,
        parameters: tuple[Fraction, Fraction, Fraction],
    ) -> tuple[Fraction, Fraction, Fraction]: ...


@dataclass(frozen=True, slots=True)
class WyckoffPosition:
    space_group_key: SpaceGroupKey
    letter: str
    multiplicity: int
    site_symmetry_symbol: str
    representatives: tuple[AffineCoordinateMap, ...]

    @property
    def degrees_of_freedom(self) -> int:
        return _matrix_rank_over_rationals(self.representatives[0].parameter_matrix)


@dataclass(frozen=True, slots=True)
class SpaceGroupRecord:
    key: SpaceGroupKey
    number: int
    hm_short: str
    hm_full: str
    point_group: str
    centering: str
    crystal_system: str
    operations: tuple[AffineOperation, ...]
    wyckoff_positions: tuple[WyckoffPosition, ...]

    def definition(self, *, provenance: SymmetryProvenance) -> SpaceGroupDefinition: ...
```

`SpaceGroupRecord.definition()` stores `key.choice or None` in the existing
`SpaceGroupDefinition.setting` field so axis, setting and origin-choice codes
survive conversion. It does not invent a separate origin value when the
upstream choice is composite.

`AffineCoordinateMap.from_xyz()` must use a tokenizer/parser and never `eval`.
Unlike `AffineOperation`, its parameter matrix need not be unimodular and may
contain coefficients such as `2x` or `x-y`. Normalize evaluated coordinates
modulo one only at the comparison boundary, not inside `evaluate()`.

Implement private helper
`_matrix_rank_over_rationals(matrix: tuple[tuple[Fraction, ...], ...]) -> int`
in `wyckoff.py` using exact Gaussian elimination; the property calls this
helper. Do not use floating-point `numpy.linalg.matrix_rank` for catalog DOF.

Validation includes:

- Hall number in `1..530`;
- IT number in `1..230`;
- nonempty symbols and operations;
- Wyckoff letter `a..z`;
- positive multiplicity;
- representative count equals multiplicity;
- all positions in one record refer to the same `SpaceGroupKey`.

- [ ] **Step 4: Run model tests**

Run:

```bash
pytest tests/crystallography/test_models.py -q
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/cristma/crystallography tests/crystallography/test_models.py
git commit -m "feat: add crystallographic catalog models"
```

---

### Task 4: Validated runtime catalog and ReferenceData integration

**Files:**
- Create: `src/cristma/crystallography/catalog.py`
- Create: `tests/crystallography/test_catalog.py`
- Modify: `src/cristma/crystallography/__init__.py`
- Modify: `src/cristma/reference_data/facade.py`
- Modify: `src/cristma/reference_data/__init__.py`

**Interfaces:**
- Consumes: Task 2 resources and Task 3 typed records.
- Produces: `SpaceGroupCatalog.default()`, `by_hall_number()`, `by_hall_symbol()`, `by_number()`, `by_hm_symbol()`, and `ReferenceData.crystallography`.

- [ ] **Step 1: Write failing catalog tests**

```python
def test_default_catalog_has_all_hall_settings() -> None:
    catalog = SpaceGroupCatalog.default()
    assert len(catalog) == 530
    assert catalog.by_hall_number(390).number == 113
    assert catalog.by_hall_number(390).key.hall_symbol == "P -4 2ab"


def test_number_lookup_preserves_setting_ambiguity() -> None:
    records = SpaceGroupCatalog.default().by_number(5)
    assert len(records) > 1
    assert len({record.key.hall_number for record in records}) == len(records)


def test_hall_symbol_lookup_is_exact_after_whitespace_normalization() -> None:
    record = SpaceGroupCatalog.default().by_hall_symbol("  P   -4   2ab ")
    assert record.key.hall_number == 390


def test_reference_data_exposes_same_cached_catalog() -> None:
    first = ReferenceData.default().crystallography
    second = SpaceGroupCatalog.default()
    assert first is second
```

Also test unknown Hall number/symbol (`KeyError`), unsupported schema
(`ValueError`), mismatched resource metadata (`ValueError`), duplicate Hall
number (`ValueError`) and a Wyckoff representative-count mismatch
(`ValueError`).

- [ ] **Step 2: Run tests and verify failure**

Run:

```bash
pytest tests/crystallography/test_catalog.py -q
```

Expected: import/attribute failures for `SpaceGroupCatalog`.

- [ ] **Step 3: Implement immutable loading and indexes**

Use:

```python
@dataclass(frozen=True, slots=True)
class SpaceGroupCatalog:
    records: tuple[SpaceGroupRecord, ...]
    dataset_id: str
    schema_version: str
    resource_sha256: tuple[str, str]

    @classmethod
    @lru_cache(maxsize=1)
    def default(cls) -> "SpaceGroupCatalog": ...

    def by_hall_number(self, hall_number: int) -> SpaceGroupRecord: ...
    def by_hall_symbol(self, hall_symbol: str) -> SpaceGroupRecord: ...
    def by_number(self, number: int) -> tuple[SpaceGroupRecord, ...]: ...
    def by_hm_symbol(self, symbol: str) -> tuple[SpaceGroupRecord, ...]: ...
```

Lookup normalization may collapse whitespace and case only. It must not remove
minus signs, screw-axis subscripts, axis choices or origin-choice information.
`by_hall_symbol()` raises a dedicated `LookupError` if normalization still
matches more than one entry; it never selects the first record.

Load resource bytes with `importlib.resources.files`, compute SHA-256 before
JSON decoding, compare metadata between both resources, build all typed
records, then validate:

- Hall numbers are exactly `1..530` with no duplicates;
- each Wyckoff map key has one space-group record;
- every operation has an exact 3x3 integral rotation and rational translation;
- every group contains identity and its operation IDs are
  `hall:{hall_number}:op:{one_based_index}`;
- the greatest Wyckoff multiplicity equals the operation count;
- letters are unique within one Hall setting;
- multiplicities are positive and non-increasing in source order.

- [ ] **Step 4: Add the catalog to `ReferenceData` without a dependency cycle**

```python
@dataclass(frozen=True, slots=True)
class ReferenceData:
    elements: ElementCatalog
    covalent_radii: CovalentRadii
    chemical: ChemicalReference
    crystallography: SpaceGroupCatalog
```

`reference_data.facade` may import `SpaceGroupCatalog`; catalog modules must not
import `ReferenceData`.

- [ ] **Step 5: Run focused tests**

Run:

```bash
pytest tests/crystallography/test_models.py tests/crystallography/test_catalog.py tests/reference_data -q
```

Expected: all selected tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/cristma/crystallography src/cristma/reference_data tests/crystallography/test_catalog.py
git commit -m "feat: add space group and Wyckoff catalog"
```

---

### Task 5: Orbit multiplicity, stabilizer and site symmetry

**Files:**
- Create: `src/cristma/crystallography/orbit.py`
- Create: `tests/crystallography/test_orbit_analyzer.py`
- Modify: `src/cristma/crystallography/__init__.py`

**Interfaces:**
- Consumes: `IndependentSite`, `UnitCell`, `SpaceGroupDefinition`, optional exact `SpaceGroupRecord`, and existing `expand_orbit()`.
- Produces: `SiteSymmetry`, `CrystallographicOrbit`, `OrbitAnalyzer.analyze()`.

- [ ] **Step 1: Write failing analytic orbit tests**

```python
def test_p_minus_one_general_position_has_multiplicity_two() -> None:
    record = SpaceGroupCatalog.default().by_hall_number(2)
    result = OrbitAnalyzer().analyze(
        site_at(0.1, 0.2, 0.3),
        record.definition(provenance="derived"),
        catalog_record=record,
        cell=cubic_cell(),
    )
    assert result.calculated_multiplicity == 2
    assert len(result.stabilizer) == 1
    assert result.site_symmetry.order == 1


def test_p_minus_one_origin_is_special_position() -> None:
    record = SpaceGroupCatalog.default().by_hall_number(2)
    result = OrbitAnalyzer().analyze(
        site_at(0, 0, 0, wyckoff="a", reported_multiplicity=1),
        record.definition(provenance="derived"),
        catalog_record=record,
        cell=cubic_cell(),
    )
    assert result.calculated_multiplicity == 1
    assert len(result.stabilizer) == 2
    assert result.site_symmetry.order == 2
    assert result.site_symmetry.symbol == "-1"


def test_reported_multiplicity_mismatch_is_diagnostic() -> None:
    result = analyze_p_minus_one(site_at(0.1, 0.2, 0.3, reported_multiplicity=1))
    assert "crystallography.orbit.reported_multiplicity_mismatch" in {
        item.code for item in result.diagnostics
    }
```

- [ ] **Step 2: Run tests and verify failure**

Run:

```bash
pytest tests/crystallography/test_orbit_analyzer.py -q
```

Expected: import failures for `OrbitAnalyzer` and result types.

- [ ] **Step 3: Implement result types and analyzer configuration**

```python
@dataclass(frozen=True, slots=True)
class SiteSymmetry:
    symbol: str | None
    stabilizer: tuple[SymmetryImageProvenance, ...]

    @property
    def order(self) -> int:
        return len(self.stabilizer)


@dataclass(frozen=True, slots=True)
class CrystallographicOrbit:
    source_site_id: str
    atoms: tuple[ExpandedAtom, ...]
    calculated_multiplicity: int
    stabilizer: tuple[SymmetryImageProvenance, ...]
    site_symmetry: SiteSymmetry
    wyckoff_position: WyckoffPosition | None
    diagnostics: tuple[Diagnostic, ...]


@dataclass(frozen=True, slots=True)
class OrbitAnalyzer:
    tolerance: float = 1e-6

    def get_config(self) -> dict[str, float]: ...
    def clone(self, **changes: float) -> "OrbitAnalyzer": ...
    def analyze(
        self,
        site: IndependentSite,
        space_group: SpaceGroupDefinition,
        *,
        cell: UnitCell,
        catalog_record: SpaceGroupRecord | None = None,
        structure_id: str | None = None,
    ) -> CrystallographicOrbit: ...
```

`analyze()` calls `expand_orbit()` once. The stabilizer is the complete
`equivalent_images` tuple of the expanded atom periodically equal to the source
coordinates. This works because `expand_orbit()` already merges operations
that map a special position onto one expanded atom while retaining every
operation and normalization translation.

Validate the orbit-stabilizer theorem explicitly:

```python
len(space_group.operations) == calculated_multiplicity * len(stabilizer)
```

If it fails, raise `ValueError("inconsistent orbit and stabilizer sizes")`;
this indicates an invalid group/tolerance, not uncertain source reporting.

Before Task 6, `site_symmetry.symbol` is populated only when exactly one
catalog position with matching multiplicity and stabilizer order has a unique
symbol; otherwise it is `None`. Full coordinate matching follows in Task 6.

- [ ] **Step 4: Emit reported-multiplicity diagnostics**

When `site.reported_multiplicity` is non-null and differs from the calculated
value, append:

```python
Diagnostic(
    Severity.WARNING,
    "crystallography.orbit.reported_multiplicity_mismatch",
    f"Site {site.label!r} reports multiplicity {reported}, "
    f"but its symmetry orbit contains {calculated} positions",
)
```

- [ ] **Step 5: Run focused tests**

Run:

```bash
pytest tests/crystallography/test_orbit_analyzer.py tests/symmetry/test_orbit.py -q
```

Expected: all selected tests pass; low-level expansion behavior is unchanged.

- [ ] **Step 6: Commit**

```bash
git add src/cristma/crystallography tests/crystallography/test_orbit_analyzer.py
git commit -m "feat: calculate crystallographic orbits and stabilizers"
```

---

### Task 6: Exact Wyckoff matching and diagnostics

**Files:**
- Modify: `src/cristma/crystallography/orbit.py`
- Modify: `tests/crystallography/test_orbit_analyzer.py`

**Interfaces:**
- Consumes: Task 5 orbit/stabilizer result and Task 3 exact Wyckoff representatives.
- Produces: deterministic `wyckoff_position`, validated `SiteSymmetry.symbol`, unresolved/ambiguous/mismatch diagnostics.

- [ ] **Step 1: Add failing Wyckoff fixtures**

```python
@pytest.mark.parametrize(
    ("fractional", "letter", "multiplicity", "site_symmetry"),
    [
        ((0.0, 0.0, 0.0), "a", 2, "-4.."),
        ((0.0, 0.0, 0.5), "b", 2, "-4.."),
        ((0.0, 0.5, 0.23), "c", 2, "2.mm"),
        ((0.0, 0.0, 0.23), "d", 4, "2.."),
        ((0.17, 0.67, 0.23), "e", 4, "..m"),
        ((0.17, 0.29, 0.23), "f", 8, "1"),
    ],
)
def test_p421m_wyckoff_positions_are_identified(
    fractional: tuple[float, float, float],
    letter: str,
    multiplicity: int,
    site_symmetry: str,
) -> None:
    record = SpaceGroupCatalog.default().by_hall_number(390)
    result = OrbitAnalyzer(tolerance=1e-6).analyze(
        site_at(*fractional),
        record.definition(provenance="derived"),
        catalog_record=record,
        cell=tetragonal_cell(),
    )
    assert result.wyckoff_position is not None
    assert result.wyckoff_position.letter == letter
    assert result.calculated_multiplicity == multiplicity
    assert result.site_symmetry.symbol == site_symmetry


def test_rounded_special_coordinate_respects_explicit_tolerance() -> None:
    record = SpaceGroupCatalog.default().by_hall_number(390)
    close = analyze(record, (0.0, 0.5, 0.2300004), tolerance=1e-6)
    assert close.wyckoff_position.letter == "c"


def test_reported_wyckoff_mismatch_is_not_silently_rewritten() -> None:
    result = analyze_p421m(site_at(0, 0, 0, wyckoff="b"))
    assert result.wyckoff_position.letter == "a"
    assert "crystallography.orbit.reported_wyckoff_mismatch" in {
        item.code for item in result.diagnostics
    }
```

Add synthetic tests where two candidate representatives match and ensure the
result is `None` with `crystallography.orbit.wyckoff_ambiguous`, plus a no-match
case with `crystallography.orbit.wyckoff_unresolved`.

- [ ] **Step 2: Run tests and verify failure**

Run:

```bash
pytest tests/crystallography/test_orbit_analyzer.py -q
```

Expected: P-42₁m letter assertions fail because only coarse multiplicity and
stabilizer-order filtering exists.

- [ ] **Step 3: Implement affine-subspace matching**

For each representative map `p -> A p + t`, determine whether the measured
fractional coordinate `x` belongs to its periodic affine image by minimizing:

```text
|| A p + t - x - n ||∞
```

over integer translations `n` required by the three wrapped coordinate
differences. Solve the small linear least-squares system with NumPy for each
candidate translation, then verify the reconstructed coordinate using periodic
component distance and `self.tolerance`. Do not infer a Wyckoff match from
multiplicity alone.

After a representative matches, expand the representative with the selected
group operations and compare its normalized position set with the actual
`ExpandedAtom.fractional` set using a one-to-one periodic matching. A candidate
passes only if both membership and complete orbit equality pass.

Outcomes:

```python
if len(matches) == 1:
    wyckoff_position = matches[0]
elif not matches:
    wyckoff_position = None
    diagnostics += (Diagnostic(Severity.WARNING, "crystallography.orbit.wyckoff_unresolved", ...),)
else:
    wyckoff_position = None
    diagnostics += (Diagnostic(Severity.WARNING, "crystallography.orbit.wyckoff_ambiguous", ...),)
```

When matched, set `SiteSymmetry.symbol` from the catalog record and retain the
computed stabilizer. Compare `site.wyckoff`, reported multiplicity and the
catalog multiplicity independently; never overwrite source-reported values.

If a reported Wyckoff letter exists, resolve that letter in the selected
catalog record before comparing it with the calculated match. Emit
`crystallography.orbit.wyckoff_multiplicity_mismatch` when the reported
letter's catalog multiplicity differs from the calculated orbit. Emit
`crystallography.orbit.reported_wyckoff_mismatch` when the unique calculated
letter differs. If a representative passes coordinate membership but its
expanded stabilizer fails the catalog orbit/stabilizer consistency check, emit
`crystallography.orbit.site_symmetry_mismatch` and do not accept it.

- [ ] **Step 4: Run focused tests**

Run:

```bash
pytest tests/crystallography/test_orbit_analyzer.py tests/crystallography/test_catalog.py -q
```

Expected: all selected tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/cristma/crystallography/orbit.py tests/crystallography/test_orbit_analyzer.py
git commit -m "feat: match Wyckoff positions from exact orbits"
```

---

### Task 7: Catalog-aware CIF mapping and structure analysis

**Files:**
- Modify: `src/cristma/io/cif/mapper.py`
- Modify: `tests/io/cif/test_mapper_advanced.py`
- Create: `tests/crystallography/test_cif_integration.py`
- Modify: `src/cristma/crystallography/orbit.py`
- Modify: `src/cristma/crystallography/__init__.py`

**Interfaces:**
- Consumes: canonical CIF metadata and Task 4 catalog.
- Produces: catalog-resolved `SpaceGroupDefinition` when operations are absent and `analyze_structure_orbits(crystal, *, catalog=None, tolerance=1e-6) -> tuple[CrystallographicOrbit, ...]`.

- [ ] **Step 1: Write failing CIF resolution tests**

```python
def test_cif_without_operations_resolves_exact_hall_catalog_entry() -> None:
    result = cristma.read_text(
        """
data_demo
_cell_length_a 5
_cell_length_b 5
_cell_length_c 5
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_Hall 'P -4 2ab'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0 0 0
""",
        format="cif",
    )
    crystal = result.structures[0]
    assert crystal.space_group.provenance == "derived"
    assert crystal.space_group.number == 113
    assert len(crystal.space_group.operations) == 8
    assert "cif.map.symmetry_operations_derived" in {d.code for d in result.diagnostics}


def test_ambiguous_hm_symbol_does_not_guess_setting() -> None:
    result = read_ambiguous_group_5_without_operations()
    assert result.structures[0].space_group.provenance == "identity_fallback"
    assert "cif.map.space_group_lookup_ambiguous" in {d.code for d in result.diagnostics}
```

Also test exact Hall metadata contradicting reported operations. Preserve the
reported operations and emit `cif.map.space_group_operations_mismatch`; do not
replace explicit source data silently.

- [ ] **Step 2: Run CIF tests and verify failure**

Run:

```bash
pytest tests/io/cif/test_mapper_advanced.py tests/crystallography/test_cif_integration.py -q
```

Expected: the missing-operation case still uses identity fallback.

- [ ] **Step 3: Inject an optional catalog at the mapper boundary**

Use signatures:

```python
def map_cif_structures(
    document: CifDocument,
    *,
    crystallography: SpaceGroupCatalog | None = None,
) -> tuple[tuple[CrystalStructure, ...], tuple[Diagnostic, ...]]: ...
```

`None` selects `SpaceGroupCatalog.default()`, so the existing built-in CIF
handler needs no format-specific dependency injection. Tests may pass a small
catalog explicitly. Resolution order when operations are absent is:

1. exact Hall symbol;
2. IT number plus explicit setting/origin choice;
3. uniquely resolved HM symbol;
4. identity fallback with the existing warning.

When catalog resolution succeeds, use
`record.definition(provenance="derived")` and emit an INFO diagnostic. If an
alias yields multiple settings, emit WARNING and retain identity fallback.

When operations are explicit, keep them authoritative. If catalog identity is
also unambiguous, compare normalized exact operation sets and emit a warning on
mismatch.

Replace the mapper's existing direct `expand_orbit()` multiplicity loop with
`OrbitAnalyzer`. Append its diagnostics to the CIF result, and build a new
`IndependentSite` snapshot with `calculated_multiplicity` set from the result.
This removes the second high-level multiplicity path while preserving the
immutable structure contract.

- [ ] **Step 4: Add whole-structure orbit analysis**

```python
def analyze_structure_orbits(
    crystal: CrystalStructure,
    *,
    catalog: SpaceGroupCatalog | None = None,
    tolerance: float = 1e-6,
) -> tuple[CrystallographicOrbit, ...]:
    selected = SpaceGroupCatalog.default() if catalog is None else catalog
    record = selected.resolve_definition(crystal.space_group)
    analyzer = OrbitAnalyzer(tolerance=tolerance)
    return tuple(
        analyzer.analyze(
            site,
            crystal.space_group,
            cell=crystal.cell,
            catalog_record=record,
            structure_id=crystal.id,
        )
        for site in crystal.sites
    )
```

`resolve_definition()` requires exact Hall identity or a unique fully qualified
match. It raises `LookupError` for unresolved/ambiguous definitions rather than
guessing.

- [ ] **Step 5: Test the real gehlenite fixture and geometry compatibility**

Read `/Users/artem/Desktop/-_B2_Ba_O4_-.cif` only if it is available; the
repository test must use a compact committed CIF fixture containing P-42₁m
metadata and representative gehlenite sites. Assert:

```python
crystal = cristma.read(FIXTURE).structures[0]
orbits = analyze_structure_orbits(crystal)
assert [orbit.calculated_multiplicity for orbit in orbits]
assert all(orbit.wyckoff_position is not None for orbit in orbits)

view = crystal.atomic_view()
graph = NeighborFinder(cutoff=3.0).find(view)
assert len(graph.atoms) == sum(orbit.calculated_multiplicity for orbit in orbits)
```

This verifies the existing periodic geometry consumes the same expanded
positions without adding chemical bond semantics.

- [ ] **Step 6: Run focused integration tests**

Run:

```bash
pytest tests/io/cif/test_mapper_advanced.py tests/crystallography tests/symmetry/test_orbit.py tests/geometry -q
```

Expected: all selected tests pass.

- [ ] **Step 7: Commit**

```bash
git add src/cristma/io/cif/mapper.py src/cristma/crystallography tests/io/cif/test_mapper_advanced.py tests/crystallography
git commit -m "feat: resolve CIF symmetry and validate site orbits"
```

---

### Task 8: Final package, attribution and clean-install verification

**Files:**
- Modify: `README.md`
- Modify: `docs/superpowers/specs/2026-09-01-cristma-structural-crystallography-design.md` only if implementation revealed a factual mismatch
- Test: all focused and full suites

**Interfaces:**
- Consumes: Tasks 1-7.
- Produces: installable wheel with complete catalog resources and documented public usage.

- [ ] **Step 1: Add concise README usage**

Document:

```python
import cristma
from cristma.crystallography import (
    SpaceGroupCatalog,
    analyze_structure_orbits,
)

crystal = cristma.read("sample.cif").structures[0]
group = SpaceGroupCatalog.default().by_hall_symbol("P -4 2ab")
orbits = analyze_structure_orbits(crystal)
```

State explicitly that spglib is a data-compilation dependency, not a runtime
dependency, and link to packaged `SOURCE.md`/license.

- [ ] **Step 2: Run all focused crystallography tests**

Run:

```bash
pytest tests/tools/test_compile_spglib_crystallography.py tests/crystallography tests/symmetry tests/io/cif tests/geometry tests/reference_data -q
```

Expected: all selected tests pass.

- [ ] **Step 3: Run the full suite once**

Run:

```bash
pytest -q
```

Expected: all tests pass with no failures or errors.

- [ ] **Step 4: Build the wheel**

Run:

```bash
python -m build
```

Expected: sdist and wheel are created successfully.

- [ ] **Step 5: Verify wheel contents**

Run:

```bash
python -m zipfile -l dist/cristma-0.1.0.dev0-py3-none-any.whl
```

Expected: both JSON resources, `SOURCE.md`, and `SPGLIB_LICENSE.txt` appear in
the wheel.

- [ ] **Step 6: Install and smoke-test in a clean environment**

Run:

```bash
python -m venv /private/tmp/cristma-crystallography-wheel
/private/tmp/cristma-crystallography-wheel/bin/python -m pip install --no-deps dist/cristma-0.1.0.dev0-py3-none-any.whl
/private/tmp/cristma-crystallography-wheel/bin/python -c "from cristma.crystallography import SpaceGroupCatalog; c=SpaceGroupCatalog.default(); assert len(c)==530; assert c.by_hall_number(390).number==113"
/private/tmp/cristma-crystallography-wheel/bin/python -c "import importlib.metadata as m; assert m.requires('cristma') == ['numpy>=1.26']"
```

Expected: both commands exit 0; no spglib, pymatgen or Gemmi dependency is
reported.

- [ ] **Step 7: Inspect repository state and commit**

Run:

```bash
git status --short
git diff --check
```

Do not add user-owned files. Then commit:

```bash
git add README.md docs/superpowers/specs/2026-09-01-cristma-structural-crystallography-design.md
git commit -m "docs: document crystallographic catalog and orbit API"
```

If the documentation was already committed and unchanged, skip the empty
commit.
