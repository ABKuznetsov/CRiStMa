# Native SHELX Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a native, loss-preserving, Qt-free SHELX RES/INS reader and writer that maps structural content to CRiStMa's canonical `CrystalStructure` without external crystallography packages.

**Architecture:** A line-preserving format document is parsed independently from scientific mapping. `ShelxDocument` owns source fidelity; mapper functions produce immutable canonical structures; registry handlers expose reading; explicit writer options select preserve or canonical output. No application workflow or CRAFT types enter CRiStMa.

**Tech Stack:** Python 3.11+, standard library, NumPy, pytest, existing CRiStMa core/symmetry/io contracts.

**Spec:** `docs/superpowers/specs/2026-08-30-native-shelx-and-craft-integration-design.md`

## Global Constraints

- [ ] Work in an isolated CRiStMa git worktree; inspect `git status` before every commit.
- [ ] Follow strict TDD: add one focused failing test, run it to observe the intended failure, implement the smallest scientific behavior, rerun the focused test.
- [ ] Preserve exact source text in preserve mode, including CRLF, blank lines, continuations, unknown records, and records after `END`.
- [ ] Keep CRiStMa independent of CRAFT, Sci, Qt, Gemmi, pymatgen, and SHELX executables.
- [ ] Do not model refinement/restraint instructions as canonical constraints in this slice.
- [ ] Keep version `0.1.0.dev0` and build only an internal wheel; do not publish to PyPI in this plan. Public release waits for proven use by CRAFT and Finder.
- [ ] Run only focused tests within tasks. Run the complete CRiStMa suite once at the final gate.
- [ ] Use the real fixture `tests/fixtures/shelx/zdk288.res`; do not mutate it.

---

## Task 1: Loss-preserving physical document and record assembly

**Files:**

- Create: `src/cristma/io/shelx/__init__.py`
- Create: `src/cristma/io/shelx/document.py`
- Create: `src/cristma/io/shelx/parser.py`
- Create: `tests/io/shelx/test_document.py`
- Create: `tests/io/shelx/test_parser.py`
- Create: `tests/io/shelx/test_preserve_writer.py`

- [ ] Add a failing test proving that parsing and preserve rendering return byte-for-byte-equivalent decoded text for LF and CRLF sources, including blank lines and material after `END`.

```python
source = "TITL demo\r\nREM keep  two spaces\r\n\r\nEND\r\nQ1 1 0 0 0 11 0.05\r\n"
result = parse_shelx(source, source_name="demo.res")
assert result.document.render_preserved() == source
```

- [ ] Run `pytest tests/io/shelx/test_document.py -q` and confirm failure because the module/types do not exist.
- [ ] Implement immutable source-layer types:

```python
@dataclass(frozen=True, slots=True)
class ShelxPhysicalLine:
    text: str
    newline: str
    span: SourceSpan

@dataclass(frozen=True, slots=True)
class ShelxRecord:
    keyword: str | None
    fields: tuple[str, ...]
    physical_line_indices: tuple[int, ...]
    span: SourceSpan

@dataclass(frozen=True, slots=True)
class ShelxSourceEdit:
    span: SourceSpan
    replacement: str

@dataclass(frozen=True, slots=True)
class ShelxDocument:
    raw_source: str
    physical_lines: tuple[ShelxPhysicalLine, ...]
    records: tuple[ShelxRecord, ...]
    source_name: str | None = None
    edits: tuple[ShelxSourceEdit, ...] = ()

    def render_preserved(self) -> str:
        return apply_source_edits(self.raw_source, self.edits)
```

- [ ] Assemble logical records across physical lines ending in `=` while retaining every physical line and its offsets.
- [ ] Classify blank, `REM`, known instruction, atom-like, Q-peak, and unknown records without yet mapping scientific values.
- [ ] Add immutable `ShelxSourceEdit(span, replacement)` support. Reject overlapping or out-of-range edits and render valid edits from highest to lowest source offset so untouched text remains exact.
- [ ] Add tests for inline `!` comments, case-insensitive instruction recognition, continuation spans, and records before/after `HKLF` and `END`.
- [ ] Add tests proving a single supported source edit changes only its requested span and that overlapping edits fail explicitly.
- [ ] Run `pytest tests/io/shelx/test_document.py tests/io/shelx/test_parser.py tests/io/shelx/test_preserve_writer.py -q`.
- [ ] Commit: `git add src/cristma/io/shelx tests/io/shelx && git commit -m "feat(io): add loss-preserving SHELX document parser"`

## Task 2: Typed SHELX instructions and exact symmetry construction

**Files:**

- Create: `src/cristma/io/shelx/records.py`
- Create: `src/cristma/io/shelx/symmetry.py`
- Modify: `src/cristma/io/shelx/parser.py`
- Create: `tests/io/shelx/test_records.py`
- Create: `tests/io/shelx/test_symmetry.py`

- [ ] Add failing tests for typed `CELL`, `ZERR`, `LATT`, `SYMM`, `SFAC`, `FVAR`, `PART`, `RESI`, `HKLF`, and `END` records. Assert that `CELL` retains wavelength separately from its `UnitCell`.
- [ ] Run the focused tests and confirm the expected missing-type failures.
- [ ] Implement typed records, including:

```python
@dataclass(frozen=True, slots=True)
class ShelxCellInstruction(ShelxRecord):
    wavelength: MeasuredValue
    cell: UnitCell

@dataclass(frozen=True, slots=True)
class ShelxLattInstruction(ShelxRecord):
    code: int

@dataclass(frozen=True, slots=True)
class ShelxSymmInstruction(ShelxRecord):
    operation: AffineOperation
```

- [ ] Implement exact rational SHELX symmetry expansion from `LATT` and explicit `SYMM` using existing `AffineOperation` primitives. Cover centring codes `1=P`, `2=I`, `3=R`, `4=F`, `5=A`, `6=B`, `7=C`; positive `LATT` adds inversion, negative does not.
- [ ] Deduplicate operations by exact affine coefficients and preserve deterministic order with identity first.
- [ ] Add tests for all seven centring codes, sign convention, duplicate explicit operations, rational translations, and a known `P2(1)/n` example.
- [ ] Run `pytest tests/io/shelx/test_records.py tests/io/shelx/test_symmetry.py -q`.
- [ ] Commit: `git add src/cristma/io/shelx tests/io/shelx && git commit -m "feat(symmetry): map SHELX LATT and SYMM exactly"`

## Task 3: SFAC, FVAR, and symbolic occupancy

**Files:**

- Create: `src/cristma/io/shelx/occupancy.py`
- Modify: `src/cristma/io/shelx/records.py`
- Modify: `src/cristma/io/shelx/parser.py`
- Create: `tests/io/shelx/test_occupancy.py`
- Create: `tests/io/shelx/test_sfac.py`

- [ ] Add failing table-driven tests for fixed occupancy, positive FVAR references, complemented negative references, absent FVAR indices, non-finite values, and evaluated values outside `[0, 1]`.
- [ ] Run `pytest tests/io/shelx/test_occupancy.py -q` and confirm failure.
- [ ] Implement a value object whose dependency is never discarded:

```python
@dataclass(frozen=True, slots=True)
class ShelxOccupancyExpression:
    raw: str
    free_variable_index: int | None
    multiplier: float
    complement: bool

    @classmethod
    def parse(cls, token: str) -> "ShelxOccupancyExpression":
        """Decode the control number and multiplier without evaluating it."""

    def evaluate(self, free_variables: tuple[float, ...]) -> float:
        """Return the physical occupancy or raise for an invalid dependency."""
```

- [ ] Encode SHELX occupancy semantics in one module only: values with control part `1` are fixed multipliers; higher control parts address `FVAR[control-1]`; negative controls evaluate the complement form. Keep parsing and evaluation separately testable.
- [ ] Parse both `SFAC C H O` element lists and one-element coefficient records such as `SFAC O 3.05 ...`; normalize ionic/special labels to CRiStMa `Species` without importing pymatgen.
- [ ] Add source-spanned diagnostics for invalid SFAC indices and occupancy expressions.
- [ ] Run `pytest tests/io/shelx/test_occupancy.py tests/io/shelx/test_sfac.py -q`.
- [ ] Commit: `git add src/cristma/io/shelx tests/io/shelx && git commit -m "feat(io): preserve SHELX occupancy dependencies"`

## Task 4: Scientific mapping to canonical CrystalStructure

**Files:**

- Create: `src/cristma/io/shelx/mapper.py`
- Modify: `src/cristma/io/shelx/records.py`
- Create: `tests/io/shelx/test_mapper.py`
- Create: `tests/io/shelx/test_displacement.py`
- Create: `tests/io/shelx/test_diagnostics.py`

- [ ] Add a failing mapping test for a minimal structure and assert cell, title, elements, evaluated occupancies, exact symmetry, and stable site IDs.
- [ ] Add failing tests proving that Q peaks and atom-looking records after `END` do not become sites.
- [ ] Run the mapper tests and confirm failure because mapping is absent.
- [ ] Implement:

```python
structures, diagnostics = map_shelx_structures(document)
assert isinstance(structures, StructureCollection)
assert all(isinstance(item, Diagnostic) for item in diagnostics)
```

- [ ] Map one-value displacement records to `U_iso`; map six-value records in SHELX order `U11 U22 U33 U23 U13 U12` to a symmetric CRiStMa tensor.
- [ ] Retain `ShelxOccupancyExpression` in `SiteComponent.metadata` while storing its evaluated occupancy in the canonical component.
- [ ] Track active `PART` and `RESI` context in site metadata/disorder fields. Preserve `AFIX`, `HFIX`, `EXYZ`, `EADP`, `DFIX`, `DANG`, `SADI`, `SAME`, `FLAT`, `RIGU`, `SIMU`, `DELU`, weighting commands, and unknown records only in `ShelxDocument`.
- [ ] Return errors with spans for missing/invalid `CELL`, unresolved elements, invalid SFAC references, invalid FVAR references, and unphysical occupancy; never construct an invalid canonical site.
- [ ] Add tests proving anisotropic tensors rotate through existing symmetry expansion rather than being copied numerically.
- [ ] Run `pytest tests/io/shelx/test_mapper.py tests/io/shelx/test_displacement.py tests/io/shelx/test_diagnostics.py -q`.
- [ ] Commit: `git add src/cristma/io/shelx tests/io/shelx && git commit -m "feat(structure): map SHELX records to canonical crystals"`

## Task 5: Format registry and top-level package reading API

**Files:**

- Create: `src/cristma/io/shelx/probe.py`
- Create: `src/cristma/io/shelx/handler.py`
- Modify: `src/cristma/io/formats.py`
- Modify: `src/cristma/io/shelx/__init__.py`
- Modify: `tests/io/test_builtin_formats.py`
- Create: `tests/io/shelx/test_handler.py`

- [ ] Add failing tests that `.res`, `.ins`, explicit `format="shelx"`, and content-only `read_text` select the lazy SHELX handler, while arbitrary text does not beat CIF or produce an ambiguous match.
- [ ] Run `pytest tests/io/test_builtin_formats.py tests/io/shelx/test_handler.py -q` and confirm the missing descriptor failure.
- [ ] Add `_shelx_handler()` and a descriptor with aliases `("res", "ins")`, suffixes `(".res", ".ins")`, `multiple=False`, and a conservative structural-content probe requiring a credible combination such as `CELL` plus `LATT`/`SFAC`.
- [ ] Have `ShelxFormatHandler.read_text()` return the same `ReadResult` contract as CIF: typed document, canonical `StructureCollection`, combined diagnostics, and source information supplied by `FormatRegistry`.
- [ ] Verify lazy selection does not import the SHELX parser before its descriptor is selected.
- [ ] Run `pytest tests/io/test_builtin_formats.py tests/io/shelx/test_handler.py tests/io/test_registry.py -q`.
- [ ] Commit: `git add src/cristma/io tests/io && git commit -m "feat(io): register native SHELX reader"`

## Task 6: Preserve and canonical SHELX writers

**Files:**

- Create: `src/cristma/io/shelx/writer.py`
- Modify: `src/cristma/io/shelx/__init__.py`
- Modify: `src/cristma/__init__.py`
- Create: `tests/io/shelx/test_writer.py`
- Modify: `tests/test_public_api.py`

- [ ] Add failing public-API tests for exact document preservation and canonical crystal output:

```python
cristma.write(result.document, copy_path, mode="preserve")
cristma.write(
    crystal,
    ins_path,
    format="shelx",
    mode="canonical",
    options=ShelxWriteOptions(wavelength=0.71073),
)
```

- [ ] Assert that canonical writing without wavelength fails with a precise `ValueError`, because wavelength is not structural state.
- [ ] Run the focused writer/API tests and confirm failure.
- [ ] Implement `ShelxWriteOptions` and two explicit renderers:

```python
preserved_text = write_shelx_document(document, mode="preserve")
canonical_text = write_crystal_shelx(crystal, options=options)
```

- [ ] Extend `cristma.write()` to dispatch by document type or explicit `format`, while preserving existing CIF defaults and errors. Do not create a generic mutable writer manager.
- [ ] Canonical output must contain deterministic `TITL`, `CELL`, `LATT`/`SYMM`, `SFAC`, `UNIT`, required `FVAR`, atom records, `HKLF`, and `END`. Compute `UNIT` from expanded unit-cell contents.
- [ ] Add canonical read-write-read tests comparing canonical cells, symmetry operations, sites, components, occupancies, coordinates, and displacement parameters within reported numeric precision.
- [ ] Run `pytest tests/io/shelx/test_writer.py tests/test_public_api.py tests/io/cif/test_writer.py -q`.
- [ ] Commit: `git add src/cristma tests && git commit -m "feat(io): write preserved and canonical SHELX"`

## Task 7: Real fixture through Structure Core geometry

**Files:**

- Create: `tests/io/shelx/test_real_fixture.py`
- Create: `tests/integration/test_structure_core_shelx.py`
- Modify only if a demonstrated bug requires it: `src/cristma/io/shelx/*.py`

- [ ] Add assertions against `tests/fixtures/shelx/zdk288.res`: CRLF preserved, continuation assembled, title/cell/wavelength correct, expected independent-site count, Q peaks present only in the document, anisotropic atoms mapped, and content after `END` retained.
- [ ] Add an integration test for the actual scientific path:

```python
result = cristma.read(FIXTURE)
crystal = result.structures[0]
view = crystal.atomic_view()
graph = NeighborFinder(cutoff=3.0).find(view)
coordination = CoordinationAnalyzer().analyze(view, graph)
assert result.ok
assert len(view.atoms) > len(crystal.sites)
assert coordination.environments
```

- [ ] Run `pytest tests/io/shelx/test_real_fixture.py tests/integration/test_structure_core_shelx.py -q`.
- [ ] If a failure appears, add the smallest analytic regression test beside the responsible module before changing implementation.
- [ ] Commit: `git add src/cristma tests && git commit -m "test(io): verify real SHELX structure workflow"`

## Task 8: Package, documentation, and final scientific gate

**Files:**

- Modify: `README.md`
- Modify: `docs/formats.md` if present; otherwise create `docs/shelx.md`
- Modify: `src/cristma/io/shelx/__init__.py`
- Modify: `src/cristma/__init__.py`

- [ ] Document supported SHELX semantics, preserve versus canonical writing, wavelength ownership, Q-peak exclusion, retained-but-unmapped instructions, and the absence of external crystallography dependencies.
- [ ] Review package exports so the future-facing names are intentional: `ShelxDocument`, `ShelxWriteOptions`, `ShelxOccupancyExpression`, reader/writer helpers; keep parser internals private. This is not authorization to publish the package.
- [ ] Run a placeholder scan: `rg -n "TODO|FIXME|NotImplemented|pass$" src/cristma/io/shelx tests/io/shelx tests/integration/test_structure_core_shelx.py` and resolve every hit introduced by this plan.
- [ ] Run focused SHELX and integration tests once: `pytest tests/io/shelx tests/integration/test_structure_core_shelx.py -q`.
- [ ] Run the complete CRiStMa suite once: `pytest -q`.
- [ ] Build the wheel: `python -m build`.
- [ ] Install the wheel into a fresh temporary virtual environment and run a smoke script that imports CRiStMa, reads `zdk288.res`, and preserve-writes its document.
- [ ] Confirm `pyproject.toml` still reports `0.1.0.dev0` and perform no upload or release action.
- [ ] Review the diff against the specification and confirm no CRAFT/Sci/Qt imports, no cross-structure comparison classes, and no accidental canonical constraint mapping.
- [ ] Commit: `git add README.md docs src tests && git commit -m "docs(io): complete native SHELX structure slice"`
