# Native CIF Foundation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an installable, Qt-free CrIStMa package that parses and writes CIF 1.1 natively, preserves source information, and maps structural blocks to an asymmetric-unit `Crystal` model with exact symmetry provenance.

**Architecture:** Parsing is split from scientific interpretation: source text becomes a loss-preserving `CifDocument`, then a mapper creates canonical CrIStMa objects. Exact affine symmetry operations use rational coefficients, while numerical cell and coordinate operations use NumPy. Preserve-mode writing edits source spans; canonical writing emits a normalized CIF from `Crystal`.

**Tech Stack:** Python 3.11+, NumPy, pytest, standard-library dataclasses/enums/fractions/pathlib/decimal; no mandatory Gemmi, pymatgen, PyXtal, CrysPy, GSAS-II, or Qt.

**Spec:** `docs/plans/2026-08-30-structure-formats-design.md`

## Global Constraints

- Ordinary CIF reading and writing must work without specialized crystallographic packages.
- The primary structural representation is an asymmetric unit; expanded atoms are derived and retain provenance.
- Parsing and semantic mapping remain separate public operations.
- Unknown CIF tags, loops, comments, ordering, and original numeric tokens survive preserve-mode round trips.
- Scientific recovery is never silent; diagnostics carry stable codes and source spans.
- CIF missing states `absent`, `?`, and `.` remain distinguishable.
- Exact symmetry uses rational affine operations, never canonical floating-point translations.
- No Qt type or application service enters the package.
- TDD is required; run only the focused test named by each step before the final slice suite.

---

## Planned file structure

```text
pyproject.toml                         package metadata and dependencies
README.md                              minimal public introduction
src/cristma/__init__.py                stable top-level API
src/cristma/core/values.py             measured/missing scientific values
src/cristma/core/cell.py               validated unit cell and metric
src/cristma/core/structure.py          canonical asymmetric-unit model
src/cristma/symmetry/affine.py         exact affine symmetry operations
src/cristma/symmetry/orbit.py          derived sites with provenance
src/cristma/io/diagnostics.py          source spans and diagnostics
src/cristma/io/result.py               read result contract
src/cristma/io/registry.py             format probing and dispatch
src/cristma/io/cif/tokens.py           CIF token types
src/cristma/io/cif/lexer.py            CIF 1.1 lexical analysis
src/cristma/io/cif/document.py         loss-preserving CIF AST
src/cristma/io/cif/parser.py           tokens to CIF document
src/cristma/io/cif/names.py            structural tag aliases
src/cristma/io/cif/mapper.py           CIF document to Crystal
src/cristma/io/cif/writer.py           preserve and canonical writers
tests/...                              focused unit and integration tests
```

The SHELX subsystem is deliberately excluded from this plan. It receives its
own plan after this slice passes real-file CIF round trips.

### Task 1: Installable package and diagnostic contracts

**Files:**
- Create: `pyproject.toml`
- Create: `README.md`
- Create: `src/cristma/__init__.py`
- Create: `src/cristma/io/__init__.py`
- Create: `src/cristma/io/diagnostics.py`
- Create: `src/cristma/io/result.py`
- Test: `tests/io/test_diagnostics.py`
- Test: `tests/test_public_api.py`

**Interfaces:**
- Consumes: none.
- Produces: `Severity`, `SourcePosition`, `SourceSpan`, `Diagnostic`, `SourceInfo`, and `ReadResult`; package version and importable `cristma` namespace.

- [ ] **Step 1: Write package and diagnostic contract tests**

```python
from cristma import __version__
from cristma.io.diagnostics import Diagnostic, Severity, SourcePosition, SourceSpan
from cristma.io.result import ReadResult, SourceInfo


def test_public_package_has_semantic_version():
    assert __version__ == "0.1.0"


def test_diagnostic_carries_stable_code_and_source_span():
    span = SourceSpan(
        start=SourcePosition(offset=7, line=2, column=3),
        end=SourcePosition(offset=12, line=2, column=8),
    )
    diagnostic = Diagnostic(Severity.WARNING, "cif.loop.width", "short row", span)
    result = ReadResult(document=None, diagnostics=(diagnostic,))
    assert result.ok
    assert result.diagnostics[0].code == "cif.loop.width"
    assert result.diagnostics[0].span == span


def test_error_diagnostic_makes_result_not_ok():
    diagnostic = Diagnostic(Severity.ERROR, "io.empty", "empty source")
    assert not ReadResult(document=None, diagnostics=(diagnostic,)).ok


def test_read_result_preserves_decoding_and_newline_metadata():
    info = SourceInfo(name="sample.cif", format="cif", encoding="utf-8-sig", newline="\r\n")
    assert ReadResult(document=None, source_info=info).source_info == info
```

- [ ] **Step 2: Run the tests and verify collection fails**

Run: `pytest -q tests/io/test_diagnostics.py tests/test_public_api.py`

Expected: FAIL during import because `cristma` does not exist.

- [ ] **Step 3: Add packaging and minimal contracts**

Use a `src` layout, Python `>=3.11`, runtime dependency `numpy>=1.26`, and test dependency `pytest>=8`. Implement frozen, slotted dataclasses and this result rule:

```python
class Severity(str, Enum):
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


@dataclass(frozen=True, slots=True)
class SourcePosition:
    offset: int
    line: int
    column: int


@dataclass(frozen=True, slots=True)
class SourceSpan:
    start: SourcePosition
    end: SourcePosition


@dataclass(frozen=True, slots=True)
class Diagnostic:
    severity: Severity
    code: str
    message: str
    span: SourceSpan | None = None
    recovery: str | None = None


@dataclass(frozen=True, slots=True)
class SourceInfo:
    name: str | None = None
    format: str | None = None
    encoding: str = "utf-8"
    newline: str = "\n"


@dataclass(frozen=True, slots=True)
class ReadResult:
    document: object | None
    structures: tuple[object, ...] = ()
    diagnostics: tuple[Diagnostic, ...] = ()
    source_info: SourceInfo | None = None

    @property
    def ok(self) -> bool:
        return not any(item.severity is Severity.ERROR for item in self.diagnostics)
```

- [ ] **Step 4: Run the focused tests**

Run: `pytest -q tests/io/test_diagnostics.py tests/test_public_api.py`

Expected: `4 passed`.

- [ ] **Step 5: Commit the package foundation**

```bash
git add pyproject.toml README.md src/cristma tests/io/test_diagnostics.py tests/test_public_api.py
git commit -m "feat: establish cristma package contracts"
```

### Task 2: Scientific values and canonical asymmetric-unit model

**Files:**
- Create: `src/cristma/core/__init__.py`
- Create: `src/cristma/core/values.py`
- Create: `src/cristma/core/cell.py`
- Create: `src/cristma/core/structure.py`
- Test: `tests/core/test_values.py`
- Test: `tests/core/test_structure.py`

**Interfaces:**
- Consumes: diagnostic-independent standard Python and NumPy.
- Produces: `MissingKind`, `MeasuredValue`, `parse_measured_value`, `UnitCell`, `SiteComponent`, `DisplacementParameters`, `IndependentSite`, and `Crystal`.

- [ ] **Step 1: Write failing measured-value tests**

```python
from cristma.core.values import MissingKind, parse_measured_value


def test_parses_standard_uncertainty_at_last_digits():
    value = parse_measured_value("7.6959(2)", unit="angstrom")
    assert value.value == 7.6959
    assert value.uncertainty == 0.0002
    assert value.raw == "7.6959(2)"


def test_preserves_distinct_cif_missing_states():
    assert parse_measured_value("?").missing is MissingKind.UNKNOWN
    assert parse_measured_value(".").missing is MissingKind.INAPPLICABLE
    assert parse_measured_value(None).missing is MissingKind.ABSENT
```

- [ ] **Step 2: Run measured-value tests and verify failure**

Run: `pytest -q tests/core/test_values.py`

Expected: FAIL because `cristma.core.values` is absent.

- [ ] **Step 3: Implement measured values**

Implement a regex accepting signs, decimals, exponent notation, and a final parenthesized uncertainty. Scale uncertainty by the decimal exponent of the mantissa; never remove arbitrary parenthesized text from a token.

```python
@dataclass(frozen=True, slots=True)
class MeasuredValue:
    value: float | None
    uncertainty: float | None
    raw: str | None
    unit: str | None = None
    missing: MissingKind = MissingKind.PRESENT
```

- [ ] **Step 4: Write failing structure invariants tests**

```python
import pytest
from cristma.core.cell import UnitCell
from cristma.core.structure import Crystal, IndependentSite, SiteComponent
from cristma.core.values import MeasuredValue


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value=value, uncertainty=None, raw=str(value))


def test_cell_rejects_non_positive_edge():
    with pytest.raises(ValueError, match="cell edge"):
        UnitCell(number(0), number(4), number(5), number(90), number(90), number(90))


def test_crystal_keeps_asymmetric_sites_primary():
    site = IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(SiteComponent("Si", number(1.0)),),
        fractional=(number(0), number(0), number(0)),
    )
    crystal = Crystal(name="silicon", cell=UnitCell.cubic(number(5.43)), sites=(site,))
    assert crystal.sites == (site,)
    assert crystal.expanded_sites is None
```

- [ ] **Step 5: Implement canonical core dataclasses**

`UnitCell` validates positive edges, angles strictly between 0 and 180 degrees, and positive metric determinant. It exposes `matrix`, `metric`, `volume`, and `cubic(edge)`.

`SiteComponent` contains `element`, `occupancy: MeasuredValue`, optional `oxidation_state: MeasuredValue`, and source metadata. `IndependentSite` contains stable ID, label, components, fractional coordinates, optional reported and calculated Wyckoff multiplicity, disorder assembly/group, isotropic or anisotropic displacement data, and source metadata. Component occupancies must be finite and non-negative; totals above `1.0 + 1e-6` raise `ValueError`.

`Crystal` contains name, cell, asymmetric sites, optional space-group definition, formula, metadata, and an optional derived expansion cache excluded from equality.

- [ ] **Step 6: Run core tests**

Run: `pytest -q tests/core/test_values.py tests/core/test_structure.py`

Expected: all tests PASS.

- [ ] **Step 7: Commit the canonical model**

```bash
git add src/cristma/core tests/core
git commit -m "feat: add canonical asymmetric-unit model"
```

### Task 3: Exact affine symmetry and orbit provenance

**Files:**
- Create: `src/cristma/symmetry/__init__.py`
- Create: `src/cristma/symmetry/affine.py`
- Create: `src/cristma/symmetry/orbit.py`
- Test: `tests/symmetry/test_affine.py`
- Test: `tests/symmetry/test_orbit.py`

**Interfaces:**
- Consumes: `IndependentSite`, `MeasuredValue`.
- Produces: `AffineOperation`, `parse_xyz_operation(text, operation_id=None)`, `SpaceGroupDefinition`, `ExpandedSite`, and `expand_orbit(site, operations, tolerance=1e-8)`.

- [ ] **Step 1: Write failing exact-operation tests**

```python
from fractions import Fraction
from cristma.symmetry.affine import parse_xyz_operation


def test_parses_rational_affine_operation_exactly():
    op = parse_xyz_operation("-x+1/2, y+1/3, z")
    assert op.rotation[0] == (Fraction(-1), Fraction(0), Fraction(0))
    assert op.translation == (Fraction(1, 2), Fraction(1, 3), Fraction(0))
    assert op.apply_fractional((0.1, 0.2, 0.3)) == (0.4, 0.5333333333333333, 0.3)


def test_operation_normalizes_integer_translation():
    op = parse_xyz_operation("x+1,y-1,z")
    assert op.normalized().translation == (Fraction(0),) * 3
```

- [ ] **Step 2: Verify affine tests fail**

Run: `pytest -q tests/symmetry/test_affine.py`

Expected: FAIL because the module is absent.

- [ ] **Step 3: Implement the restricted affine-expression parser**

Parse only linear expressions composed of `x`, `y`, `z`, integer signs, and rational constants. Do not call `eval`. Reject products, powers, functions, repeated variables with non-integral coefficients, and malformed triplets with `ValueError` containing the source expression.

```python
@dataclass(frozen=True, slots=True)
class AffineOperation:
    rotation: tuple[tuple[Fraction, Fraction, Fraction], ...]
    translation: tuple[Fraction, Fraction, Fraction]
    source: str | None = None
    id: str | None = None
```

- [ ] **Step 4: Write failing orbit-provenance test**

```python
def test_orbit_deduplicates_special_position_and_keeps_operation_id():
    from cristma.core.structure import IndependentSite, SiteComponent
    from cristma.core.values import MeasuredValue

    def number(value: float) -> MeasuredValue:
        return MeasuredValue(value=value, uncertainty=None, raw=str(value))

    si_site = IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(SiteComponent("Si", number(1.0)),),
        fractional=(number(0), number(0), number(0)),
    )
    operations = (
        parse_xyz_operation("x,y,z", operation_id="op:1"),
        parse_xyz_operation("-x,-y,-z", operation_id="op:2"),
    )
    expanded = expand_orbit(si_site, operations)
    assert len(expanded) == 1
    assert expanded[0].independent_site_id == si_site.id
    assert expanded[0].equivalent_operation_ids == ("op:1", "op:2")
```

- [ ] **Step 5: Implement orbit expansion and space-group container**

`ExpandedSite` stores fractional coordinates, independent-site ID, representative operation ID, all equivalent operation IDs, and integer translation. Deduplication uses wrapped fractional coordinates and a configurable tolerance. `SpaceGroupDefinition` stores reported symbols/numbers, setting/origin, exact operations, and provenance (`reported`, `derived`, or `identity_fallback`).

- [ ] **Step 6: Run symmetry tests**

Run: `pytest -q tests/symmetry/test_affine.py tests/symmetry/test_orbit.py`

Expected: all tests PASS.

- [ ] **Step 7: Commit exact symmetry**

```bash
git add src/cristma/symmetry tests/symmetry
git commit -m "feat: add exact symmetry operations and orbit provenance"
```

### Task 4: Format registry and public read dispatch

**Files:**
- Create: `src/cristma/io/registry.py`
- Modify: `src/cristma/__init__.py`
- Test: `tests/io/test_registry.py`

**Interfaces:**
- Consumes: `ReadResult`.
- Produces: `FormatHandler` protocol, `FormatRegistry`, `read(path_or_text, format=None)`, and `register_format(handler)`.

- [ ] **Step 1: Write failing registry tests**

```python
from pathlib import Path
from cristma.io.registry import FormatRegistry


class StubHandler:
    name = "stub"
    suffixes = (".stub",)

    def probe(self, source: str) -> float:
        return 1.0 if source.startswith("STUB") else 0.0

    def read_text(self, source: str, source_name: str | None = None):
        return (source_name, source)


def test_registry_uses_content_probe_without_suffix(tmp_path: Path):
    path = tmp_path / "unknown.data"
    path.write_text("STUB value", encoding="utf-8")
    registry = FormatRegistry((StubHandler(),))
    assert registry.read(path) == (str(path), "STUB value")
```

- [ ] **Step 2: Verify registry test fails**

Run: `pytest -q tests/io/test_registry.py`

Expected: FAIL because `FormatRegistry` is absent.

- [ ] **Step 3: Implement deterministic probing and dispatch**

The registry accepts paths and explicit text sources. Explicit `format=` wins; otherwise combine suffix confidence and content probe, reject ties at the highest score, and include the candidates in the error message. Decode bytes as UTF-8/UTF-8-with-BOM; on invalid UTF-8, decode Latin-1 and add `io.encoding_fallback` warning. Populate `SourceInfo` with encoding and the detected `\n`, `\r\n`, or `\r` newline style.

- [ ] **Step 4: Run registry test**

Run: `pytest -q tests/io/test_registry.py`

Expected: PASS.

- [ ] **Step 5: Commit registry infrastructure**

```bash
git add src/cristma/__init__.py src/cristma/io/registry.py tests/io/test_registry.py
git commit -m "feat: add structure format registry"
```

### Task 5: CIF 1.1 lexer with source spans

**Files:**
- Create: `src/cristma/io/cif/__init__.py`
- Create: `src/cristma/io/cif/tokens.py`
- Create: `src/cristma/io/cif/lexer.py`
- Test: `tests/io/cif/test_lexer.py`

**Interfaces:**
- Consumes: `SourcePosition`, `SourceSpan`, `Diagnostic`.
- Produces: `CifTokenKind`, `CifToken`, and `lex_cif(source: str) -> tuple[tuple[CifToken, ...], tuple[Diagnostic, ...]]`.

- [ ] **Step 1: Write failing lexer tests**

```python
from cristma.io.cif.lexer import lex_cif
from cristma.io.cif.tokens import CifTokenKind


def test_lexer_preserves_quotes_comments_and_spans():
    source = "data_a\n_tag 'two words' # note\n"
    tokens, diagnostics = lex_cif(source)
    assert not diagnostics
    assert [token.kind for token in tokens] == [
        CifTokenKind.DATA,
        CifTokenKind.TAG,
        CifTokenKind.VALUE,
        CifTokenKind.COMMENT,
    ]
    assert tokens[2].value == "two words"
    assert tokens[2].raw == "'two words'"
    assert source[tokens[2].span.start.offset:tokens[2].span.end.offset] == "'two words'"


def test_lexer_reads_semicolon_text_only_from_column_one():
    source = "data_a\n_note\n;line one\nline two\n;\n"
    tokens, diagnostics = lex_cif(source)
    assert not diagnostics
    assert tokens[-1].value == "line one\nline two"
```

- [ ] **Step 2: Run lexer tests and verify failure**

Run: `pytest -q tests/io/cif/test_lexer.py`

Expected: FAIL because CIF lexer modules are absent.

- [ ] **Step 3: Implement a state-machine lexer**

Use explicit states for whitespace, comment, unquoted value, quoted value, and semicolon text. Classify `data_...`, `loop_`, `_tag`, `save_...`, `stop_`, `global_`, and ordinary values case-insensitively while retaining raw case. Unterminated quoted/multiline values produce error diagnostics with spans and stop only the affected token.

- [ ] **Step 4: Add malformed-token tests**

```python
def test_unterminated_quote_reports_location():
    tokens, diagnostics = lex_cif("data_a\n_tag 'broken\n")
    assert diagnostics[0].code == "cif.lex.unterminated_quote"
    assert diagnostics[0].span.start.line == 2
```

- [ ] **Step 5: Run focused lexer suite**

Run: `pytest -q tests/io/cif/test_lexer.py`

Expected: all tests PASS.

- [ ] **Step 6: Commit lexer**

```bash
git add src/cristma/io/cif tests/io/cif/test_lexer.py
git commit -m "feat: add native cif lexer"
```

### Task 6: Loss-preserving CIF document parser

**Files:**
- Create: `src/cristma/io/cif/document.py`
- Create: `src/cristma/io/cif/parser.py`
- Test: `tests/io/cif/test_parser.py`

**Interfaces:**
- Consumes: `lex_cif`, tokens, diagnostics.
- Produces: `CifScalar`, `CifLoop`, `CifBlock`, `CifDocument`, and `parse_cif(source: str, source_name=None) -> ReadResult`.

- [ ] **Step 1: Write failing parser tests**

```python
from cristma.io.cif.parser import parse_cif


def test_parser_keeps_multiple_blocks_scalars_and_loops():
    source = """data_first
_cell_length_a 5.0
loop_
_atom_site_label
_atom_site_fract_x
Si1 0
data_second
_audit_note ?
"""
    result = parse_cif(source)
    assert result.ok
    assert [block.name for block in result.document.blocks] == ["first", "second"]
    assert result.document.blocks[0].scalar("_CELL_LENGTH_A").value == "5.0"
    assert result.document.blocks[0].loops[0].rows == (("Si1", "0"),)
    assert result.document.blocks[1].scalar("_audit_note").value == "?"
```

- [ ] **Step 2: Verify parser test fails**

Run: `pytest -q tests/io/cif/test_parser.py::test_parser_keeps_multiple_blocks_scalars_and_loops`

Expected: FAIL because `parse_cif` is absent.

- [ ] **Step 3: Implement immutable document nodes**

Every scalar stores tag token and value token. Every loop stores ordered tag tokens and row value tokens. `CifDocument` stores raw source, source name, ordered blocks, comments, and pending source edits. Lookup normalizes tag names with `casefold()` but never changes stored spelling.

- [ ] **Step 4: Implement parser and loop-width recovery**

The parser requires a data block before data items, associates comments without discarding them, groups loop values by tag count, and emits `cif.parse.loop_width` if the final row is short. Preserve complete rows and the incomplete raw tokens; do not pad missing scientific values.

- [ ] **Step 5: Add loop-width and unknown-content tests**

```python
def test_short_loop_row_is_reported_without_fabricated_cell():
    result = parse_cif("data_a\nloop_\n_a\n_b\n1\n")
    assert not result.ok
    assert result.diagnostics[-1].code == "cif.parse.loop_width"
    assert result.document.blocks[0].loops[0].rows == ()


def test_parser_retains_unknown_tag_and_comment():
    source = "data_a\n# instrument note\n_local_detector_mode fast\n"
    result = parse_cif(source)
    assert result.document.blocks[0].scalar("_local_detector_mode").value == "fast"
    assert "# instrument note" in result.document.raw_source
```

- [ ] **Step 6: Run parser tests**

Run: `pytest -q tests/io/cif/test_parser.py`

Expected: all tests PASS.

- [ ] **Step 7: Commit document parser**

```bash
git add src/cristma/io/cif/document.py src/cristma/io/cif/parser.py tests/io/cif/test_parser.py
git commit -m "feat: parse loss-preserving cif documents"
```

### Task 7: CIF cell, symmetry, and atom-site mapping

**Files:**
- Create: `src/cristma/chemistry/__init__.py`
- Create: `src/cristma/chemistry/elements.py`
- Create: `src/cristma/io/cif/names.py`
- Create: `src/cristma/io/cif/mapper.py`
- Test: `tests/io/cif/test_mapper_basic.py`

**Interfaces:**
- Consumes: `CifDocument`, canonical core models, `parse_xyz_operation`.
- Produces: `map_cif_structures(document: CifDocument) -> tuple[tuple[Crystal, ...], tuple[Diagnostic, ...]]` and element normalization independent of external packages.

- [ ] **Step 1: Write failing basic mapper test**

```python
from cristma.io.cif.mapper import map_cif_structures
from cristma.io.cif.parser import parse_cif


MINIMAL = """data_si
_cell_length_a 5.43
_cell_length_b 5.43
_cell_length_c 5.43
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt 'P -1'
loop_
_space_group_symop_operation_xyz
'x,y,z'
'-x,-y,-z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Si1 Si 0 0 0 1
"""


def test_maps_asymmetric_site_and_exact_reported_symmetry():
    document = parse_cif(MINIMAL).document
    structures, diagnostics = map_cif_structures(document)
    assert not [d for d in diagnostics if d.severity.value == "error"]
    crystal = structures[0]
    assert crystal.cell.volume == 5.43 ** 3
    assert crystal.sites[0].components[0].element == "Si"
    assert crystal.space_group.provenance == "reported"
    assert len(crystal.space_group.operations) == 2
```

- [ ] **Step 2: Verify mapper test fails**

Run: `pytest -q tests/io/cif/test_mapper_basic.py`

Expected: FAIL because mapper and element modules are absent.

- [ ] **Step 3: Implement tag aliases and internal element normalization**

Define ordered aliases for modern and legacy cell, space-group, symmetry-operation, atom-site, formula, oxidation, and publication tags. Include a complete immutable set of IUPAC element symbols; normalize capitalization but retain reported type symbol in site metadata. Do not infer two-letter elements from arbitrary lowercase label prefixes without checking the symbol set.

- [ ] **Step 4: Implement block-to-crystal mapping**

Require six valid cell values and complete fractional coordinates. Map reported IT number, Hermann-Mauguin symbol, Hall symbol, setting, and origin fields without conflating them. Parse reported symmetry operations exactly and derive expanded-site provenance from those operations. If operations are absent, use identity only and emit `cif.map.symmetry_operations_missing` with provenance `identity_fallback`; do not synthesize a named group yet. Map one atom-loop row to one independent site with stable ID `<block>:<label>:<row-index>`. Preserve formula, mineral/common/systematic names, authors, journal, year, volume, pages, and DOI in typed metadata fields plus the original CIF document.

- [ ] **Step 5: Add explicit failure tests**

```python
def test_missing_cell_does_not_fabricate_structure():
    document = parse_cif("data_a\n_atom_site_label Si1\n").document
    structures, diagnostics = map_cif_structures(document)
    assert structures == ()
    assert {d.code for d in diagnostics} >= {"cif.map.cell_missing"}


def test_incomplete_coordinate_row_is_error():
    source = MINIMAL.replace("Si1 Si 0 0 0 1", "Si1 Si 0 ? 0 1")
    structures, diagnostics = map_cif_structures(parse_cif(source).document)
    assert structures == ()
    assert "cif.map.coordinate_missing" in {d.code for d in diagnostics}
```

- [ ] **Step 6: Run basic mapper tests**

Run: `pytest -q tests/io/cif/test_mapper_basic.py`

Expected: all tests PASS.

- [ ] **Step 7: Commit basic CIF mapping**

```bash
git add src/cristma/chemistry src/cristma/io/cif/names.py src/cristma/io/cif/mapper.py tests/io/cif/test_mapper_basic.py
git commit -m "feat: map cif cells symmetry and sites"
```

### Task 8: Occupancy, disorder, Wyckoff, and displacement mapping

**Files:**
- Modify: `src/cristma/io/cif/mapper.py`
- Modify: `src/cristma/core/structure.py`
- Test: `tests/io/cif/test_mapper_advanced.py`
- Create: `tests/io/cif/conftest.py`
- Create: `tests/fixtures/cif/mixed_disorder.cif`
- Create: `tests/fixtures/cif/anisotropic.cif`

**Interfaces:**
- Consumes: basic CIF mapper and canonical site model.
- Produces: safe mixed-position grouping, reported Wyckoff/multiplicity, isotropic U/B, anisotropic U tensors, and disorder provenance.

Create a local fixture loader in `tests/io/cif/conftest.py` that parses a named
fixture, maps its structures, combines parser and mapper diagnostics, and
returns `ReadResult(document, structures, diagnostics)`:

```python
@pytest.fixture
def read_fixture():
    def read(name: str) -> ReadResult:
        parsed = parse_cif((FIXTURE_DIR / name).read_text(encoding="utf-8"))
        structures, mapped = map_cif_structures(parsed.document)
        return ReadResult(parsed.document, structures, parsed.diagnostics + mapped)
    return read
```

- [ ] **Step 1: Add a mixed/disorder fixture and failing test**

Create a CIF fixture with two rows at the same coordinates, distinct elements,
occupancies `0.6` and `0.4`, identical disorder assembly/group, and reported
Wyckoff label. Then assert:

```python
def test_complementary_disorder_rows_form_one_mixed_site(read_fixture):
    result = read_fixture("mixed_disorder.cif")
    site = result.structures[0].sites[0]
    assert [(c.element, c.occupancy.value) for c in site.components] == [
        ("La", 0.6), ("Zr", 0.4)
    ]
    assert site.disorder_assembly == "A"
    assert site.disorder_group == "1"
    assert site.wyckoff == "4a"
```

- [ ] **Step 2: Run the mixed-position test and verify failure**

Run: `pytest -q tests/io/cif/test_mapper_advanced.py::test_complementary_disorder_rows_form_one_mixed_site`

Expected: FAIL because rows remain separate.

- [ ] **Step 3: Implement conservative mixed-position grouping**

Group coincident rows only when disorder assembly/group and label identity are compatible, elements differ, every occupancy is explicitly partial, and the sum is at most `1.0 + 1e-6`. Emit `cif.map.coincident_sites_unmerged` when coincident rows are intentionally left separate. Never merge solely by rounded coordinates.

- [ ] **Step 4: Add anisotropic displacement fixture and failing test**

```python
def test_maps_anisotropic_u_tensor_by_atom_label(read_fixture):
    site = read_fixture("anisotropic.cif").structures[0].sites[0]
    assert site.displacement.kind == "U_aniso"
    assert site.displacement.tensor[0][1].value == 0.0012
    assert site.displacement.tensor[1][0].value == 0.0012


def test_maps_reported_oxidation_and_checks_orbit_multiplicity(read_fixture):
    result = read_fixture("anisotropic.cif")
    site = result.structures[0].sites[0]
    assert site.components[0].oxidation_state.value == 4
    assert "cif.map.multiplicity_mismatch" not in {d.code for d in result.diagnostics}
```

- [ ] **Step 5: Implement U/B and anisotropic-loop mapping**

Join `_atom_site_aniso_*` rows by label, preserve uncertainties, construct a symmetric tensor, and reject conflicting `U_12/U_21` values. Convert Biso to Uiso only as a derived property using `U = B/(8*pi**2)` while retaining the reported representation. Validate positive semidefiniteness with NumPy eigenvalues and emit `cif.map.adp_not_positive_semidefinite` without discarding the reported tensor. Map reported oxidation values to site components. Expand each independent site with the exact reported operations and emit `cif.map.multiplicity_mismatch` when a reported multiplicity disagrees with the deduplicated orbit; retain both reported and calculated values.

- [ ] **Step 6: Run advanced mapper tests**

Run: `pytest -q tests/io/cif/test_mapper_advanced.py`

Expected: all tests PASS.

- [ ] **Step 7: Commit advanced CIF semantics**

```bash
git add src/cristma/core/structure.py src/cristma/io/cif/mapper.py tests/io/cif/test_mapper_advanced.py tests/fixtures/cif
git commit -m "feat: map cif disorder and displacement data"
```

### Task 9: Preserve-mode and canonical CIF writing

**Files:**
- Modify: `src/cristma/io/cif/document.py`
- Create: `src/cristma/io/cif/writer.py`
- Test: `tests/io/cif/test_writer.py`

**Interfaces:**
- Consumes: source spans, `CifDocument`, and `Crystal`.
- Produces: `SourceEdit`, `replace_scalar(document, block, tag, raw_value)`, `write_cif_document(document, mode="preserve")`, and `write_crystal_cif(crystal, block_name=None)`.

- [ ] **Step 1: Write failing preserve-mode tests**

```python
from cristma.io.cif.document import replace_scalar
from cristma.io.cif.parser import parse_cif
from cristma.io.cif.writer import write_cif_document


def test_unchanged_document_round_trips_byte_for_byte():
    source = "data_a\r\n# keep me\r\n_local_unknown 'A B'\r\n_cell_length_a 5.0\r\n"
    document = parse_cif(source).document
    assert write_cif_document(document, mode="preserve") == source


def test_scalar_edit_preserves_unknown_content_and_newlines():
    source = "data_a\r\n# keep me\r\n_local_unknown 'A B'\r\n_cell_length_a 5.0\r\n"
    document = replace_scalar(parse_cif(source).document, "a", "_cell_length_a", "5.1(2)")
    rendered = write_cif_document(document, mode="preserve")
    assert "# keep me\r\n_local_unknown 'A B'" in rendered
    assert "_cell_length_a 5.1(2)" in rendered
```

- [ ] **Step 2: Verify preserve tests fail**

Run: `pytest -q tests/io/cif/test_writer.py -k preserve`

Expected: FAIL because writer/edit APIs are absent.

- [ ] **Step 3: Implement non-overlapping source edits**

`SourceEdit` stores start offset, end offset, and replacement. `replace_scalar` targets only the value token span. Rendering sorts edits in descending offset order and rejects overlapping edits with `ValueError`. With no edits, preserve mode returns the original source exactly.

- [ ] **Step 4: Write failing canonical-writer test**

```python
def test_canonical_writer_emits_parseable_asymmetric_structure():
    source = """data_si
_cell_length_a 5.43
_cell_length_b 5.43
_cell_length_c 5.43
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_space_group_symop_operation_xyz
'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Si1 Si 0 0 0 1
"""
    parsed_source = parse_cif(source)
    original, original_diagnostics = map_cif_structures(parsed_source.document)
    assert not [d for d in original_diagnostics if d.severity.value == "error"]
    text = write_crystal_cif(original[0], block_name="silicon")
    reparsed = parse_cif(text)
    structures, diagnostics = map_cif_structures(reparsed.document)
    assert reparsed.ok
    assert not [d for d in diagnostics if d.severity.value == "error"]
    assert structures[0].sites[0].id.endswith("Si1:0")
    assert "_space_group_symop_operation_xyz" in text
```

- [ ] **Step 5: Implement deterministic canonical writing**

Emit UTF-8 CIF 1.1 with cell values, reported space-group metadata, explicit symmetry-operation loop, atom-site loop, occupancy, Wyckoff/multiplicity, isotropic or anisotropic displacement loops, disorder fields, and formula when available. Quote values according to CIF 1.1 lexical rules and use semicolon text for multiline values. Order sites by asymmetric-unit order, never by expanded coordinates.

- [ ] **Step 6: Run writer tests**

Run: `pytest -q tests/io/cif/test_writer.py`

Expected: all tests PASS.

- [ ] **Step 7: Commit CIF writers**

```bash
git add src/cristma/io/cif/document.py src/cristma/io/cif/writer.py tests/io/cif/test_writer.py
git commit -m "feat: write preserved and canonical cif"
```

### Task 10: Public CIF handler, real fixtures, and slice verification

**Files:**
- Create: `src/cristma/io/cif/handler.py`
- Modify: `src/cristma/io/cif/__init__.py`
- Modify: `src/cristma/io/registry.py`
- Modify: `src/cristma/__init__.py`
- Modify: `README.md`
- Create: `tests/io/cif/test_end_to_end.py`
- Create: `tests/fixtures/cif/fixture_provenance.md`
- Create: `tests/fixtures/cif/lithium_triborate.cif`
- Create: `tests/fixtures/cif/mixed_occupancy_positions.cif`

**Interfaces:**
- Consumes: parser, mapper, writers, registry.
- Produces: public `cristma.read(...)`, `cristma.read_text(...)`, `cristma.write(...)`, and `CifFormatHandler`.

- [ ] **Step 1: Add provenance-recorded real fixtures**

Add repository fixtures derived from these existing application fixtures:

```text
/Users/artem/Yandex.Disk.localized/Python/XRD/XRD_Analysis_Toolkit/XRD_Craft/tests/data/structures/lithium_triborate.cif
/Users/artem/Yandex.Disk.localized/Python/XRD/XRD_Analysis_Toolkit/XRD_Craft/tests/data/mixed_occupancy_positions.cif
```

Record source path, SHA-256, date copied, and purpose in
`tests/fixtures/cif/fixture_provenance.md`. Do not import Craft code at runtime.

- [ ] **Step 2: Write failing end-to-end tests**

```python
from pathlib import Path
import cristma


def test_public_read_maps_real_inorganic_cif():
    path = Path("tests/fixtures/cif/lithium_triborate.cif")
    result = cristma.read(path)
    assert result.ok
    assert result.structures
    assert {component.element for site in result.structures[0].sites for component in site.components} >= {"Li", "B", "O"}


def test_public_preserve_round_trip_keeps_mixed_occupancy_fixture(tmp_path):
    source = Path("tests/fixtures/cif/mixed_occupancy_positions.cif")
    result = cristma.read(source)
    target = tmp_path / "copy.cif"
    cristma.write(result.document, target, mode="preserve")
    assert target.read_bytes() == source.read_bytes()
```

- [ ] **Step 3: Verify end-to-end tests fail**

Run: `pytest -q tests/io/cif/test_end_to_end.py`

Expected: FAIL because the handler is not registered.

- [ ] **Step 4: Implement and register the CIF handler**

Probe confidence is `1.0` for text whose first non-comment reserved token is `data_`, `0.8` for `.cif` suffix with recognizable tags, and `0.0` otherwise. `read_text` chains lexer/parser and mapper diagnostics into one `ReadResult`. Public `write` accepts `CifDocument` in preserve mode and `Crystal` in canonical mode; reject incompatible combinations explicitly.

- [ ] **Step 5: Document the public API**

Add concise examples showing `read`, diagnostic inspection, asymmetric sites, preserve writing, and canonical writing. State explicitly that Gemmi/pymatgen/PyXtal are not required.

- [ ] **Step 6: Run the complete CIF slice, not unrelated application tests**

Run: `pytest -q tests/core tests/symmetry tests/io`

Expected: all CrIStMa slice tests PASS.

- [ ] **Step 7: Verify package installation and imports in a clean environment**

Run:

```bash
python -m venv .venv-check
.venv-check/bin/pip install -e '.[test]'
.venv-check/bin/python -c "import cristma; print(cristma.__version__)"
```

Expected output ends with `0.1.0` and installs no specialized crystallographic package.

- [ ] **Step 8: Commit the completed CIF vertical slice**

```bash
git add src tests README.md
git commit -m "feat: expose native cif structure io"
```

## Completion gate

Before declaring this plan complete:

- run `pytest -q tests/core tests/symmetry tests/io`;
- run `python -m build` and inspect wheel contents;
- confirm `python -c "import cristma"` works in the clean check environment;
- confirm a real CIF preserve round-trip is byte-identical;
- confirm a canonical write can be read back into an equivalent asymmetric unit;
- confirm dependency metadata contains no Gemmi, pymatgen, PyXtal, CrysPy, GSAS-II, or Qt package;
- run `git status --short` and verify only intentional files remain.
