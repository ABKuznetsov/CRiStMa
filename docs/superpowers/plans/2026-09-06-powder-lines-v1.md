# Powder Diffraction Lines v1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build deterministic intrinsic powder lines from a `StructureFactorSet` and an explicit multi-component X-ray spectrum, including a reproducible Cu K-alpha1/K-alpha2 preset.

**Architecture:** Immutable radiation and powder result models live in `powder_models.py`; `powder.py` performs only Friedel grouping, multiplicity weighting, and Bragg-angle calculation. The Cu K-alpha preset is a packaged, checksummed build-time extraction from pinned xraylib data with no runtime xraylib dependency.

**Tech Stack:** Python 3.11+, standard library, NumPy only through existing diffraction dependencies, pytest for temporary development tests, setuptools build.

**Spec:** `docs/superpowers/specs/2026-09-06-powder-lines-v1-design.md`

## Global Constraints

- Public flow is exactly `StructureFactorSet + RadiationSpectrum -> PowderLineCalculator -> PowderLineSet`.
- The calculator does not accept a structure, space-group setting, experiment, profile, or correction model.
- The public field is named `intrinsic_line_intensity`, never `intensity`.
- LP, preferred orientation, absorption, profile broadening, background, matching, and refinement are out of scope.
- A multi-component spectrum never uses an effective or averaged wavelength.
- xraylib is build-time reference provenance only and must not become a project dependency.
- Development tests remain under `/tmp/cristma-powder-line-tests` and are neither committed nor packaged.
- Repository `docs` and `tests` remain excluded from wheel and sdist artifacts.
- Do not publish, tag, push, or change the local beta version during this plan.

---

### Task 1: Immutable radiation models

**Files:**
- Create: `src/cristma/diffraction/powder_models.py`
- Test: `/tmp/cristma-powder-line-tests/test_radiation_models.py`

**Interfaces:**
- Consumes: no new project interfaces.
- Produces: `RadiationComponent`, `RadiationSpectrumProvenance`, and `RadiationSpectrum`.

- [ ] **Step 1: Write failing validation and normalization tests**

Create tests that require the following interface and behavior:

```python
component_1 = RadiationComponent("cu-ka1", "Cu Kalpha1", 1.5406, 2.0)
component_2 = RadiationComponent("cu-ka2", "Cu Kalpha2", 1.5444, 1.0)
provenance = RadiationSpectrumProvenance.user_supplied("test-spectrum")
spectrum = RadiationSpectrum((component_2, component_1), provenance)

assert spectrum.normalized_weights == pytest.approx((1 / 3, 2 / 3))
assert tuple(item.component_id for item in spectrum.components) == (
    "cu-ka2",
    "cu-ka1",
)
```

Also reject bool/non-real/non-finite/non-positive wavelength or weight, empty
IDs and labels, empty spectra, duplicate component IDs, and invalid SHA-256
values when preset provenance supplies a checksum.

- [ ] **Step 2: Run the tests and verify RED**

Run:

```bash
PYTHONPATH=src pytest -q /tmp/cristma-powder-line-tests/test_radiation_models.py
```

Expected: collection fails because `cristma.diffraction.powder_models` does not
exist.

- [ ] **Step 3: Implement the minimal immutable models**

Use frozen slotted dataclasses with these signatures:

```python
@dataclass(frozen=True, slots=True)
class RadiationComponent:
    component_id: str
    label: str
    wavelength_angstrom: float
    relative_weight: float

@dataclass(frozen=True, slots=True)
class RadiationSpectrumProvenance:
    dataset_id: str
    dataset_version: str
    source: str
    energy_source: str | None = None
    radiative_rate_source: str | None = None
    energy_to_wavelength_formula: str | None = None
    hc_value: float | None = None
    hc_units: str | None = None
    xraylib_version: str | None = None
    xraylib_commit: str | None = None
    resource_checksum: str | None = None
    generator: str | None = None

    @classmethod
    def user_supplied(cls, dataset_id: str) -> "RadiationSpectrumProvenance": ...

@dataclass(frozen=True, slots=True)
class RadiationSpectrum:
    components: tuple[RadiationComponent, ...]
    provenance: RadiationSpectrumProvenance

    @property
    def normalized_weights(self) -> tuple[float, ...]: ...
```

Normalize with `math.fsum`; preserve component order and original weights.
Require either all Cu-preset provenance fields or none of them, preventing
partially reproducible provenance.

- [ ] **Step 4: Run the model tests and verify GREEN**

Run the Step 2 command. Expected: all tests pass.

- [ ] **Step 5: Commit the radiation models**

```bash
git add src/cristma/diffraction/powder_models.py
git commit -m "Add powder radiation models"
```

---

### Task 2: Reproducible Cu K-alpha preset

**Files:**
- Create: `tools/build_xray_radiation.py`
- Create: `src/cristma/reference_data/resources/xray/xray_radiation.json`
- Modify: `src/cristma/reference_data/resources/xray/SOURCE.md`
- Modify: `src/cristma/reference_data/resources/xray/XRAYLIB_LICENSE.txt` only if the pinned source license text differs.
- Modify: `src/cristma/diffraction/powder_models.py`
- Test: `/tmp/cristma-powder-line-tests/test_radiation_preset.py`

**Interfaces:**
- Consumes: Task 1 models.
- Produces: `RadiationSpectrum.copper_k_alpha() -> RadiationSpectrum`.

- [ ] **Step 1: Write failing preset and provenance tests**

Require two ordered components named `cu-ka1` and `cu-ka2`, distinct positive
wavelengths with K-alpha1 shorter than K-alpha2, normalized weights summing to
one, and complete non-null preset provenance:

```python
spectrum = RadiationSpectrum.copper_k_alpha()
assert tuple(item.component_id for item in spectrum.components) == (
    "cu-ka1",
    "cu-ka2",
)
assert spectrum.components[0].wavelength_angstrom < spectrum.components[1].wavelength_angstrom
assert math.fsum(spectrum.normalized_weights) == pytest.approx(1.0)
assert spectrum.provenance.energy_to_wavelength_formula == "lambda_angstrom = hc_keV_angstrom / energy_keV"
assert spectrum.provenance.hc_units == "keV angstrom"
assert len(spectrum.provenance.resource_checksum or "") == 64
```

Test that changing a resource component makes loading fail with
`ValueError("X-ray radiation resource checksum mismatch")`.

- [ ] **Step 2: Run the preset tests and verify RED**

Run:

```bash
PYTHONPATH=src pytest -q /tmp/cristma-powder-line-tests/test_radiation_preset.py
```

Expected: `RadiationSpectrum.copper_k_alpha` is missing.

- [ ] **Step 3: Implement the pinned build tool and resource**

Pin the same xraylib release and commit already used by `build_xray_f0.py`:

```python
XRAYLIB_VERSION = "4.3.0"
XRAYLIB_COMMIT = "f94a3f5008dfd1c882b88ff26cd5052559423c83"
HC_KEV_ANGSTROM = 12.398419843320026
```

At build time call `LineEnergy(29, KA1_LINE)`, `LineEnergy(29, KA2_LINE)`,
`RadRate(29, KA1_LINE)`, and `RadRate(29, KA2_LINE)`. Convert each energy with
`wavelength = HC_KEV_ANGSTROM / energy_keV`. Write canonical sorted-key JSON
containing components plus every provenance field required by the spec. Store
the SHA-256 of canonical component data as `resource_checksum`.

- [ ] **Step 4: Load and validate the packaged preset**

Implement a cached `RadiationSpectrum.copper_k_alpha()` that reads the resource
with `importlib.resources`, recomputes the component checksum before creating
models, and imports no xraylib module. Add a subprocess test that blocks
`xraylib` import and still loads the preset.

- [ ] **Step 5: Run the preset and model tests**

```bash
PYTHONPATH=src pytest -q \
  /tmp/cristma-powder-line-tests/test_radiation_models.py \
  /tmp/cristma-powder-line-tests/test_radiation_preset.py
```

Expected: all tests pass.

- [ ] **Step 6: Commit the preset**

```bash
git add tools/build_xray_radiation.py \
  src/cristma/reference_data/resources/xray/xray_radiation.json \
  src/cristma/reference_data/resources/xray/SOURCE.md \
  src/cristma/diffraction/powder_models.py
git commit -m "Bundle Cu K-alpha radiation preset"
```

---

### Task 3: Powder result models and diagnostic codes

**Files:**
- Modify: `src/cristma/diffraction/powder_models.py`
- Modify: `src/cristma/diffraction/diagnostics.py`
- Test: `/tmp/cristma-powder-line-tests/test_powder_models.py`

**Interfaces:**
- Consumes: Task 1 radiation models and existing `MillerIndex`, `ReflectionSetStatus`, `StructureFactorSet`.
- Produces: `PowderLine`, `PowderReflectionFamily`, `PowderLineProvenance`, and `PowderLineSet`.

- [ ] **Step 1: Write failing result-invariant tests**

Require the exact public fields below. Test non-empty IDs, positive finite
wavelength and d-spacing, `0 <= two_theta_deg <= 180`, non-negative finite
strengths, unique sorted member reflection IDs, line-to-family identity,
declared component order inside a family, family order by
`(family_sort_angle, family_id)`, and source status inheritance.

- [ ] **Step 2: Run the result tests and verify RED**

```bash
PYTHONPATH=src pytest -q /tmp/cristma-powder-line-tests/test_powder_models.py
```

Expected: result classes are missing.

- [ ] **Step 3: Add stable powder diagnostic constants**

```python
POWDER_MISSING_FRIEDEL_MATE = "diffraction.powder.missing_friedel_mate"
POWDER_NONRECIPROCAL_FRIEDEL_LINK = "diffraction.powder.nonreciprocal_friedel_link"
POWDER_FRIEDEL_D_SPACING_MISMATCH = "diffraction.powder.friedel_d_spacing_mismatch"
```

Export them from `diagnostics.py`.

- [ ] **Step 4: Implement result models**

Use frozen slotted dataclasses with these fields:

```python
@dataclass(frozen=True, slots=True)
class PowderLine:
    line_id: str
    family_id: str
    radiation_component_id: str
    wavelength_angstrom: float
    normalized_radiation_weight: float
    two_theta_deg: float
    intrinsic_line_intensity: float

@dataclass(frozen=True, slots=True)
class PowderReflectionFamily:
    family_id: str
    reflection_ids: tuple[str, ...]
    representative_hkls: tuple[MillerIndex, ...]
    d_spacing: float
    multiplicity_crystallographic: int
    family_strength: float
    lines: tuple[PowderLine, ...]

    @property
    def family_sort_angle(self) -> float: ...

@dataclass(frozen=True, slots=True)
class PowderLineProvenance:
    method: str
    version: str
    reflections_considered: int
    families_emitted: int
    radiation_components_skipped: int

@dataclass(frozen=True, slots=True)
class PowderLineSet:
    families: tuple[PowderReflectionFamily, ...]
    structure_factors: StructureFactorSet
    spectrum: RadiationSpectrum
    provenance: PowderLineProvenance
    diagnostics: tuple[Diagnostic, ...] = ()

    @property
    def status(self) -> ReflectionSetStatus: ...

    @property
    def lines(self) -> tuple[PowderLine, ...]: ...

    @property
    def lines_by_angle(self) -> tuple[PowderLine, ...]: ...
```

The stored flat `lines` property follows family order and spectrum order;
`lines_by_angle` returns a derived `(two_theta_deg, line_id)` ordering.

- [ ] **Step 5: Run model tests and verify GREEN**

Run the Step 2 command. Expected: all tests pass.

- [ ] **Step 6: Commit result models and diagnostics**

```bash
git add src/cristma/diffraction/powder_models.py src/cristma/diffraction/diagnostics.py
git commit -m "Add powder line result models"
```

---

### Task 4: Powder line calculator

**Files:**
- Create: `src/cristma/diffraction/powder.py`
- Test: `/tmp/cristma-powder-line-tests/test_powder_calculator.py`
- Test: `/tmp/cristma-powder-line-tests/test_powder_invariants.py`

**Interfaces:**
- Consumes: `StructureFactorSet`, `RadiationSpectrum`, Task 3 result models and diagnostic codes.
- Produces: `PowderLineCalculator.calculate(structure_factors, spectrum) -> PowderLineSet`.

- [ ] **Step 1: Write failing analytical calculation tests**

Build small in-memory reflection and structure-factor fixtures. Require:

```python
result = PowderLineCalculator().calculate(structure_factors, spectrum)
line = result.families[0].lines[0]
assert line.two_theta_deg == pytest.approx(
    math.degrees(2 * math.asin(wavelength / (2 * d_spacing)))
)
assert line.intrinsic_line_intensity == pytest.approx(
    normalized_weight * multiplicity * f_squared
)
```

Cover one singleton family, one reciprocal Friedel pair, an extinct reflection
omitted entirely, an allowed `F_squared == 0` family retained, one unreachable
component skipped, and a fully unreachable family omitted.

- [ ] **Step 2: Run calculation tests and verify RED**

```bash
PYTHONPATH=src pytest -q /tmp/cristma-powder-line-tests/test_powder_calculator.py
```

Expected: `PowderLineCalculator` is missing.

- [ ] **Step 3: Implement deterministic family grouping**

Index reflections and structure factors by reflection ID. Iterate the source
reflection order, skipping already consumed IDs and systematic extinctions.
For a non-null mate link, require existence and reciprocity before creating a
sorted two-member tuple. Define `family_id` from the sorted member IDs with a
length-prefixed SHA-256 serialization so separators inside IDs cannot collide.

Compute:

```python
family_strength = math.fsum(
    reflection.multiplicity_crystallographic * factor.f_squared
    for reflection, factor in members
)
family_multiplicity = sum(
    reflection.multiplicity_crystallographic for reflection, _ in members
)
```

- [ ] **Step 4: Implement component lines and physical exclusions**

For each spectrum component in declared order, use its normalized weight and
calculate the first-order Bragg angle. If `wavelength / (2*d) > 1`, increment
`radiation_components_skipped` and emit no line. Use deterministic line IDs
derived from `(family_id, component_id)`. Sort completed families by
`(min(line.two_theta_deg), family_id)`.

- [ ] **Step 5: Write and run strict invariant tests**

Use deliberately corrupted frozen objects via `object.__setattr__` only inside
tests to prove the defensive calculator checks. Require exact error codes and
evidence for missing mate, nonreciprocal mate, and Friedel d-spacing mismatch.
Use the existing reciprocal boundary tolerance `rel_tol=1e-10`, matching
`Reflection`'s d-spacing/reciprocal-norm validation scale.

```bash
PYTHONPATH=src pytest -q \
  /tmp/cristma-powder-line-tests/test_powder_calculator.py \
  /tmp/cristma-powder-line-tests/test_powder_invariants.py
```

Expected: all tests pass.

- [ ] **Step 6: Add the high-angle Cu doublet regression**

For one family with a reachable high-angle d-spacing, assert K-alpha2 has a
larger `two_theta_deg` than K-alpha1 and that their separation is larger than
for a second low-angle family. Reorder a custom spectrum as K-alpha2,
K-alpha1 and assert family order is unchanged while line order follows the
custom declaration.

- [ ] **Step 7: Run all powder tests**

```bash
PYTHONPATH=src pytest -q /tmp/cristma-powder-line-tests
```

Expected: all tests pass.

- [ ] **Step 8: Commit the calculator**

```bash
git add src/cristma/diffraction/powder.py
git commit -m "Calculate intrinsic powder diffraction lines"
```

---

### Task 5: Public API, regression, and package verification

**Files:**
- Modify: `src/cristma/diffraction/__init__.py`
- Modify: `README.md`
- Test: `/tmp/cristma-powder-line-tests/test_powder_public_api.py`
- Test: `/tmp/cristma-powder-line-tests/test_powder_real_structure.py`

**Interfaces:**
- Consumes: Tasks 1-4 public types.
- Produces: the intentional `cristma.diffraction` powder-lines v1 API.

- [ ] **Step 1: Write failing public API tests**

Require direct imports from `cristma.diffraction` for:

```python
PowderLineCalculator
PowderLine
PowderLineProvenance
PowderLineSet
PowderReflectionFamily
RadiationComponent
RadiationSpectrum
RadiationSpectrumProvenance
```

Verify every public name appears in `cristma.diffraction.__all__` and no build
helper becomes public.

- [ ] **Step 2: Export the intentional API and document the boundary**

Add only the listed names to `src/cristma/diffraction/__init__.py`. Add a short
README example that calls `PowderLineCalculator` and labels the result
"intrinsic line intensity" while stating that profile and experimental
comparison belong to later/external layers.

- [ ] **Step 3: Add a real-structure regression**

Read an existing CIF fixture, resolve its catalog `SpaceGroupSetting`, generate
a bounded `ReflectionSet`, calculate neutral-X-ray structure factors, and then
calculate Cu K-alpha powder lines. Assert:

```python
assert result.status is structure_factors.status
assert result.families
assert all(family.lines for family in result.families)
assert any(len(family.lines) == 2 for family in result.families)
assert all(line.intrinsic_line_intensity >= 0 for line in result.lines)
```

- [ ] **Step 4: Run focused and baseline verification**

```bash
PYTHONPATH=src pytest -q /tmp/cristma-powder-line-tests
PYTHONPATH=src pytest -q --import-mode=importlib \
  /tmp/cristma-beta-baseline.UWss30/tests \
  --ignore=/tmp/cristma-beta-baseline.UWss30/tests/tools
python -m compileall -q src
```

Expected: all commands exit zero.

- [ ] **Step 5: Build and inspect local artifacts**

Build wheel and sdist into a new `/tmp/cristma-powder-dist.*` directory. Run
`twine check`, inspect both archives, and verify:

- `xray_radiation.json`, `SOURCE.md`, and `XRAYLIB_LICENSE.txt` are present;
- no repository `tests` or `docs` directory is present;
- metadata has no xraylib runtime dependency;
- version remains `0.1.0b4`.

- [ ] **Step 6: Install the wheel and smoke-test without repository imports**

Create a new `/tmp` virtual environment, install the wheel with `--no-deps`,
run from `/tmp`, and assert the imported module path is inside that environment.
Load `RadiationSpectrum.copper_k_alpha()` and calculate at least one two-line
family from the installed package.

- [ ] **Step 7: Commit the API and documentation**

```bash
git add src/cristma/diffraction/__init__.py README.md
git commit -m "Expose powder diffraction lines v1"
```

- [ ] **Step 8: Final local-only audit**

Run `git status --short --branch` and `git log -5 --oneline`. Report test
counts, artifact paths, checksums, and any diagnostics. Do not push, tag, or
upload artifacts.
