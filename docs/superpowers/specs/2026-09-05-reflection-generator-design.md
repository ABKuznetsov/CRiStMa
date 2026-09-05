# CrIStMa Diffraction Reflection Generator Design

**Date:** 2026-09-05  
**Status:** approved in discussion; written specification awaiting final review  
**Implementation scope:** reciprocal metric, bounded reflection generation,
reciprocal symmetry orbits, Friedel relations, and systematic absences

## 1. Purpose

This milestone establishes `cristma.diffraction` as an independent scientific
layer. It generates crystallographic reflections from a unit cell and one
unambiguously identified catalog space-group setting. It does not require an
atomic structure.

The milestone implements this complete calculation:

```text
UnitCell
    +
SpaceGroupSetting
    +
d_min
    |
    v
ReciprocalMetric
    |
    v
bounded integer ellipsoid enumeration
    |
    v
reciprocal point-group orbits
    |
    v
systematic-extinction analysis
    |
    v
Friedel linking
    |
    v
ReflectionSet
```

The package is stateless and consumer-neutral. Its public types, diagnostics,
provenance, and implementation contain no application-specific concepts.

## 2. Scope

The first milestone includes:

- the crystallographic reciprocal basis without a `2*pi` factor;
- reciprocal metric and `d_hkl` calculations;
- complete `hkl` generation down to an explicit `d_min`;
- a bounded integer-lattice enumeration over a positive-definite quadratic
  form;
- reciprocal point-group orbits and deterministic representatives;
- crystallographic multiplicity;
- explicit Friedel relations without unconditional inversion merging;
- exact symmetry-only systematic-extinction decisions;
- structured extinction evidence and secondary cause classification;
- immutable results, diagnostics, and reproducibility provenance.

The milestone excludes:

- atomic coordinates, species, occupancy, and displacement parameters;
- scattering factors, structure factors, `F`, and `F^2`;
- radiation, wavelength, `2theta`, and diffraction geometry;
- powder-line merging and multiplicity corrections;
- profile functions, broadening, backgrounds, grids, and experimental data;
- indexing, phase identification, optimization, and refinement;
- resolution of an arbitrary `SpaceGroupDefinition` to a catalog setting.

No empty modules are created for later milestones.

## 3. Package boundary

The first package slice is:

```text
cristma/diffraction/
|- __init__.py
|- reciprocal.py
|- reflections.py
|- extinction.py
|- models.py
`- diagnostics.py
```

Responsibilities are intentionally narrow:

- `models.py` owns immutable public value and result types;
- `reciprocal.py` owns reciprocal metric calculations and integer ellipsoid
  enumeration;
- `extinction.py` owns exact phase-bucket analysis and explanatory cause
  classification;
- `reflections.py` owns `ReflectionGenerator`, orbit construction, stable IDs,
  and Friedel linking;
- `diagnostics.py` owns diffraction diagnostic codes and invariant errors;
- `__init__.py` exports only the intentional public API.

The layer depends only on existing CrIStMa types, NumPy, and the Python standard
library. It introduces no new runtime dependency.

## 4. Public API and accepted symmetry input

The primary API is:

```python
from cristma.diffraction import ReflectionGenerator

result = ReflectionGenerator().generate(
    cell=cell,
    space_group=setting,
    d_min=0.8,
)
```

Its contract is:

```text
ReflectionGenerator.generate(
    cell: UnitCell,
    space_group: SpaceGroupSetting,
    d_min: float,
) -> ReflectionSet
```

Only `SpaceGroupSetting` is accepted. It fixes the Hall setting, `setting_id`,
direct and reciprocal index basis, origin/choice, and exact symmetry-operation
set. `ReflectionGenerator` never guesses or resolves a setting.

Every input operation must have a non-empty ID unique within the setting. The
packaged catalog already satisfies this requirement; rejecting missing or
duplicate IDs keeps exact evidence stable for manually constructed invalid
inputs.

`SpaceGroupDefinition` is not accepted directly. A future independent adapter
may resolve a definition through a catalog and return an explicit resolved,
ambiguous, or unresolved result. Such resolution is outside this milestone.

`ReflectionGenerator` contains only immutable, inspectable configuration:

```text
max_candidates
boundary_tolerance
metric_compatibility_tolerance
```

It stores no current cell, setting, reflection set, or hidden cache. It exposes
the same `get_config()` and `clone(...)` conventions as other configured CrIStMa
tools.

## 5. Public models

### 5.1 Miller index

```python
@dataclass(frozen=True, slots=True, order=True)
class MillerIndex:
    h: int
    k: int
    l: int
```

Each component must be an integer and must not be `bool`. `MillerIndex(0, 0,
0)` is a valid general value object, but the generator never emits it and a
`Reflection` rejects any orbit containing it. Dataclass ordering defines the
canonical lexicographic order used throughout the milestone.

### 5.2 Reflection

```text
Reflection
|- reflection_id: str
|- representative_hkl: MillerIndex
|- equivalent_hkls: tuple[MillerIndex, ...]
|- d_spacing: float
|- reciprocal_norm: float
|- multiplicity_crystallographic: int
|- friedel_mate_id: str | None
|- extinction: ExtinctionResult
`- provenance: ReflectionProvenance
```

One `Reflection` represents one reciprocal orbit under the actual point-group
rotations of the selected setting. It enforces:

```text
equivalent_hkls is non-empty, unique, and sorted
representative_hkl == equivalent_hkls[0]
multiplicity_crystallographic == len(equivalent_hkls)
(0, 0, 0) is absent from equivalent_hkls
reciprocal_norm > 0
d_spacing > 0
```

The stable ID is:

```text
reflection:{space_group_setting_id}:{h}:{k}:{l}
```

where `h`, `k`, and `l` belong to `representative_hkl`. Cell parameters,
`d_spacing`, enumeration order, and calculation limits do not enter the ID.

### 5.3 Reflection set

```text
ReflectionSet
|- reflections: tuple[Reflection, ...]
|- space_group_setting_id: int
|- d_min: float
|- status: ReflectionSetStatus
|- diagnostics: tuple[Diagnostic, ...]
`- provenance: ReflectionGenerationProvenance
```

`ReflectionSetStatus` has exactly two values:

```text
COMPLETE
INCOMPLETE
```

`COMPLETE` means the entire `d_min` region was covered. `INCOMPLETE` means a
correct calculation was deliberately truncated only because `max_candidates`
was reached. An incomplete result must contain a corresponding diagnostic.

`reflections` is unique and sorted by `representative_hkl`. Derived immutable
views expose allowed and systematically extinct reflections. Sorting by
`d_spacing`, reciprocal norm, or another presentation criterion is a derived
operation and never changes scientific identity.

### 5.4 Friedel contract

For a reflection represented by `h`:

```text
friedel_mate_id is None
iff -h belongs to the same crystallographic reciprocal orbit
```

If inversion is not in the point-group orbit, the two reciprocal orbits remain
separate and link to one another:

```text
A.friedel_mate_id == B.reflection_id
B.friedel_mate_id == A.reflection_id
```

Even an `INCOMPLETE` set contains no dangling or ambiguous Friedel link. Orbit
construction therefore closes every accepted seed atomically under both the
point group and its Friedel partner before materializing results. A completely
untested orbit may be absent from an incomplete set, but a present orbit is
internally complete and Friedel-consistent.

## 6. Reciprocal-space convention

`UnitCell.matrix` contains row-wise direct Cartesian basis vectors. The
reciprocal basis is defined without `2*pi`, and the reciprocal metric is the
inverse of the direct metric:

```text
G  = A A^T
G* = G^-1
```

For a column Miller index `h = (h, k, l)^T`:

```text
reciprocal_norm_squared = h^T G* h
reciprocal_norm         = sqrt(h^T G* h)
d_hkl                   = 1 / reciprocal_norm
```

The physical scattering-vector magnitude `Q = 2*pi/d` is a separate future
quantity and is not exposed as `reciprocal_norm`.

`ReciprocalMetric` is an immutable public calculation value containing a
read-only reciprocal basis and metric plus methods for norm and spacing. It
rejects the zero Miller index for norm/spacing calculations.

## 7. Cell and point-group compatibility

Before reflection generation, each unique exact point-group rotation is checked
against the numerical reciprocal metric. Under the convention used by this
milestone:

```text
h' = R^T h
```

and the rotation must preserve `h^T G* h` within one documented metric
compatibility tolerance. A clear mismatch raises:

```text
DiffractionInvariantError
code = diffraction.reflections.incompatible_cell_and_setting
```

Evidence contains the setting ID, six numerical cell values, operation ID or
rotation identity, and measured metric residual.

This is only a consistency guard. The diffraction layer does not infer a
crystal system, symmetrize or repair a cell, analyze pseudosymmetry, or interpret
experimental cell uncertainty.

## 8. Bounded integer ellipsoid enumeration

Reflection generation is governed by the physical inequality:

```text
h^T G* h <= 1 / d_min^2
```

There is no fixed `max_index`. Integer points are enumerated directly inside
the positive-definite quadratic form using a Cholesky- or LDL-based recursive
branch-and-bound algorithm. At each level, the remaining quadratic budget
produces the finite integer interval for the next index.

Every complete `(h, k, l)` reaching the leaf is checked again against the
original quadratic form. A single documented relative `boundary_tolerance` is
applied only to this floating-point ellipsoid boundary and to conservative
branch bounds. It never enters reciprocal symmetry, orbit identity, phase
relations, or extinction decisions.

The zero index is discarded. Accepted indices are sorted canonically before
orbit construction, so traversal details cannot alter complete scientific
results or reflection IDs.

`max_candidates` counts complete integer triples that reach the final original
quadratic-form check. It does not count recursive branches. When unvisited
branches remain after that limit, enumeration stops and the result becomes
`INCOMPLETE` with:

```text
diffraction.reflections.search_limit_reached
```

The partial seed set is still closed atomically under point-group symmetry and
Friedel relations before the result is returned.

## 9. Reciprocal symmetry orbits

The existing direct-space convention is:

```text
x' = R x + t
```

This milestone consistently uses the reciprocal action:

```text
h' = R^T h
```

Rotation matrices remain exact rational/integer values. Applying a valid
catalog rotation to a Miller index must produce another exact integer Miller
index; failure raises `DiffractionInvariantError`.

Translations do not determine reciprocal orbit membership. Duplicate
rotations are removed before orbit construction. For every accepted seed, the
generator calculates the complete unique point-group orbit, sorts it, and uses:

```text
representative_hkl = min(equivalent_hkls)
multiplicity_crystallographic = len(equivalent_hkls)
```

Global inversion is never inserted automatically. It contributes to the same
orbit only when it is generated by the setting's actual rotations.

The final quadratic form is evaluated for generated orbit members as a
defensive invariant. Metric compatibility requires every member to have the
same reciprocal norm and remain within the same `d_min` boundary.

## 10. Exact systematic extinctions

### 10.1 Scientific decision

Systematic absence is determined only by the complete exact symmetry-operation
set. For one representative Miller index `h`, every operation contributes:

```text
k   = R^T h
phi = h dot t mod 1
```

`k` is an exact Miller index and `phi` is a `Fraction` normalized into `[0, 1)`.
Operations are grouped into buckets by equal `k`.

Each bucket is a coset of the stabilizer of the reciprocal vector. Relative
phase factors define an exact character of that stabilizer:

```text
reference_phase = phase of a deterministic reference operation
relative_phase  = phase - reference_phase mod 1

all relative_phase == 0  -> trivial character; bucket survives
any relative_phase != 0  -> non-trivial character; bucket cancels exactly
```

This uses the group property of a complete space-group operation set and does
not require numerical complex sums or general cyclotomic arithmetic.

All independent buckets must give the same cancellation verdict. Conflicting
bucket verdicts violate the symmetry contract and raise:

```text
DiffractionInvariantError
code = diffraction.extinction.inconsistent_phase_buckets
```

The exception evidence contains the Miller index, setting ID, bucket reciprocal
vectors, operation IDs, exact translations and phases, and bucket verdicts.

For a consistent operation set:

```text
ExtinctionResult.absent = bucket cancellation verdict
```

No atom, scattering factor, expected reflection list, or hard-coded table of
space-group extinction conditions participates in the decision.

### 10.2 Exact evidence

```text
PhaseBucketEvidence
|- transformed_hkl: MillerIndex
|- operation_ids: tuple[str, ...]
|- translation_parts: tuple[tuple[Fraction, Fraction, Fraction], ...]
|- exact_phases: tuple[Fraction, ...]
|- relative_phases: tuple[Fraction, ...]
`- cancels: bool
```

Evidence is immutable, deterministically ordered, and exact. Extinction uses no
floating-point tolerance.

```text
ExtinctionResult
|- absent: bool
|- causes: tuple[ExtinctionCause, ...]
`- evidence: tuple[PhaseBucketEvidence, ...]
```

`absent` is always a strict boolean. There is no indeterminate extinction
state; mathematically inconsistent input produces an exception.

### 10.3 Secondary cause classification

Cause classification is explanatory and cannot affect `absent`. It examines
relative operations `(R_rel, t_rel)` responsible for non-trivial phase
characters.

```text
R_rel = identity
+ non-lattice fractional translation
-> CENTERING

det(R_rel) = +1
+ one-dimensional invariant subspace
+ non-zero intrinsic translation along that subspace
-> SCREW_AXIS

det(R_rel) = -1
+ two-dimensional invariant plane
+ non-zero intrinsic translation in that plane
-> GLIDE_PLANE

multiple mechanisms, multiple valid interpretations, or non-unique geometry
-> COMBINED
```

Classification uses the origin-invariant component of relative translation
inside the invariant subspace of the relative rotation, not the raw translation
alone. If a unique classification cannot be proved, `COMBINED` is mandatory.

```text
ExtinctionCause
|- kind: ExtinctionCauseKind
|- operation_ids: tuple[str, ...]
|- evidence: tuple[PhaseBucketEvidence, ...]
`- condition: str
```

`condition` is a deterministic human-readable modular relation derived from
the exact relative translation, such as `h dot delta_t is an integer`; it may
be simplified to an equivalent conventional form such as `00l: l = 2n` when
that simplification is unambiguous. It is never parsed or used as the source of
the extinction decision. A `COMBINED` cause retains the full operation IDs,
phase relations, translation parts, affected reciprocal orbit, and derived
condition through its referenced evidence.

Extinction is calculated from the representative index. Tests require every
member of its reciprocal orbit to produce the same `absent` verdict. Evidence
is canonicalized relative to the representative and need not have identical
orientation-dependent presentation for another orbit member.

## 11. Friedel linking

After crystallographic orbits are built, the generator negates the canonical
representative and constructs its actual point-group orbit.

- If the negated index canonicalizes to the same orbit, the reflection is
  self-Friedel and `friedel_mate_id` is `None`.
- Otherwise the mate remains a separate `Reflection`, and the two stable IDs
  are linked bidirectionally.

The generator materializes an accepted orbit and its distinct Friedel orbit as
one atomic unit. This preserves the Friedel contract in both complete and
incomplete results without globally merging the two orbits.

No assumption about anomalous scattering or powder-line equivalence is made in
this milestone.

## 12. Diagnostics, invariant errors, and provenance

Expected controlled truncation uses ordinary CrIStMa diagnostics:

```text
diffraction.reflections.search_limit_reached
```

Invalid scalar configuration or API argument types raise `ValueError` or
`TypeError`. Mathematical contradictions raise `DiffractionInvariantError`,
which carries a stable code and immutable structured evidence. Initial codes
include:

```text
diffraction.reflections.incompatible_cell_and_setting
diffraction.reflections.non_integral_reciprocal_action
diffraction.reflections.orbit_metric_mismatch
diffraction.extinction.inconsistent_phase_buckets
```

These exceptions are not converted to `INCOMPLETE` results.

`ReflectionGenerationProvenance` contains:

```text
generator method and version
space_group_setting_id
cell_fingerprint
d_min
reciprocal convention
boundary_tolerance
metric_compatibility_tolerance
max_candidates
integer_points_tested
reflections_within_d_min
orbits_created
completeness status
```

`reflections_within_d_min` counts tested non-zero integer triples accepted by
the original quadratic-form check before symmetry/Friedel closure.
`orbits_created` counts final `Reflection` objects.

The cell fingerprint is a SHA-256 digest over a canonical numeric serialization
of the six finite values actually used by the calculation:

```text
a, b, c, alpha, beta, gamma
```

Reported text and uncertainties do not enter the first-milestone fingerprint,
because they do not enter the reciprocal calculation.

## 13. Validation rules

Public immutable models reject inconsistent construction. In addition to the
orbit rules above:

- IDs and required provenance fields must be non-empty;
- symmetry-operation IDs must be non-empty and unique within the setting;
- `d_min`, reciprocal norms, spacings, and tolerances must be finite and
  positive where mathematically required;
- candidate limits must be positive integers and reject `bool`;
- reflection and orbit IDs must be unique within a set;
- a complete result must not contain a search-limit diagnostic;
- an incomplete result must contain one;
- all non-null Friedel IDs must exist in the same set and be reciprocal;
- every reflection belongs to the declared setting ID;
- all exact phases and relative phases are normalized into `[0, 1)`;
- an allowed reflection has no extinction causes;
- an absent reflection has at least one complete structured cause.

## 14. Verification strategy

Tests begin with analytic cells and operation sets rather than expected lists
for named materials.

### 14.1 Reciprocal metric and enumeration

- cubic, orthorhombic, and triclinic reciprocal metrics;
- `reciprocal_norm == 1/d_hkl` without `2*pi`;
- comparison against exhaustive small-box enumeration for randomly selected
  positive-definite test metrics;
- boundary inclusion and exclusion under the documented tolerance;
- a highly skew cell that would make a rectangular search inefficient;
- a large real-space cell producing valid indices above 12;
- exact `max_candidates` counting and explicit incomplete diagnostics;
- exclusion of `(0, 0, 0)`.

### 14.2 Reciprocal orbits and identity

- primitive, body-centered, and face-centered settings;
- centrosymmetric and non-centrosymmetric point groups;
- exact `R^T h` action and integer validation;
- unique sorted orbit members and canonical representative;
- multiplicity equal to orbit length;
- stable reflection IDs when cell parameters change;
- distinct IDs for different settings;
- cell/setting metric incompatibility raises the invariant error.

### 14.3 Friedel relations

- inversion already inside the reciprocal orbit produces `None`;
- a distinct Friedel orbit remains separate and links bidirectionally;
- an incomplete result contains no dangling Friedel references.

### 14.4 Systematic extinctions

- primitive settings without translational extinction;
- I and F centering;
- analytic screw-axis and glide-plane examples;
- combined or ambiguous evidence classified as `COMBINED`;
- exact `Fraction` phases with no numerical tolerance;
- every member of a reciprocal orbit has the same absence verdict;
- no production branch uses a setting name, expected reflection list, or
  hard-coded extinction table;
- an intentionally damaged operation set raises
  `diffraction.extinction.inconsistent_phase_buckets`.

The packaged catalog audit runs every non-zero Miller index in `[-3, 3]^3`
across all 530 `SpaceGroupSetting` records and requires that none produces
conflicting phase-bucket verdicts. This domain exercises axes, planes, general
directions, and every residue class needed by the packaged translation
denominators 1, 2, 3, 4, and 6. It is a finite regression audit, not a claim to
enumerate the infinite reciprocal lattice; the general result follows from the
exact group-character algorithm.

### 14.5 Distribution and regressions

- public imports work from a built wheel;
- packaged catalog data remain available in the installed wheel;
- the full existing CrIStMa suite remains green;
- imports introduce no GUI, specialized crystallographic, or new runtime
  dependency.

## 15. Acceptance criteria

The milestone is complete when:

- `ReflectionGenerator` accepts exactly `UnitCell`, `SpaceGroupSetting`, and
  positive `d_min` as scientific inputs;
- generation is complete by the physical `d_min` ellipsoid, never a fixed
  Miller-index cap;
- every emitted reflection is one complete point-group orbit with deterministic
  identity and crystallographic multiplicity;
- Friedel-related orbits are self-contained or linked bidirectionally without
  unconditional sign merging;
- systematic absence is decided exactly from symmetry-operation phases;
- explanatory extinction classification cannot alter the absence verdict;
- search truncation is explicit and mathematically inconsistent states raise
  invariant errors;
- results and evidence are immutable, deterministic, and reproducible;
- analytic, catalog-wide, integration, full-suite, and installed-wheel tests
  pass;
- no excluded diffraction layer or application-specific behavior enters the
  implementation.

## 16. Subsequent independent milestones

Later work composes new result types over this layer without changing the
meaning of `ReflectionSet`:

```text
CrystalStructure
    +
ReflectionSet
    +
ScatteringContext
    |
    v
StructureFactorSet
```

Then:

```text
StructureFactorSet + Radiation + DiffractionGeometry
    -> PowderLineSet
```

Then:

```text
PowderLineSet + ProfileModel + CalculationGrid
    -> CalculatedProfile
```

`ReflectionSet`, `PowderLineSet`, and `CalculatedProfile` remain distinct
scientific levels and are never collapsed into one object.
