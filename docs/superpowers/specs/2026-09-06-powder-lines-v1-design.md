# Powder diffraction lines v1

## Scope

This milestone converts an existing `StructureFactorSet` into intrinsic powder
lines for an explicitly supplied X-ray spectrum. It belongs to CRiStMa's
forward crystallographic physics layer.

The milestone includes:

- one or more monochromatic radiation components;
- a built-in Cu K-alpha doublet preset;
- Bragg angles for every observable component;
- crystallographic multiplicity;
- Friedel-pair grouping;
- deterministic immutable results.

It does not include peak profiles, peak widths, instrumental broadening,
background, preferred orientation, absorption, Lorentz-polarization
corrections, experimental peak matching, similarity metrics, or refinement.

## Public flow

```text
StructureFactorSet + RadiationSpectrum
        -> PowderLineCalculator
        -> PowderLineSet
```

`PowderLineCalculator` does not receive a structure or a space-group setting.
Those inputs have already been validated and recorded by `StructureFactorSet`.
It does not know which application will consume the result.

## Radiation model

`RadiationComponent` is an immutable value containing:

- `component_id`: a stable non-empty identifier;
- `label`: a display label such as `Cu Kalpha1`;
- `wavelength_angstrom`: a finite positive wavelength;
- `relative_weight`: a finite positive relative intensity.

`RadiationSpectrum` contains a non-empty ordered tuple of components and
provenance. Component identifiers must be unique. Public calculations use
normalized weights `relative_weight / sum(relative_weight)`, so multiplying
all input weights by the same constant does not change a result. The original
relative weights remain available for provenance.

The built-in Cu K-alpha preset contains separate K-alpha1 and K-alpha2
components. Its immutable packaged values are generated from a pinned xraylib
release and commit using line energies and radiative rates, then stored with
the source version, conversion convention, checksum, and license notice. Its
provenance explicitly records `energy_source`, `radiative_rate_source`,
`energy_to_wavelength_formula`, `hc_value`, `hc_units`, `xraylib_version`,
`xraylib_commit`, and `resource_checksum`. xraylib is a reference-data build
dependency only and is not imported at runtime.

## Powder families and Friedel grouping

A powder family contains either one crystallographic reflection orbit or the
two orbits connected by `friedel_mate_id`. Each pair is consumed once, using a
deterministic family identifier derived from the sorted reflection IDs.

The calculator defensively validates Friedel links even though `ReflectionSet`
already enforces them at model construction. A non-null mate ID must identify
a member of the source `StructureFactorSet`, and the two links must be
reciprocal. A missing mate raises `DiffractionInvariantError` with code
`diffraction.powder.missing_friedel_mate`; a non-reciprocal link raises the
corresponding `diffraction.powder.nonreciprocal_friedel_link` invariant error.
Neither case may degrade to a singleton family.

For family members `i`, the intrinsic family strength is

```text
family_strength = sum(multiplicity_i * F_squared_i)
```

This sum remains valid if a later scattering model makes Friedel mates
unequal. The family multiplicity is the sum of member multiplicities. Members
must have the same d-spacing within the documented reciprocal numerical
tolerance; disagreement is a `DiffractionInvariantError`.

Systematically extinct reflections do not produce powder families. An allowed
reflection whose particular atomic structure gives zero `F_squared` remains a
valid zero-strength family; it is not reclassified as an extinction.

## Radiation components and Bragg angles

Each family produces one `PowderLine` per radiation component for which the
first-order Bragg condition is reachable:

```text
two_theta_deg = degrees(2 * asin(wavelength_angstrom / (2 * d_spacing)))
intrinsic_line_intensity = normalized_component_weight * family_strength
```

If `wavelength_angstrom / (2 * d_spacing) > 1`, that component produces no
line. This is a physical exclusion, not an incomplete calculation. A family
is omitted only when none of the supplied radiation components can reach it.
The provenance records considered reflections, emitted families and skipped
components.

No effective or averaged wavelength is used for a multi-component spectrum.
K-alpha1 and K-alpha2 therefore separate increasingly at high angles.

## Results and order

`PowderReflectionFamily` stores:

- `family_id`;
- member `reflection_ids` and representative Miller indices;
- `d_spacing`;
- total crystallographic multiplicity;
- `family_strength`;
- the component `PowderLine` values.

`PowderLine` stores:

- deterministic `line_id`;
- parent `family_id` and `radiation_component_id`;
- wavelength and normalized radiation weight;
- `two_theta_deg`;
- `intrinsic_line_intensity`.

`PowderLineSet` retains the source `StructureFactorSet`, spectrum, diagnostics,
provenance, and inherited `COMPLETE` or `INCOMPLETE` reflection-search status.
For every family, `family_sort_angle` is the minimum `two_theta_deg` among its
emitted lines. Families are ordered by `(family_sort_angle, family_id)`, so
reordering otherwise identical spectrum components cannot change family
order. Lines within a family follow the declared spectrum component order.
Derived flat and angle-sorted views may be provided without changing stored
identity.

The public scientific API consistently names this value
`intrinsic_line_intensity`, never simply `intensity`, because it excludes
Lorentz-polarization, preferred orientation, absorption, instrument response,
and profile integration.

## Errors and diagnostics

Invalid scalar values and malformed standalone models raise `TypeError` or
`ValueError`. Cross-object scientific contradictions raise
`DiffractionInvariantError` with `diffraction.powder.*` codes and evidence.
Normal Bragg exclusions and zero-strength allowed families are not errors.

## Files

```text
src/cristma/diffraction/powder_models.py
    radiation and powder result models

src/cristma/diffraction/powder.py
    family grouping and line calculation

src/cristma/diffraction/diagnostics.py
    powder invariant codes

src/cristma/diffraction/__init__.py
    intentional public exports

src/cristma/reference_data/resources/xray/
    immutable radiation preset data and provenance

tools/
    reproducible preset-data build tool
```

## Verification

Tests cover model validation, spectrum weight normalization, Cu K-alpha
provenance, an analytical Bragg-angle case, high-angle doublet separation,
multiplicity, Friedel grouping, extinction omission, allowed zero-strength
families, unreachable Bragg conditions, inherited incomplete status,
missing and non-reciprocal Friedel-link failures, spectrum-order-independent
family ordering, deterministic line ordering, package contents, and absence
of a runtime xraylib dependency.
