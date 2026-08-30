# Structure Model and Reader Infrastructure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish the stable CRiStMa structure, species, multi-model, provenance, source-resolution, and lazy format-registry contracts required by every native reader.

**Architecture:** Periodic `CrystalStructure` and non-periodic `MolecularStructure` remain distinct immutable models and expose a shared `AtomicView`. Finite structures use `StructureCollection`; large trajectories use an indexed lazy `StructureSequence`. I/O handlers map preserved format documents into those public types through a declarative, lazily loaded registry.

**Tech Stack:** Python 3.11+, NumPy, standard-library dataclasses/protocols/pathlib/threading/importlib.metadata/gzip/bz2/lzma; pytest; no mandatory specialized crystallographic, molecular, DFT, Qt, or application package.

**Spec:** `docs/superpowers/specs/2026-08-30-native-structure-readers-design.md`

## Global Constraints

- Public structure imports live under `cristma.structure`.
- Existing `Crystal`, tuple-like structure access, and CIF behavior remain compatible during migration.
- Canonical scientific objects are immutable snapshots.
- Unknown species and unknown properties remain explicit and retain provenance.
- Reader-created structures, independent sites, and expanded atoms have stable semantic IDs; labels are never used as identity by downstream topology.
- Source readers never scan neighboring files or resolve includes without an explicit `SourceResolver`.
- Format selection order is explicit format, content confidence, special basename, then suffix.
- Large frame sources are indexed without eager structure materialization.
- Use TDD and run only the focused tests named by each task.
- Do not import or run tests from XRD Manager, Crystal Blocks, Craft, Finder, Organic, or other applications.

---

## Planned file structure

```text
src/cristma/chemistry/species.py       typed known and unknown species
src/cristma/structure/__init__.py      stable public structure namespace
src/cristma/structure/identity.py      stable source and expanded-atom identity
src/cristma/structure/crystal.py       crystal model and compatibility alias
src/cristma/structure/molecular.py     molecule, atoms, bonds, groups
src/cristma/structure/properties.py    immutable typed atomic properties
src/cristma/structure/view.py          shared numerical AtomicView
src/cristma/structure/collection.py    finite and lazy multi-structure APIs
src/cristma/io/source.py               decoded source and safe source resolver
src/cristma/io/formats.py              lazy format descriptors
src/cristma/io/registry.py             selection and dispatch
src/cristma/io/result.py               new collection/sequence result contract
src/cristma/io/cif/probe.py            lightweight CIF content probe
```

### Task 1: Typed chemical species

**Files:**
- Create: `src/cristma/chemistry/species.py`
- Modify: `src/cristma/chemistry/__init__.py`
- Test: `tests/chemistry/test_species.py`

**Interfaces:**
- Consumes: `normalize_element(symbol: str) -> str` from `cristma.chemistry.elements`.
- Produces: `ChemicalSpecies`, `ElementSpecies`, `IsotopeSpecies`, `ChargedSpecies`, `UnknownSpecies`, and `as_species(value)`.

- [ ] **Step 1: Write failing species tests**

```python
import pytest

from cristma.chemistry.species import (
    ChargedSpecies,
    ElementSpecies,
    IsotopeSpecies,
    UnknownSpecies,
    as_species,
)


def test_string_becomes_normalized_element_species():
    species = as_species("si")
    assert species == ElementSpecies("Si")
    assert species.element == "Si"


def test_isotope_and_charge_keep_element_identity():
    assert IsotopeSpecies("C", 13).element == "C"
    assert ChargedSpecies("Fe", 3).label == "Fe3+"


def test_unknown_species_is_explicit_and_not_an_element():
    species = UnknownSpecies("species:1", source_label="type 1")
    assert species.element is None
    with pytest.raises(ValueError, match="known element"):
        species.require_element()
```

- [ ] **Step 2: Verify the new module is absent**

Run: `pytest -q tests/chemistry/test_species.py`

Expected: collection FAIL because `cristma.chemistry.species` does not exist.

- [ ] **Step 3: Implement immutable species values**

```python
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from .elements import normalize_element


@runtime_checkable
class ChemicalSpecies(Protocol):
    @property
    def label(self) -> str: ...

    @property
    def element(self) -> str | None: ...

    def require_element(self) -> str: ...


@dataclass(frozen=True, slots=True)
class ElementSpecies:
    symbol: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "symbol", normalize_element(self.symbol))

    @property
    def label(self) -> str:
        return self.symbol

    @property
    def element(self) -> str:
        return self.symbol

    def require_element(self) -> str:
        return self.symbol
```

Implement `IsotopeSpecies(element_symbol, mass_number)` with a positive integer
mass, `ChargedSpecies(element_symbol, oxidation)` with labels such as `Fe3+`
and `O2-`, and `UnknownSpecies(id, source_label=None)` whose
`require_element()` raises `ValueError("species has no known element")`.
`as_species()` accepts an existing protocol instance or converts a string to
`ElementSpecies`.

- [ ] **Step 4: Run the focused species tests**

Run: `pytest -q tests/chemistry/test_species.py`

Expected: `3 passed`.

- [ ] **Step 5: Commit typed species**

```bash
git add src/cristma/chemistry tests/chemistry/test_species.py
git commit -m "feat: add typed chemical species"
```

### Task 2: Public CrystalStructure and compatibility layer

**Files:**
- Create: `src/cristma/structure/__init__.py`
- Create: `src/cristma/structure/identity.py`
- Create: `src/cristma/structure/crystal.py`
- Modify: `src/cristma/core/structure.py`
- Modify: `src/cristma/symmetry/orbit.py`
- Test: `tests/structure/test_crystal.py`
- Modify: `tests/core/test_structure.py`

**Interfaces:**
- Consumes: Task 1 `ChemicalSpecies` and `as_species`; existing `UnitCell`, `MeasuredValue`, `SpaceGroupDefinition`, and `ExpandedSite`.
- Produces: `CrystalStructure`, compatibility alias `Crystal`, species-aware `SiteComponent`, `SourceReference`, `StructureProvenance`, `ExpandedAtomRef`, compatibility alias `ExpandedSite`, and per-axis periodicity.

- [ ] **Step 1: Write failing public crystal tests**

```python
from cristma.chemistry.species import ElementSpecies
from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
from cristma.structure import Crystal, CrystalStructure, IndependentSite, SiteComponent
from cristma.symmetry.affine import parse_xyz_operation
from cristma.symmetry.orbit import expand_orbit


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def test_public_crystal_name_and_compatibility_alias():
    site = IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(SiteComponent("Si", number(1)),),
        fractional=(number(0), number(0), number(0)),
    )
    structure = CrystalStructure("Si", UnitCell.cubic(number(5.43)), (site,))
    assert isinstance(structure, Crystal)
    assert structure.periodic == (True, True, True)
    assert structure.sites[0].components[0].species == ElementSpecies("Si")
    assert structure.sites[0].components[0].element == "Si"


def test_identity_for_explicit_dft_sites_is_not_reported_p1():
    structure = CrystalStructure.explicit(
        name="simulation",
        cell=UnitCell.cubic(number(4)),
        sites=(),
    )
    assert structure.space_group.provenance == "unreported_identity"
    assert structure.space_group.number is None


def test_expanded_atom_identity_resolves_to_independent_site():
    site = IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(SiteComponent("Si", number(1)),),
        fractional=(number(0), number(0), number(0)),
    )
    atom = expand_orbit(
        site,
        (parse_xyz_operation("x,y,z", operation_id="op:1"),),
        structure_id="structure:Si",
    )[0]
    assert atom.structure_id == "structure:Si"
    assert atom.source_site_id == "site:Si1"
    assert atom.independent_site_id == "site:Si1"
    assert atom.id == "expanded:structure:Si:site:Si1:op:1:0,0,0"
```

- [ ] **Step 2: Verify public imports fail**

Run: `pytest -q tests/structure/test_crystal.py`

Expected: collection FAIL because `cristma.structure` does not exist.

- [ ] **Step 3: Move the canonical definitions behind the public namespace**

Define in `cristma.structure.crystal`:

```python
@dataclass(frozen=True, slots=True)
class SourceReference:
    source_name: str | None = None
    format: str | None = None
    record_id: str | None = None
    start_offset: int | None = None
    end_offset: int | None = None


@dataclass(frozen=True, slots=True)
class StructureProvenance:
    source: SourceReference | None = None
    parent_structure_id: str | None = None


@dataclass(frozen=True, slots=True)
class ExpandedAtomRef:
    id: str
    structure_id: str | None
    fractional: tuple[float, float, float]
    source_site_id: str
    representative_operation_id: str
    equivalent_operation_ids: tuple[str, ...]
    cell_translation: tuple[int, int, int]

    @property
    def independent_site_id(self) -> str:
        return self.source_site_id


@dataclass(frozen=True, slots=True)
class SiteComponent:
    species: ChemicalSpecies | str
    occupancy: MeasuredValue
    oxidation_state: MeasuredValue | None = None
    metadata: Mapping[str, object] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "species", as_species(self.species))
        # Retain the existing finite, non-negative occupancy validation.

    @property
    def element(self) -> str | None:
        return self.species.element


@dataclass(frozen=True, slots=True)
class CrystalStructure:
    name: str
    cell: UnitCell
    sites: tuple[IndependentSite, ...]
    id: str | None = None
    space_group: SpaceGroupDefinition | None = None
    formula: str | None = None
    periodic: tuple[bool, bool, bool] = (True, True, True)
    provenance: StructureProvenance = field(default_factory=StructureProvenance)
    metadata: Mapping[str, object] = field(default_factory=dict, compare=False)
    expanded_sites: tuple[ExpandedSite, ...] | None = field(default=None, compare=False)


Crystal = CrystalStructure
ExpandedSite = ExpandedAtomRef
```

`CrystalStructure.explicit()` constructs an exact identity operation with no
space-group number or symbol and provenance `unreported_identity`. Extend
`SymmetryProvenance` with that literal. Reject a crystal with no periodic axis.
Keep `id` optional for compatibility with manually constructed objects, but
require every format mapper to assign a stable document-derived structure ID.
Extend `expand_orbit(site, operations, *, structure_id=None)` to create the
documented deterministic expanded ID; preserve all equivalent operation IDs
when special-position images merge.
Make `cristma.core.structure` re-export the public definitions so current CIF
imports and application code keep working.

- [ ] **Step 4: Run crystal and existing core/symmetry tests**

Run: `pytest -q tests/structure/test_crystal.py tests/core/test_structure.py tests/symmetry`

Expected: all selected tests PASS.

- [ ] **Step 5: Commit the public crystal model**

```bash
git add src/cristma/structure src/cristma/core/structure.py src/cristma/symmetry/orbit.py tests/structure tests/core/test_structure.py
git commit -m "refactor: expose public crystal structure model"
```

### Task 3: MolecularStructure, atoms, bonds, and groups

**Files:**
- Create: `src/cristma/structure/molecular.py`
- Modify: `src/cristma/structure/__init__.py`
- Test: `tests/structure/test_molecular.py`

**Interfaces:**
- Consumes: Task 1 species and Task 2 `StructureProvenance`.
- Produces: `MolecularAtom`, `MolecularBond`, `MolecularGroup`, `MolecularStructure`, and runtime-checkable `Structure` protocol.

- [ ] **Step 1: Write failing molecular model tests**

```python
import pytest

from cristma.structure import MolecularAtom, MolecularBond, MolecularGroup, MolecularStructure, Structure


def test_molecule_has_no_artificial_periodic_cell():
    atoms = (
        MolecularAtom("atom:C1", "C1", "C", (0.0, 0.0, 0.0)),
        MolecularAtom("atom:O1", "O1", "O", (1.2, 0.0, 0.0)),
    )
    molecule = MolecularStructure(
        name="CO",
        atoms=atoms,
        bonds=(MolecularBond("bond:1", "atom:C1", "atom:O1", order=2.0),),
    )
    assert isinstance(molecule, Structure)
    assert molecule.cell is None
    assert molecule.periodic == (False, False, False)


def test_bond_rejects_missing_atom_identity():
    atom = MolecularAtom("atom:C1", "C1", "C", (0.0, 0.0, 0.0))
    with pytest.raises(ValueError, match="unknown atom"):
        MolecularStructure(
            "bad",
            atoms=(atom,),
            bonds=(MolecularBond("bond:1", "atom:C1", "atom:X", 1.0),),
        )
```

- [ ] **Step 2: Verify molecular imports fail**

Run: `pytest -q tests/structure/test_molecular.py`

Expected: collection FAIL for missing public molecular names.

- [ ] **Step 3: Implement immutable molecular records**

```python
@runtime_checkable
class Structure(Protocol):
    name: str
    periodic: tuple[bool, bool, bool]
    provenance: StructureProvenance


@dataclass(frozen=True, slots=True)
class MolecularAtom:
    id: str
    label: str
    species: ChemicalSpecies | str
    cartesian: tuple[float, float, float]
    occupancy: float = 1.0
    metadata: Mapping[str, object] = field(default_factory=dict, compare=False)


@dataclass(frozen=True, slots=True)
class MolecularBond:
    id: str
    atom1_id: str
    atom2_id: str
    order: float | str
    aromatic: bool = False
    stereo: str | None = None


@dataclass(frozen=True, slots=True)
class MolecularStructure:
    name: str
    atoms: tuple[MolecularAtom, ...]
    bonds: tuple[MolecularBond, ...] = ()
    groups: tuple[MolecularGroup, ...] = ()
    periodic: tuple[bool, bool, bool] = (False, False, False)
    provenance: StructureProvenance = field(default_factory=StructureProvenance)
```

Validate unique IDs, finite coordinates, occupancies in `[0, 1]`, bonds and
groups referencing existing atom IDs, and all-false periodicity. Provide
`cell = None` as a property.

- [ ] **Step 4: Run molecular tests**

Run: `pytest -q tests/structure/test_molecular.py`

Expected: both tests PASS.

- [ ] **Step 5: Commit molecular structures**

```bash
git add src/cristma/structure tests/structure/test_molecular.py
git commit -m "feat: add canonical molecular structures"
```

### Task 4: Typed atomic properties and shared AtomicView

**Files:**
- Create: `src/cristma/structure/properties.py`
- Create: `src/cristma/structure/view.py`
- Modify: `src/cristma/structure/crystal.py`
- Modify: `src/cristma/structure/molecular.py`
- Modify: `src/cristma/structure/__init__.py`
- Test: `tests/structure/test_properties.py`
- Test: `tests/structure/test_atomic_view.py`

**Interfaces:**
- Consumes: Task 2 and 3 structure/site types and NumPy.
- Produces: `PropertyProvenance`, `AtomicProperty`, `AtomicPropertyTable`, `AtomicView`, and both structures' `atomic_view()`.

- [ ] **Step 1: Write failing property-table tests**

```python
import numpy as np
import pytest

from cristma.structure import AtomicProperty, AtomicPropertyTable


def test_property_values_are_immutable_and_typed():
    prop = AtomicProperty("magnetic_moment", np.array([1.0, -1.0]), unit="mu_B")
    table = AtomicPropertyTable(2, (prop,))
    assert table["magnetic_moment"].unit == "mu_B"
    with pytest.raises(ValueError):
        table["magnetic_moment"].values[0] = 0


def test_property_length_must_match_atoms():
    with pytest.raises(ValueError, match="leading dimension"):
        AtomicPropertyTable(2, (AtomicProperty("charge", np.array([0.0])),))
```

- [ ] **Step 2: Write failing atomic-view tests**

```python
import numpy as np

from cristma.structure import MolecularAtom, MolecularStructure


def test_molecular_atomic_view_has_cartesian_coordinates_and_no_cell():
    molecule = MolecularStructure(
        "water",
        atoms=(MolecularAtom("atom:O", "O", "O", (1.0, 2.0, 3.0)),),
    )
    view = molecule.atomic_view()
    assert np.array_equal(view.cartesian, [[1.0, 2.0, 3.0]])
    assert view.fractional is None
    assert view.cell is None
```

- [ ] **Step 3: Verify both modules are absent**

Run: `pytest -q tests/structure/test_properties.py tests/structure/test_atomic_view.py`

Expected: collection FAIL for missing names.

- [ ] **Step 4: Implement typed immutable arrays and views**

```python
@dataclass(frozen=True, slots=True)
class AtomicProperty:
    name: str
    values: np.ndarray
    unit: str | None = None
    missing: np.ndarray | None = None
    source_name: str | None = None
    provenance: PropertyProvenance = field(default_factory=PropertyProvenance)

    def __post_init__(self) -> None:
        values = np.array(self.values, copy=True)
        values.flags.writeable = False
        object.__setattr__(self, "values", values)


@dataclass(frozen=True, slots=True)
class AtomicView:
    ids: tuple[str, ...]
    species: tuple[ChemicalSpecies, ...]
    cartesian: np.ndarray
    fractional: np.ndarray | None
    cell: np.ndarray | None
    periodic: tuple[bool, bool, bool]
    properties: AtomicPropertyTable
    source_site_ids: tuple[str | None, ...]
```

`AtomicPropertyTable(atom_count, properties=())` rejects duplicate names and
arrays whose leading dimension differs from `atom_count`. It implements
`Mapping[str, AtomicProperty]`.

`MolecularStructure.atomic_view()` maps atoms directly. For
`CrystalStructure.atomic_view(expanded=False)`, use independent fractional
positions; for `expanded=True`, use `ExpandedSite` and its
`independent_site_id`, map each expanded site back to its components, and
calculate row-vector Cartesian positions as `fractional @ cell.matrix`.
Mixed sites remain one geometric row with their component tuple recorded in a
reserved `site_components` object property; do not duplicate them into fully
occupied atoms.

Now that `AtomicView` exists, extend the public runtime `Structure` protocol
with `atomic_view(self, *, expanded: bool = True) -> AtomicView`. The molecular
implementation accepts the keyword and ignores it because molecules have no
symmetry expansion.

- [ ] **Step 5: Run focused view/property and existing orbit tests**

Run: `pytest -q tests/structure/test_properties.py tests/structure/test_atomic_view.py tests/symmetry/test_orbit.py`

Expected: all selected tests PASS.

- [ ] **Step 6: Commit atomic views**

```bash
git add src/cristma/structure tests/structure
git commit -m "feat: add typed atomic views and properties"
```

### Task 5: Finite collections and lazy structure sequences

**Files:**
- Create: `src/cristma/structure/collection.py`
- Modify: `src/cristma/structure/__init__.py`
- Test: `tests/structure/test_collection.py`
- Test: `tests/structure/test_sequence.py`

**Interfaces:**
- Consumes: Task 3 runtime `Structure` protocol and Task 2 `SourceReference`.
- Produces: `StructureRole`, `StructureEntry`, `StructureCollection`, `FrameReference`, `StructureSequence`, and `StructureSeries` protocol.

- [ ] **Step 1: Write failing finite-collection tests**

```python
from cristma.structure import MolecularAtom, MolecularStructure, StructureCollection


def molecule(name: str) -> MolecularStructure:
    return MolecularStructure(name, (MolecularAtom(f"atom:{name}", name, "H", (0, 0, 0)),))


def test_collection_is_sequence_with_primary_and_final_roles():
    first, last = molecule("first"), molecule("last")
    collection = StructureCollection.from_structures(
        (first, last), primary_index=0, final_index=1
    )
    assert tuple(collection) == (first, last)
    assert collection.primary is first
    assert collection.final is last
```

- [ ] **Step 2: Write failing lazy-sequence tests**

```python
from cristma.structure import FrameReference, StructureSequence


def test_sequence_does_not_materialize_frames_until_requested():
    loaded = []
    refs = (FrameReference(0), FrameReference(1, role="final"))

    def load(reference):
        loaded.append(reference.index)
        return molecule(str(reference.index))

    sequence = StructureSequence(refs, load)
    assert loaded == []
    assert sequence.final.name == "1"
    assert loaded == [1]
    assert sequence.final is sequence.final
    assert loaded == [1]
```

- [ ] **Step 3: Verify collection imports fail**

Run: `pytest -q tests/structure/test_collection.py tests/structure/test_sequence.py`

Expected: collection FAIL for missing names.

- [ ] **Step 4: Implement sequence-compatible immutable containers**

```python
StructureRole = Literal["model", "primary", "intermediate", "final"]


@dataclass(frozen=True, slots=True)
class StructureEntry:
    structure: Structure
    role: StructureRole = "model"
    source_index: int | None = None
    source: SourceReference | None = None


@dataclass(frozen=True, slots=True)
class FrameReference:
    index: int
    role: StructureRole = "intermediate"
    source: SourceReference | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)
```

Define `StructureSeries` as a runtime-checkable protocol extending
`Sequence[Structure]` with read-only `primary: Structure | None` and
`final: Structure | None` properties.

`StructureCollection` stores entries, implements `Sequence[Structure]`, and
validates at most one primary and one final role. `from_structures()` assigns
roles by explicit indexes. `StructureSequence` stores references and a loader,
implements integer/slice indexing, and uses a private dictionary plus
`threading.RLock` to cache each successfully loaded frame exactly once. Slices
return `StructureCollection`. Loader failures are not cached.

- [ ] **Step 5: Run collection/sequence tests**

Run: `pytest -q tests/structure/test_collection.py tests/structure/test_sequence.py`

Expected: both focused files PASS.

- [ ] **Step 6: Commit collections and sequences**

```bash
git add src/cristma/structure tests/structure/test_collection.py tests/structure/test_sequence.py
git commit -m "feat: add structure collections and lazy sequences"
```

### Task 6: ReadResult migration and CIF compatibility

**Files:**
- Modify: `src/cristma/io/result.py`
- Modify: `src/cristma/io/cif/handler.py`
- Modify: `src/cristma/io/cif/mapper.py`
- Modify: `src/cristma/io/cif/writer.py`
- Modify: `src/cristma/__init__.py`
- Test: `tests/io/test_result.py`
- Modify: `tests/io/cif/conftest.py`
- Modify: `tests/io/cif/test_end_to_end.py`

**Interfaces:**
- Consumes: Task 5 `StructureCollection`, `StructureSequence`, and `StructureSeries`; Task 2 `CrystalStructure`.
- Produces: `ReadResult.structures: StructureCollection | StructureSequence` while preserving iteration/index behavior and current `cristma.read/write` behavior.

- [ ] **Step 1: Write failing result compatibility tests**

```python
from cristma.io.result import ReadResult
from cristma.structure import StructureCollection, StructureSequence


def test_result_converts_legacy_structure_tuple_to_collection():
    result = ReadResult(document=None, structures=())
    assert isinstance(result.structures, StructureCollection)
    assert len(result.structures) == 0


def test_result_keeps_explicit_lazy_sequence():
    sequence = StructureSequence((), lambda reference: None)
    assert ReadResult(document=None, structures=sequence).structures is sequence
```

- [ ] **Step 2: Verify the result still exposes tuples**

Run: `pytest -q tests/io/test_result.py`

Expected: FAIL because empty structures are currently a tuple.

- [ ] **Step 3: Implement result normalization and migrate CIF types**

```python
@dataclass(frozen=True, slots=True)
class ReadResult:
    document: object | None
    structures: StructureCollection | StructureSequence | tuple[Structure, ...] = ()
    diagnostics: tuple[Diagnostic, ...] = ()
    source_info: SourceInfo | None = None

    def __post_init__(self) -> None:
        if isinstance(self.structures, tuple):
            object.__setattr__(
                self,
                "structures",
                StructureCollection.from_structures(self.structures),
            )
```

Change CIF mapper annotations to `tuple[CrystalStructure, ...]`, assign each
mapped structure `id=f"cif:{block.name}"`, pass that ID into every
`expand_orbit()` call, wrap mapped
tuples in `StructureCollection.from_structures()` in the handler, and make the
canonical writer accept `CrystalStructure`. Existing indexing, truth testing,
and loops remain unchanged.

- [ ] **Step 4: Run result and complete CIF tests**

Run: `pytest -q tests/io/test_result.py tests/io/cif tests/core tests/symmetry`

Expected: all selected CRiStMa tests PASS.

- [ ] **Step 5: Commit result integration**

```bash
git add src/cristma tests/io/test_result.py tests/io/cif
git commit -m "refactor: return canonical structure collections"
```

### Task 7: Decoded sources, compression, and explicit SourceResolver

**Files:**
- Create: `src/cristma/io/source.py`
- Modify: `src/cristma/io/__init__.py`
- Test: `tests/io/test_source.py`

**Interfaces:**
- Consumes: `Diagnostic`, `Severity`, and `SourceInfo`.
- Produces: `DecodedSource`, `ResolvedSource`, `SourceResolver`, `MappingSourceResolver`, `decode_source(path)`, and `decode_bytes(raw, source_name)`.

- [ ] **Step 1: Write failing decoding and resolver tests**

```python
import gzip

from cristma.io.source import MappingSourceResolver, decode_source


def test_gzip_source_is_decoded_with_inner_suffix(tmp_path):
    path = tmp_path / "POSCAR.gz"
    path.write_bytes(gzip.compress(b"title\n1.0\n"))
    source = decode_source(path)
    assert source.text == "title\n1.0\n"
    assert source.logical_name.endswith("POSCAR")
    assert source.compression == "gzip"


def test_mapping_resolver_is_explicit_and_blocks_parent_escape():
    resolver = MappingSourceResolver({"POTCAR": b"TITEL = PAW_PBE Si"})
    assert resolver.resolve("POTCAR", from_source="POSCAR").raw.startswith(b"TITEL")
    assert resolver.resolve("../POTCAR", from_source="POSCAR") is None
```

- [ ] **Step 2: Verify source infrastructure is absent**

Run: `pytest -q tests/io/test_source.py`

Expected: collection FAIL because `cristma.io.source` does not exist.

- [ ] **Step 3: Implement deterministic decoding and safe resolution**

```python
@dataclass(frozen=True, slots=True)
class DecodedSource:
    raw: bytes
    text: str
    logical_name: str | None
    encoding: str
    newline: str
    compression: Literal["gzip", "bzip2", "xz"] | None
    diagnostics: tuple[Diagnostic, ...] = ()


@dataclass(frozen=True, slots=True)
class ResolvedSource:
    reference: str
    raw: bytes
    source_name: str | None = None


@runtime_checkable
class SourceResolver(Protocol):
    def resolve(
        self,
        reference: str,
        *,
        from_source: str | None,
    ) -> ResolvedSource | None: ...
```

Detect gzip/bzip2/xz by magic bytes, strip only the compression suffix from
`logical_name`, decode UTF-8 with BOM before Latin-1 fallback, and reuse the
existing stable encoding diagnostic. `MappingSourceResolver` accepts an
explicit immutable name-to-bytes mapping, normalizes POSIX relative names, and
rejects absolute paths and any `..` component.

- [ ] **Step 4: Run source tests**

Run: `pytest -q tests/io/test_source.py`

Expected: both tests PASS.

- [ ] **Step 5: Commit source infrastructure**

```bash
git add src/cristma/io/source.py src/cristma/io/__init__.py tests/io/test_source.py
git commit -m "feat: add explicit reader source infrastructure"
```

### Task 8: Declarative lazy format descriptors

**Files:**
- Create: `src/cristma/io/formats.py`
- Modify: `src/cristma/io/registry.py`
- Create: `src/cristma/io/cif/probe.py`
- Modify: `src/cristma/io/cif/handler.py`
- Modify: `src/cristma/__init__.py`
- Modify: `tests/io/test_registry.py`
- Test: `tests/io/test_builtin_formats.py`

**Interfaces:**
- Consumes: Task 7 `DecodedSource` and `decode_source`; current handler `read_text` contract.
- Produces: `FormatDescriptor`, `FormatCapabilities`, `builtin_format_descriptors()`, lazy `FormatRegistry`, basename-aware selection, and a lazily registered built-in CIF descriptor.

- [ ] **Step 1: Write failing lazy-selection tests in `tests/io/test_registry.py`**

```python
from cristma.io.formats import FormatCapabilities, FormatDescriptor
from cristma.io.registry import FormatRegistry


def test_special_basename_beats_uninformative_suffix_without_loading_handler():
    loaded = []
    descriptor = FormatDescriptor(
        name="vasp",
        aliases=("poscar",),
        suffixes=(),
        basenames=("POSCAR", "CONTCAR"),
        probe=lambda source: 0.0,
        factory=lambda: loaded.append(True),
        capabilities=FormatCapabilities(text=True),
    )
    registry = FormatRegistry((descriptor,))
    selected = registry.select("title", basename="POSCAR")
    assert selected.name == "vasp"
    assert loaded == []


def test_handler_factory_runs_only_when_selected_reader_is_used():
    loaded = []
    handler = StubHandler()
    descriptor = FormatDescriptor(
        name="stub",
        aliases=(),
        suffixes=(".stub",),
        basenames=(),
        probe=handler.probe,
        factory=lambda: loaded.append(handler) or handler,
        capabilities=FormatCapabilities(text=True),
    )
    registry = FormatRegistry((descriptor,))
    assert loaded == []
    assert registry.read_text("STUB value") == (None, "STUB value")
    assert loaded == [handler]
```

In `tests/io/test_builtin_formats.py`, add:

```python
from cristma.io.formats import builtin_format_descriptors
from cristma.io.registry import FormatRegistry


def test_builtin_cif_descriptor_is_content_aware_and_multiple():
    registry = FormatRegistry(builtin_format_descriptors())
    descriptor = registry.select("data_demo\n_tag value\n")
    assert descriptor.name == "cif"
    assert descriptor.capabilities.multiple
```

- [ ] **Step 2: Verify descriptor imports fail**

Run: `pytest -q tests/io/test_registry.py tests/io/test_builtin_formats.py`

Expected: collection FAIL because `cristma.io.formats` does not exist.

- [ ] **Step 3: Implement descriptor-first selection**

```python
@dataclass(frozen=True, slots=True)
class FormatCapabilities:
    text: bool = True
    binary: bool = False
    multiple: bool = False
    lazy_frames: bool = False


@dataclass(frozen=True, slots=True)
class FormatDescriptor:
    name: str
    aliases: tuple[str, ...]
    suffixes: tuple[str, ...]
    basenames: tuple[str, ...]
    probe: Callable[[str], float]
    factory: Callable[[], FormatHandler]
    capabilities: FormatCapabilities
```

Store the loaded handler in a registry-owned cache rather than mutating the
frozen descriptor. Selection returns a descriptor and validates probe scores.
Explicit names match canonical names or aliases. Score candidates as content
confidence first, then basename confidence `0.7`, then suffix confidence `0.6`;
equal best candidates raise the existing ambiguity error. Calling `read` or
`read_text` loads only the selected handler.

Retain compatibility by accepting old handler instances in `FormatRegistry`
and converting them with `descriptor_for()`. Define the built-in CIF probe as a
lightweight function in `cif/probe.py`; it may import the lexer but must not
import the parser, mapper, writer, or handler. Register the handler factory in
`builtin_format_descriptors()` and construct the public registry from that
function in `cristma.__init__`.
Move path decoding to Task 7 `decode_source()` and propagate its diagnostics and
`SourceInfo`.

- [ ] **Step 4: Run registry, source, and CIF public API tests**

Run: `pytest -q tests/io/test_registry.py tests/io/test_builtin_formats.py tests/io/test_source.py tests/io/cif/test_end_to_end.py`

Expected: all selected tests PASS.

- [ ] **Step 5: Commit the lazy registry**

```bash
git add src/cristma/io src/cristma/__init__.py tests/io
git commit -m "refactor: add lazy native format registry"
```

### Task 9: Foundation integration and package gate

**Files:**
- Modify: `README.md`
- Modify: `tests/test_public_api.py`
- Create: `tests/structure/test_public_api.py`
- Modify: `pyproject.toml` only if test discovery needs no source changes

**Interfaces:**
- Consumes: all prior tasks.
- Produces: documented stable foundation ready for every format plan in the reader roadmap.

- [ ] **Step 1: Add public import and installed-package tests**

```python
def test_stable_structure_namespace_is_importable():
    from cristma.structure import (
        AtomicView,
        CrystalStructure,
        MolecularStructure,
        StructureCollection,
        StructureSequence,
    )

    assert all(
        value is not None
        for value in (
            AtomicView,
            CrystalStructure,
            MolecularStructure,
            StructureCollection,
            StructureSequence,
        )
    )
```

Add a README example importing `CrystalStructure`, reading CIF into a
`StructureCollection`, inspecting `.primary`, and explaining that trajectories
will use `StructureSequence`.

- [ ] **Step 2: Run the foundation slice only**

Run: `pytest -q tests/chemistry tests/structure tests/core tests/symmetry tests/io`

Expected: all foundation and existing CIF tests PASS.

- [ ] **Step 3: Build and install the wheel in a fresh environment**

Run:

```bash
python3 -m venv .venv-foundation-check
.venv-foundation-check/bin/pip install '.[test]'
.venv-foundation-check/bin/python -m build
.venv-foundation-check/bin/pip install --force-reinstall --no-deps dist/cristma-0.1.0-py3-none-any.whl
.venv-foundation-check/bin/python -c "from cristma.structure import CrystalStructure, StructureSequence; print(CrystalStructure.__name__, StructureSequence.__name__)"
```

Expected output ends with `CrystalStructure StructureSequence`.

- [ ] **Step 4: Audit the built wheel and repository state**

Run:

```bash
unzip -l dist/cristma-0.1.0-py3-none-any.whl
unzip -p dist/cristma-0.1.0-py3-none-any.whl cristma-0.1.0.dist-info/METADATA
git diff --check
git status --short
```

Expected: wheel contains the new `cristma/structure` and I/O modules; runtime
metadata contains only NumPy; Git status contains only the intentional README,
test, and metadata changes for this task.

- [ ] **Step 5: Commit the verified reader foundation**

```bash
git add README.md tests/test_public_api.py tests/structure pyproject.toml
git commit -m "docs: publish structure reader foundation"
```

## Foundation completion gate

Before starting the SHELX implementation plan:

- run `pytest -q tests/chemistry tests/structure tests/core tests/symmetry tests/io`;
- confirm existing real CIF preserve and canonical round trips still pass;
- confirm lazy sequence tests prove no eager frame loading;
- confirm registry tests prove unselected handlers are not imported;
- confirm old `Crystal` imports and new `CrystalStructure` imports refer to the
  same class;
- confirm clean wheel installation and dependency metadata;
- run `git status --short` and resolve only intentional changes.
