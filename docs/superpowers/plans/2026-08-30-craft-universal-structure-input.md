# CRAFT Universal CRiStMa Structure Input Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make CRiStMa the only structure-format boundary used by CRAFT, so every registered structure format is auto-detected, normalized to a canonical scientific object, and made available without adding parser or dispatch code to the application.

**Architecture:** Format-specific reader/document/mapper implementations live only in CRiStMa and converge on `CrystalStructure` or `MolecularStructure`. CRAFT calls `cristma.read(path)` once, adapts the returned canonical object for presentation, and retains the opaque source document and diagnostics. CRiStMa's registry is the single source of truth for supported suffixes and basenames. XPFF remains a separate application-project route.

**Tech Stack:** Python 3.11+, CRiStMa package, CRAFT's existing Python/PySide6 stack, pytest.

**Spec:** `/Users/artem/Yandex.Disk.localized/Python/XRD/CRiStMa/docs/superpowers/specs/2026-08-30-native-shelx-and-craft-integration-design.md`

## Global Constraints

- [ ] Execute the application cutover only after CRiStMa readers for every format currently advertised by CRAFT pass their focused gates: CIF/mmCIF, RES/INS, POSCAR/CONTCAR/VASP, PDB, and XYZ.
- [ ] CRAFT is currently not a git repository: do not initialize git; record explicit changed-file and verification checkpoints.
- [ ] CRAFT may special-case application project containers such as XPFF, but never a scientific structure format.
- [ ] When XPFF contains embedded structural text, CRAFT may extract it but must pass it to `cristma.read_text`; it must not parse that payload itself.
- [ ] Do not add format-specific parsing, occupancy, symmetry, chemistry, or normalization logic to CRAFT.
- [ ] A newly registered CRiStMa structure format must become loadable without editing CRAFT's loader dispatch or supported-extension constants.
- [ ] Run focused tests per task and the complete CRAFT suite once at the final gate.

---

## Task 1: Expose the CRiStMa format registry as package data

**Files:**

- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/CRiStMa/src/cristma/io/registry.py`
- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/CRiStMa/src/cristma/__init__.py`
- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/CRiStMa/tests/io/test_registry.py`
- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/CRiStMa/tests/test_public_api.py`

- [ ] Add failing tests for a read-only package API returning registered format information, supported suffixes, and supported basenames without exposing factories or instantiating lazy handlers.

```python
formats = cristma.structure_formats()
assert {item.name for item in formats} >= {"cif", "shelx", "vasp", "pdb", "xyz"}
assert ".res" in cristma.structure_suffixes()
assert "poscar" in cristma.structure_basenames()
```

- [ ] Run `pytest tests/io/test_registry.py tests/test_public_api.py -q` and confirm the API is missing.
- [ ] Add an immutable `StructureFormatInfo(name, aliases, suffixes, basenames, capabilities)` projection and tuple-returning query functions over the existing descriptor registry. Do not expose probes, factories, registry mutation, or handler instances.
- [ ] Add `can_read_structure(path)` that examines basename/suffix for UI prefiltering but leaves authoritative content detection to `read(path)`.
- [ ] Verify these calls do not import individual reader/parser modules.
- [ ] Run `pytest tests/io/test_registry.py tests/test_public_api.py -q`.
- [ ] Commit in CRiStMa: `git add src/cristma tests && git commit -m "feat(io): expose registered structure format capabilities"`

## Task 2: Generalize CRAFT's canonical presentation adapter

**Files:**

- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/вивер/src/crystal_viewer/core/model.py`
- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/вивер/src/crystal_viewer/adapters/cristma.py`
- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/вивер/tests/test_cristma_adapter.py`
- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/вивер/tests/test_document.py`

- [ ] Add failing tests proving that an arbitrary CRiStMa format document survives projection by identity and that CIF `CifSourceData` extraction still works unchanged.
- [ ] Add a failing molecular test: a CRiStMa `MolecularStructure` becomes a non-periodic CRAFT display model with a padded display cell, while the canonical object remains molecular.
- [ ] Run `pytest tests/test_cristma_adapter.py tests/test_document.py -q` and confirm both missing capabilities.
- [ ] Add opaque provenance/canonical fields to the temporary display model:

```python
canonical_structure: CrystalStructure | MolecularStructure | None
source_document: object | None = field(default=None, repr=False, compare=False)
```

- [ ] Split presentation projection by canonical type inside the adapter, not by source format. `CrystalStructure` uses its real cell/symmetry; `MolecularStructure` receives only a CRAFT display cell and retains non-periodic canonical identity.
- [ ] Call `source_data_from_cif_document()` only for a `CifDocument`; all other documents remain opaque and yield empty CIF-specific source data.
- [ ] Run `pytest tests/test_cristma_adapter.py tests/test_document.py -q`.
- [ ] Record the exact modified files because CRAFT has no git history.

## Task 3: Introduce one CRiStMa loader for every structure file

**Files:**

- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/вивер/src/crystal_viewer/core/structure_io.py`
- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/вивер/tests/test_structure_io.py`
- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/вивер/tests/test_progressive_load.py`

- [ ] Add parameterized failing tests for `.cif`, `.res`, `.ins`, `.vasp`, `POSCAR`, `CONTCAR`, `.pdb`, and `.xyz`. Spy on `cristma.read` and assert every case calls it as `cristma.read(path)` with no explicit format argument.
- [ ] Assert each returned display model retains the canonical CRiStMa object, source document, and diagnostics; Q peaks remain excluded from SHELX structures.
- [ ] Monkeypatch CRiStMa's generic `read` result with a canonical structure from a novel source name and prove the CRAFT loader needs no format-specific branch to accept it.
- [ ] Run `pytest tests/test_structure_io.py tests/test_progressive_load.py -q` and confirm the current suffix switch fails.
- [ ] Replace all structure-format branches with one boundary:

```python
def _load_cristma_structures(path: Path) -> list[CraftCrystalStructure]:
    result = cristma.read(path)
    return [
        from_cristma(
            structure,
            source_path=path,
            diagnostics=tuple(result.diagnostics),
            source_document=result.document,
        )
        for structure in result.structures
    ]
```

- [ ] Keep only the explicit XPFF project-container branch before the generic CRiStMa call. Audit the XPFF loader so any embedded CIF/other structural payload is delegated to `cristma.read_text` rather than parsed by the application.
- [ ] Keep the narrow small-CIF-occupancy recovery as a temporary adapter policy selected by the returned `SourceInfo.format == "cif"`, not by CRAFT's filename dispatch.
- [ ] Remove `_load_cristma_cif`, `_load_shelx`, `_load_vasp`, `_load_pdb`, `_load_xyz`, and their parser-specific helpers only after parity tests pass.
- [ ] Run `pytest tests/test_structure_io.py tests/test_progressive_load.py tests/test_cristma_adapter.py -q`.
- [ ] Search for remaining application parsers: `rg -n "def _load_(cif|shelx|vasp|pdb|xyz)|_shelx_sfac|FVAR|LATT|CRYST1" src/crystal_viewer/core src/crystal_viewer/adapters` and resolve every format-parsing hit.

## Task 4: Derive CRAFT file support and dialogs from CRiStMa

**Files:**

- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/вивер/src/crystal_viewer/core/structure_io.py`
- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/вивер/src/crystal_viewer/ui/main_window.py`
- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/вивер/tests/test_structure_io.py`
- Modify: `/Users/artem/Yandex.Disk.localized/Python/XRD/вивер/tests/test_main_window_collection.py`

- [ ] Add failing tests proving `is_supported_structure_path`, supported suffixes/basenames, and the open-dialog filter are derived from CRiStMa registry data plus the explicit XPFF project suffix.
- [ ] Monkeypatch the CRiStMa format-information query with a novel suffix and assert it appears without editing CRAFT constants.
- [ ] Run `pytest tests/test_structure_io.py tests/test_main_window_collection.py -q` and confirm the hard-coded constants fail.
- [ ] Replace duplicated format constants with immutable values produced by `cristma.structure_suffixes()` and `cristma.structure_basenames()`. Add XPFF only in the application-owned project-filter layer.
- [ ] Build the user-facing dialog filter from descriptor labels/suffixes while keeping a separate `XRD Finder projects (*.xpff)` entry.
- [ ] Run `pytest tests/test_structure_io.py tests/test_main_window_collection.py -q`.
- [ ] Record the exact modified files and focused test result.

## Task 5: Shared-environment and final application gate

**Files:**

- Modify only if a demonstrated packaging failure requires it: `/Users/artem/Yandex.Disk.localized/Python/XRD/вивер/run_viewer.command`
- Modify only if a demonstrated packaging failure requires it: `/Users/artem/Yandex.Disk.localized/Python/XRD/вивер/run_viewer.bat`

- [ ] Install the completed internal CRiStMa wheel into the same Sci environment used by CRAFT launchers; do not introduce a CRAFT-local environment or a Sci runtime API dependency.
- [ ] Smoke-load one real file for each registered CRAFT format through `cristma.read(path)` and then through CRAFT's generic loader.
- [ ] Run CRAFT's complete test suite once: `pytest -q`.
- [ ] Launch CRAFT through `run_viewer.command`; open the real RES fixture and at least one CIF, POSCAR/PDB, and XYZ fixture. Verify names, cells/display cells, atoms, symmetry, and progressive analysis.
- [ ] Use a test-only CRiStMa format-information/read provider for a novel suffix and verify the loader path recognizes it without a CRAFT code change; this is the acceptance proof for future formats.
- [ ] Perform the boundary review against both hard rules: applications own project formats, CRiStMa owns structural formats; no application-specific structural parser, mapper, writer, or registry exists outside CRiStMa.
- [ ] Produce a final changed-file list and test transcript because CRAFT has no git history.
