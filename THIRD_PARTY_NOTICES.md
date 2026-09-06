# Third-party data and notices

CrIStMa is distributed under the BSD-3-Clause license in [`LICENSE`](LICENSE).
That license covers original CrIStMa code, documentation, and project-authored
test data. It does not replace the licenses or provenance of bundled
third-party resources.

CrIStMa has no runtime dependency on the projects named below. Selected source
data are compiled into versioned local resources so calculations remain
reproducible. Upstream versions, commits, file hashes, scientific references,
and license texts are retained with the package.

## Crystallographic catalog

The packaged space-group and Wyckoff catalog is normalized from spglib 2.7.0:

- upstream: [spglib/spglib](https://github.com/spglib/spglib), tag `v2.7.0`;
- commit: `12355c77fb7c505a55f52cae36341d73b781a065`;
- source files: `database/spg.csv` and `database/Wyckoff.csv`;
- license: BSD-3-Clause;
- generated content: 530 Hall settings covering all 230 three-dimensional
  space-group types and 3,467 Wyckoff records.

No International Tables pages or Bilbao Crystallographic Server records are
copied into CrIStMa. Exact hashes and rebuild instructions are recorded in
[`src/cristma/reference_data/resources/crystallography/SOURCE.md`](src/cristma/reference_data/resources/crystallography/SOURCE.md).
The upstream license is reproduced in
[`SPGLIB_LICENSE.txt`](src/cristma/reference_data/resources/crystallography/SPGLIB_LICENSE.txt).

## Cordero covalent radii

`covalent_radii.json` is compiled from a pinned QCElemental artifact:

- upstream: [MolSSI/QCElemental](https://github.com/MolSSI/QCElemental);
- commit: `c4eb31cff9c7041f4767804a0076e35343df8177`;
- source path: `qcelemental/data/alvarez_2008_covalent_radii.py`;
- source license: BSD-3-Clause;
- source SHA-256:
  `9ac22bedfc04ead3567ebf0484fe09583e959d679aec99b69a0aef13388cb63e`.

Scientific reference: B. Cordero et al., “Covalent radii revisited”,
*Dalton Transactions* (2008), 2832–2838,
[doi:10.1039/B801115J](https://doi.org/10.1039/B801115J).

The complete provenance and selection policy are documented in
[`COVALENT_RADII_SOURCE.md`](src/cristma/reference_data/resources/COVALENT_RADII_SOURCE.md).
The upstream license is reproduced in
[`QCELEMENTAL_LICENSE.txt`](src/cristma/reference_data/resources/QCELEMENTAL_LICENSE.txt).

## Shannon ionic and crystal radii

`shannon_radii.json` is compiled from a pinned pymatgen artifact:

- upstream: [materialsproject/pymatgen](https://github.com/materialsproject/pymatgen);
- commit: `0428f232a569ffe6b16fa030d38ea35a56d70fd6`;
- source path: `dev_scripts/periodic_table_resources/Shannon_Radii.csv`;
- source license: MIT;
- source SHA-256:
  `d71d42ef465b7ab48bc9fec7e60c7b2fe500b5787b107d1c813410eb3581b52e`.

Scientific reference: R. D. Shannon, “Revised effective ionic radii and
systematic studies of interatomic distances in halides and chalcogenides”,
*Acta Crystallographica A* **32** (1976), 751–767,
[doi:10.1107/S0567739476001551](https://doi.org/10.1107/S0567739476001551).

The pinned pymatgen generator states that the earlier provenance of its CSV is
unknown. CrIStMa therefore attributes the numerical artifact to the exact
pymatgen commit and cites the Shannon paper as its scientific reference; it
does not claim a stronger undocumented lineage. See
[`SHANNON_SOURCE.md`](src/cristma/reference_data/resources/SHANNON_SOURCE.md)
and the reproduced [`PYMATGEN_LICENSE.txt`](src/cristma/reference_data/resources/PYMATGEN_LICENSE.txt).

## Neutral-atom X-ray form factors

`xray_f0.json` is normalized from xraylib 4.3.0 `data/FF.dat`:

- upstream: [tschoonj/xraylib](https://github.com/tschoonj/xraylib), tag
  `xraylib-4.3.0`;
- commit: `f94a3f5008dfd1c882b88ff26cd5052559423c83`;
- source SHA-256:
  `9aca1801042adee51aac62ab32c9d9445e37ce5c947a7e685b42311f520c530a`;
- normalized-data SHA-256:
  `02a187ebaf5a66d599d9cf1df781c01b42f1144e32870da349f3079c82ee1d32`;
- source license: xraylib BSD-style license.

The xraylib documentation attributes its elastic-scattering data to D. E.
Cullen, J. H. Hubbell, and J. H. Kissel, *EPDL97: The Evaluated Photon Data
Library, '97 Version*, LLNL report UCRL-50400, Vol. 6, Rev. 5. CrIStMa stores
the native tabulation and second derivatives in its own documented convention
`s = 1/(2d)` and does not copy or import xraylib code at runtime.

Complete build provenance is recorded in
[`SOURCE.md`](src/cristma/reference_data/resources/xray/SOURCE.md), and the
upstream license is reproduced in
[`XRAYLIB_LICENSE.txt`](src/cristma/reference_data/resources/xray/XRAYLIB_LICENSE.txt).

## Curated chemical reference knowledge

`chemical_reference_v3.json` and `chemical_reference_v3_1.json` are original
CrIStMa compilations of machine-readable classification rules. Their embedded
`sources` sections cite 22 books, papers, reviews, and nomenclature references
that support the scientific concepts. CrIStMa does not reproduce the source
publications or copy their prose or tables. The JSON records the project's
curated interpretation and links each rule to its supporting bibliography.

## Crystallographic Open Database fixtures

Selected CIF files used only as tests come from the Crystallography Open
Database, whose data are dedicated to the public domain under CC0. Each file
retains publication metadata, source identity, and a recorded SHA-256 digest.
The fixture inventory is in
[`tests/fixtures/crystal_chemistry/PROVENANCE.md`](tests/fixtures/crystal_chemistry/PROVENANCE.md).

## Project-authored format fixtures

Small VASP, XYZ/extXYZ, PDB, and analytic crystallographic fixtures are
hand-authored or hand-reduced test inputs. They exercise public file-format
grammars and are not copied calculation outputs, pseudopotentials, or upstream
scientific datasets. Their provenance is documented beside the corresponding
fixtures where applicable.

## Downstream redistribution

Commercial and closed-source software may use CrIStMa under BSD-3-Clause.
Distributors must retain the CrIStMa license and the applicable notices and
license texts for bundled third-party resources. Materials with incompatible,
copyleft, non-commercial, or unclear redistribution terms must not be added to
CrIStMa without a separate compatibility review.
