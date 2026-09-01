# CRiStMa crystallography reference data

Dataset ID: `cristma.crystallography.spglib`  
Schema version: `1.0.0`  
Compiled: `2026-09-01`

The normalized JSON is generated from spglib 2.7.0 database files under the
BSD-3-Clause license. No International Tables pages or Bilbao database records
are copied into this package.

## Upstream

- Repository: `https://github.com/spglib/spglib`
- Tag: `v2.7.0`
- Commit: `12355c77fb7c505a55f52cae36341d73b781a065`
- License: `SPGLIB_LICENSE.txt`
- `database/spg.csv`:
  `4457df1042b14a65ea62af0bad7b5b609a4fc33592df245802cfe005b221f95e`
- `database/Wyckoff.csv`:
  `d3d786a1f0187e5c6d69a3ade35648ffab34fd1b977d61ad84d8b0434b8b7ca0`

## Generated resources

- `space_groups.json`:
  `d5a63e2af0143e7b5004fcc5fdd5ef12616467f1bedfa251f681e9b31c3c000d`
- `wyckoff_positions.json`:
  `750b892ae0e22c82b634643aa53967115d34e1c6549c986aca7f9a65676f3744`

The output contains 530 Hall settings representing all 230 three-dimensional
space-group types and 3467 Wyckoff records.

## Rebuild

```bash
python tools/compile_spglib_crystallography.py \
  --spg /private/tmp/cristma-spg-v270.csv \
  --wyckoff /private/tmp/cristma-Wyckoff-v270.csv \
  --output src/cristma/reference_data/resources/crystallography \
  --upstream-commit 12355c77fb7c505a55f52cae36341d73b781a065 \
  --compiled-date 2026-09-01
```

The compiler requires the optional development dependency `spglib==2.7.0`.
The installed CRiStMa runtime does not require spglib.
