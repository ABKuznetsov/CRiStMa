# Inorganic crystal-chemistry acceptance fixtures

Retrieved/compiled on 2026-09-02. The COD fixtures are redistributed under
the Crystallography Open Database CC0/public-domain policy and retain their
original publication metadata inside each CIF. Stable download URLs follow
`https://www.crystallography.net/cod/{COD_ID}.cif`.

| File | Formula | Space group | Source | SHA-256 |
|---|---|---|---|---|
| `NaF_9007457.cif` | NaF | Fm-3m | COD 9007457 | `677881bbde65d3cf34361a7e5ceda7388f7e4e76ee1e7e5c1d65435d12fb1f8c` |
| `SiC_9008856.cif` | SiC | F-43m | COD 9008856 | `d6d737b7905398a944b54abbf464b41db6045a83c402cb73d72e73adce9a3d6f` |
| `Si3N4_9013139.cif` | Si3N4 | P31c | COD 9013139 | `ac73a14e751b0f3227f23a988572f3326e9df1d5a3a91acb8b2f141aec816ec1` |
| `FeS2_9000594.cif` | FeS2 | Pa-3 | COD 9000594 | `0caef0969b6bcc5697310ab6f4316cbd01daf5fbcb7cfb1b5d6036647548a9d8` |
| `Na3P_1010294.cif` | Na3P | P63/mmc | COD 1010294 | `3af85993da9ab52d6a50f8071937719d14f1ccd9c0ec54c2d0bc14c442878d04` |
| `Bi2Te3_9011962.cif` | Bi2Te3 | R-3m:H | COD 9011962 | `442636ecaeea4c4bceb21c622b2e8e3369b860a943aad4aca5332fb9c9e97e27` |
| `CaMoO4_9009632.cif` | CaMoO4 | I41/a:2 | COD 9009632 | `ad4223a1d2812cf8e8301fda2f76172b4002f48f93653f98e89cf1d505771039` |
| `LiB3O5_3000122.cif` | LiB3O5 | Pna21 | COD/AMCSD-derived entry 3000122, previously retained by CRAFT | `6876e2e2436d610ae23ca2bfc62da30f97c62a060d3bdb39a4aa07b225b47476` |
| `anorthite_9000361.cif` | CaAl2Si2O8 | P-1 | COD 9000361 | `f4cc64d582ff84523232a99acdee8fe6eda446d7e7b15704139cb671ab9743e6` |

## Analytic fixtures

`CaN2_analytic.cif` was reconstructed from the reported I4/mmm model with
Ca on 2a and N on 4e. The published cell is `a = 3.5747 Å`, `c = 5.9844 Å`
and the published N-N distance is `1.202 Å`; the 4e coordinate was derived as
`z = 0.5 - 1.202/(2c) = 0.399572`. Scientific source: S. B. Schneider,
R. Frankovsky and W. Schnick, Inorganic Chemistry **51** (2012), 2366–2373,
DOI `10.1021/ic2023677`. SHA-256:
`00078f31c1d699e166d99afa92efe95489aa62457136aa9e8cc78a09c2bb8547`.

`NiAl_B2_analytic.cif` is the ideal B2/CsCl model in Pm-3m with Ni on 1a,
Al on 1b and the commonly reported room-temperature lattice parameter
`a = 2.887 Å`. It is generated test data rather than a copied database file.
SHA-256:
`e6ee133c7b72bb5c0be0a0cadbc844be5c7df5e8dc888d89e55fc66cf6f47c75`.

The analytic files are released as part of CrIStMa under the repository
license. They are acceptance fixtures, not independent structure
determinations.
