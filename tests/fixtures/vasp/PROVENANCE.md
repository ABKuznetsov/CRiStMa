# VASP reference fixture provenance

Copy/reduction date: 2026-08-31.

These are deliberately small, hand-reduced structural reference fixtures, not
outputs claimed to come from a completed physical calculation. Their grammar,
field names, units, and ordering follow the official VASP Wiki examples. They
contain only the records needed to verify CRiStMa's structure-I/O boundary.

| File | Official format reference | SHA-256 | Expected content |
|---|---|---|---|
| `POSCAR` | https://vasp.at/wiki/POSCAR | `9e5bf9dd9ce954c642a88d63a8f30a5f41390390d018e751ef7cbb99a35b3324` | 1 Si atom, 1 structure |
| `XDATCAR` | https://vasp.at/wiki/XDATCAR | `8b57ce7f8e02cbf87b5994e695115a73841af8771eecc20804666264a4bb6315` | 1 Si atom, 3 complete configurations |
| `OUTCAR` | https://vasp.at/wiki/OUTCAR | `5d2fd33266320bf48861fd408e0d2e96aa7e96b855e159e708fab72985310f2b` | 1 Si atom, 2 complete ionic POSITION/TOTAL-FORCE blocks |
| `vasprun.xml` | https://vasp.at/wiki/Vasprun.xml | `602508ef2cdc198272ae5818d43b76540d37e7575bb80e0cfdeac9c0bf3cb2a5` | 1 Si atom, 2 complete calculation blocks |

No fixture is intentionally truncated. Truncation and malformed-record
recovery are covered by reduced inline parser tests, where the expected
diagnostic is visible beside the source mutation.
