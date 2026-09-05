# Neutral-atom X-ray form factors

`xray_f0.json` is a deterministic normalization of `data/FF.dat` from
xraylib 4.3.0, tag `xraylib-4.3.0`, peeled commit
`f94a3f5008dfd1c882b88ff26cd5052559423c83`.

- Source SHA-256: `9aca1801042adee51aac62ab32c9d9445e37ce5c947a7e685b42311f520c530a`
- Normalized data SHA-256: `02a187ebaf5a66d599d9cf1df781c01b42f1144e32870da349f3079c82ee1d32`
- Supported elements: atomic numbers 1 through 98
- Source variable: `q = sin(theta) / wavelength`
- CrIStMa variable: `s = 1 / (2d)`, in inverse angstroms
- Interpolation: cubic spline using the source table's second derivatives
- Generator: `tools/build_xray_f0.py`, version 1

The corresponding xraylib documentation attributes elastic-scattering data
to D. E. Cullen, J. H. Hubbell, and J. H. Kissel, *EPDL97: The Evaluated
Photon Data Library, '97 Version*, Lawrence Livermore National Laboratory
Report UCRL-50400, Vol. 6, Rev. 5.

The installed CrIStMa package does not import or require xraylib. The upstream
BSD-style license is reproduced in `XRAYLIB_LICENSE.txt`.
