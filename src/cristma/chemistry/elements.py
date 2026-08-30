"""Dependency-free chemical element symbol normalization."""

from __future__ import annotations

import re


ELEMENT_SYMBOLS = frozenset(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe "
    "Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In "
    "Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf "
    "Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm "
    "Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split()
)


def normalize_element(value: str) -> str:
    """Extract and validate an IUPAC symbol from a CIF type or site label."""

    match = re.match(r"[A-Za-z]+", value.strip())
    if match is None:
        raise ValueError(f"Cannot determine chemical element from {value!r}")
    letters = match.group(0)
    for width in (2, 1):
        if len(letters) < width:
            continue
        candidate = letters[:width].capitalize()
        if candidate in ELEMENT_SYMBOLS:
            return candidate
    raise ValueError(f"Unknown chemical element in {value!r}")
