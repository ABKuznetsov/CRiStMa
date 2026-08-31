"""Native VASP structure I/O with lazy public exports."""

from __future__ import annotations

from importlib import import_module


_EXPORTS = {
    "CoordinateMode": ("document", "CoordinateMode"),
    "OutcarDocument": ("document", "OutcarDocument"),
    "PoscarDocument": ("document", "PoscarDocument"),
    "VaspAtomRow": ("document", "VaspAtomRow"),
    "VaspFormatHandler": ("handler", "VaspFormatHandler"),
    "VaspFrameSpan": ("document", "VaspFrameSpan"),
    "VaspHeader": ("document", "VaspHeader"),
    "VaspScale": ("document", "VaspScale"),
    "VaspSnapshot": ("document", "VaspSnapshot"),
    "VasprunDocument": ("document", "VasprunDocument"),
    "XdatcarDocument": ("document", "XdatcarDocument"),
    "fractional_from_cartesian": ("numeric", "fractional_from_cartesian"),
    "map_vasp_snapshot": ("mapper", "map_vasp_snapshot"),
    "parse_poscar": ("poscar", "parse_poscar"),
    "poscar_snapshot": ("poscar", "poscar_snapshot"),
    "probe_vasp": ("probe", "probe_vasp"),
    "scaled_cartesian": ("numeric", "scaled_cartesian"),
    "scaled_lattice": ("numeric", "scaled_lattice"),
    "load_xdatcar_snapshot": ("xdatcar", "load_xdatcar_snapshot"),
    "parse_xdatcar": ("xdatcar", "parse_xdatcar"),
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    try:
        module_name, attribute_name = _EXPORTS[name]
    except KeyError as error:
        raise AttributeError(name) from error
    value = getattr(import_module(f"{__name__}.{module_name}"), attribute_name)
    globals()[name] = value
    return value
