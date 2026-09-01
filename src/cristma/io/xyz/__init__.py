"""Native XYZ API with implementation modules imported only on demand."""

from __future__ import annotations

from importlib import import_module


_EXPORTS = {
    "XyzDocument": (".document", "XyzDocument"),
    "XyzFrame": (".document", "XyzFrame"),
    "XyzFrameSpan": (".document", "XyzFrameSpan"),
    "XyzPropertySpec": (".document", "XyzPropertySpec"),
    "XyzMetadata": (".metadata", "XyzMetadata"),
    "index_xyz": (".index", "index_xyz"),
    "load_xyz_frame": (".parser", "load_xyz_frame"),
    "map_xyz_frame": (".mapper", "map_xyz_frame"),
    "parse_property_schema": (".metadata", "parse_property_schema"),
    "parse_xyz": (".parser", "parse_xyz"),
    "parse_xyz_metadata": (".metadata", "parse_xyz_metadata"),
    "probe_xyz": (".probe", "probe_xyz"),
    "validate_xyz_frame": (".parser", "validate_xyz_frame"),
    "XyzFormatHandler": (".handler", "XyzFormatHandler"),
}


def __getattr__(name: str):
    try:
        module_name, attribute = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(name) from exc
    value = getattr(import_module(module_name, __name__), attribute)
    globals()[name] = value
    return value


__all__ = sorted(_EXPORTS)
