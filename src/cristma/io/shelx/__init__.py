"""Native SHELX source documents and structure I/O with lazy public exports."""

from __future__ import annotations

from importlib import import_module


_EXPORTS = {
    "ShelxAtomRecord": ("document", "ShelxAtomRecord"),
    "ShelxBlankRecord": ("document", "ShelxBlankRecord"),
    "ShelxCellInstruction": ("records", "ShelxCellInstruction"),
    "ShelxCommentRecord": ("document", "ShelxCommentRecord"),
    "ShelxDocument": ("document", "ShelxDocument"),
    "ShelxEndInstruction": ("records", "ShelxEndInstruction"),
    "ShelxFvarInstruction": ("records", "ShelxFvarInstruction"),
    "ShelxHklfInstruction": ("records", "ShelxHklfInstruction"),
    "ShelxInstructionRecord": ("document", "ShelxInstructionRecord"),
    "ShelxLattInstruction": ("records", "ShelxLattInstruction"),
    "ShelxOccupancyExpression": ("occupancy", "ShelxOccupancyExpression"),
    "ShelxPartInstruction": ("records", "ShelxPartInstruction"),
    "ShelxPhysicalLine": ("document", "ShelxPhysicalLine"),
    "ShelxQPeakRecord": ("document", "ShelxQPeakRecord"),
    "ShelxRecord": ("document", "ShelxRecord"),
    "ShelxResiInstruction": ("records", "ShelxResiInstruction"),
    "ShelxScatteringEntry": ("sfac", "ShelxScatteringEntry"),
    "ShelxSfacInstruction": ("records", "ShelxSfacInstruction"),
    "ShelxSourceEdit": ("document", "ShelxSourceEdit"),
    "ShelxSymmInstruction": ("records", "ShelxSymmInstruction"),
    "ShelxUnknownRecord": ("document", "ShelxUnknownRecord"),
    "ShelxZerrInstruction": ("records", "ShelxZerrInstruction"),
    "build_shelx_operations": ("symmetry", "build_shelx_operations"),
    "extract_sfac_entries": ("sfac", "extract_sfac_entries"),
    "format_shelx_symmetry": ("symmetry", "format_shelx_symmetry"),
    "map_shelx_structures": ("mapper", "map_shelx_structures"),
    "parse_shelx": ("parser", "parse_shelx"),
    "parse_shelx_symmetry": ("symmetry", "parse_shelx_symmetry"),
    "probe_shelx": ("probe", "probe_shelx"),
    "write_shelx_document": ("writer", "write_shelx_document"),
    "ShelxWriteOptions": ("writer", "ShelxWriteOptions"),
    "write_crystal_shelx": ("writer", "write_crystal_shelx"),
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
