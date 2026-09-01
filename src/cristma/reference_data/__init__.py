"""Versioned scientific reference data used by CRiStMa tools."""

from .elements import ElementCatalog, ElementCategory, ElementRecord
from .chemical_reference import (
    ChemicalReference,
    ChemicalReferenceIntegrityReport,
    load_chemical_reference,
    validate_reference_integrity,
)
from .facade import ReferenceData
from .radii import CovalentRadii, CovalentRadiusRecord

__all__ = [
    "ChemicalReference",
    "ChemicalReferenceIntegrityReport",
    "CovalentRadii",
    "CovalentRadiusRecord",
    "ElementCatalog",
    "ElementCategory",
    "ElementRecord",
    "ReferenceData",
    "load_chemical_reference",
    "validate_reference_integrity",
]
