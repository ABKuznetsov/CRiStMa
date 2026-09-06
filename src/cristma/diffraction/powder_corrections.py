"""Explicit Lorentz-polarization corrections for intrinsic powder lines."""

from __future__ import annotations

import math

from .diagnostics import DiffractionInvariantError, POWDER_CORRECTION_SINGULAR_ANGLE
from .powder_correction_models import (
    BraggBrentanoGeometry,
    CorrectedPowderLine,
    CorrectedPowderLineSet,
    PowderCorrectionProvenance,
)
from .powder_models import PowderLineSet, RadiationProbe


CALCULATOR_METHOD = "bragg_brentano_lorentz_polarization"
CALCULATOR_VERSION = "1"
LORENTZ_FORMULA = "1 / (sin(theta)^2 * cos(theta))"
POLARIZATION_FORMULA = "p_perpendicular + (1 - p_perpendicular) * cos(2theta)^2"


class PowderCorrectionCalculator:
    """Apply explicit geometry factors without profiles or sample corrections."""

    def calculate(
        self,
        powder_lines: PowderLineSet,
        geometry: BraggBrentanoGeometry,
    ) -> CorrectedPowderLineSet:
        if not isinstance(powder_lines, PowderLineSet):
            raise TypeError("powder_lines must be PowderLineSet")
        if not isinstance(geometry, BraggBrentanoGeometry):
            raise TypeError("geometry must be BraggBrentanoGeometry")
        if powder_lines.spectrum.probe is not RadiationProbe.XRAY:
            raise TypeError(
                "Bragg-Brentano X-ray polarization requires an X-ray spectrum"
            )

        corrected: list[CorrectedPowderLine] = []
        perpendicular = geometry.perpendicular_polarization_fraction
        for line in powder_lines.lines:
            if not 0.0 < line.two_theta_deg < 180.0:
                raise DiffractionInvariantError(
                    POWDER_CORRECTION_SINGULAR_ANGLE,
                    "Lorentz correction is singular at this diffraction angle",
                    {
                        "powder_line_id": line.line_id,
                        "two_theta_deg": line.two_theta_deg,
                        "geometry_id": geometry.geometry_id,
                    },
                )
            theta = math.radians(line.two_theta_deg / 2.0)
            lorentz = 1.0 / (math.sin(theta) ** 2 * math.cos(theta))
            polarization = perpendicular + (1.0 - perpendicular) * math.cos(
                2.0 * theta
            ) ** 2
            intensity = (
                line.intrinsic_line_intensity * lorentz * polarization
            )
            corrected.append(
                CorrectedPowderLine(
                    powder_line_id=line.line_id,
                    two_theta_deg=line.two_theta_deg,
                    intrinsic_line_intensity=line.intrinsic_line_intensity,
                    lorentz_factor=lorentz,
                    polarization_factor=polarization,
                    corrected_line_intensity=intensity,
                )
            )

        provenance = PowderCorrectionProvenance(
            method=CALCULATOR_METHOD,
            version=CALCULATOR_VERSION,
            geometry_id=geometry.geometry_id,
            radiation_probe=powder_lines.spectrum.probe,
            radiation_source_id=powder_lines.spectrum.source_id,
            lorentz_formula=LORENTZ_FORMULA,
            polarization_formula=POLARIZATION_FORMULA,
            lines_corrected=len(corrected),
        )
        return CorrectedPowderLineSet(
            lines=tuple(corrected),
            powder_lines=powder_lines,
            geometry=geometry,
            provenance=provenance,
        )


__all__ = ["PowderCorrectionCalculator"]
