"""Resource-bounded calculation of sampled powder profiles."""

from __future__ import annotations

import math

import numpy as np

from .powder_correction_models import CorrectedPowderLineSet
from .powder_models import PowderLineSet
from .profile_models import (
    CalculatedProfile,
    ConstantWidthProfile,
    PowderProfileProvenance,
    ProfileIntensityBasis,
    TchProfile,
    UniformTwoThetaGrid,
)


_LOCAL_WINDOW_FWHM = 25.0
_GAUSSIAN_NORMALIZATION = math.sqrt(math.pi / (4.0 * math.log(2.0)))


class PowderProfileCalculator:
    """Apply an instrument-only broadening function on an explicit grid."""

    def calculate(
        self,
        lines: PowderLineSet | CorrectedPowderLineSet,
        grid: UniformTwoThetaGrid,
        broadening: ConstantWidthProfile | TchProfile,
        *,
        zero_shift_deg: float = 0.0,
    ) -> CalculatedProfile:
        if not isinstance(lines, (PowderLineSet, CorrectedPowderLineSet)):
            raise TypeError("lines must be PowderLineSet or CorrectedPowderLineSet")
        if not isinstance(grid, UniformTwoThetaGrid):
            raise TypeError("grid must be UniformTwoThetaGrid")
        if not isinstance(broadening, (ConstantWidthProfile, TchProfile)):
            raise TypeError("broadening must be ConstantWidthProfile or TchProfile")
        if isinstance(zero_shift_deg, bool) or not isinstance(
            zero_shift_deg, (int, float)
        ):
            raise TypeError("zero_shift_deg must be a real number")
        shift = float(zero_shift_deg)
        if not math.isfinite(shift):
            raise ValueError("zero_shift_deg must be finite")

        angles = np.asarray(grid.values, dtype=float)
        profile = np.zeros(angles.size, dtype=float)
        if isinstance(lines, CorrectedPowderLineSet):
            source_lines = lines.lines
            basis = ProfileIntensityBasis.CORRECTED
            intensities = tuple(item.corrected_line_intensity for item in source_lines)
        else:
            source_lines = lines.lines
            basis = ProfileIntensityBasis.INTRINSIC
            intensities = tuple(item.intrinsic_line_intensity for item in source_lines)

        contributed = 0
        for line, line_intensity in zip(source_lines, intensities, strict=True):
            if line_intensity == 0.0:
                continue
            original_center = line.two_theta_deg
            center = original_center + shift
            if isinstance(broadening, ConstantWidthProfile):
                fwhm = broadening.fwhm_deg
                mixing = 0.0
            else:
                fwhm, mixing = broadening._width_and_mixing_at(original_center)
            radius = _LOCAL_WINDOW_FWHM * fwhm
            raw_left = math.ceil(
                (center - radius - grid.start_deg) / grid.step_deg
            )
            raw_right = (
                math.floor((center + radius - grid.start_deg) / grid.step_deg) + 1
            )
            left = max(0, raw_left)
            right = min(angles.size, raw_right)
            if left >= right:
                continue
            offsets = angles[left:right] - center
            peak = _pseudo_voigt_density(offsets, fwhm, mixing)
            full_offsets = (
                grid.start_deg
                + np.arange(raw_left, raw_right, dtype=float) * grid.step_deg
                - center
            )
            sampled_window_area = float(
                np.sum(_pseudo_voigt_density(full_offsets, fwhm, mixing))
                * grid.step_deg
            )
            if sampled_window_area <= 0.0 or not math.isfinite(sampled_window_area):
                raise RuntimeError("profile kernel has no finite sampled area")
            profile[left:right] += peak * (line_intensity / sampled_window_area)
            contributed += 1

        model_id = (
            "constant_fwhm_gaussian"
            if isinstance(broadening, ConstantWidthProfile)
            else "tch_pseudo_voigt"
        )
        provenance = PowderProfileProvenance(
            method="cristma.diffraction.PowderProfileCalculator",
            version="1",
            broadening_model=model_id,
            broadening_scope="instrument_only",
            instrument_profile=broadening,
            intensity_basis=basis,
            zero_shift_deg=shift,
            local_window_fwhm=_LOCAL_WINDOW_FWHM,
            lines_considered=len(source_lines),
            lines_contributed=contributed,
        )
        return CalculatedProfile(
            two_theta_deg=tuple(float(value) for value in angles),
            intensity=tuple(float(value) for value in profile),
            intensity_basis=basis,
            source_lines=lines,
            provenance=provenance,
        )


def _pseudo_voigt_density(
    offsets: np.ndarray,
    fwhm: float,
    mixing: float,
) -> np.ndarray:
    gaussian = np.exp(-4.0 * math.log(2.0) * (offsets / fwhm) ** 2)
    gaussian /= fwhm * _GAUSSIAN_NORMALIZATION
    half_width = fwhm / 2.0
    lorentzian = half_width / (math.pi * (offsets**2 + half_width**2))
    return (1.0 - mixing) * gaussian + mixing * lorentzian


__all__ = ["PowderProfileCalculator"]
