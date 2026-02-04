"""Package-wide constants.

This module centralizes calibration limits and other shared numeric constants.
"""

from __future__ import annotations

# Calibration limits from Salcido et al. (2023)
CALIBRATED_Z_MAX: float = 3.0
CALIBRATED_Z_MIN: float = 0.125
CALIBRATED_K_MAX: float = 12.0
K_NYQUIST: float = 8.0

SUPPORTED_SOS: tuple[int, int] = (200, 500)
