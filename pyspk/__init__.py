"""Public package interface for pySPk.

The primary entrypoint is `sup_model`, which returns the scale-dependent power spectrum suppression.
"""

from __future__ import annotations

from .model import get_limits, optimal_mass, sup_model

__all__ = [
    "get_limits",
    "optimal_mass",
    "sup_model",
]

__version__ = "1.8.0"
