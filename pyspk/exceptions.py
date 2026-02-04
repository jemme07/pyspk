"""Exceptions for pySPk."""

from __future__ import annotations


class PySPkError(Exception):
    """Base exception for pySPk."""


class InputValidationError(PySPkError):
    """Raised when user inputs are invalid."""


class CalibrationRangeError(PySPkError):
    """Raised when inputs are outside the calibrated range."""


class MissingResourceError(PySPkError):
    """Raised when a required packaged data resource is missing."""


class MissingOptionalDependencyError(PySPkError, ImportError):
    """Raised when an optional dependency is required but not installed."""
