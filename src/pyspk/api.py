"""User-facing Pydantic models for structured API usage.

These are always available because Pydantic is a required dependency of `pyspk`.
"""

from __future__ import annotations

from .schema import (
    BinnedRelation,
    CosmoPowerLawRelation,
    DoublePowerLawRelation,
    PowerLawRelation,
    SupModelRequest,
)

__all__ = [
    "BinnedRelation",
    "CosmoPowerLawRelation",
    "DoublePowerLawRelation",
    "PowerLawRelation",
    "SupModelRequest",
]
