"""Pydantic schemas for validating the user-facing API.

These models are used internally to validate inputs to the public functions, and can
also be imported by users who want structured inputs.
"""

from __future__ import annotations

from typing import Any, Literal, Optional, Union

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator
from typing_extensions import TypeAlias

from .constants import CALIBRATED_K_MAX, CALIBRATED_Z_MAX, SUPPORTED_SOS


class _Base(BaseModel):
    model_config = ConfigDict(extra="forbid")


class PowerLawRelation(_Base):
    """Power-law baryon fraction relation."""

    kind: Literal["power_law"] = "power_law"
    fb_a: float
    fb_pow: float
    fb_pivot: float = 1.0

    @field_validator("fb_pivot")
    @classmethod
    def _pivot_positive(cls, v: float) -> float:
        if not np.isfinite(v) or v <= 0:
            raise ValueError("fb_pivot must be finite and > 0")
        return v


class BinnedRelation(_Base):
    """Binned baryon fraction relation."""

    kind: Literal["binned"] = "binned"
    M_halo: list[float]
    fb: list[float]
    extrapolate: bool = False

    @model_validator(mode="after")
    def _validate_binned(self) -> BinnedRelation:
        if len(self.M_halo) != len(self.fb):
            raise ValueError("M_halo and fb must have the same length")
        if len(self.M_halo) < 2:
            raise ValueError("M_halo and fb must have at least 2 points")

        m = np.asarray(self.M_halo, dtype=float)
        f = np.asarray(self.fb, dtype=float)
        if not np.all(np.isfinite(m)) or not np.all(np.isfinite(f)):
            raise ValueError("M_halo and fb must be finite")
        if np.any(m <= 0):
            raise ValueError("M_halo must be > 0")
        if np.any(f <= 0):
            raise ValueError("fb must be > 0")
        if np.any(np.diff(m) <= 0):
            raise ValueError("M_halo must be strictly increasing")
        return self


class CosmoPowerLawRelation(_Base):
    """Cosmology-based redshift-dependent power-law relation.

    Notes:
        This corresponds to the functional form motivated by Akino et al. (2022).
    """

    kind: Literal["cosmo_power_law"] = "cosmo_power_law"
    alpha: float
    beta: float
    gamma: float


class DoublePowerLawRelation(_Base):
    """Redshift-dependent double power-law relation."""

    kind: Literal["double_power_law"] = "double_power_law"
    epsilon: float
    alpha: float
    beta: float
    gamma: float
    m_pivot: float

    @field_validator("m_pivot")
    @classmethod
    def _pivot_positive(cls, v: float) -> float:
        if not np.isfinite(v) or v <= 0:
            raise ValueError("m_pivot must be finite and > 0")
        return v


Relation: TypeAlias = Union[
    PowerLawRelation,
    BinnedRelation,
    CosmoPowerLawRelation,
    DoublePowerLawRelation,
]


class SupModelRequest(_Base):
    """Validated request object for `pyspk.model.sup_model`.

    Notes:
        For `cosmo_power_law` and `double_power_law` relations, the caller must provide a compatible
        `astropy` cosmology when calling `sup_model`.
    """

    model_config = ConfigDict(extra="forbid", arbitrary_types_allowed=True)

    SO: int
    z: float
    relation: Relation

    cosmo: Optional[Any] = None

    k_array: Optional[list[float]] = None
    k_min: float = 0.1
    k_max: float = 8.0
    n: int = Field(default=100, ge=2)

    errors: bool = False
    verbose: bool = False

    @field_validator("SO")
    @classmethod
    def _so_supported(cls, v: int) -> int:
        if v not in SUPPORTED_SOS:
            raise ValueError(f"SO must be one of {SUPPORTED_SOS}")
        return v

    @field_validator("z")
    @classmethod
    def _z_range(cls, v: float) -> float:
        if not np.isfinite(v) or v < 0:
            raise ValueError("z must be finite and >= 0")
        if v > CALIBRATED_Z_MAX:
            raise ValueError(f"z must be <= {CALIBRATED_Z_MAX}")
        return float(v)

    @model_validator(mode="after")
    def _k_validation(self) -> SupModelRequest:
        if self.k_array is not None:
            k = np.asarray(self.k_array, dtype=float)
            if k.size < 2:
                raise ValueError("k_array must contain at least 2 values")
            if not np.all(np.isfinite(k)):
                raise ValueError("k_array must be finite")
            if np.any(k <= 0):
                raise ValueError("k_array must be > 0")
            if float(np.max(k)) > CALIBRATED_K_MAX:
                raise ValueError(f"k_array max must be <= {CALIBRATED_K_MAX}")
        else:
            if not np.isfinite(self.k_min) or not np.isfinite(self.k_max):
                raise ValueError("k_min and k_max must be finite")
            if self.k_min <= 0:
                raise ValueError("k_min must be > 0")
            if self.k_max <= 0:
                raise ValueError("k_max must be > 0")
            if self.k_min >= self.k_max:
                raise ValueError("k_min must be < k_max")
            if self.k_max > CALIBRATED_K_MAX:
                raise ValueError(f"k_max must be <= {CALIBRATED_K_MAX}")

        if self.relation.kind in {"cosmo_power_law", "double_power_law"} and self.cosmo is None:
            raise ValueError("cosmo is required for cosmo_power_law and double_power_law relations")

        return self
