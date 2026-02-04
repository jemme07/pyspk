"""Optional Pydantic models for the user-facing API.

This module is optional and is only functional when Pydantic is installed.

Notes:
    - `import pyspk.api` is always safe.
    - Accessing any of the model classes requires Pydantic.
    - Install with `pip install pyspk[api]`.
"""

from __future__ import annotations

from typing import Any, Literal

from .exceptions import MissingOptionalDependencyError

_PYDANTIC_AVAILABLE = False

try:  # pragma: no cover
    from pydantic import BaseModel  # type: ignore

    _PYDANTIC_AVAILABLE = True
except Exception:  # pragma: no cover
    BaseModel = object  # type: ignore[misc,assignment]


if _PYDANTIC_AVAILABLE:  # pragma: no cover

    class PowerLawRelation(BaseModel):
        """Power-law baryon fraction relation."""

        kind: Literal["power_law"] = "power_law"
        fb_a: float
        fb_pow: float
        fb_pivot: float = 1.0

    class BinnedRelation(BaseModel):
        """Binned baryon fraction relation."""

        kind: Literal["binned"] = "binned"
        M_halo: list[float]
        fb: list[float]
        extrapolate: bool = False

    class AkinoRelation(BaseModel):
        """Akino et al. (2022) redshift-dependent power-law relation."""

        kind: Literal["akino"] = "akino"
        alpha: float
        beta: float
        gamma: float

    class DoublePowerLawRelation(BaseModel):
        """Redshift-dependent double power-law relation."""

        kind: Literal["double_power_law"] = "double_power_law"
        epsilon: float
        alpha: float
        beta: float
        gamma: float
        m_pivot: float

    class SupModelRequest(BaseModel):
        """Validated request object for `pyspk.model.sup_model`.

        Notes:
            For `akino` and `double_power_law` relations, the caller must provide a compatible
            `astropy` cosmology when calling `sup_model`.
        """

        SO: int
        z: float
        relation: PowerLawRelation | BinnedRelation | AkinoRelation | DoublePowerLawRelation
        k_array: list[float] | None = None
        k_min: float = 0.1
        k_max: float = 8.0
        n: int = 100
        errors: bool = False
        verbose: bool = False


def __getattr__(name: str) -> Any:  # pragma: no cover
    if name in {
        "PowerLawRelation",
        "BinnedRelation",
        "AkinoRelation",
        "DoublePowerLawRelation",
        "SupModelRequest",
    }:
        raise MissingOptionalDependencyError(
            "Pydantic is required for pyspk.api models. Install via `pip install pyspk[api]`."
        )
    raise AttributeError(name)
