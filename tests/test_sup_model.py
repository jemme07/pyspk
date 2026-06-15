"""Regression tests for the public pyspk API."""

import numpy as np
import pytest

import pyspk
from pyspk.exceptions import InputValidationError


class _DummyCosmo:
    def __init__(self, omega_m: float):
        self.omega_m = float(omega_m)

    def efunc(self, z: float) -> float:
        # Simple flat LCDM-like E(z) with Omega_L = 1 - Omega_m.
        zf = float(z)
        return float(np.sqrt(self.omega_m * (1 + zf) ** 3 + (1 - self.omega_m)))


def test_sup_model_power_law_shapes() -> None:
    """Returns finite suppression array for a basic power-law relation."""
    k_in = np.array([0.1, 0.5, 1.0, 5.0], dtype=float)
    k, sup = pyspk.sup_model(SO=200, z=0.5, fb_a=0.75, fb_pow=0.1, fb_pivot=1e14, k_array=k_in)

    assert np.allclose(k, k_in)
    assert sup.shape == k_in.shape
    assert np.isfinite(sup).all()


def test_sup_model_errors_returns_intervals() -> None:
    """When errors=True, returns confidence interval arrays."""
    k_in = np.array([0.1, 0.5, 1.0, 5.0], dtype=float)
    out = pyspk.sup_model(
        SO=200,
        z=0.5,
        fb_a=0.75,
        fb_pow=0.1,
        fb_pivot=1e14,
        k_array=k_in,
        errors=True,
    )

    assert len(out) == 6
    k, sup, err68_m, err68_p, err95_m, err95_p = out

    assert np.allclose(k, k_in)
    assert sup.shape == k_in.shape
    assert err68_m.shape == k_in.shape
    assert err68_p.shape == k_in.shape
    assert err95_m.shape == k_in.shape
    assert err95_p.shape == k_in.shape


def test_invalid_so_raises() -> None:
    """Invalid SO values are rejected by input validation."""
    with pytest.raises(InputValidationError):
        pyspk.sup_model(SO=123, z=0.5, fb_a=0.75, fb_pow=0.1)


def test_evaluator_matches_sup_model_power_law() -> None:
    """Evaluator matches sup_model for a power-law relation."""
    evaluator = pyspk.build_sup_model_evaluator(SO=200, relation_kind="power_law", k_max=2.0, n=32)

    z = 0.3
    fb_a = 0.5
    fb_pow = 0.2
    fb_pivot = 10**13.5

    k1, sup1 = pyspk.sup_model(
        SO=200,
        z=z,
        fb_a=fb_a,
        fb_pow=fb_pow,
        fb_pivot=fb_pivot,
        k_max=2.0,
        n=32,
        errors=False,
    )
    k2, sup2 = evaluator(z=z, fb_a=fb_a, fb_pow=fb_pow, fb_pivot=fb_pivot)

    assert np.allclose(k1, k2)
    assert np.allclose(sup1, sup2, equal_nan=True)


def test_evaluator_matches_sup_model_cosmo_power_law() -> None:
    """Evaluator matches sup_model for the cosmology-based power-law relation."""
    cosmo = _DummyCosmo(omega_m=0.3)
    evaluator = pyspk.build_sup_model_evaluator(
        SO=500,
        relation_kind="cosmo_power_law",
        k_max=2.0,
        n=32,
    )

    z = 0.7
    alpha = 4.16
    beta = 1.2
    gamma = 0.39

    k1, sup1 = pyspk.sup_model(
        SO=500,
        z=z,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        cosmo=cosmo,
        k_max=2.0,
        n=32,
        errors=False,
    )
    k2, sup2 = evaluator(z=z, alpha=alpha, beta=beta, gamma=gamma, cosmo=cosmo)

    assert np.allclose(k1, k2)
    assert np.allclose(sup1, sup2, equal_nan=True)


def test_evaluator_raises_for_z_above_calibrated_max() -> None:
    """Fast evaluator defaults to strict z validation (z <= 3)."""
    evaluator = pyspk.build_sup_model_evaluator(SO=200, relation_kind="power_law", k_max=2.0, n=32)

    with pytest.raises(InputValidationError, match=r"z must be <= 3\.0"):
        evaluator(z=3.1, fb_a=0.5, fb_pow=0.2, fb_pivot=10**13.5)


def test_evaluator_nan_policy_for_z_above_calibrated_max() -> None:
    """Optional NaN policy returns finite k and NaN suppression for z > 3."""
    evaluator = pyspk.build_sup_model_evaluator(
        SO=200,
        relation_kind="power_law",
        k_max=2.0,
        n=32,
        z_out_of_range="nan",
    )

    k, sup = evaluator(z=3.1, fb_a=0.5, fb_pow=0.2, fb_pivot=10**13.5)

    assert np.all(np.isfinite(k))
    assert np.isnan(sup).all()
