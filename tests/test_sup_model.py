"""Regression tests for the public pyspk API."""

import numpy as np
import pytest

import pyspk
from pyspk.exceptions import InputValidationError


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
