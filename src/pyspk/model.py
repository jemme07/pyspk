"""SP(k) model implementation.

This module provides the public API: `sup_model`, `optimal_mass`, and `get_limits`.

All user-facing inputs are validated using Pydantic models in `pyspk.schema`.
"""

from __future__ import annotations

if __name__ == "__main__" and __package__ is None:  # pragma: no cover
    raise SystemExit(
        "This module is meant to be imported as part of the 'pyspk' package. "
        "Run it in module mode: 'python -m pyspk.model'."
    )

import warnings as _warnings
from collections.abc import Mapping
from typing import Any, Literal, Optional, Protocol, overload

import numpy as _np
from pydantic import ValidationError
from scipy.interpolate import Akima1DInterpolator as _Akima1DInterpolator
from scipy.interpolate import LinearNDInterpolator as _LinearNDInterpolator
from typing_extensions import TypeAlias

from .constants import (
    CALIBRATED_K_MAX,
    CALIBRATED_Z_MIN,
    K_NYQUIST,
)
from .exceptions import (
    CalibrationRangeError,
    InputValidationError,
    MissingResourceError,
    PySPkError,
)
from .fit_vals import best_fit_vals as _best_fit_vals
from .fit_vals import limits as _limits
from .resources import read_text as _read_text
from .schema import (
    AkinoRelation,
    BinnedRelation,
    DoublePowerLawRelation,
    PowerLawRelation,
    SupModelRequest,
)

__all__ = [
    "PySPkError",
    "get_limits",
    "optimal_mass",
    "sup_model",
]

ArrayLike: TypeAlias = Any


class CosmologyLike(Protocol):
    """Protocol for cosmology objects used by the Akino relations."""

    def efunc(self, z: float) -> float:  # pragma: no cover
        """Dimensionless Hubble parameter E(z)."""
        ...


def _power_law(m_halo: ArrayLike, fb_a: float, fb_pow: float, fb_pivot: float = 1.0) -> _np.ndarray:
    """Compute a simple power-law baryon fraction relation.

    Args:
        m_halo: Halo mass in M_sun units.
        fb_a: Power-law normalization.
        fb_pow: Power-law exponent.
        fb_pivot: Pivot mass in M_sun units.

    Returns:
        Baryon fraction normalized by the universal baryon fraction: fb / (Omega_b / Omega_m).
    """
    m = _np.asarray(m_halo, dtype=float)
    fb = float(fb_a) * _np.power(m / float(fb_pivot), float(fb_pow))
    return fb


def _poly_2(x: ArrayLike, vals: ArrayLike) -> _np.ndarray:
    """Evaluate a second-order polynomial.

    Args:
        x: Independent variable.
        vals: Polynomial coefficients as `(c0, c1, c2)`.

    Returns:
        Dependent variable $y = c_2 x^2 + c_1 x + c_0$.
    """
    x_arr = _np.asarray(x, dtype=float)
    coeffs = _np.asarray(vals, dtype=float)
    y = coeffs[2] * _np.power(x_arr, 2) + coeffs[1] * x_arr + coeffs[0]
    return y


def _optimal_mass_funct(k: ArrayLike, params: Mapping[str, float]) -> _np.ndarray:
    """Compute the optimal mass (log10) as in eq. (2) of Salcido et al. (2023).

    Args:
        k: Co-moving wavenumber in units of [h/Mpc].
        params: Best-fit SP(k) parameters at the desired redshift.

    Returns:
        log10 of the optimal mass in M_sun units.
    """
    k_arr = _np.asarray(k, dtype=float)
    output = params["alpha"] - (params["alpha"] - params["beta"]) * _np.power(
        k_arr, params["gamma"]
    )
    return output


def optimal_mass(SO: int, z: float, k: ArrayLike, verbose: bool = False) -> _np.ndarray:
    """Compute the optimal halo mass as a function of scale and redshift.

    Args:
        SO: Spherical over-density. Supported values: 200 or 500.
        z: Redshift. Calibrated for z <= 3.
        k: Co-moving wavenumber in units of [h/Mpc].
        verbose: Whether to run in verbose mode.

    Returns:
        Optimal mass in M_sun units.
    """
    # Validate via request model (relation is unused here).
    try:
        SupModelRequest.model_validate(
            {
                "SO": SO,
                "z": z,
                "relation": {
                    "kind": "power_law",
                    "fb_a": 1.0,
                    "fb_pow": 0.0,
                    "fb_pivot": 1.0,
                },
                "k_array": _np.asarray(k, dtype=float).tolist(),
            }
        )
    except ValidationError as exc:
        raise InputValidationError(str(exc)) from exc

    k_arr = _np.asarray(k, dtype=float)
    k_max = float(k_arr.max())
    if k_max > CALIBRATED_K_MAX:
        raise CalibrationRangeError(
            f"pyspk was calibrated up to k_max = {CALIBRATED_K_MAX} [h/Mpc]. "
            f"Please specify k <= {CALIBRATED_K_MAX} [h/Mpc]."
        )
    if k_max > K_NYQUIST:
        _warnings.warn(
            (
                f"Scales with k_max > k_ny = {K_NYQUIST} [h/Mpc] "
                "may not be accurately reproduced by the model."
            ),
            stacklevel=2,
        )

    params = _get_params(SO, z)
    output = params["alpha"] - (params["alpha"] - params["beta"]) * _np.power(
        k_arr, params["gamma"]
    )
    return _np.power(10, output)


def _lambda_funct(x: _np.ndarray, params: Mapping[str, float]) -> _np.ndarray:
    """Lambda function at a specific redshift (eq. 6 of Salcido et al. 2023).

    Args:
        x: log10(k) values.
        params: Best-fit parameter mapping for a given redshift.

    Returns:
        Lambda values.
    """
    return 1 + params["lambda_a"] * _np.exp(params["lambda_b"] * x)


def _mu_funct(x: _np.ndarray, params: Mapping[str, float]) -> _np.ndarray:
    """Mu function at a specific redshift (eq. 7 of Salcido et al. 2023).

    Args:
        x: log10(k) values.
        params: Best-fit parameter mapping for a given redshift.

    Returns:
        Mu values.
    """
    A = params["mu_a"]
    B = 1 - params["mu_a"]
    C = 1 + _np.exp(params["mu_b"] * x + params["mu_c"])
    return A + (B / C)


def _nu_func(x: _np.ndarray, params: Mapping[str, float]) -> _np.ndarray:
    """Nu function at a specific redshift (eq. 8 of Salcido et al. 2023).

    Args:
        x: log10(k) values.
        params: Best-fit parameter mapping for a given redshift.

    Returns:
        Nu values.
    """
    A = params["nu_a"]
    B = _np.exp(-0.5 * ((x - params["nu_b"]) / params["nu_c"]) ** 2)
    return A * B


def get_limits(SO: int, z: float, m_halo: ArrayLike) -> tuple[_np.ndarray, _np.ndarray]:
    """Return baryon-fraction fitting limits as a function of mass and redshift.

    Args:
        SO: Spherical over-density. Supported values: 200 or 500.
        z: Redshift.
        m_halo: Halo mass in M_sun units.

    Returns:
        Tuple of `(min_fb, max_fb)` where each is the fitting limit for the baryon fraction
        normalized by the universal baryon fraction.
    """
    # Validate SO/z via request model.
    try:
        SupModelRequest.model_validate(
            {
                "SO": SO,
                "z": z,
                "relation": {
                    "kind": "power_law",
                    "fb_a": 1.0,
                    "fb_pow": 0.0,
                    "fb_pivot": 1.0,
                },
            }
        )
    except ValidationError as exc:
        raise InputValidationError(str(exc)) from exc

    m_arr = _np.asarray(m_halo, dtype=float)

    inter_min_x0 = _Akima1DInterpolator(_limits[str(SO)]["z"], _limits[str(SO)]["min_x0"])
    inter_min_x1 = _Akima1DInterpolator(_limits[str(SO)]["z"], _limits[str(SO)]["min_x1"])
    inter_min_x2 = _Akima1DInterpolator(_limits[str(SO)]["z"], _limits[str(SO)]["min_x2"])

    inter_max_x0 = _Akima1DInterpolator(_limits[str(SO)]["z"], _limits[str(SO)]["max_x0"])
    inter_max_x1 = _Akima1DInterpolator(_limits[str(SO)]["z"], _limits[str(SO)]["max_x1"])
    inter_max_x2 = _Akima1DInterpolator(_limits[str(SO)]["z"], _limits[str(SO)]["max_x2"])

    logm = _np.log10(m_arr)
    min_fb = 10 ** (inter_min_x0(z) + inter_min_x1(z) * logm + inter_min_x2(z) * logm**2)
    max_fb = 10 ** (inter_max_x0(z) + inter_max_x1(z) * logm + inter_max_x2(z) * logm**2)

    return min_fb * 0.8, max_fb * 1.2


def _get_params(SO: int, z: float) -> dict[str, float]:
    """Compute best-fit SP(k) parameters at a specific redshift.

    Args:
        SO: Spherical over-density. Supported values: 200 or 500.
        z: Redshift.

    Returns:
        Dictionary of best-fit parameters for the specified `SO` and `z`.
    """
    params: dict[str, float] = {}
    try:
        for param_i in _best_fit_vals[str(SO)]:
            params[param_i] = float(_poly_2(1 + z, _best_fit_vals[str(SO)][param_i]))
    except Exception as exc:
        raise PySPkError("Failed to compute best-fit parameters for the given SO/z.") from exc
    return params


def _akino(rel: AkinoRelation, m_halo: _np.ndarray, z: float, cosmo: CosmologyLike) -> _np.ndarray:
    """Evaluate the Akino relation."""
    A = _np.exp(rel.alpha) / 100
    B = _np.power(m_halo / 1e14, rel.beta - 1)
    C = _np.power(cosmo.efunc(z) / cosmo.efunc(0.3), rel.gamma)
    return A * B * C


def _double_power_law(
    rel: DoublePowerLawRelation, m_halo: _np.ndarray, z: float, cosmo: CosmologyLike
) -> _np.ndarray:
    """Evaluate the double power-law relation."""
    A = 0.5 * rel.epsilon * _np.power(cosmo.efunc(z) / cosmo.efunc(0.3), rel.gamma)
    B = _np.power(m_halo / rel.m_pivot, rel.alpha)
    C = _np.power(m_halo / rel.m_pivot, rel.beta)
    return A * (B + C)


@overload
def sup_model(
    SO: int,
    z: float,
    fb_a: Optional[float] = ...,
    fb_pow: Optional[float] = ...,
    fb_pivot: float = ...,
    M_halo: Optional[ArrayLike] = ...,
    fb: Optional[ArrayLike] = ...,
    extrapolate: bool = ...,
    epsilon: Optional[float] = ...,
    alpha: Optional[float] = ...,
    beta: Optional[float] = ...,
    gamma: Optional[float] = ...,
    m_pivot: Optional[float] = ...,
    cosmo: Optional[CosmologyLike] = ...,
    k_array: Optional[ArrayLike] = ...,
    k_min: float = ...,
    k_max: float = ...,
    n: int = ...,
    *,
    errors: Literal[False] = ...,
    verbose: bool = ...,
) -> tuple[_np.ndarray, _np.ndarray]: ...


@overload
def sup_model(
    SO: int,
    z: float,
    fb_a: Optional[float] = ...,
    fb_pow: Optional[float] = ...,
    fb_pivot: float = ...,
    M_halo: Optional[ArrayLike] = ...,
    fb: Optional[ArrayLike] = ...,
    extrapolate: bool = ...,
    epsilon: Optional[float] = ...,
    alpha: Optional[float] = ...,
    beta: Optional[float] = ...,
    gamma: Optional[float] = ...,
    m_pivot: Optional[float] = ...,
    cosmo: Optional[CosmologyLike] = ...,
    k_array: Optional[ArrayLike] = ...,
    k_min: float = ...,
    k_max: float = ...,
    n: int = ...,
    *,
    errors: Literal[True],
    verbose: bool = ...,
) -> tuple[
    _np.ndarray,
    _np.ndarray,
    _np.ndarray,
    _np.ndarray,
    _np.ndarray,
    _np.ndarray,
]: ...


def sup_model(
    SO: int,
    z: float,
    fb_a: Optional[float] = None,
    fb_pow: Optional[float] = None,
    fb_pivot: float = 1,
    M_halo: Optional[ArrayLike] = None,
    fb: Optional[ArrayLike] = None,
    extrapolate: bool = False,
    epsilon: Optional[float] = None,
    alpha: Optional[float] = None,
    beta: Optional[float] = None,
    gamma: Optional[float] = None,
    m_pivot: Optional[float] = None,
    cosmo: Optional[CosmologyLike] = None,
    k_array: Optional[ArrayLike] = None,
    k_min: float = 0.1,
    k_max: float = 8,
    n: int = 100,
    errors: bool = False,
    verbose: bool = False,
):
    """Compute power spectrum suppression using the SP(k) model.

    The model requires a baryon fraction - halo mass relation, provided either parametrically
    (power-law, Akino et al. 2022, double power-law) or as binned arrays.

    Args:
        SO: Spherical over-density. Supported values: 200 or 500.
        z: Redshift (calibrated for z <= 3).
        fb_a: Power-law normalization.
        fb_pow: Power-law exponent.
        fb_pivot: Power-law pivot mass in M_sun units.
        M_halo: Binned halo mass array for a non-parametric relation (M_sun).
        fb: Binned baryon fraction array normalized by the universal baryon fraction.
        extrapolate: Whether to extrapolate binned relations beyond provided bounds.
        epsilon: Double power-law normalization parameter.
        alpha: Akino normalization parameter OR low-mass slope for double power-law.
        beta: Akino slope parameter OR high-mass slope for double power-law.
        gamma: Redshift dependence parameter.
        m_pivot: Double power-law pivot mass in M_sun units.
        cosmo: Astropy cosmology object (required for Akino and double power-law modes).
        k_array: Explicit k array in [h/Mpc]. If provided, `k_min`, `k_max`, and `n` are ignored.
        k_min: Minimum k in [h/Mpc] for generated grid.
        k_max: Maximum k in [h/Mpc] for generated grid.
        n: Number of log-spaced k samples.
        errors: Whether to return statistical confidence intervals.
        verbose: Whether to print mode information.

    Returns:
        If `errors` is False: `(k, sup)`.
        If `errors` is True: `(k, sup, err68_minus, err68_plus, err95_minus, err95_plus)`.

    Raises:
        InputValidationError: If inputs are invalid.
        CalibrationRangeError: If k or z exceed calibrated ranges.
        MissingResourceError: If packaged statistical error tables are missing.
    """
    relation: dict
    if (fb_a is not None) or (fb_pow is not None):
        relation = {"kind": "power_law", "fb_a": fb_a, "fb_pow": fb_pow, "fb_pivot": fb_pivot}
    elif (M_halo is not None) or (fb is not None):
        relation = {
            "kind": "binned",
            "M_halo": _np.asarray(M_halo, dtype=float).tolist() if M_halo is not None else None,
            "fb": _np.asarray(fb, dtype=float).tolist() if fb is not None else None,
            "extrapolate": extrapolate,
        }
    elif (epsilon is not None) or (m_pivot is not None):
        relation = {
            "kind": "double_power_law",
            "epsilon": epsilon,
            "alpha": alpha,
            "beta": beta,
            "gamma": gamma,
            "m_pivot": m_pivot,
        }
    else:
        relation = {"kind": "akino", "alpha": alpha, "beta": beta, "gamma": gamma}

    req_payload = {
        "SO": SO,
        "z": z,
        "relation": relation,
        "cosmo": cosmo,
        "k_array": _np.asarray(k_array, dtype=float).tolist() if k_array is not None else None,
        "k_min": k_min,
        "k_max": k_max,
        "n": n,
        "errors": errors,
        "verbose": verbose,
    }

    try:
        req = SupModelRequest.model_validate(req_payload)
    except ValidationError as exc:
        raise InputValidationError(str(exc)) from exc

    if req.z < CALIBRATED_Z_MIN:
        _warnings.warn(
            (
                f"pyspk was calibrated down to z = {CALIBRATED_Z_MIN}. "
                f"Redshifts z < {CALIBRATED_Z_MIN} may not be accurately reproduced by the model."
            ),
            stacklevel=2,
        )

    if req.k_array is not None:
        k = _np.asarray(req.k_array, dtype=float)
        k_max_val = float(k.max())
    else:
        k = _np.round(_np.logspace(_np.log10(req.k_min), _np.log10(req.k_max), req.n), 6)
        k_max_val = float(k.max())

    logk = _np.log10(k)

    if k_max_val > K_NYQUIST:
        _warnings.warn(
            (
                f"Scales with k_max > k_ny = {K_NYQUIST} [h/Mpc] "
                "may not be accurately reproduced by the model."
            ),
            stacklevel=2,
        )

    params = _get_params(req.SO, req.z)
    best_mass = _optimal_mass_funct(k, params)

    if isinstance(req.relation, PowerLawRelation):
        if verbose:
            _warnings.warn(
                (
                    "Using power-law fit for fb - M_halo at "
                    f"z={req.z:.3f}, normalised at M_halo = {req.relation.fb_pivot:.2e} [M_sun]."
                ),
                stacklevel=2,
            )
        f_b = _power_law(
            10**best_mass,
            req.relation.fb_a,
            req.relation.fb_pow,
            req.relation.fb_pivot,
        )

    elif isinstance(req.relation, BinnedRelation):
        if verbose:
            _warnings.warn(f"Using binned data for fb - M_halo at z={req.z:.3f}.", stacklevel=2)
        m_arr = _np.asarray(req.relation.M_halo, dtype=float)
        fb_arr = _np.asarray(req.relation.fb, dtype=float)
        fb_inter = _Akima1DInterpolator(_np.log10(m_arr), _np.log10(fb_arr))
        fb_inter.extrapolate = bool(req.relation.extrapolate)
        f_b = 10 ** fb_inter(best_mass)

    elif isinstance(req.relation, DoublePowerLawRelation):
        if verbose:
            _warnings.warn(
                f"Using double power law for fb - M_halo at z={req.z:.3f}.", stacklevel=2
            )
        if req.cosmo is None:
            raise InputValidationError(
                "A cosmology object with an `efunc(z)` method is required for the "
                "double power-law relation."
            )
        f_b = _double_power_law(req.relation, 10**best_mass, float(req.z), req.cosmo)

    else:
        if verbose:
            _warnings.warn(
                f"Using an Akino et al. 2022 power-law fit for fb - M_halo at z={req.z:.3f}.",
                stacklevel=2,
            )
        if req.cosmo is None:
            raise InputValidationError(
                "A cosmology object with an `efunc(z)` method is required for the Akino relation."
            )
        f_b = _akino(req.relation, 10**best_mass, float(req.z), req.cosmo)

    min_fb, max_fb = get_limits(req.SO, req.z, 10**best_mass)
    out_min = f_b < min_fb
    out_max = f_b > max_fb

    if _np.any(out_min):
        mass_out_min = best_mass[out_min]
        _warnings.warn(
            (
                "Found baryon fraction values outside fitting limits. "
                f"fb < lower_limit between {10 ** mass_out_min.min():.1e} <= M_halo [M_sun] <= "
                f"{10 ** mass_out_min.max():.1e}. sup_model() will return NaNs within those limits."
            ),
            stacklevel=2,
        )

    if _np.any(out_max):
        mass_out_max = best_mass[out_max]
        _warnings.warn(
            (
                "Found baryon fraction values outside fitting limits. "
                f"fb > upper_limit between {10 ** mass_out_max.min():.1e} <= M_halo [M_sun] <= "
                f"{10 ** mass_out_max.max():.1e}. sup_model() will return NaNs within those limits."
            ),
            stacklevel=2,
        )

    mask = _np.logical_or(out_min, out_max)

    x0 = _lambda_funct(logk, params)
    x1 = _mu_funct(logk, params)
    x2 = _nu_func(logk, params)
    sup = x0 - (x0 - x1) * _np.exp(-x2 * f_b)
    sup[mask] = _np.nan

    if not req.errors:
        return k, sup

    resource_name = f"stat_errors_{req.SO}.csv"
    try:
        raw_text = _read_text("pyspk", resource_name, encoding="utf-8")
    except FileNotFoundError as exc:
        raise MissingResourceError(f"Missing packaged resource: {resource_name}") from exc

    table = _np.loadtxt(raw_text.splitlines(), delimiter=",", skiprows=1)
    coords = table[:, [0, 1, 2]]
    interp_68_m = _LinearNDInterpolator(coords, table[:, 4], rescale=True)
    interp_68_p = _LinearNDInterpolator(coords, table[:, 5], rescale=True)
    interp_95_m = _LinearNDInterpolator(coords, table[:, 6], rescale=True)
    interp_95_p = _LinearNDInterpolator(coords, table[:, 7], rescale=True)

    z_array = _np.full_like(k, req.z)
    data = _np.column_stack([k, f_b, z_array])
    error_68_m = interp_68_m(data)
    error_68_p = interp_68_p(data)
    error_95_m = interp_95_m(data)
    error_95_p = interp_95_p(data)

    return k, sup, error_68_m, error_68_p, error_95_m, error_95_p
