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
from dataclasses import dataclass
from typing import Any, Callable, Literal, Optional, Protocol, overload

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
    "SupModelEvaluator",
    "build_sup_model_evaluator",
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


def _get_efunc(cosmo: Any) -> Callable[[float], Any]:
    """Get an `efunc(z)` callable from a cosmology-like object.

    Args:
        cosmo: Cosmology-like object.

    Returns:
        Callable returning the dimensionless Hubble parameter $E(z)$.

    Raises:
        InputValidationError: If the object does not provide a callable `efunc`.
    """
    efunc = getattr(cosmo, "efunc", None)
    if not callable(efunc):
        raise InputValidationError(
            "A cosmology object with a callable `efunc(z)` method is required for this relation."
        )
    return efunc


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


def _akino(rel: AkinoRelation, m_halo: _np.ndarray, z: float, cosmo: Any) -> _np.ndarray:
    """Evaluate the Akino relation."""
    A = _np.exp(rel.alpha) / 100
    B = _np.power(m_halo / 1e14, rel.beta - 1)
    efunc = _get_efunc(cosmo)
    C = _np.power(efunc(z) / efunc(0.3), rel.gamma)
    return A * B * C


def _double_power_law(
    rel: DoublePowerLawRelation, m_halo: _np.ndarray, z: float, cosmo: Any
) -> _np.ndarray:
    """Evaluate the double power-law relation."""
    efunc = _get_efunc(cosmo)
    A = 0.5 * rel.epsilon * _np.power(efunc(z) / efunc(0.3), rel.gamma)
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
    cosmo: Optional[Any] = ...,
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
    cosmo: Optional[Any] = ...,
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
    cosmo: Optional[Any] = None,
    k_array: Optional[ArrayLike] = None,
    k_min: float = 0.1,
    k_max: float = 8,
    n: int = 100,
    errors: bool = False,
    verbose: bool = False,
):
    r"""Compute power spectrum suppression using the SP(k) model.

    The model requires a baryon fraction - halo mass relation, provided either parametrically
    (power-law, Akino et al. 2022, double power-law) or as binned arrays.

    Relation specification
        Provide exactly one of the relation kinds below by passing its required parameters.
        (If you pass parameters from multiple kinds, the function will pick one using an
        internal precedence rule; to avoid surprises, provide only one kind.)

        - Power-law relation
          Required: `fb_a`, `fb_pow`.
          Optional: `fb_pivot` (pivot mass in M_sun; default is 1).
          e.g. `sup_model(SO=200, z=0.125, fb_a=0.4, fb_pow=0.3, fb_pivot=10**13.5)`

        - Binned relation
          Required: `M_halo`, `fb` (same length).
          Optional: `extrapolate` (if True, extrapolates in log10-space beyond the provided range).
          e.g. `sup_model(SO=200, z=0.5, M_halo=masses, fb=fractions, extrapolate=True)`

        - Akino (redshift-dependent) relation
          Required: `alpha`, `beta`, `gamma`, `cosmo`.
          `cosmo` must provide a callable `efunc(z)`.
          e.g. `sup_model(SO=500, z=0.7, alpha=4.16, beta=1.2, gamma=0.39, cosmo=cosmo)`

        - Double power-law relation
          Required: `epsilon`, `alpha`, `beta`, `gamma`, `m_pivot`, `cosmo`.
          `alpha` and `beta` are the low- and high-mass slopes; `m_pivot` is in M_sun.
          `cosmo` must provide a callable `efunc(z)`.
          e.g. `sup_model(SO=500, z=0.7, epsilon=0.3, alpha=1.1, beta=0.2, gamma=0.5, `
          `m_pivot=1e14, cosmo=cosmo)`

    Args:
        SO: Spherical over-density. Supported values: 200 or 500.
        z: Redshift (calibrated for z <= 3).
        fb_a: Power-law normalization (required for power-law relation).
        fb_pow: Power-law exponent (required for power-law relation).
        fb_pivot: Power-law pivot mass in M_sun units (power-law relation).
        M_halo: Binned halo mass array for a binned relation (M_sun).
        fb:
            Binned baryon fraction array normalized by the universal baryon fraction
            (binned relation).
        extrapolate:
            Whether to extrapolate binned relations beyond provided bounds (binned relation).
        epsilon: Double power-law normalization parameter (double power-law relation).
        alpha:
            Akino normalization parameter (Akino relation) OR low-mass slope (double power-law
            relation).
        beta: Akino slope parameter (Akino relation) OR high-mass slope (double power-law relation).
        gamma: Redshift dependence parameter (Akino and double power-law relations).
        m_pivot: Double power-law pivot mass in M_sun units (double power-law relation).
        cosmo: Astropy cosmology object (required for Akino and double power-law relations).
        k_array: Explicit k array in [h/Mpc]. If provided, `k_min`, `k_max`, and `n` are ignored.
        k_min: Minimum k in [h/Mpc] for generated grid.
        k_max: Maximum k in [h/Mpc] for generated grid.
        n: Number of log-spaced k samples.
        errors: Whether to return statistical confidence intervals.
        verbose: Whether to print mode information.

    Returns:
        When `errors=False`, only `k` and `sup` are returned.

        k:
            1D array of co-moving wavenumbers in units of [h/Mpc].
        sup:
            1D array of suppression values (the SP(k) prediction).
        err68_minus:
            1D array of the lower 68% statistical interval (only if `errors=True`).
        err68_plus:
            1D array of the upper 68% statistical interval (only if `errors=True`).
        err95_minus:
            1D array of the lower 95% statistical interval (only if `errors=True`).
        err95_plus:
            1D array of the upper 95% statistical interval (only if `errors=True`).

    Raises:
        InputValidationError: If inputs are invalid.
        CalibrationRangeError: If k or z exceed calibrated ranges.
        MissingResourceError: If packaged statistical error tables are missing.
    """
    relation: dict
    if (fb_a is not None) or (fb_pow is not None):
        relation = {
            "kind": "power_law",
            "fb_a": fb_a,
            "fb_pow": fb_pow,
            "fb_pivot": fb_pivot,
        }
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
        relation = {
            "kind": "akino",
            "alpha": alpha,
            "beta": beta,
            "gamma": gamma,
        }

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


def _build_limits_interpolators(
    SO: int,
) -> tuple[
    _Akima1DInterpolator,
    _Akima1DInterpolator,
    _Akima1DInterpolator,
    _Akima1DInterpolator,
    _Akima1DInterpolator,
    _Akima1DInterpolator,
]:
    """Build Akima interpolators for the fitting-limit polynomials.

    This mirrors the work done in `get_limits`, but is designed to be done once and reused.
    """
    limits_so = _limits[str(SO)]
    inter_min_x0 = _Akima1DInterpolator(limits_so["z"], limits_so["min_x0"])
    inter_min_x1 = _Akima1DInterpolator(limits_so["z"], limits_so["min_x1"])
    inter_min_x2 = _Akima1DInterpolator(limits_so["z"], limits_so["min_x2"])

    inter_max_x0 = _Akima1DInterpolator(limits_so["z"], limits_so["max_x0"])
    inter_max_x1 = _Akima1DInterpolator(limits_so["z"], limits_so["max_x1"])
    inter_max_x2 = _Akima1DInterpolator(limits_so["z"], limits_so["max_x2"])

    return (
        inter_min_x0,
        inter_min_x1,
        inter_min_x2,
        inter_max_x0,
        inter_max_x1,
        inter_max_x2,
    )


def _get_limits_fast(
    *,
    z: float,
    m_halo: _np.ndarray,
    inter_min_x0: _Akima1DInterpolator,
    inter_min_x1: _Akima1DInterpolator,
    inter_min_x2: _Akima1DInterpolator,
    inter_max_x0: _Akima1DInterpolator,
    inter_max_x1: _Akima1DInterpolator,
    inter_max_x2: _Akima1DInterpolator,
) -> tuple[_np.ndarray, _np.ndarray]:
    """Fast fitting limits computation with cached interpolators."""
    logm = _np.log10(m_halo)
    min_fb = 10 ** (inter_min_x0(z) + inter_min_x1(z) * logm + inter_min_x2(z) * logm**2)
    max_fb = 10 ** (inter_max_x0(z) + inter_max_x1(z) * logm + inter_max_x2(z) * logm**2)
    return min_fb * 0.8, max_fb * 1.2


def _resolve_efunc(
    *,
    cosmo: Any,
    efunc: Optional[Callable[[float], Any]],
) -> Callable[[float], Any]:
    """Resolve an `efunc(z)` callable from either a direct callable or a cosmology-like object."""
    if efunc is not None:
        if not callable(efunc):
            raise InputValidationError("`efunc` must be callable.")
        return efunc
    return _get_efunc(cosmo)


@dataclass(frozen=True)
class SupModelEvaluator:
    """Fast evaluator for repeated SP(k) calls (e.g., MCMC inner loops).

    This is an additive performance API:

    - It keeps `sup_model(...)` unchanged.
    - It avoids per-call Pydantic validation and avoids list<->array roundtrips.
    - It caches the k-grid and the fitting-limit interpolators.

    Notes:
        This evaluator only supports `errors=False`.
    """

    SO: int
    relation_kind: Literal["power_law", "binned", "akino", "double_power_law"]
    k: _np.ndarray
    logk: _np.ndarray
    inter_min_x0: _Akima1DInterpolator
    inter_min_x1: _Akima1DInterpolator
    inter_min_x2: _Akima1DInterpolator
    inter_max_x0: _Akima1DInterpolator
    inter_max_x1: _Akima1DInterpolator
    inter_max_x2: _Akima1DInterpolator

    def __call__(
        self,
        *,
        z: float,
        fb_a: Optional[float] = None,
        fb_pow: Optional[float] = None,
        fb_pivot: float = 1.0,
        M_halo: Optional[ArrayLike] = None,
        fb: Optional[ArrayLike] = None,
        extrapolate: bool = False,
        epsilon: Optional[float] = None,
        alpha: Optional[float] = None,
        beta: Optional[float] = None,
        gamma: Optional[float] = None,
        m_pivot: Optional[float] = None,
        cosmo: Optional[Any] = None,
        efunc: Optional[Callable[[float], Any]] = None,
        verbose: bool = False,
    ) -> tuple[_np.ndarray, _np.ndarray]:
        """Evaluate suppression for a given redshift and relation parameters.

        Args:
            z: Redshift.
            fb_a: Power-law normalization (power-law relation).
            fb_pow: Power-law exponent (power-law relation).
            fb_pivot: Power-law pivot mass in M_sun units (power-law relation).
            M_halo: Binned halo mass array (binned relation).
            fb: Binned baryon fraction array (binned relation).
            extrapolate: Extrapolate binned relations beyond bounds (binned relation).
            epsilon: Double power-law normalization parameter (double power-law relation).
            alpha: Akino normalization OR low-mass slope for double power-law.
            beta: Akino slope OR high-mass slope for double power-law.
            gamma: Redshift dependence parameter (Akino/double power-law).
            m_pivot: Double power-law pivot mass in M_sun units.
            cosmo: Cosmology-like object providing `efunc(z)` (Akino/double power-law).
            efunc: Optional direct callable for `E(z)`; if provided, `cosmo` is not used.
            verbose: Whether to emit informational warnings.

        Returns:
            Tuple `(k, sup)`.

        Raises:
            InputValidationError: If required parameters are missing for the selected relation kind.
        """
        if z < CALIBRATED_Z_MIN:
            _warnings.warn(
                (
                    f"pyspk was calibrated down to z = {CALIBRATED_Z_MIN}. "
                    f"Redshifts z < {CALIBRATED_Z_MIN} may not be accurately "
                    "reproduced by the model."
                ),
                stacklevel=2,
            )

        params = _get_params(self.SO, float(z))
        best_mass = _optimal_mass_funct(self.k, params)

        if self.relation_kind == "power_law":
            if fb_a is None or fb_pow is None:
                raise InputValidationError("Power-law relation requires `fb_a` and `fb_pow`.")
            if verbose:
                _warnings.warn(
                    (
                        "Using power-law fit for fb - M_halo at "
                        f"z={z:.3f}, normalised at M_halo = {fb_pivot:.2e} [M_sun]."
                    ),
                    stacklevel=2,
                )
            f_b = _power_law(10**best_mass, float(fb_a), float(fb_pow), float(fb_pivot))

        elif self.relation_kind == "binned":
            if M_halo is None or fb is None:
                raise InputValidationError("Binned relation requires both `M_halo` and `fb`.")
            if verbose:
                _warnings.warn(f"Using binned data for fb - M_halo at z={z:.3f}.", stacklevel=2)
            m_arr = _np.asarray(M_halo, dtype=float)
            fb_arr = _np.asarray(fb, dtype=float)
            fb_inter = _Akima1DInterpolator(_np.log10(m_arr), _np.log10(fb_arr))
            fb_inter.extrapolate = bool(extrapolate)
            f_b = 10 ** fb_inter(best_mass)

        elif self.relation_kind == "double_power_law":
            if epsilon is None:
                raise InputValidationError("Double power-law relation requires `epsilon`.")
            if alpha is None:
                raise InputValidationError("Double power-law relation requires `alpha`.")
            if beta is None:
                raise InputValidationError("Double power-law relation requires `beta`.")
            if gamma is None:
                raise InputValidationError("Double power-law relation requires `gamma`.")
            if m_pivot is None:
                raise InputValidationError("Double power-law relation requires `m_pivot`.")
            if cosmo is None and efunc is None:
                raise InputValidationError(
                    "Double power-law relation requires either `cosmo` (with `efunc(z)`) "
                    "or `efunc`."
                )
            if verbose:
                _warnings.warn(
                    f"Using double power law for fb - M_halo at z={z:.3f}.",
                    stacklevel=2,
                )

            epsilon_f = float(epsilon)
            alpha_f = float(alpha)
            beta_f = float(beta)
            gamma_f = float(gamma)
            m_pivot_f = float(m_pivot)

            efunc_callable = _resolve_efunc(cosmo=cosmo, efunc=efunc)
            A = 0.5 * epsilon_f * _np.power(efunc_callable(z) / efunc_callable(0.3), gamma_f)
            B = _np.power(10**best_mass / m_pivot_f, alpha_f)
            C = _np.power(10**best_mass / m_pivot_f, beta_f)
            f_b = A * (B + C)

        else:  # akino
            if alpha is None:
                raise InputValidationError("Akino relation requires `alpha`.")
            if beta is None:
                raise InputValidationError("Akino relation requires `beta`.")
            if gamma is None:
                raise InputValidationError("Akino relation requires `gamma`.")
            if cosmo is None and efunc is None:
                raise InputValidationError("Akino relation requires either `cosmo` or `efunc`.")
            if verbose:
                _warnings.warn(
                    f"Using an Akino et al. 2022 power-law fit for fb - M_halo at z={z:.3f}.",
                    stacklevel=2,
                )

            alpha_f = float(alpha)
            beta_f = float(beta)
            gamma_f = float(gamma)

            efunc_callable = _resolve_efunc(cosmo=cosmo, efunc=efunc)
            A = _np.exp(alpha_f) / 100
            B = _np.power(10**best_mass / 1e14, beta_f - 1)
            C = _np.power(efunc_callable(z) / efunc_callable(0.3), gamma_f)
            f_b = A * B * C

        min_fb, max_fb = _get_limits_fast(
            z=float(z),
            m_halo=10**best_mass,
            inter_min_x0=self.inter_min_x0,
            inter_min_x1=self.inter_min_x1,
            inter_min_x2=self.inter_min_x2,
            inter_max_x0=self.inter_max_x0,
            inter_max_x1=self.inter_max_x1,
            inter_max_x2=self.inter_max_x2,
        )
        out_min = f_b < min_fb
        out_max = f_b > max_fb
        mask = _np.logical_or(out_min, out_max)

        x0 = _lambda_funct(self.logk, params)
        x1 = _mu_funct(self.logk, params)
        x2 = _nu_func(self.logk, params)
        sup = x0 - (x0 - x1) * _np.exp(-x2 * f_b)
        sup = _np.asarray(sup)
        sup[mask] = _np.nan

        return self.k, sup


def build_sup_model_evaluator(
    *,
    SO: int,
    relation_kind: Literal["power_law", "binned", "akino", "double_power_law"],
    k_array: Optional[ArrayLike] = None,
    k_min: float = 0.1,
    k_max: float = 8,
    n: int = 100,
) -> SupModelEvaluator:
    """Build a fast SP(k) evaluator for repeated calls.

    Args:
        SO: Spherical over-density. Supported values: 200 or 500.
        relation_kind: Relation kind. This fixes which parameters are required at call time.
        k_array: Explicit k array in [h/Mpc]. If provided, `k_min`, `k_max`, and `n` are ignored.
        k_min: Minimum k in [h/Mpc] for generated grid.
        k_max: Maximum k in [h/Mpc] for generated grid.
        n: Number of log-spaced k samples.

    Returns:
        A callable `SupModelEvaluator` instance.
    """
    if SO not in (200, 500):
        raise InputValidationError("SO must be 200 or 500.")

    if k_array is not None:
        k = _np.asarray(k_array, dtype=float)
        k_max_val = float(k.max())
    else:
        k = _np.round(_np.logspace(_np.log10(k_min), _np.log10(k_max), int(n)), 6)
        k_max_val = float(k.max())

    if k_max_val > CALIBRATED_K_MAX:
        raise CalibrationRangeError(
            f"pyspk was calibrated up to k_max = {CALIBRATED_K_MAX} [h/Mpc]. "
            f"Please specify k <= {CALIBRATED_K_MAX} [h/Mpc]."
        )
    if k_max_val > K_NYQUIST:
        _warnings.warn(
            (
                f"Scales with k_max > k_ny = {K_NYQUIST} [h/Mpc] "
                "may not be accurately reproduced by the model."
            ),
            stacklevel=2,
        )

    logk = _np.log10(k)
    (
        inter_min_x0,
        inter_min_x1,
        inter_min_x2,
        inter_max_x0,
        inter_max_x1,
        inter_max_x2,
    ) = _build_limits_interpolators(SO)

    return SupModelEvaluator(
        SO=SO,
        relation_kind=relation_kind,
        k=k,
        logk=logk,
        inter_min_x0=inter_min_x0,
        inter_min_x1=inter_min_x1,
        inter_min_x2=inter_min_x2,
        inter_max_x0=inter_max_x0,
        inter_max_x1=inter_max_x1,
        inter_max_x2=inter_max_x2,
    )
