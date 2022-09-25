import numpy as _np
from scipy.interpolate import Akima1DInterpolator as _Akima1DInterpolator
from scipy.interpolate import LinearNDInterpolator as _LinearNDInterpolator
from .fit_vals import best_fit_vals as _best_fit_vals
from .fit_vals import limits as _limits
import warnings as _warnings
import pkgutil as _pkgutil

_warnings.filterwarnings("always", category=UserWarning)

def _power_law(m_halo, fb_a, fb_pow, fb_pivot=1):
    """
    Simple power-law function fit to the baryon fraction - halo mass relation

    Parameters
    ----------
    m_halo: array of float
        halo mass in M_sun units
    fb_a : float
        power law constant
    fb_pow : float
        power law exponent
    fb_pivot : float, default 1
        power law pivot point. 

    Returns
    -------
    output: array of float
        baryon fraction normalised by the Universal baryon fraction: f_b / (Omega_b / Omega_m)
    """
    fb = fb_a * _np.power(m_halo / fb_pivot, fb_pow)
    return fb


def _poly_2(x, vals):
    """
    Simple polynomial of order 2. 

    Parameters
    ----------
    x : array of float
        independent variable 
    vals : array of float
        array containing the polynomial coefficients

    Returns
    -------
    y: array of float
        dependent variable 
    """
    y = vals[2] * _np.power(x, 2) + vals[1] * x + vals[0] 
    return y


def _optimal_mass_funct(k, params):
    """
    Optimal mass function defined in eq.(2) in Salcido et al. (2022).

    Parameters
    ----------
    k : array of float
        co-moving wavenumber in units [h/Mpc]
    prams : dict
        dictionary containing the best fitting parameters of the SP(k) model at a specific redshift. 

    Returns
    -------
    output: array of float
        log_10 of the optimal mass in M_sun units
    """
    output = params['alpha'] - (params['alpha'] - params['beta']) * _np.power(k, params['gamma'])
    return output


def optimal_mass(SO, z, k, verbose=False):
    """
    Optimal mass function as a function of scale and redshift for a specific 
    spherical over-density. Defined in eq.(2) in Salcido et al. (2022).

    Parameters
    ----------
    SO : int
        spherical over-density. Only accepts 200 or 500
    z : float
        redshift z. Only accepts values of z <= 3
    k : array of float
        co-moving wavenumber in units [h/Mpc]
    verbose : boolean, default True
        Run in verbose mode

    Returns
    -------
    output: array of float
        optimal mass in M_sun units
    """
    k_max = _np.max(k)
    if k_max > 12:
        raise Exception('\033[91mpy-spk was calibrated up to k_max = 12 [h/Mpc] '
                        'Please specify values of k <= 12 [h/Mpc] \033[0m')
    if k_max > 8:
        _warnings.warn('\033[33mScales with k_max > k_ny = 8 [h/Mpc] '
                       'may not be accurately reproduced by the model. \033[0m', stacklevel=2)

    params = _get_params(SO, z)
    output = params['alpha'] - (params['alpha'] - params['beta']) * _np.power(k, params['gamma'])
    return _np.power(10, output)


def _lambda_funct(x, params):
    """
    lambda function for a specific redshift. Defined in eq.(6) in Salcido et al. (2022).

    Parameters
    ----------
    k : array of float
        co-moving wavenumber in units [h/Mpc]
    prams : dict
        dictionary containing the best fitting parameters of the SP(k) model at a specific redshift. 

    Returns
    -------
    output: array of float
        lambda
    """
    output = 1 + params['lambda_a'] * _np.exp(params['lambda_b'] * x)
    return output


def _mu_funct(x, params):
    """
    mu function for a specific redshift. Defined in eq.(7) in Salcido et al. (2022).

    Parameters
    ----------
    k : array of float
        co-moving wavenumber in units [h/Mpc]
    prams : dict
        dictionary containing the best fitting parameters of the SP(k) model at a specific redshift. 

    Returns
    -------
    output: array of float
        mu
    """
    A = params['mu_a']
    B = 1 - params['mu_a']
    C = 1 + _np.exp(params['mu_b'] * x + params['mu_c'])
    output = A + (B / C)
    return output


def _nu_func(x, params):
    """
    nu function for a specific redshift. Defined in eq.(8) in Salcido et al. (2022).

    Parameters
    ----------
    k : array of float
        co-moving wavenumber in units [h/Mpc]
    prams : dict
        dictionary containing the best fitting parameters of the SP(k) model at a specific redshift. 

    Returns
    -------
    output: array of float
        nu
    """
    A = params['nu_a']
    B = _np.exp(-.5 * ((x - params['nu_b']) / params['nu_c']) ** 2)
    output = A * B
    return output


def get_limits(SO, z, m_halo):
    """
    Function that returns the baryon fraction fitting limits as a function of halo mass
    and redshift.

    Parameters
    ----------
    SO : int
        spherical over-density. Only accepts 200 or 500
    z : float
        redshift z. Only accepts values of z <= 3
    m_halo: array of float
        halo mass in M_sun units

    Returns
    -------
    min_fb: array of float
        lower fitting limit for the baryon fraction normalised by the universal baryon fraction
    max_fb: array of float
        upper fitting limit for the baryon fraction normalised by the universal baryon fraction
    """

    inter_min_x0 = _Akima1DInterpolator(_limits[str(SO)]['z'], _limits[str(SO)]['min_x0'])
    inter_min_x1 = _Akima1DInterpolator(_limits[str(SO)]['z'], _limits[str(SO)]['min_x1'])
    inter_min_x2 = _Akima1DInterpolator(_limits[str(SO)]['z'], _limits[str(SO)]['min_x2'])

    inter_max_x0 = _Akima1DInterpolator(_limits[str(SO)]['z'], _limits[str(SO)]['max_x0'])
    inter_max_x1 = _Akima1DInterpolator(_limits[str(SO)]['z'], _limits[str(SO)]['max_x1'])
    inter_max_x2 = _Akima1DInterpolator(_limits[str(SO)]['z'], _limits[str(SO)]['max_x2'])

    min_fb = 10 ** (inter_min_x0(z) + inter_min_x1(z) * _np.log10(m_halo) + inter_min_x2(z) * _np.log10(m_halo) ** 2)
    max_fb = 10 ** (inter_max_x0(z) + inter_max_x1(z) * _np.log10(m_halo) + inter_max_x2(z) * _np.log10(m_halo) ** 2)

    return min_fb, max_fb


def _get_params(SO, z):
    """
    Computes the best fit parameters for the SP(k) model at a specific redshift. 

    Parameters
    ----------
    SO : int
        spherical over-density. Only accepts 200 or 500
    z : float
        redshift z. Only accepts values of z <= 3 

    Returns
    -------
    params: dict
        best fit parameters for the specified SO and z.
    """
    params = {}
    try:
        for param_i in _best_fit_vals[str(SO)]:
            params[param_i] = _poly_2(1 + z, _best_fit_vals[str(SO)][param_i])
    except:
        raise Exception('''\033[91m
                        Spherical overdensity should be specified. 
                        Please use 200 or 500 
                        \033[0m''')
    return params


def sup_model(SO, z, fb_a=None, fb_pow=None, fb_pivot=1, M_halo=None, fb=None, extrapolate=False,
              k_min=0.1, k_max=8, n=100, errors=False, verbose=False):
    """
    Returns the suppression of the total matter power spectrum as a function of scale 'k' using the SP(k) model.
    Automatically selects the required optimal mass as a function of scale and redshift. Requires the baryon 
    fraction - halo mass relation, either by specifying power-law fitting parameters, or non-parametric 
    proving arrays with fb and M_halo at a specific redshift. If a non-parametric relation is given, an 
    interpolator is used to compute fb at the optimal mass. 

    Parameters
    ----------
    SO : int
        spherical over-density. Only accepts 200 or 500
    z : float
        redshift z. Only accepts values of z <= 3
    fb_a : float, optional
        power law constant
    fb_pow : float, option
        power law exponent
    fb_pivot : float, optional, default 1
        power law pivot point. 
    fb : array of float, optional
        array containing the (binned) baryon fraction normalised by the universal baryon fraction: 
        f_b / (Omega_b / Omega_m) for the fb - M-halo relation.
    M_halo : array of float, optional
        array containing the (binned) halo mass for the fb - M-halo relation in M_sun units.
        For interpolation, out-of-bounds points return NaNs.
    extrapolate: boolean, default False
        Whether to extrapolate to out-of-bounds points based on first and last intervals, or to return NaNs.
    k_min : float, default 0.1
        minimum co-moving wavenumber in units [h/Mpc]
    k_max : float, default 8, max 12
        maximum co-moving wavenumber in units [h/Mpc]. Default is set to the Nyquist frequency of the Antilles
        simulations (see Salcido et al. 2022). 
    n : int, default 100
        number of equally spaced co-moving wavenumber in log-spaced between k_min and k_max.
    errors : boolean, default False
        enables additional output with the bootstrapped 68% and 95% confidence intervals from statistical errors.
    verbose : boolean, default True
        run in verbose mode
        
    Returns
    -------
    k: array of float
        array with the co-moving wavenumber in units [h/Mpc]
    sup: array of float
        array with the suppression of the total matter power spectrum as a function of scale
    error_68_m : array of float, optional
        array with the -1 sigma confidence interval 
    error_68_p : array of float, optional
        array with the +1 sigma confidence interval 
    error_95_m : array of float, optional
        array with the -2 sigma confidence interval 
    error_95_p : array of float, optional
        array with the +2 sigma confidence interval 

    Notes
    ----------
    The maximum co-moving wavenumber in units [h/Mpc] is set to the Nyquist frequency of the Antilles simulations 
    (see Salcido et al. 2022). Although SP(k) was fitted up to k = 12 [h/Mpc], we caution the user that setting 
    𝑘 > 𝑘Ny (8 [h/Mpc]) might not be representative of the true uncertainties in the data. 
    """

    if z < 0:
        raise Exception('\033[91mIncorrect redshift.\033[0m') from None
        
    if z > 3:
        raise Exception('\033[91mpy-spk was calibrated up to z = 3.0. '
                        'Please specify z <= 3.0 \033[0m') from None

    if z < 0.125:
        _warnings.warn('\033[33mpy-spk was calibrated down to z = 0.125. Redshifts '
                       'z < 0.125 may not be accurately reproduced by the model. \033[0m', stacklevel=2)

    if k_max > 12:
        raise Exception('\033[91mpy-spk was calibrated up to k_max = 12 [h/Mpc] '
                        'Please specify k_max <= 12 [h/Mpc] \033[0m') from None
        
    if k_max > 8:
        _warnings.warn('\033[33mScales with k_max > k_ny = 8 [h/Mpc] '
                       'may not be accurately reproduced by the model. \033[0m', stacklevel=2)

    params = _get_params(SO, z)
    k = _np.round(_np.logspace(_np.log10(k_min), _np.log10(k_max), n), 6)
    logk = _np.log10(k)

    best_mass = _optimal_mass_funct(k, params)

    if fb_a or fb_pow:
        if verbose:
            print('\033[36mUsing power-law fit for fb - M_halo at z=%.3f, ' 
                  'normalised at M_halo = %.2e [M_sun] \033[0m' % (z, fb_pivot))
        try:
            f_b = _power_law(10 ** best_mass, fb_a, fb_pow, fb_pivot)
        except:
            raise Exception('\033[91mWhen using a power-law, both parameters should be given. '
                            'Please specify: fb_a and fb_pow \033[0m') from None
    else:
        if verbose:
            print('\033[36mUsing binned data for fb - M_halo at z=%.3f \033[0m' % z)
        try:
            fb_inter = _Akima1DInterpolator(_np.log10(M_halo), _np.log10(fb))
            if extrapolate:
                fb_inter.extrapolate = True
            f_b = 10 ** fb_inter(best_mass)
        except:
            raise Exception('\033[91mWhen using binned data, both halo mass and baryon '
                            'fraction should be given as monotonically increasing arrays. '
                            'Please specify: M_halo (array) and fb (array). \033[0m')

    min_fb, max_fb = get_limits(SO, z, 10 ** best_mass)

    out_min = f_b < min_fb
    out_max = f_b > max_fb

    mass_out_min = best_mass[out_min]

    if any(out_min):
        _warnings.warn('\033[91mFound baryon fraction values outside fitting limits. '
                       'fb < lower_limit between %.1e <= M_halo [M_sun] <= %.1e. ' 
                       'sup_model() will return NaNs within those limits. \033[0m' 
                       % (10 ** mass_out_min.min(), 10 ** mass_out_min.max()), stacklevel=2)

    mass_out_max = best_mass[out_max]

    if any(out_max):
        _warnings.warn('\033[91mFound baryon fraction values outside fitting limits. '
                       'fb > upper_limit between %.1e <= M_halo [M_sun] <= %.1e. ' 
                       'sup_model() will return NaNs within those limits. \033[0m' 
                       % (10 ** mass_out_max.min(), 10 ** mass_out_max.max()), stacklevel=2)

    mask = _np.logical_or(out_min,out_max)

    x0 = _lambda_funct(logk, params)
    x1 = _mu_funct(logk, params)
    x2 = _nu_func(logk, params)

    sup = x0 - (x0 - x1) * _np.exp(-x2 * f_b)

    sup[mask] = _np.NAN

    if errors:
        raw_table = _pkgutil.get_data(__name__, 'stat_errors_' + str (SO) + '.csv').decode('utf-8').splitlines()
        table = _np.loadtxt(raw_table, delimiter=",", skiprows=1)
        coords = table[:, [0, 1, 2]]
        interp_68_m = _LinearNDInterpolator(coords, table[:, 4], rescale=True)
        interp_68_p = _LinearNDInterpolator(coords, table[:, 5], rescale=True)
        interp_95_m = _LinearNDInterpolator(coords, table[:, 6], rescale=True)
        interp_95_p = _LinearNDInterpolator(coords, table[:, 7], rescale=True)

        z_array = _np.full_like(k, z)
        data = _np.column_stack([k, f_b, z_array])
        error_68_m = interp_68_m(data)
        error_68_p = interp_68_p(data)
        error_95_m = interp_95_m(data)
        error_95_p = interp_95_p(data)

        return k, sup, error_68_m, error_68_p, error_95_m, error_95_p

    else:
        return k, sup
