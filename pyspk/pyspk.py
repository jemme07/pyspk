import numpy as _np
from scipy.interpolate import PchipInterpolator as _PchipInterpolator
from .fit_vals import best_fit_vals as _best_fit_vals


def _power_law(m_halo, fb_a, fb_pow, fb_norm=1):
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
    fb_norm : float, default 1
        power law normalisation. 

    Returns
    -------
    output: array of float
        baryon fraction normalised by the Universal baryon fraction: f_b / (Omega_b / Omega_m)
    """
    fb = fb_a * _np.power(m_halo / fb_norm, fb_pow)
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
        co-moving wavenumber in units [k/Mpc]
    prams : dict
        dictionary containing the best fitting parameters of the SP(k) model at a specific redshift. 

    Returns
    -------
    output: array of float
        log_10 of the optimal mass in M_sun units
    """
    output = params['alpha'] - (params['alpha'] - params['beta']) * _np.power(k, params['gamma'])
    return output


def optimal_mass(SO, z, k):
    """
    Optimal mass function as a function of scale and redshift for a specific 
    spherical over-density. Defined in eq.(2) in Salcido et al. (2022).

    Parameters
    ----------
    SO : int
        spherical over-density. (Only accepts 200 or 500)
    z : float
        redshift z
    k : array of float
        co-moving wavenumber in units [k/Mpc]

    Returns
    -------
    output: array of float
        optimal mass in M_sun units
    """
    params = _get_params(SO, z)
    output = params['alpha'] - (params['alpha'] - params['beta']) * _np.power(k, params['gamma'])
    return _np.power(10, output)


def _lambda_funct(x, params):
    """
    lambda function for a specific redshift. Defined in eq.(6) in Salcido et al. (2022).

    Parameters
    ----------
    k : array of float
        co-moving wavenumber in units [k/Mpc]
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
        co-moving wavenumber in units [k/Mpc]
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
        co-moving wavenumber in units [k/Mpc]
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


def _get_params(SO, z):
    """
    Computes the best fit parameters for the SP(k) model at a specific redshift. 

    Parameters
    ----------
    SO : int
        spherical over-density. (Only accepts 200 or 500)
    z : float
        redshift z 

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


def sup_model(SO, z, fb_a=None, fb_pow=None, fb_norm=1, M_halo=None, fb=None, k_min=0.1, k_max=8, n=100):
    """
    Returns the suppression of the total matter power spectrum as a function of scale 'k' using the SP(k) model.
    Automatically selects the required optimal mass as a function of scale and redshift. Requires the baryon 
    fraction - halo mass relation, either by specifying power-law fitting parameters, or non-parametric 
    proving arrays with fb and M_halo at a specific redshift. If a non-parametric relation is given, an 
    interpolator is used to compute fb at the optimal mass. 

    Parameters
    ----------
    SO : int
        spherical over-density. (Only accepts 200 or 500)
    z : float
        redshift z
    fb_a : float, optional
        power law constant
    fb_pow : float, option
        power law exponent
    fb_norm : float, optional, default 1
        power law normalisation. 
    fb : array of float, optional
        array containing the (binned) baryon fraction normalised by the universal baryon fraction: 
        f_b / (Omega_b / Omega_m) for the fb - M-halo relation.
    M_halo : array of float, optional
        array containing the (binned) halo mass for the fb - M-halo relation in M_sun units
    k_min : float, default 0.1
        minimum co-moving wavenumber in units [k/Mpc]
    k_max : float, default 8, max 12
        maximum co-moving wavenumber in units [k/Mpc]. Default is set to the Nyquist frequency of the Antilles
        simulations (see Salcido et al. 2022). 
    n : int, default 100
        number of equally spaced co-moving wavenumber in log-spaced between k_min and k_max.
    Returns
    -------
    k: array of float
        array with the co-moving wavenumber in units [k/Mpc]
    sup: array of float
        array with the suppression of the total matter power spectrum as a function of scale

    Notes
    ----------
    The maximum co-moving wavenumber in units [k/Mpc] is set to the Nyquist frequency of the Antilles simulations 
    (see Salcido et al. 2022). Although SP(k) was fitted up to k = 12 [k/Mpc], we caution the user that setting 
    𝑘 > 𝑘Ny (8 [k/Mpc]) might not be representative of the true uncertainties in the data. 
    """

    params = _get_params(SO, z)
    k = _np.logspace(_np.log10(k_min), _np.log10(k_max), n)
    logk = _np.log10(k)

    best_mass = _optimal_mass_funct(k, params)

    if fb_a or fb_pow:
        print('''\033[33m 
              Using power-law fit for the fb - M_halo at z=%.2f
              \033[0m''' % z)
        try:
            f_b = _power_law(10 ** best_mass, fb_a, fb_pow, fb_norm)
        except:
            raise Exception('''\033[91m 
                            When using power-law, both parameters should be given. 
                            Please specify: fb_a and fb_pow 
                            \033[0m''')
    else:
        print('''\033[33m
              Using binned data for fb - M_halo at z=%.2f
              \033[0m''' % z)
        try:
            fb_inter = _PchipInterpolator(M_halo, fb)
            f_b = fb_inter(10 ** best_mass)
        except:
            raise Exception('''\033[91m 
                            When using binned data, both halo mass and baryon fraction arrays should be given. 
                            Please specify: M_halo(array) and fb(array) 
                            \033[0m''')

    x0 = _lambda_funct(logk, params)
    x1 = _mu_funct(logk, params)
    x2 = _nu_func(logk, params)

    sup = x0 - (x0 - x1) * _np.exp(-x2 * f_b)

    return k, sup
