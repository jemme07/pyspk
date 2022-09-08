import numpy as _np
from scipy.interpolate import PchipInterpolator as _PchipInterpolator


_fit_vals = {}
_fit_vals['200'] = {'alpha': [15.274930840056541,-1.245047539257722,0.14540677973676638],
                   'beta': [14.930101598751966,-1.062673017386771,0.12120874422958225],
                   'gamma': [0.973317059802909,-0.3341132933728971,0.147971236001855],
                   'lambda_a': [0.03257049282030497,-0.012994181725415443,0.0016012935859615331],
                   'lambda_b': [2.936213543501669,0.4647835283419599,0.015469588056462733],
                   'mu_a': [0.7185313695191637,-0.17419453040304703,0.04344834913957274],
                   'mu_b': [7.061473347061119,-2.6451585009317475,0.693874112309198],
                   'mu_c': [9.72782645402674,-6.692331015886623,1.3184537648431263],
                   'nu_a': [659.6062166482205,190.65546760275708,25.12315520958482],
                   'nu_b': [-10.627984535501893,0.38934321021631624,0.2461189904914769],
                   'nu_c': [3.153654283048472,-0.13183283346451197,-0.07021924767085688]}

_fit_vals['500'] = {'alpha': [14.76469953943553,-0.9861664024147414,0.11918239669227224],
                   'beta': [14.61129247278941,-0.9086389425569597,0.10810207298747473],
                   'gamma': [0.8419196246836991,0.06967440823217347,0.014730497492187217],
                   'lambda_a': [0.019909493214134613,-0.006723408341682784,0.0007220653696607669],
                   'lambda_b': [3.04279503307814,0.5256217945904988,0.00660016011437811],
                   'mu_a': [0.6875723754738109,-0.15501931108472428,0.04106325613876646],
                   'mu_b': [4.372736593576084,0.4034162144260609,-0.02740227359038959],
                   'mu_c': [5.514096648948588,-2.8918746793434114,0.5961162742741664],
                   'nu_a': [585.1128565518618,867.0756296967672,38.90230731032294],
                   'nu_b': [-11.243766644267883,-0.09618413993416547,0.3767515346494691],
                   'nu_c': [3.3546510607155344,-0.12716883213945795,-0.08194814244409948]}


def _power_law(x, vals):
    output = vals[0] * _np.power(x, vals[1])
    return output


def _poly_2(x, vals):
    output = vals[2] * _np.power(x, 2) + vals[1] * x + vals[0] 
    return output


def _optimal_mass_funct(k, params):
    output = params['alpha'] - (params['alpha'] - params['beta']) * _np.power(k, params['gamma'])
    return output


def optimal_mass(SO, z, k):
    params = _get_params(SO, z)
    output = params['alpha'] - (params['alpha'] - params['beta']) * _np.power(k, params['gamma'])
    return _np.power(10, output)


def _lambda_funct(x, params):
    output = 1 + params['lambda_a'] * _np.exp(params['lambda_b'] * x)
    return output


def _mu_funct(x, params):
    A = params['mu_a']
    B = 1 - params['mu_a']
    C = 1 + _np.exp(params['mu_b'] * x + params['mu_c'])
    output = A + (B / C)
    return output


def _nu_func(x, params):
    A = params['nu_a']
    B = _np.exp(-.5 * ((x - params['nu_b']) / params['nu_c']) ** 2)
    output = A * B
    return output


def _get_params(SO, z):
    params = {}
    try:
        for param_i in _fit_vals[str(SO)]:
            params[param_i] = _poly_2(1 + z, _fit_vals[str(SO)][param_i])
    except:
        raise Exception('''\033[91m
                        Spherical overdensity should be specified. 
                        Please use 200 or 500 
                        \033[0m''')
    return params


def sup_model(SO, z, fb_a=None, fb_pow=None, M_halo=None, fb=None, k_min=0.1, k_max=8, n=100):

    params = _get_params(SO, z)
    k = _np.logspace(_np.log10(k_min), _np.log10(k_max), n)
    logk = _np.log10(k)

    best_mass = _optimal_mass_funct(k, params)

    if fb_a or fb_pow:
        print('''\033[33m 
              Using power-law fit for the fb - M_halo at z=%.2f
              \033[0m''' % z)
        try:
            f_b = _power_law(10 ** best_mass, [fb_a, fb_pow])
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

    output = x0 - (x0 - x1) * _np.exp(-x2 * f_b)

    return k, output
