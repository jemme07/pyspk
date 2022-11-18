# py-SP(k) - A hydrodynamical simulation-based model for the impact of baryon physics on the non-linear matter power spectrum
                          _____ ____   ____   _ 
        ____  __  __     / ___// __ \_/_/ /__| |
       / __ \/ / / /_____\__ \/ /_/ / // //_// /
      / /_/ / /_/ /_____/__/ / ____/ // ,<  / / 
     / .___/\__, /     /____/_/   / //_/|_|/_/  
    /_/    /____/                 |_|    /_/    

py-SP(k) is a python package aimed at predicting the suppression of the total matter power spectrum due to baryonic physics as a function of the baryon fraction of haloes and redshift.

## Requirements

The module requires the following:

- numpy
- scipy

## Installation

The easiest way to install py-SP(k) is using pip:

```
pip install pyspk [--user]
```

The --user flag may be required if you do not have root privileges.

## Usage

py-SP(k) is not restrictive to a particular shape of the baryon fraction – halo mass relation. In order to provide flexibility to the user, we have implemented 3 different methods to provide py-SP(k) with the required $f_b$ - $M_\mathrm{halo}$ relation. In the following sections we describe these implementations. A jupyter notebook with more detailed examples can be found within this [repository](https://github.com/jemme07/pyspk/blob/main/examples/pySPk_Examples.ipynb). 

### Method 1: Using a power-law fit to the $f_b$ - $M_\mathrm{halo}$ relation

py-SP(k) can be provided with power-law fitted parameters to the $f_b$ - $M_\mathrm{halo}$ relation using the functional form:

$$f_b/(\Omega_b/\Omega_m)=a\left(\frac{M_{SO}}{M_{\mathrm{pivot}}}\right)^{b},$$

where $M_{SO}$ could be either $M_{200c}$ or $M_{500c}$ in $\mathrm{M}_ \odot$, $a$ is the normalisation of the $f_b$ - $M_\mathrm{halo}$ relation at $M_\mathrm{pivot}$, and $b$ is the power-law slope. The power-law can be normalised at any pivot point in units of $\mathrm{M}_ {\odot}$. If a pivot point is not given, `spk.sup_model()` uses a default pivot point of $M_{\mathrm{pivot}} = 1 \mathrm{M}_ \odot$. $a$, $b$ and $M_\mathrm{pivot}$ can be specified at each redshift independently.  

Next, we show a simple example using power-law fit parameters:

```
import pyspk as spk

z = 0.125
fb_a = 0.4
fb_pow = 0.3
fb_pivot = 10 ** 13.5

k, sup = spk.sup_model(SO=200, z=z, fb_a=fb_a, fb_pow=fb_pow, fb_pivot=fb_pivot)
```

### Method 2: Redshift-dependent power-law fit to the $f_b$ - $M_\mathrm{halo}$ relation. 

For the mass range that can be relatively well probed in current X-ray and Sunyaev-Zel'dovich effect observations (roughly $10^{13} \leq M_{500c} [\mathrm{M}_ \odot] \leq 10^{15}$), the total baryon fraction of haloes can be roughly approximated by a power-law with constant slope (e.g. Mulroy et al. 2019; Akino et al. 2022). Akino et al. 2022 determined the of the baryon budget for X-ray-selected galaxy groups and clusters using weak-lensing mass measurements. They provide a parametric redshift-dependent power-law fit to the gas mass - halo mass and stellar mass - halo mass relations, finding very little redshift evolution. 

We implemented a modified version of the functional form presented in Akino et al. 2022, to fit the total $f_b$ - $M_\mathrm{halo}$ relation as follows:

$$f_b/(\Omega_b/\Omega_m)= \left(\frac{0.1658}{\Omega_b/\Omega_m}\right) \left(\frac{e^\alpha}{100}\right) \left(\frac{M_{500c}}{10^{14} \mathrm{M}_ \odot}\right)^{\beta - 1} \left(\frac{E(z)}{E(0.3)}\right)^{\gamma},$$

where $\alpha$ sets the power-law normalisation, $\beta$ sets power-law slope, $\gamma$ provides the redshift dependence and $E(z)$ is the usual dimensionless Hubble parameter. For simplicity, we use the cosmology implementation of `astropy` to specify the cosmological parameters in py-SP(k).

Note that this power-law has a normalisation that is redshift dependent, while the the slop is constant in redshift. While this provides a less flexible approach compared with Methods 1 (simple power-law) and Method 3 (binned data), we find that this parametrisation agrees well with our simulations up to redshift $z=1$, which is the redshift range proved by Akino et al. 2022. For higher redshifts, we find that simulations require a mass-dependent slope, especially at the lower mass range required to predict the suppression of the total matter power spectrum at such redshifts. 

In the following example we use the redshift-dependent power-law fit parameters with a flat LambdaCDM cosmology. Note that any `astropy` cosmology could be used instead.

```
import pyspk.model as spk
from astropy.cosmology import FlatLambdaCDM

H0 = 70 
Omega_b = 0.0463
Omega_m = 0.2793

cosmo = FlatLambdaCDM(H0=H0, Om0=Omega_m, Ob0=Omega_b) 

alpha = 4.189
beta = 1.273
gamma = 0.298
z = 0.5

k, sup = spk.sup_model(SO=500, z=z, alpha=alpha, beta=beta, gamma=gamma, cosmo=cosmo)
```

### Method 3: Binned data for the $f_b$ - $M_\mathrm{halo}$ relation. 

The final, and most flexible method is to provide py-SP(k) with the baryon fraction binned in bins of halo mass. This could be, for example, obtained from observational constraints, measured directly form simulations, or sampled from a predefined distribution or functional form. For an example using data obtained from the BAHAMAS simulations ([McCarthy et al. 2017](https://academic.oup.com/mnras/article/465/3/2936/2417021)), please refer to the [examples](https://github.com/jemme07/pyspk/blob/main/examples/pySPk_Examples.ipynb) provided. 


## Priors

While py-SP(k) was calibrated using a wide range of sub-grid feedback parameters, some applications may require a more limited range of baryon fractions that encompass current observational constraints. For such applications, we used the gas mass - halo mass and stellar mass - halo mass constraints from Table 5 in Akino et al. 2022, and find the subset of simulations from our 400 models that agree with the inferred baryon budget at redshift $z=0.1$. We note that we constrained our simulations to within a normalisation of $\pm 3 \times \sigma$ at $M_{500c} = 10^{14} \mathrm{M}_ \odot$. 

Using the simulations that fall within these constraints, we can impose observational priors for the redshift-dependent power-law fitting parameters for the $f_b$ - $M_\mathrm{halo}$ relation in Method 3 as follows:

| Parameter   | Description        | Prior           |
| ----------- | ------------------ | --------------- |
| $\alpha$    | Normaliasation     | G(4.189, 0.066) |
| $\beta$     | Slope              | G(1.273, 0.044) |
| $\gamma$    | Redshift evolution | G(0.298, 0.063) |

where G(x, y) is a Gaussian distribution with center x and width y.

A less conservative approach could be to use a flat priors over the entire range of parameters fitted to simulations that fall within Akino et al. 2022 these constraints:

| Parameter   | Description        | Prior           |
| ----------- | ------------------ | --------------- |
| $\alpha$    | Normaliasation     | U(4.060, 4.306) |
| $\beta$     | Slope              | U(1.199, 1.347) |
| $\gamma$    | Redshift evolution | U(0.159, 0.414) |

where U(x, y) is a uniform distribution over [x, y]. 

Finally, the full range of fitted parameters spanned by our simulations, regardless of whether or not they agree with observational constraints is:

| Parameter   | Description        | Prior           |
| ----------- | ------------------ | --------------- |
| $\alpha$    | Normaliasation     | U(3.060, 4.508) |
| $\beta$     | Slope              | U(0.989, 1.620) |
| $\gamma$    | Redshift evolution | U(0.046, 0.631) |

where U(x, y) is a uniform distribution over [x, y]. 

## Acknowledging the code

Please cite py-SP(k) using:

```
@ARTICLE{spk,
       author = {},
        title = "{}",
      journal = {\mnras},
     keywords = {},
         year = ,
        month = ,
       volume = {},
       number = {},
        pages = {},
          doi = {},
archivePrefix = {arXiv},
       eprint = {},
 primaryClass = {astro-ph.CO},
       adsurl = {},
      adsnote = {}
}
```
For any questions and enquires please contact me via email at *j.salcidonegrete@ljmu.ac.uk*


