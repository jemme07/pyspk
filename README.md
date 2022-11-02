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


### Example 1: Using a power-law fit to the $f_b$ - $M_\mathrm{halo}$ relation

In our first example, we will provide py-SP(k) with the (approximate)power-law fitted parameters to the fiducial BAHAMAS simulations ([McCarthy et al. 2017](https://academic.oup.com/mnras/article/465/3/2936/2417021)) at redshift $z=0.125$. We will use a spherical-overdensity of 200 times the critical density of the Universe. The functional form used is:

$$f_b / (\Omega_b/\Omega_m) = a \left(\frac{M_{200c}}{10^{13.5}\mathrm{M}_\odot} \right)^{b}$$

For this example, the power-law has been normalised at pivot point of $M_\mathrm{halo} = 10^{13.5} \mathrm{M}_\odot$. If a pivot point is not given `spk.sup_model()` a default pivot point of $M_\mathrm{halo} = 1 \mathrm{M}_\odot$ will be assumed.  
A simple example:

```
import pyspk as spk

z = 0.125
fb_a = 8.44e-05
fb_pow = 0.275
k, sup = spk.sup_model(SO=200, z=z, fb_a=fb_a, fb_pow=fb_pow, k_max=10)

```

A jupyter notebook with more examples can be found within this [repository](https://github.com/jemme07/pyspk/blob/main/examples/pySPk_Examples.ipynb). 


## Acknowledging the code

Please cite spk using:


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


