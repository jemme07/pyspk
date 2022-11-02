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


py-SP(k) is not restrictive to a particular shape of the baryon fraction – halo mass relation, and in order to provide flexibility to the user, we have implemented 3 different method to provide the $f_b$ - $M_\mathrm{halo}$ relation to the mode.

### Method 1: Using a power-law fit to the $f_b$ - $M_\mathrm{halo}$ relation

py-SP(k) can be provided with power-law fitted parameters to the $f_b$ - $M_\mathrm{halo}$ relation using the functional form:

$$f_b/(\Omega_b/\Omega_m)=a\left(\frac{M_{200c}}{f_{b,\mathrm{pivot}}\mathrm{M}_\odot}\right)^{b}$$

The power-law can be normalised at any pivot point in units of $\mathrm{M}_ {\odot}$. If a pivot point is not given, `spk.sup_model()` uses a default pivot point of $M_\mathrm{halo} = 1 \mathrm{M}_\odot$.  


### Method 2: Redshift-dependent power-law fit to the $f_b$ - $M_\mathrm{halo}$ relation. 

We have implemented 


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


