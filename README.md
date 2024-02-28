# EUCLID photometric window functions 

This repository builts EUCLID photometric window functions

## Installation

```
git clone https://github.com/paganol/euclid_windows.git 
cd euclid_windows
pip install -e .
```

## Usage
```
import euclid_windows as EW

#This initializes the container 
Win = EW.Windows()

#This computes the windows and other useful variables
Win.get_distributions()

#This returns the list of window functions ready for camb
sources = Win.get_camb_distributions()
```

## Parameters

- ``zmin``: minimum redshift, replaced if you provide bin ranges in the variable bintype.

- ``zmax``: maximum redshift, replaced if you provide bin ranges in the variable bintype.

- ``zminsampled``: minimum redshift sampled, if not provided ``zmin`` is used.

- ``zmaxsampled``: maximum redshift sampled, if not provided ``zmax`` is used. 

- ``nbin``: number of bins, replaced if you provide bin ranges in the variable bintype.

-  ``use_true_galactic_dist``: use the true galactic distribution for the windows, no 
   convolution with photo z distribution.

- ``dz``: integration step in redshift.

- ``cb``, ``zb``, ``sigmab``, ``c0``, ``z0``, ``sigma0``, ``fout``: parameters of the 
    photo z distribution, see equation 115 and table 5 of 1910.09273.

- ``bintype``: three options here, "equipopulated", "equispaced", numpy array or list with 
  bin edges. 

- ``normalize``: normalization of the windows.

- ``biastype``: several options here:
   - "stepwise" with a different constant value for each bin, using $f(z)=\sqrt{1+z}$;
   - "continuous" which implements a continuous function $f(z)=\sqrt{1+z}$;
   - "tutusaus_Flag1" which implements Tutusaus bias (Flagship1);
   - "tutusaus_Flag2" which implements Tutusaus bias (Flagship2);
   - numpy array (or list) with bias provided by the user for each bin.
 
- ``errortype``: the default option is "gauss_err", because we expect a gaussian error and 
    then we can compute the galaxy selection functions via an erf function; if the error is not 
    gaussian, we need to compute the integral of the probability distribution function to determine 
    the galaxy selection functions.
