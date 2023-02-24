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

- ``zmin``: minimum redshift, replaced if you provide bin ranges in the variable bintype

- ``zmax``: maximum redshift, replaced if you provide bin ranges in the variable bintype

- ``zmaxsampled``: maximum redshift sampled, if not provied ``zmax`` is used 

- ``nbin``: number of bins, replaced if you provide bin ranges in the variable bintype

-  ``use_true_galactic_dist``: use the true galactic distribution for the windows, no 
   convolution with photo z distribution

- ``dz``: integration step in redshift

- ``cb``, ``zb``, ``sigmab``, ``c0``, ``z0``, ``sigma0``, ``fout``: parameters of the 
    photo z distribution, see equation 115 and table 5 of 1910.09273

- ``bintype``: three options here, "equipopulated", "equispaced", numpy array or list with 
  bin edges 

- ``normalize``: normalization of the windows

- ``biastype``: three options here: "stepwise" with a different constant value for each bin,
  "continuous" which implements a continuous function (in both cases $\sqrt{1+z}$ is used), or
  a numpy array (or list) with bias for each bin.
 
