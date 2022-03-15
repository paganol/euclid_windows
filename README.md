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

    - ``nbin``: number of bins, replaced if you provide bin ranges in the variable bintype

    - ``dz``: integration step in redshift

    - ``cb``, ``zb``, ``sigmab``, ``c0``, ``z0``, ``sigma0``, ``fout``: parameters of the 
    photo z distribution, see equation 115 and table 5 of [1910.09273](https://arxiv.org/abs/1910.09273)

    - ``bintype``: three options here, "equipopulated", "equispaced", numpy array with bin edges 

    - ``normalize``: normalization of the windows 
