# The Transient Sky Simulator

This is the Transient Sky Simulator, which provides the user with a way to estimate transient numbers and brightnesses for synoptic sky surveys.

## Prerequisites

The following Python (2.7) packages are required:

* numpy
* scipy
* matplotlib
* time
* glob
* h5py
* csv
* astropy
* astLib
* itertools
* [Cosmolopy](https://roban.github.io/CosmoloPy/)
* [dustmaps](https://github.com/gregreen/dustmaps#dustmaps)
* [astroplan](https://github.com/astropy/astroplan) (recommended for observation plan creation)
* [sncosmo](https://github.com/sncosmo/sncosmo/tree/v1.6.x)   (Only necessary if you want to add extragalactic transients to the preexisting ones)

For the creation of observation plans you will also want to install Jupyter Notebooks.

It is strongly recommended to use the online version of the dustmaps package. For the purpose of TSS this is much faster. In case you don't have an internet connection while running TSS, you can also use the off-line version of dustmaps. Before first time off-line usage one needs to download the dust maps 'bayestar' and 'sfd' for the dustmaps package as described [here](https://github.com/gregreen/dustmaps#dustmaps)
