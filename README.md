# POLARIS (POLArized RadIation Simulator)

[![arXiv](https://img.shields.io/badge/arXiv-1604.05305-red)](https://arxiv.org/abs/1604.05305)
[![Bibcode](https://img.shields.io/badge/Bibcode-2016A%26A...593A..87R-blue)](https://ui.adsabs.harvard.edu/abs/2016A&A...593A..87R)
[![doi](https://img.shields.io/badge/doi-10.1051%2F0004--6361%2F201424930-yellow)](https://doi.org/10.1051/0004-6361/201424930)

is a 3D Monte Carlo radiative transfer code that

- allows to simulate intensity and polarization of light emerging from analytical astrophysical models as well as complex magneto-hydrodynamic simulations on various grids

- is capable to perform dust heating, -emission, -scattering, -grain alignment, line radiative transfer, and synchrotron simulations

- calculates synthetic intensity and polarization maps

- makes use of a full set of physical quantities (density, temperature, velocity, magnetic field distribution, and dust grain properties as well as different sources of radiation) as input


## Requirements

The following packages are required for the installation:

- gcc (preferred), icc, or clang++

- cmake (preferred), or ninja

- python3 (packages: *numpy*, *setuptools*)


## Installation

Open a terminal/console and move into the POLARIS directory:
```bash
cd /YOUR/POLARIS/PATH/
```

To install POLARIS, run the installation script:
```bash
./compile.sh -f
```
For the first installation, the option `-f` is required to install the cfitsio and CCfits libraries.
For more information, type:
```bash
./compile.sh -h
```
POLARIS can now be executed from any newly opened terminal/console.
However, to use it in already open terminals/consoles, execute the following command to update the environmental paths:
```bash
source ~/.bashrc
```


## Use POLARIS

For a guide how to run first simulations, please take a look in our [quickstart](quickstart.pdf).

For more information about POLARIS and its capabilities, please take a look in our [manual](manual.pdf).

If you use results from POLARIS in a publication, please cite [Reissl et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016A%26A...593A..87R) and the [website](http://www1.astrophysik.uni-kiel.de/~polaris).
If line radiative transfer and/or Zeeman simulations are used, please cite [Brauer et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017A%26A...601A..90B) as well.


## Copyright

The code is free of charge for any scientific purpose. This software is provided in the hope that it will
be useful but without any warranty of ability or fitness of a particular purpose. We also reject any
responsibility for incorrect result that may be result from this code.

For a list of authors, please see [COPYING](COPYING.md).

**contact**: polaris@astrophysik.uni-kiel.de
