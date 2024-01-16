# POLARIS (POLArized RadIation Simulator)

is a 3D Monte Carlo radiative transfer code that numerically solves the radiative transfer equation. It is particularly capable of simulating the Mie scattering of laser light in optically thick nanodusty plasmas.

## Project structure

```
bin             <- folder to store executable after compilation
ext             <- folder containing scripts for testing
input           <- folder containing input files for dust properties
lib             <- folder containing the cfitsio and CCfits libraries
projects        <- folder for outputs containing examples
src             <- main folder containing the project source files
tools           <- folder containing python scripts for creating the grid file

COPYING.md      <- information on authors and license
QUICKSTART.md   <- a short introduction on how to install and use polaris (markdown)
README.md       <- this file
compile.sh      <- the compile script
manual.pdf      <- a detailed documentation of the software
quickstart.pdf  <- a short introduction on how to install and use polaris (pdf)

```

## Requirements

The following packages are required for the installation:

- gcc (preferred), icc, or clang++

- cmake (preferred), or ninja

- python3 (packages: *numpy*, *setuptools*)


## Installation (Linux)

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

**HINT**: Please refer to the [manual](manual.pdf) for installation on **macOS**. An installer to use POLARIS with Windows is not available yet.


## Use POLARIS

For a guide how to run first simulations, please take a look in our [quickstart](quickstart.pdf).

For more information about POLARIS and its capabilities, please take a look in our [manual](manual.pdf).

Pre-calculated exemplary simulation results can be found in `projects/constantCylinder/example1/data/` and `projects/constantCylinder/example1/data/`. To re-run the corresponding exemplary `example1.cmd` command file in `projects`, move into the POLARIS directory and execute `polaris` followed by the command file:
```bash
cd /YOUR/POLARIS/PATH/
polaris projects/example1.cmd
```

If you use results from POLARIS in a publication, please cite [Reissl et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016A%26A...593A..87R) or [Reissl et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018ascl.soft07001R).
If line radiative transfer and/or Zeeman simulations are used, please cite [Brauer et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017A%26A...601A..90B) as well.

## Copyright

Copyright (C) 2024  Stefan Reißl

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

For a list of authors and a copy of the GNU General Public License, please see [COPYING](COPYING.md).

**contact**: polaris@astrophysik.uni-kiel.de

[![arXiv](https://img.shields.io/badge/arXiv-1604.05305-b31b1b)](https://arxiv.org/abs/1604.05305)
[![ascl](https://img.shields.io/badge/ascl-1807.001-262255)](https://ascl.net/1807.001)
[![bibcode](https://img.shields.io/badge/bibcode-2016A%26A...593A..87R-1c459b)](https://ui.adsabs.harvard.edu/abs/2016A&A...593A..87R)
[![doi](https://img.shields.io/badge/doi-10.1051%2F0004--6361%2F201424930-fab70c)](https://doi.org/10.1051/0004-6361/201424930)
[![License](https://img.shields.io/badge/License-GPLv3-blue)](https://www.gnu.org/licenses/gpl-3.0)
