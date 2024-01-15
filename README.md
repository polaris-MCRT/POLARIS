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
