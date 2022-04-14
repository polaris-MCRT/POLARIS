# POLARIS (POLArized RadIation Simulator)

is a 3D Monte Carlo continuum radiative transfer code.

- Simulate intensity and polarization of light emerging from analytical astrophysical models as well as complex magneto-hydrodynamic simulations on various grids

- Perform dust heating, -emission, -scattering, -grain alignment, line radiative transfer, and synchrotron simulations

- Calculate synthetic intensity and polarization maps

- Make use of a full set of physical quantities (density, temperature, velocity, magnetic field distribution, and dust grain properties as well as different sources of radiation) as input

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

## Use POLARIS

For a guide how to run first simulations, please take a look in our [quickstart](quickstart.pdf).

For more information about POLARIS and its capabilities, please take a look in our [manual](manual.pdf).

contact: polaris@astrophysik.uni-kiel.de