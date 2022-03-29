<!-- create PDF file with 'pandoc --pdf-engine=pdflatex -V colorlinks -V --highlight-style tango quickstart.md -o quickstart.pdf' -->
# POLARIS Quickstart Guide


## Download

Download zip package from the [homepage](http://www1.astrophysik.uni-kiel.de/~polaris/) or clone the [github repository](https://github.com/polaris-MCRT/POLARIS).


## Requirements

The following packages are required for the installation:

- gcc (preferred), icc, or clang++

- cmake (preferred), or ninja

- python3 and numpy


## Installation

- Extract the zip file

- Open a terminal/console and move into the POLARIS directory:
```bash
cd /YOUR/POLARIS/PATH/
```

- Run the installation script:
```bash
./compile.sh -f
```
For the first installation, the option `-f` is required.
For more information, type:
```bash
./compile.sh -h
```


## Start a simulation

POLARIS is shipped with a prebuild grid and command files to perform example simulations.
The grid (binary) file `grid.dat` of an example disk can be found in `projects/disk/`.
The command files `.cmd` of the temperature, thermal emission, and scattered stellar emission can be found in

- `projects/disk/example/temp/`,
- `projects/disk/example/dust/`, and
- `projects/disk/example/dust_mc/`, respectively.

Before starting the simulation, change `/YOUR/POLARIS/PATH/` in the command file at `<dust_component>`, `<path_grid>`, and `<path_out>` to your POLARIS path.

To start a temperature simulation, type:
```bash
polaris projects/disk/example/temp/temp.cmd
```
The results are stored at `projects/disk/example/temp/data/` as `.fits` files. These files can be opened with, for example, [SAOImageDS9](https://sites.google.com/cfa.harvard.edu/saoimageds9/home).

Similar, the simulation for thermal emission and scattered stellar emission are performed.
For available options in the command file, please read the manual.


## Create a grid


### Predefined models

The grid (binary) file will be created with the command `polaris-gen`.
There are already three models available:

- disk: A circumstellar disk with a [Shakura & Sunyaev](https://ui.adsabs.harvard.edu/abs/1973A&A....24..337S) density distribution
([Lynden-Bell & Pringle 1974](https://ui.adsabs.harvard.edu/abs/1974MNRAS.168..603L); [Hartmann et al. 1998](https://ui.adsabs.harvard.edu/abs/1998ApJ...495..385H))
$$ \rho(r, z) = \rho_0 \left( \frac{r}{r_0} \right)^{-\alpha} \times \exp\left[ -\frac{1}{2} \left( \frac{z}{h(r)} \right)^2 \right] $$
$$ h(r) = h_0 \left( \frac{r}{r_0} \right)^\beta $$
Here, the default values are $r_0 = 100\,\mathrm{AU}$, $h_0 = 10\,\mathrm{AU}$, $\alpha = 0.9$, and $\beta = \frac{2 \alpha + 3}{6} = 0.8$.

- globule: A Bok globule with a [Bonnor](https://ui.adsabs.harvard.edu/abs/1956MNRAS.116..351B)-[Ebert](https://ui.adsabs.harvard.edu/abs/1955ZA.....37..217E) sphere density distribution
([Harvey et al. 2001](https://ui.adsabs.harvard.edu/abs/2001ApJ...563..903H); [Kaminski et al. 2014](https://ui.adsabs.harvard.edu/abs/2014ApJ...790...70K))
$$ \rho(r) = \rho_0 \begin{cases}
r_\mathrm{t}^{-2} & \text{if}\ r \leq r_\mathrm{t}\\
r^{-2} & \text{if}\ r_\mathrm{t} < r
\end{cases}$$
Here, the default value is $r_\mathrm{t} = 10^3\,\mathrm{AU}$.

- sphere: A sphere with a constant density distribution
$$ \rho(r) = \rho_0 $$

By default, the density distribution is normalized to the given total mass.
To create a grid file with a globule model, type:
```bash
polaris-gen globule grid.dat
```
The grid file will be stored at `projects/globule/`.
It is also possible to modify some grid parameters with the command `polaris-gen`.
For more information, type:
```bash
polaris-gen -h
```


### Extra parameter

To modify further parameter values, the user can parse a list of parameter values using the option `--extra` followed by a list of values (int, float, or str).
These additional parameter values can be used in the function `update_parameter` in the file `model.py` to vary the model.
For example, the user can parse

- 4 values for the `disk` model: reference radius $r_0$, reference scale height $h_0$, $\alpha$, and $\beta$,

- 1 value for the `globule` model: truncation radius $r_\mathrm{t}$, and

- 1 value for the `sphere` model: the geometry of the magnetic field (toroidal, vertical, or radial).

**Hint**: For any changes in the files, the user has to recompile with:
```bash
./compile.sh -u
```


### Custom model

For a more complex model modification, it is recommended that users define their own models in `tools/polaris_tools_custom/model.py`.
Therein, each model is defined as a class with a corresponding entry in the dictionary at the top of `model.py`.
**Hint**: For any changes in the files, the user has to recompile with:
```bash
./compile.sh -u
```
Similar, to create the grid file `grid.dat` with the model named *custom*, type:
```bash
polaris-gen custom grid.dat
```


### Write a grid file

It is also possible, to write their own grid file.
For the general structure and available options in the grid file, please read the manual.
For this purpose, the command `polaris-gen` has an ascii to binary converter (and vice versa) for the grid files.
To convert an existing ascii grid file of the `disk` model to a binary grid file, type:
```bash
polaris-gen --convert ascii2binary disk grid.txt
```
The ascii file has to be located in `projects/disk/` and the new binary grid file will be stored at `projects/disk/`.