# POLARIS Quickstart Guide


## Download

Download zip package from the [homepage](http://www1.astrophysik.uni-kiel.de/~polaris/)
or clone the [github repository](https://github.com/polaris-MCRT/POLARIS).


## Requirements

The following packages are required for the installation:

- gcc (preferred), icc, or clang++

- cmake (preferred), or ninja

- python3 and numpy


## Installation

- Extract the zip file

- Open a terminal / console and move into the POLARIS directory:\
`cd /YOUR/POLARIS/PATH/`

- Run the installation script:\
`./compile.sh -f`\
For the first installation, the option `-f` is required.
For more information, type:\
`./compile.sh -h`


## Start a simulation

POLARIS is shipped with a prebuild grid and command files to perform example simulations.

The grid (binary) file `grid.dat` of an example disk can be found in `projects/disk/`.

The command files of the temperature `temp.cmd`, thermal emission `dust.cmd`, and scattered stellar emission `dust_mc.cmd` can be found in `projects/disk/example/temp/`, `projects/disk/example/dust/`, and `projects/disk/example/dust_mc/`, respectively.

Before starting the simulation, the POLARIS path in these command files has to be adjusted:\
Change `/YOUR/POLARIS/PATH/` at `<dust_component>`, `<path_grid>`, and `<path_out>` to your POLARIS path.

To start a temperature simulation, type:\
`polaris /YOUR/POLARIS/PATH/projects/disk/example/temp/temp.cmd`\
The results are stored at `projects/disk/example/temp/data/` as `.fits` files. These files can be opened with, for example, [SAOImageDS9](https://sites.google.com/cfa.harvard.edu/saoimageds9/home).

Similar, the simulation for thermal emission and scattered stellar emission are performed.
For available options in the command file, please read the manual.


## Create a grid

The grid (binary) file will be created with the command `polaris-gen`.

There are already three models avaiable in `tools/polaris_tools_modules/model.py`:

- disk: A circumstellar disk with a Shakura & Sunyaev density distribution

- globule: A Bok globule with a Bonnor-Ebert sphere density distribution

- sphere: A sphere with a constant density distribution

To create a grid file with a globule model, type:\
`polaris-gen globule grid.dat`\
The grid file will be stored at `projects/globule/`.
It is also possible to modify some grid parameters with the command `polaris-gen`.
For more information, type:\
`polaris-gen -h`

To modify further parameter values such as the density distribution, it is recommended that the user defines their own models in `tools/polaris_tools_custom/model.py`.
Therein, each model is defined as a class with a corresponding entry in the dictionary at the top of `model.py`.
For any changes, the user has to recompile with:\
`./compile.sh -u`\
Similar, to create the grid file with a model named *custom*, type:\
`polaris-gen custom grid.dat`

It is also possible, to write their own grid file.
For the general structure and available options in the grid file, please read the manual.
For this purpose, the command `polaris-gen` has an ascii to binary converter (and vice versa) for the grid files.
To convert an existing ascii grid file of the `disk` model to a binary grid file, type:\
`polaris-gen --convert ascii2binary disk grid.txt`\
The ascii file has to be located in `projects/disk/` and the new binary grid file will be stored at `projects/disk/`.