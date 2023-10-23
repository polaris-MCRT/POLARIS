#!/bin/bash

# script to run POLARIS tests on CI
# parameters:
#  - POLARIS_BINARY = command to run polaris, default: polaris
#  - PROJECTS_PATH = location of projects folder, default: projects

set -e -x

POLARIS_BINARY=${1:-polaris}
PROJECTS_PATH=${2:-projects}

# set path to CCfits and cfitsio libraries
export POLARIS_FITS_PATH=./lib/lib
if [[ :$LD_LIBRARY_PATH: != *:$POLARIS_FITS_PATH:* ]]; then
    export LD_LIBRARY_PATH=$POLARIS_FITS_PATH:$LD_LIBRARY_PATH
fi

# run disk temp example
${POLARIS_BINARY} ${PROJECTS_PATH}/test/stellar_sed/POLARIS.cmd

# validate output
python3 ${PROJECTS_PATH}/test/stellar_sed/compare.py
