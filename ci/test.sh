#!/bin/bash

# script to run POLARIS tests on CI
# parameters:
#  - POLARIS_BINARY = command to run polaris, default: polaris
#  - PROJECTS_PATH = location of projects folder, default: projects

set -e -x

POLARIS_BINARY=${1:-polaris}
PROJECTS_PATH=${2:-projects}

# set paths in cmd files
sed -i.bak 's|/PATH/TO/POLARIS|'`pwd`'|' ${PROJECTS_PATH}/CommandList.cmd
sed -i.bak 's|/YOUR/POLARIS/PATH|'`pwd`'|' ${PROJECTS_PATH}/disk/example/temp/POLARIS.cmd

# run disk temp example
${POLARIS_BINARY} ${PROJECTS_PATH}/disk/example/temp/POLARIS.cmd

# todo: validate output, for now just display a list of generated files
ls ${PROJECTS_PATH}/disk/example/temp/*/*
