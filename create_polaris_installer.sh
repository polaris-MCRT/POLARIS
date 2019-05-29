#!/usr/bin/env bash

function missing_file () {
    echo $1" is missing to create the POLARIS package"
    exit
}

for files in 'CCfits.tar.gz' \
            'cfitsio.tar.gz' \
            'examples.tar.xz' \
            'AUTHORS' \
            'COPYING' \
            'ChangeLog' \
            'INSTALL' \
            'README' \
            'CMakeLists.txt' \
            'NEWS' \
            'manual.pdf' \
            'install_polaris.sh' \
            'src/CMakeLists.txt' \
            'tools/COPYING' \
            'tools/AUTHORS' \
            'tools/README' \
            'tools/setup.py' \
            'tools/polaris-run.in' \
            'tools/polaris-plot.in' \
            'tools/polaris-gen.in' \
            'tools/polaris-remote.in' \
            'tools/polaris-extra.in' \
            'tools/polaris-test.in' \
            'bin/polaris'
do
    if [ ! -f $files ]
    then
        missing_file $files
    fi
done

mkdir /tmp/polaris

cp -rv --parents \
    CCfits.tar.gz \
    cfitsio.tar.gz \
    examples.tar.xz \
    AUTHORS \
    COPYING \
    ChangeLog \
    INSTALL \
    README \
    CMakeLists.txt \
    NEWS \
    manual.pdf \
    install_polaris.sh \
    input/*.dat \
    input/gas \
    input/dust/silicate* \
    input/dust/graphite* \
    input/dust/aOlM5* \
    input/dust/aPyM5* \
    input/dust/CM20* \
    src/*.cpp \
    src/*.h \
    src/*.hh \
    src/*.cc \
    src/CMakeLists.txt \
    tools/COPYING \
    tools/AUTHORS \
    tools/README \
    tools/ChangeLog \
    tools/setup.py.txt \
    tools/polaris-run.in \
    tools/polaris-plot.in \
    tools/polaris-gen.in \
    tools/polaris-remote.in \
    tools/polaris-extra.in \
    tools/polaris-test.in \
    tools/polaris_tools_modules/*.py \
    tools/polaris_tools_custom/*.py.empty \
    bin/polaris \
    /tmp/polaris/
    
for f in /tmp/polaris/tools/src/custom/*py.empty; do 
    mv -- "$f" "${f%.py.empty}.py"
done

makeself --notemp /tmp/polaris ./polaris.run "the radiative transfer code POLARIS" ./install_polaris.sh
