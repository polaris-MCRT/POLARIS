#!/usr/bin/env bash

function missing_file () {
    echo $1" is missing to create the POLARIS package"
    exit
}

for files in 'autogen.sh' \
                'CCfits.tar.gz' \
                'cfitsio.tar.gz' \
                'examples.tar.xz' \
                'configure.ac' \
                'AUTHORS' \
                'COPYING' \
                'ChangeLog' \
                'INSTALL' \
                'README' \
                'Makefile.am' \
                'NEWS' \
                'manual.pdf' \
                'install_polaris.sh' \
                'src/Makefile.am' \
                'tools/configure.ac' \
                'tools/Makefile.am' \
                'tools/COPYING' \
                'tools/AUTHORS' \
                'tools/INSTALL' \
                'tools/README' \
                'tools/ChangeLog' \
                'tools/NEWS' \
                'tools/autogen.sh' \
                'tools/src/Makefile.am' \
                'tools/src/polaris-run.in' \
                'tools/src/polaris-plot.in' \
                'tools/src/polaris-gen.in' \
                'tools/src/polaris-remote.in' \
                'tools/src/polaris-extra.in' \
                'tools/src/polaris-test.in' \
                'tools/src/modules/Makefile.am' \
                'tools/src/custom/Makefile.am'
do
    if [ ! -f $files ]
    then
        missing_file $files
    fi
done

mkdir /tmp/polaris

cp -rv --parents \
    autogen.sh \
    CCfits.tar.gz \
    cfitsio.tar.gz \
    examples.tar.xz \
    configure.ac \
    AUTHORS \
    COPYING \
    ChangeLog \
    INSTALL \
    README \
    Makefile.am \
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
    src/Makefile.am \
    tools/configure.ac \
    tools/Makefile.am \
    tools/COPYING \
    tools/AUTHORS \
    tools/INSTALL \
    tools/README \
    tools/ChangeLog \
    tools/NEWS \
    tools/autogen.sh \
    tools/src/Makefile.am \
    tools/src/polaris-run.in \
    tools/src/polaris-plot.in \
    tools/src/polaris-gen.in \
    tools/src/polaris-remote.in \
    tools/src/polaris-extra.in \
    tools/src/polaris-test.in \
    tools/src/modules/*.py \
    tools/src/modules/Makefile.am \
    tools/src/custom/*.py.empty \
    tools/src/custom/Makefile.am \
    /tmp/polaris/
    
for f in /tmp/polaris/tools/src/custom/*py.empty; do 
    mv -- "$f" "${f%.py.empty}.py"
done

makeself --notemp /tmp/polaris ./polaris.run "the radiative transfer code POLARIS" ./install_polaris.sh
