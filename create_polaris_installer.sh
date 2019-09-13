#!/usr/bin/env bash

function missing_file () {
    echo $1" is missing to create the POLARIS package"
    exit
}

rm -rv /tmp/polaris

for files in 'lib/CCfits.tar.gz' \
            'lib/cfitsio.tar.gz' \
            'examples.tar.xz' \
            'AUTHORS' \
            'COPYING' \
            'ChangeLog' \
            'INSTALL' \
            'README' \
            'NEWS' \
            'manual.pdf' \
            'install_polaris.sh' \
            'src/CMakeLists.txt' \
            'bin/polaris' \
            'tools/COPYING' \
            'tools/AUTHORS' \
            'tools/README' \
            'tools/setup.py' \
            'tools/polaris-run.in' \
            'tools/polaris-plot.in' \
            'tools/polaris-gen.in' \
            'tools/polaris-remote.in' \
            'tools/polaris-extra.in' \
            'tools/polaris-test.in'

do
    if [ ! -f $files ]
    then
        missing_file $files
    fi
done

mkdir /tmp/polaris

cp -rv --parents \
    lib/CCfits.tar.gz \
    lib/cfitsio.tar.gz \
    examples.tar.xz \
    AUTHORS \
    COPYING \
    ChangeLog \
    INSTALL \
    README \
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
    bin/polaris \
    tools/COPYING \
    tools/AUTHORS \
    tools/README \
    tools/ChangeLog \
    tools/setup.py \
    tools/polaris-run.in \
    tools/polaris-plot.in \
    tools/polaris-gen.in \
    tools/polaris-remote.in \
    tools/polaris-extra.in \
    tools/polaris-test.in \
    tools/polaris_tools_modules/__init__.py \
    tools/polaris_tools_modules/*.py \
    tools/polaris_tools_custom/__init__.py \
    tools/polaris_tools_custom/*.py.empty \
    /tmp/polaris/
    
for f in /tmp/polaris/tools/polaris_tools_custom/*py.empty; do 
    mv -- "$f" "${f%.py.empty}.py"
done

makeself --notemp /tmp/polaris ./polaris.run "the radiative transfer code POLARIS" ./install_polaris.sh
