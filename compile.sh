#!/bin/bash

# Define colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'


CXX=${1:-icc}

type $CXX >/dev/null 2>&1 || { echo >&2 "I require $CXX but it's not installed.  Aborting."; exit 1; }


if [ -z "$CO" ]
then
      CO="fast"
fi

if [ -z "$TARGET" ]
then
      TARGET="host"
fi


if [ $CXX = "icc" ]; then
    if [ $TARGET = "host" ]; then
        TARGET=$(hostname)
    fi
    case $TARGET in
        prometheus)
            CXXFLAGS="-march=skylake-avx512"
            ;;
        nesh*)
            CXXFLAGS="-xCORE-AVX512"
            ;;
        fenrir)
            CXXFLAGS="-I/usr/include/x86_64-linux-gnu/c++/8/ -march=skylake"
            ;;
        oberon)
            CXXFLAGS="-march=broadwell"
            ;;
        hera|rhea|hydra)
            CXXFLAGS="-march=ivybridge"
            ;;
        atlas)
            CXXFLAGS="-march=sandybridge"
            ;;
        *)
            CXXFLAGS="-xHost"
            ;;
    esac

    if [ $CO = "debug" ]; then
	   CXXFLAGS="-O0 -g -debug inline-debug-info"

    elif [ $CO = "fast" ]; then
        CXXFLAGS="-O3 -parallel -ip -ipo $CXXFLAGS"

    elif [ $CO = "fast-multi" ]; then
        CXXFLAGS="-O3 -parallel -ip -ipo -xAVX -axCORE-AVX512,CORE-AVX2"

    elif [ $CO = "fast-debug" ]; then
        CXXFLAGS="-O3 -parallel -ip -ipo -g -debug inline-debug-info $CXXFLAGS"

    elif [ $CO = "fast-debug-multi" ]; then
        CXXFLAGS="-O3 -parallel -ip -ipo -g -debug inline-debug-info -xAVX -axCORE-AVX512,CORE-AVX2"
    fi

    if [ $TARGET = "fenrir" ]; then
            CXXFLAGS="$CXXFLAGS -I/usr/include/x86_64-linux-gnu/c++/8/"
    fi

elif [ $CXX = "g++" ]; then
    if [ $CO = "debug" ]; then
	   CXXFLAGS="-O0 -g3 -Wall -Wno-unused-function -Wno-unused-but-set-variable -Wno-comment -Wno-unknown-pragmas"
    elif [ $CO = "fast" ]; then
        CXXFLAGS="-march=native -O3 -funroll-loops -ffast-math -fno-finite-math-only -flto"
    fi
fi

CXXFLAGS="$CXXFLAGS -fopenmp"


# remove old bin
rm bin/polaris

# remove old build dir
cd src/
rm -r build

mkdir build
cd build

# configure POLARIS with appropriate flags
cmake --parallel .. -DBUILD_SHARED_LIBS=ON -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_CXX_FLAGS="$CXXFLAGS" \
    && echo -e "Configuring POLARIS [${GREEN}done${NC}]" \
    || { echo -e  "Configuring POLARIS [${RED}Error${NC}]"; exit; }

# compile and install POLARIS
make \
    && echo -e "Compiling POLARIS [${GREEN}done${NC}]" \
    || { echo -e  "Compiling POLARIS [${RED}Error${NC}]"; exit; }
make install \
    && echo -e "Installing POLARIS [${GREEN}done${NC}]" \
    || { echo -e  "Installing POLARIS [${RED}Error${NC}]"; exit; }

# go to PolarisTools dir
cd ../../tools/

# compile and install PolarisTools
python setup.py install --user &>/dev/null \
    && echo -e "Compiling and installing PolarisTools [${GREEN}done${NC}]" \
    || { echo -e  "Compiling and installing PolarisTools [${RED}Error${NC}]"; exit; }