#!/bin/bash

# Define colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

# set pipefail so "$?" gives 1 for piping of "false | true" for instance
# necessary because of piping to sed
set -o pipefail


# ================================================================================ #
# ================================ Default values ================================ #
# ================================================================================ #

# default compiler
CC="gcc"
CXX="g++"

# default cmake generator
CMAKE_GENERATOR="Ninja"

# default compiler flag
CO="fast"

# install target
TARGET="host"

# do not install fits libraries each time
DO_FITS=false

# do install PolarisTools
DO_PYTHON=true


# ================================================================================ #
# ============================ install FITS libraries ============================ #
# ================================================================================ #

function install_fits_libraries()
{
    # Check if lib directory exists
    echo -e "--- Install required libraries for FITS support ---"
    if [ ! -d "lib" ]; then
        echo -e "--- ${RED}Error:${NC} lib directory not found (incomplete polaris package?)"
        exit 1
    fi
    cd "lib"

    echo ""

    # cfitsio
    rm -rf "cfitsio/"
    tar -xf cfitsio.tar.gz
    cd "cfitsio"
    mkdir -p "build"

    echo -e "Configure cfitsio ..."
    cmake \
            -S . -B build \
            -G "$CMAKE_GENERATOR" \
            -DBUILD_SHARED_LIBS="ON" \
            -DCMAKE_C_COMPILER="$CC" \
            -DCMAKE_C_FLAGS="-w -O3" \
        | sed 's/^/        /' \
        && { echo -e "Configure cfitsio [${GREEN}done${NC}]"; echo ""; } \
        || { echo -e "Configure cfitsio [${RED}Error${NC}]"; exit 1; }

    echo -e "Compile cfitsio ..."
    cmake --build build \
        | sed 's/^/        /' \
        && { echo -e "Compile cfitsio [${GREEN}done${NC}]"; echo ""; } \
        || { echo -e "Compile cfitsio [${RED}Error${NC}]"; exit 1; }

    cd ..
    echo ""

    # CCfits
    rm -rf "CCfits/"
    tar -xf CCfits.tar.gz
    cd "CCfits"
    mkdir -p "build"

    echo -e "Configure CCfits ..."
    cmake \
            -S . -B build \
            -DCMAKE_PREFIX_PATH="../cfitsio" \
            -G "$CMAKE_GENERATOR" \
            -DBUILD_SHARED_LIBS="ON" \
            -DCMAKE_C_COMPILER="$CC" \
            -DCMAKE_CXX_COMPILER="$CXX" \
            -DCMAKE_CXX_FLAGS="-w -O3" \
        | sed 's/^/        /' \
        && { echo -e "Configure CCfits [${GREEN}done${NC}]"; echo ""; } \
        || { echo -e "Configure CCfits [${RED}Error${NC}]"; exit 1; }

    echo -e "Compile CCfits ..."
    cmake --build build \
        | sed 's/^/        /' \
        && { echo -e "Compile CCfits [${GREEN}done${NC}]"; echo ""; } \
        || { echo -e "Compile CCfits [${RED}Error${NC}]"; exit 1; }

    cd ../..
    echo ""
}


# ================================================================================ #
# =============================== install POLARIS ================================ #
# ================================================================================ #

function install_polaris()
{
    # remove old bin
    rm -f bin/polaris

    # remove and re-make build dir
    cd src/
    rm -rf build
    mkdir -p build

    # configure POLARIS with appropriate flags
    echo -e "Configure POLARIS ..."
    cmake \
            -S . -B build \
            -G "$CMAKE_GENERATOR" \
            -DBUILD_SHARED_LIBS="ON" \
            -DCMAKE_C_COMPILER="$CC" \
            -DCMAKE_CXX_COMPILER="$CXX" \
            -DCMAKE_CXX_FLAGS="$CXXFLAGS" \
        | sed 's/^/        /' \
        && { echo -e "Configure POLARIS [${GREEN}done${NC}]"; echo ""; } \
        || { echo -e "Configure POLARIS [${RED}Error${NC}]"; exit 1; }

    echo ""

    # compile and install POLARIS
    echo -e "Compile POLARIS ..."
    cmake --build build \
        | sed 's/^/        /' \
        && { echo -e "Compile POLARIS [${GREEN}done${NC}]"; echo ""; } \
        || { echo -e "Compile POLARIS [${RED}Error${NC}]"; exit 1; }

    echo ""

    echo -e "Installing ..."
    cmake --build build --target install \
        | sed 's/^/        /' \
        && { echo -e "Installing POLARIS [${GREEN}done${NC}]"; echo ""; } \
        || { echo -e "Installing POLARIS [${RED}Error${NC}]"; exit 1; }

    cd ..
    echo ""
}


# ================================================================================ #
# =========================== Get user input and verify ========================== #
# ================================================================================ #

# d -> debug mode
# u -> update (no configure)
# c -> choose c++ compiler
# g -> choose cmake generator
# f -> install fits libraries as well
# t -> test installation, eg. for the gitlab-runner
while getopts "duc:g:ft" opt; do
    case $opt in
    d)
        CO="debug"
        ;;
    u)
        cd "src"
        # compile and install POLARIS
        cmake --build build \
            | sed 's/^/        /' \
            && { echo -e "Compile POLARIS [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Compile POLARIS [${RED}Error${NC}]"; exit 1; }
        cmake --build build --target install \
            | sed 's/^/        /' \
            && { echo -e "Installing POLARIS [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Installing POLARIS [${RED}Error${NC}]"; exit 1; }
        # compile and install PolarisTools
        cd "../tools/"
        python setup.py install --user &>/dev/null \
            && { echo -e "Compile and installing PolarisTools [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Compile and installing PolarisTools [${RED}Error${NC}]"; exit 1; }
        cd ".."
        exit
        ;;
    c)
        CXX=${OPTARG}
        ;;
    g)
        if [[ ${OPTARG,,} = "ninja" ]]; then
            CMAKE_GENERATOR="Ninja"
        elif [[ ${OPTARG,,} = "make" || ${OPTARG,,} = "unix makefiles" ]]; then
            CMAKE_GENERATOR="Unix Makefiles"
        fi
        ;;
    f)
        DO_FITS=true
        ;;
    t)
        DO_FITS=true
        DO_PYTHON=false
        CO="debug"
        ;;
    esac
done


# is chosen compiler available?
type $CXX >/dev/null 2>&1 || { echo >&2 "I cannot use $CXX as compiler: it's not installed. Aborting."; exit 1; }


# is chosen cmake generator available?
if [[ $CMAKE_GENERATOR = "Ninja" ]]; then
    type "ninja" >/dev/null 2>&1 \
        || { echo >&2 "I cannot use $CMAKE_GENERATOR as CMake generator: it's not installed. Aborting."; exit 1; }
elif [[ $CMAKE_GENERATOR = "Unix Makefiles" ]]; then
    type "make" >/dev/null 2>&1 \
        || { echo >&2 "I cannot use $CMAKE_GENERATOR as CMake generator: it's not installed. Aborting."; exit 1; }
fi


# ================================================================================ #
# ============================= Set compiler options ============================= #
# ================================================================================ #

# intel
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
	   CXXFLAGS="-O1 -g3 -debug inline-debug-info"

    elif [ $CO = "fast" ]; then
        CXXFLAGS="-O3 -parallel -ip -ipo -g $CXXFLAGS"

    elif [ $CO = "fast-multi" ]; then
        CXXFLAGS="-O3 -parallel -ip -ipo -g -xAVX -axCORE-AVX512,CORE-AVX2"

    elif [ $CO = "fast-debug" ]; then
        CXXFLAGS="-O3 -parallel -ip -ipo -g3 -debug inline-debug-info $CXXFLAGS"

    elif [ $CO = "fast-debug-multi" ]; then
        CXXFLAGS="-O3 -parallel -ip -ipo -g3 -debug inline-debug-info -xAVX -axCORE-AVX512,CORE-AVX2"

    else
        echo "No valid option chosen"
        exit 1
    fi

    CXXFLAGS="$CXXFLAGS -I/usr/include/x86_64-linux-gnu/c++/8/"

# gnu
elif [ $CXX = "g++" ]; then
    if [ $CO = "debug" ]; then
	    CXXFLAGS="-O1 -g3 -Wall"
            CXXFLAGS="$CXXFLAGS -Wno-unused-function -Wno-sign-compare -Wno-unused-but-set-variable -Wno-comment -Wno-unknown-pragmas"

    elif [ $CO = "fast-debug" ]; then
	    CXXFLAGS="-march=native -O3 -funroll-loops -ffast-math -fno-finite-math-only -flto"
            CXXFLAGS="$CXXFLAGS -g3 -Wall"
            CXXFLAGS="$CXXFLAGS -Wno-unused-function -Wno-sign-compare -Wno-unused-but-set-variable -Wno-comment -Wno-unknown-pragmas"

    elif [ $CO = "fast" ]; then
        CXXFLAGS="-march=native -O3 -funroll-loops -ffast-math -fno-finite-math-only -flto"

    else
        echo "No valid option chosen"
        exit 1
    fi
    CXXFLAGS="$CXXFLAGS -fuse-linker-plugin -fuse-ld=gold"

# clang
elif [ $CXX = "clang++" ]; then
    if [ $CO = "debug" ]; then
	    CXXFLAGS="-O1 -g3 -Wall"
            CXXFLAGS="$CXXFLAGS -Wno-unused-function -Wno-sign-compare -Wno-unused-private-field -Wno-unknown-pragmas"

    elif [ $CO = "fast-debug" ]; then
	    CXXFLAGS="-march=native -O3 -ffast-math -fno-finite-math-only -flto"
            CXXFLAGS="$CXXFLAGS -g3 -Wall"
            CXXFLAGS="$CXXFLAGS -Wno-unused-function -Wno-sign-compare -Wno-unused-private-field -Wno-unknown-pragmas"

    elif [ $CO = "fast" ]; then
        CXXFLAGS="-march=native -O3 -ffast-math -fno-finite-math-only -flto"

    else
        echo "No valid option chosen"
        exit 1
    fi
    CXXFLAGS="$CXXFLAGS -Wno-comment"

    CC="clang"
else
    echo "No valid compiler chosen"
    exit 1
fi

# always add OpenMP support and C++11 standard
CXXFLAGS="$CXXFLAGS -fopenmp -std=c++11"


# ================================================================================ #
# ============================= Show set-up and do it ============================ #
# ================================================================================ #

# show the chosen options
echo ""
echo "CXX:               $CXX"
echo "CO:                $CO"
echo "CMake Generator:   $CMAKE_GENERATOR"
echo ""
echo ""
echo "Compiler flags:    $CXXFLAGS"
echo ""
echo ""


if [ "$DO_FITS" = true ]; then
    install_fits_libraries
fi


install_polaris


if [ "$DO_PYTHON" = true ]; then
    # go to PolarisTools dir
    cd tools/

    # compile and install PolarisTools
    python setup.py install --user &>/dev/null \
        && { echo -e "Compile and installing PolarisTools [${GREEN}done${NC}]"; echo ""; } \
        || { echo -e "Compile and installing PolarisTools [${RED}Error${NC}]"; exit 1; }
fi
