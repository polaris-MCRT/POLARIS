#!/bin/bash

# Define colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'


# default compiler
CC="gcc"
CXX="g++"

# default cmake generator
CMAKE_GENERATOR="Ninja"

# default compiler flag
CO="fast"


# read user input
while getopts "duc:g:" opt; do
    case $opt in
    d)
        CO="debug"
        ;;
    u)
        cd "src"
        # compile and install POLARIS
        cmake --build build \
            | sed 's/^/        /' \
            && { echo -e "Compiling POLARIS [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Compiling POLARIS [${RED}Error${NC}]"; exit; }
        cmake --build build --target install \
            | sed 's/^/        /' \
            && { echo -e "Installing POLARIS [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Installing POLARIS [${RED}Error${NC}]"; exit; }
        # compile and install PolarisTools
        cd "../tools/"
        python setup.py install --user &>/dev/null \
            && { echo -e "Compiling and installing PolarisTools [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Compiling and installing PolarisTools [${RED}Error${NC}]"; exit; }
        cd ".."
        exit
        ;;
    c)
        CXX=${OPTARG}
        ;;
    g)
        if [[ ${OPTARG,,} = "ninja" ]]; then
            CMAKE_GENERATOR="Ninja"
        elif [[ ${OPTARG,,} = "make" ]]; then
            CMAKE_GENERATOR="Unix Makefiles"
        elif [[ ${OPTARG,,} = "unix makefiles" ]]; then
            CMAKE_GENERATOR="Unix Makefiles"
        fi
        ;;
    esac
done

# is chosen compiler available?
type $CXX >/dev/null 2>&1 || { echo >&2 "I cannot use $CXX as compiler: it's not installed. Aborting."; exit 1; }


# is chosen cmake generator available?
if [ $CMAKE_GENERATOR = "Ninja" ]; then
    type "ninja" >/dev/null 2>&1 \
        || { echo >&2 "I cannot use $CMAKE_GENERATOR as CMake generator: it's not installed. Aborting."; exit 1; }
elif [ $CMAKE_GENERATOR = "Unix Makefiles" ]; then
    type "make" >/dev/null 2>&1 \
        || { echo >&2 "I cannot use $CMAKE_GENERATOR as CMake generator: it's not installed. Aborting."; exit 1; }
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

elif [ $CXX = "g++" ]; then
    if [ $CO = "debug" ]; then
	    CXXFLAGS="-O1 -g3 -Wall \
                    -Wno-unused-function -Wno-sign-compare -Wno-unused-but-set-variable -Wno-comment -Wno-unknown-pragmas"

    elif [ $CO = "fast-debug" ]; then
	    CXXFLAGS="-march=native -O3 -funroll-loops -ffast-math -fno-finite-math-only -flto \
                    -g3 -Wall \
                    -Wno-unused-function -Wno-sign-compare -Wno-unused-but-set-variable -Wno-comment -Wno-unknown-pragmas"

    elif [ $CO = "fast" ]; then
        CXXFLAGS="-march=native -O3 -funroll-loops -ffast-math -fno-finite-math-only -flto"

    else
        echo "No valid option chosen"
        exit 1
    fi
    CXXFLAGS="$CXXFLAGS -fuse-linker-plugin -fuse-ld=gold"

elif [ $CXX = "clang++" ]; then
    if [ $CO = "debug" ]; then
	    CXXFLAGS="-O1 -g3 -Wall \
                    -Wno-unused-function -Wno-sign-compare -Wno-unused-private-field -Wno-unknown-pragmas"

    elif [ $CO = "fast-debug" ]; then
	    CXXFLAGS="-march=native -O3 -ffast-math -fno-finite-math-only -flto \
                    -g3 -Wall \
                    -Wno-unused-function -Wno-sign-compare -Wno-unused-private-field -Wno-unknown-pragmas"

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

CXXFLAGS="$CXXFLAGS -fopenmp -std=c++11"

echo ""
echo "CXX:               $CXX"
echo "CO:                $CO"
echo "CMake Generator:   $CMAKE_GENERATOR"
echo ""
echo ""
echo "Compiler flags:    $CXXFLAGS"
echo ""
echo ""


# remove old bin
rm -f bin/polaris

# remove build dir
cd src/
rm -rf build
mkdir -p build


#        -DCMAKE_CXX_INCLUDE_WHAT_YOU_USE="/home/rbrunngraeber/iwyu/build/bin/include-what-you-use;-Xiwyu;any;-Xiwyu;iwyu;-Xiwyu;args" \
# configure POLARIS with appropriate flags
echo -e "Configure ..."
cmake \
        -S . -B build \
        -G "$CMAKE_GENERATOR" \
        -DBUILD_SHARED_LIBS=ON \
        -DCMAKE_C_COMPILER="$CC" \
        -DCMAKE_CXX_COMPILER="$CXX" \
        -DCMAKE_CXX_FLAGS="$CXXFLAGS" \
    | sed 's/^/        /' \
    && { echo -e "Configuring POLARIS [${GREEN}done${NC}]"; echo ""; } \
    || { echo -e "Configuring POLARIS [${RED}Error${NC}]"; exit; }

echo ""

# compile and install POLARIS
echo -e "Compiling ..."
cmake --build build \
    | sed 's/^/        /' \
    && { echo -e "Compiling POLARIS [${GREEN}done${NC}]"; echo ""; } \
    || { echo -e "Compiling POLARIS [${RED}Error${NC}]"; exit; }

echo ""

echo -e "Installing ..."
cmake --build build --target install \
    | sed 's/^/        /' \
    && { echo -e "Installing POLARIS [${GREEN}done${NC}]"; echo ""; } \
    || { echo -e "Installing POLARIS [${RED}Error${NC}]"; exit; }

echo ""

# go to PolarisTools dir
cd ../tools/

# compile and install PolarisTools
python setup.py install --user &>/dev/null \
    && { echo -e "Compiling and installing PolarisTools [${GREEN}done${NC}]"; echo ""; } \
    || { echo -e "Compiling and installing PolarisTools [${RED}Error${NC}]"; exit; }
