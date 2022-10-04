#!/bin/bash

# Define colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m'

# get current directory
current_path=$(pwd)

# set pipefail so "$?" gives 1 for piping of "false | true" for instance
# necessary because of piping to sed
set -o pipefail

# check version of cmake
# from 3.13 on, cmake can use -B and -S to specify build and source dirs
currentver="$(cmake --version | head -n 1 | awk -F "version " {'print $NF'})"
requiredver="3.13"
 if [ "$(printf '%s\n' "$requiredver" "$currentver" | sort -V | head -n1)" = "$requiredver" ]; then
        old_cmake=0
 else
        old_cmake=1
 fi


# ================================================================================ #
# ================================ Default values ================================ #
# ================================================================================ #


# default compiler
CC="gcc"
CXX="g++"

# default cmake generator
CMAKE_GENERATOR="Unix Makefiles"

# default compiler flag
CO="release"

# install target
TARGET="host"

# do not install fits libraries each time
DO_FITS=false


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
    if [ "$old_cmake" -eq "1" ]; then
        cd "build"
        cmake .. \
                -G "$CMAKE_GENERATOR" \
                -DBUILD_SHARED_LIBS="ON" \
                -DCMAKE_INSTALL_PREFIX="${current_path}/lib" \
                -DCMAKE_C_COMPILER="$CC" \
                -DCMAKE_C_FLAGS="-w -O3" \
                -DCMAKE_CXX_COMPILER="$CXX" \
                -DCMAKE_CXX_FLAGS="-w -O3" \
            | sed 's/^/        /' \
            && { echo -e "Configure cfitsio [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Configure cfitsio [${RED}Error${NC}]"; exit 1; }
    else
        cmake \
                -S . -B build \
                -G "$CMAKE_GENERATOR" \
                -DBUILD_SHARED_LIBS="ON" \
                -DCMAKE_INSTALL_PREFIX="${current_path}/lib" \
                -DCMAKE_C_COMPILER="$CC" \
                -DCMAKE_C_FLAGS="-w -O3" \
                -DCMAKE_CXX_COMPILER="$CXX" \
                -DCMAKE_CXX_FLAGS="-w -O3" \
            | sed 's/^/        /' \
            && { echo -e "Configure cfitsio [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Configure cfitsio [${RED}Error${NC}]"; exit 1; }
    fi

    echo -e "Compile cfitsio ..."
    if [ "$old_cmake" -eq "1" ]; then
        cmake --build . && cmake --build . --target install \
            | sed 's/^/        /' \
            && { echo -e "Compile cfitsio [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Compile cfitsio [${RED}Error${NC}]"; exit 1; }
        cd ..
    else
        cmake --build build && cmake --build build --target install\
            | sed 's/^/        /' \
            && { echo -e "Compile cfitsio [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Compile cfitsio [${RED}Error${NC}]"; exit 1; }
    fi

    cd ..
    echo ""

    # CCfits
    rm -rf "CCfits/"
    tar -xf CCfits.tar.gz
    cd "CCfits"
    mkdir -p "build"

    echo -e "Configure CCfits ..."
    if [ "$old_cmake" -eq "1" ]; then
        cd "build"
        cmake .. \
                -DCMAKE_PREFIX_PATH="../cfitsio" \
                -G "$CMAKE_GENERATOR" \
                -DCMAKE_PREFIX_PATH="${current_path}/lib" \
                -DCMAKE_INSTALL_PREFIX="${current_path}/lib" \
                -DBUILD_SHARED_LIBS="ON" \
                -DCMAKE_C_COMPILER="$CC" \
                -DCMAKE_C_FLAGS="-w -O3" \
                -DCMAKE_CXX_COMPILER="$CXX" \
                -DCMAKE_CXX_FLAGS="-w -O3" \
            | sed 's/^/        /' \
            && { echo -e "Configure CCfits [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Configure CCfits [${RED}Error${NC}]"; exit 1; }
    else
        cmake \
                -S . -B build \
                -DCMAKE_PREFIX_PATH="../cfitsio" \
                -G "$CMAKE_GENERATOR" \
                -DCMAKE_PREFIX_PATH="${current_path}/lib" \
                -DCMAKE_INSTALL_PREFIX="${current_path}/lib" \
                -DBUILD_SHARED_LIBS="ON" \
                -DCMAKE_C_COMPILER="$CC" \
                -DCMAKE_C_FLAGS="-w -O3" \
                -DCMAKE_CXX_COMPILER="$CXX" \
                -DCMAKE_CXX_FLAGS="-w -O3" \
            | sed 's/^/        /' \
            && { echo -e "Configure CCfits [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Configure CCfits [${RED}Error${NC}]"; exit 1; }
    fi

    echo -e "Compile CCfits ..."
    if [ "$old_cmake" -eq "1" ]; then
        cmake --build . && cmake --build . --target install \
            | sed 's/^/        /' \
            && { echo -e "Compile CCfits [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Compile CCfits [${RED}Error${NC}]"; exit 1; }
        cd ".."
    else
        cmake --build build && cmake --build build --target install \
            | sed 's/^/        /' \
            && { echo -e "Compile CCfits [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Compile CCfits [${RED}Error${NC}]"; exit 1; }
    fi

    cd ../..
    echo ""
}


# ================================================================================ #
# =============================== install POLARIS ================================ #
# ================================================================================ #


function install_polaris()
{
    # remove old bins
    rm -f bin/polaris
    rm -f src/polaris

    # remove and re-make build dir
    cd src/
    rm -rf build
    mkdir -p build

    # configure POLARIS with appropriate flags
    echo -e "Configure POLARIS ..."
    if [ "$old_cmake" -eq "1" ]; then
        cd "build"
        cmake .. \
                -G "$CMAKE_GENERATOR" \
                -DCMAKE_PREFIX_PATH="${current_path}/lib" \
                -DCMAKE_INSTALL_PREFIX="${current_path}" \
                -DBUILD_SHARED_LIBS="ON" \
                -DCMAKE_C_COMPILER="$CC" \
                -DCMAKE_CXX_COMPILER="$CXX" \
                -DCMAKE_CXX_FLAGS="$CXXFLAGS" \
            | sed 's/^/        /' \
            && { echo -e "Configure POLARIS [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Configure POLARIS [${RED}Error${NC}]"; exit 1; }
    else
        cmake \
                -S . -B build \
                -G "$CMAKE_GENERATOR" \
                -DCMAKE_PREFIX_PATH="${current_path}/lib" \
                -DCMAKE_INSTALL_PREFIX="${current_path}" \
                -DBUILD_SHARED_LIBS="ON" \
                -DCMAKE_C_COMPILER="$CC" \
                -DCMAKE_CXX_COMPILER="$CXX" \
                -DCMAKE_CXX_FLAGS="$CXXFLAGS" \
            | sed 's/^/        /' \
            && { echo -e "Configure POLARIS [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Configure POLARIS [${RED}Error${NC}]"; exit 1; }
    fi

    echo ""

    # compile and install POLARIS
    echo -e "Compile POLARIS ..."
    if [ "$old_cmake" -eq "1" ]; then
        cmake --build . \
            | sed 's/^/        /' \
            && { echo -e "Compile POLARIS [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Compile POLARIS [${RED}Error${NC}]"; exit 1; }
    else
        cmake --build build \
            | sed 's/^/        /' \
            && { echo -e "Compile POLARIS [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Compile POLARIS [${RED}Error${NC}]"; exit 1; }
    fi

    echo ""

    echo -e "Installing ..."
    if [ "$old_cmake" -eq "1" ]; then
        cmake --build . --target install \
            | sed 's/^/        /' \
            && { echo -e "Installing POLARIS [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Installing POLARIS [${RED}Error${NC}]"; exit 1; }
        cd ..
    else
        cmake --build build --target install \
            | sed 's/^/        /' \
            && { echo -e "Installing POLARIS [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Installing POLARIS [${RED}Error${NC}]"; exit 1; }
    fi

    cd ..
    echo ""
}


# ================================================================================ #
# ================================ Update POLARIS ================================ #
# ================================================================================ #


function update_installation()
{
    cd "src"
    # compile and install POLARIS
    if [ "$old_cmake" -eq "1" ]; then
        cd "build"
        cmake --build . \
            | sed 's/^/        /' \
            && { echo -e "Compile POLARIS [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Compile POLARIS [${RED}Error${NC}]"; exit 1; }
        cmake --build . --target install \
            | sed 's/^/        /' \
            && { echo -e "Installing POLARIS [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Installing POLARIS [${RED}Error${NC}]"; exit 1; }
        cd ..
    else
        cmake --build build \
            | sed 's/^/        /' \
            && { echo -e "Compile POLARIS [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Compile POLARIS [${RED}Error${NC}]"; exit 1; }
        cmake --build build --target install \
            | sed 's/^/        /' \
            && { echo -e "Installing POLARIS [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Installing POLARIS [${RED}Error${NC}]"; exit 1; }
    fi
    cd ".."

    # compile and install PolarisTools if available
    # PolarisTools is no longer maintained and is not shipped with POLARIS anymore
    if [[ -d "tools" ]]; then
        cd "tools/"
        python3 setup.py install --user &>/dev/null \
            && { echo -e "Compile and installing PolarisTools [${GREEN}done${NC}]"; echo ""; } \
            || { echo -e "Compile and installing PolarisTools [${RED}Error${NC}]"; exit 1; }
        cd ".."
    fi

    source ${HOME}/.bashrc

    exit
}


# ================================================================================ #
# ================================ Delete POLARIS ================================ #
# ================================================================================ #


function delete_installation()
{
    printf "%s\n" "Do you really want to delete your POLARIS installation [y/N]?"
    read really_delete
    case ${really_delete:=n} in
    [yY]*)
        echo "------ delete POLARIS ------"

        export_str="export POLARIS_PATH=\"${current_path}\""
        if grep -q "${export_str}" ${HOME}/.bashrc; then
            sed -i.bak "/${export_str//\//\\/}/,+4d" ${HOME}/.bashrc
        fi

        export_str="export POLARIS_FITS_PATH=\""'${POLARIS_PATH}'"/lib/CCfits/build:"'${POLARIS_PATH}'"/lib/cfitsio/build\""
        if grep -q "${export_str}" ${HOME}/.bashrc; then
            sed -i.bak "/${export_str//\//\\/}/,+4d" ${HOME}/.bashrc
        fi

        export_str="export POLARISTOOLS_PATH=\"$HOME/.local/bin\""
        if grep -q "${export_str}" ${HOME}/.bashrc; then
            sed -i.bak "/${export_str//\//\\/}/,+4d" ${HOME}/.bashrc
        fi

        cd ${current_path}/../
        rm -rfv ${current_path}
        pip3 uninstall PolarisTools
        exit
        ;;
    *)
        exit
        ;;
    esac
}


# ================================================================================ #
# =========================== Get user input and verify ========================== #
# ================================================================================ #


function usage() {
    echo ""
    echo "usage: compile.sh [-h] [-frdu] [-c CXX_COMPILER] [-g CMAKE_GENERATOR] [-D]"
    echo ""
    echo "Install and compile POLARIS"
    echo -e "${YELLOW}HINT:${NC} For first installation, use option -f"
    echo ""
    echo "optional arguments:"
    echo "-h      show this help message and exit"
    echo "-f      first installation (compile POLARIS and install the cfitsio and CCfits libraries)"
    echo "-r      clear and compile POLARIS (release mode) (default)"
    echo "-d      clear and compile POLARIS (debug mode)"
    echo "-u      re-compile POLARIS if necessary with last configuration (update)"
    echo "-c CXX_compiler"
    echo "        choose the c++ compiler you want to use:"
    echo "          - gcc (default)"
    echo "          - icc"
    echo "          - clang++"
    echo "-g CMAKE_GENERATOR"
    echo "        choose the generator for cmake:"
    echo "          - make (default)"
    echo "          - ninja"
    echo "-D delete POLARIS from your computer"
    echo ""
    exit
}

release=false
debug=false

while getopts "hfrduc:g:D" opt; do
    case $opt in
    h)
        usage
        ;;
    r)
        # this is the default -> do nothing
        # except if -d is also chosen
        # release + debug -> fast debug (nice for profiling)
        release=true
        ;;
    d)
        CO="debug"
        # release + debug -> fast debug (nice for profiling)
        debug=true
        ;;
    u)
        update_installation
        ;;
    f)
        DO_FITS=true
        ;;
    c)
        if [[ ${OPTARG,,} = "gcc" ]]; then
            CXX="g++"
        else
            CXX=${OPTARG}
        fi
        ;;
    g)
        if [[ ${OPTARG,,} = "ninja" ]]; then
            CMAKE_GENERATOR="Ninja"
        elif [[ ${OPTARG,,} = "make" || ${OPTARG,,} = "unix makefiles" ]]; then
            CMAKE_GENERATOR="Unix Makefiles"
        fi
        ;;
    D)
        delete_installation
        ;;
    *)
        usage
        ;;
    esac
done

# release + debug -> fast debug (nice for profiling)
if $release && $debug; then CO="fast-debug"; fi


# Print hint for first installation
echo ""
echo -e "${YELLOW}HINT:${NC} For first installation, use option -f"
echo ""


# Search for required packages
required_packages=true
echo -e "Search for required packages ..."
if type $CXX >/dev/null 2>&1; then
    echo -e "        -- Required package $CXX [${GREEN}found${NC}]"
else
    echo -e "        -- ${RED}Error:${NC} Required package $CXX not found!"
    required_packages=false
fi

if [[ $CMAKE_GENERATOR = "Ninja" ]]; then
    if type "ninja" >/dev/null 2>&1; then
        echo -e "        -- Required package $CMAKE_GENERATOR [${GREEN}found${NC}]"
    else
        echo -e "        -- ${RED}Error:${NC} Required package $CMAKE_GENERATOR not found!"
        required_packages=false
    fi
elif [[ $CMAKE_GENERATOR = "Unix Makefiles" ]]; then
    if type "make" >/dev/null 2>&1; then
        echo -e "        -- Required package $CMAKE_GENERATOR [${GREEN}found${NC}]"
    else
        echo -e "        -- ${RED}Error:${NC} Required package $CMAKE_GENERATOR not found!"
        required_packages=false
    fi
fi


if type "python3" >/dev/null 2>&1; then
    echo -e "        -- Required package python3 [${GREEN}found${NC}]"
    for package_name in argparse io numpy os setuptools shutil struct sys; do
        if python3 -c "import ${package_name}" &>/dev/null; then
            echo -e "        -- Required python package ${package_name} [${GREEN}found${NC}]"
        else
            if type "pip3" >/dev/null 2>&1; then
                echo "Pip installation detected. Install ${package_name}!"
                pip3 install ${package_name}
                if python3 -c "import ${package_name}" &>/dev/null; then
                    echo -e "        -- Installation succesfull. Required python package ${package_name} [${GREEN}found${NC}]"
                else
                    echo -e "        -- ${RED}Error:${NC} Installation of ${package_name} not succesfull!"
                    required_packages=false
                fi
            else
                echo -e "        -- ${RED}Error:${NC} Required python package ${package_name} not found!"
                required_packages=false
            fi
        fi
    done
else
    echo -e "        -- ${RED}Error:${NC} Required package python3 not found!"
    required_packages=false
fi

if ! $required_packages; then
    echo -e "${RED}Error:${NC} Installation aborted. Please contact your system admin to install the required packages."
    exit 1
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

    elif [ $CO = "release" ]; then
        CXXFLAGS="-O3 -parallel -ip -ipo -g $CXXFLAGS"

    elif [ $CO = "release-multi" ]; then
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
            CXXFLAGS="$CXXFLAGS -Wno-unused-function -Wno-sign-compare -Wno-unused-but-set-variable -Wno-comment -Wno-unknown-pragmas -Wno-maybe-uninitialized"

    elif [ $CO = "fast-debug" ]; then
	    CXXFLAGS="-march=native -O3 -funroll-loops -flto"
            CXXFLAGS="$CXXFLAGS -g3 -Wall"
            CXXFLAGS="$CXXFLAGS -Wno-unused-function -Wno-sign-compare -Wno-unused-but-set-variable -Wno-comment -Wno-unknown-pragmas"

    elif [ $CO = "release" ]; then
        CXXFLAGS="-march=native -O3 -funroll-loops -flto -Wno-deprecated"

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
	    CXXFLAGS="-march=native -O3 -flto"
            CXXFLAGS="$CXXFLAGS -g3 -Wall"
            CXXFLAGS="$CXXFLAGS -Wno-unused-function -Wno-sign-compare -Wno-unused-private-field -Wno-unknown-pragmas"

    elif [ $CO = "release" ]; then
        CXXFLAGS="-march=native -O3 -flto"

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


if [ ! -d "projects" ]; then
    mkdir "projects"
fi


# set POLARIS_PATH as environment variable
set_envvar_POLARIS_str="export POLARIS_PATH=\"${current_path}\""
# check if this variable is already in PATH to prevent multiple entries
if_POLARIS_str="if [[ \":"'$PATH'":\" != *\":"'${POLARIS_PATH}'"/bin:\"* ]]; then"
# prepend this variable to PATH
export_POLARIS_str="    export PATH=\""'${POLARIS_PATH}'"/bin:"'$PATH'"\""
# check if these lines are already in .bashrc
if ! grep -q "${set_envvar_POLARIS_str}" ${HOME}/.bashrc; then
    echo -e "Add POLARIS install path to ~/.bashrc ..."
    echo "${set_envvar_POLARIS_str}" >>${HOME}/.bashrc
    echo "${if_POLARIS_str}" >>${HOME}/.bashrc
    echo "${export_POLARIS_str}" >>${HOME}/.bashrc
    echo "fi" >>${HOME}/.bashrc
    echo "" >>${HOME}/.bashrc

    source ${HOME}/.bashrc
    echo -e "... [${GREEN}done${NC}]"
    echo ""
    echo ""
fi

# set POLARIS_FITS_PATH as environment variable
set_envvar_FITS_str="export POLARIS_FITS_PATH=\""'${POLARIS_PATH}'"/lib/lib\""
# check if this variable is already in LD_LIBRARY_PATH to prevent multiple entries
if_FITS_str="if [[ \":"'$LD_LIBRARY_PATH'":\" != *\":"'$POLARIS_FITS_PATH'":\"* ]]; then"
# prepend this variable to LD_LIBRARY_PATH
export_FITS_str="    export LD_LIBRARY_PATH=\""'$POLARIS_FITS_PATH'":"'$LD_LIBRARY_PATH'"\""
if ! grep -q "${set_envvar_FITS_str}" ${HOME}/.bashrc; then
    echo -e "Add FITS libraries install paths to ~/.bashrc ..."
    echo "${set_envvar_FITS_str}" >>${HOME}/.bashrc
    echo "${if_FITS_str}" >>${HOME}/.bashrc
    echo "${export_FITS_str}" >>${HOME}/.bashrc
    echo "fi" >>${HOME}/.bashrc
    echo "" >>${HOME}/.bashrc

    source ${HOME}/.bashrc
    echo -e "... [${GREEN}done${NC}]"
    echo ""
    echo ""
fi


if [ "$DO_FITS" = true ]; then
    install_fits_libraries
fi


install_polaris


# compile and install PolarisTools if available
if [[ -d "tools" ]]; then
    # go to PolarisTools dir
    cd "tools/"

    # compile and install PolarisTools
    python3 setup.py install --user &>/dev/null \
        && { echo -e "Compile and installing PolarisTools [${GREEN}done${NC}]"; echo ""; } \
        || { echo -e "Compile and installing PolarisTools [${RED}Error${NC}]"; exit 1; }

    # set POLARISTOOLS_PATH as environment variable
    set_envvar_TOOLS_str="export POLARISTOOLS_PATH=\"$HOME/.local/bin\""
    # check if this variable is already in PATH to prevent multiple entries
    if_TOOLS_str="if [[ \":"'$PATH'":\" != *\":"'$POLARISTOOLS_PATH'":\"* ]]; then"
    # prepend this variable to PATH
    export_TOOLS_str="    export PATH=\""'$POLARISTOOLS_PATH'":"'$PATH'"\""
    if ! grep -q "${set_envvar_TOOLS_str}" ${HOME}/.bashrc; then
        echo -e "Add PolarisTools install path to ~/.bashrc ..."
        echo "${set_envvar_TOOLS_str}" >>${HOME}/.bashrc
        echo "${if_TOOLS_str}" >>${HOME}/.bashrc
        echo "${export_TOOLS_str}" >>${HOME}/.bashrc
        echo "fi" >>${HOME}/.bashrc
        echo "" >>${HOME}/.bashrc

        source ${HOME}/.bashrc
        echo -e "... [${GREEN}done${NC}]"
        echo ""
        echo ""
    fi
fi

source ${HOME}/.bashrc
