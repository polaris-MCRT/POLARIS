# - Try to find CCFITS.
# Once executed, this module will define:
# Variables defined by this module:
#  CCFITS_FOUND        - system has CCFITS
#  CCFITS_INCLUDE_DIR  - the CCFITS include directory (cached)
#  CCFITS_INCLUDE_DIRS - the CCFITS include directories
#                         (identical to CCFITS_INCLUDE_DIR)
#  CCFITS_LIBRARY      - the CCFITS library (cached)
#  CCFITS_LIBRARIES    - the CCFITS libraries
#                         (identical to CCFITS_LIBRARY)

include(FindPackageHandleStandardArgs)

find_path(CCFITS_INCLUDE_DIR CCfits/CCfits.h)

find_library(CCFITS_LIBRARY CCfits)
        
find_package(CFITSIO REQUIRED)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(CCFITS DEFAULT_MSG CCFITS_LIBRARY CCFITS_INCLUDE_DIR)