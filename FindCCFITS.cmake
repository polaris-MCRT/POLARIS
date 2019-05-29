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

if(NOT CCFITS_FOUND)

    find_path(CCFITS_INCLUDE_DIR 
        NAMES CCfits.h)
        
    find_library(CCFITS_LIBRARY 
        NAMES libCCfits.a libCCfits.so libCCfits.la libCCfits.dylib)
        
    find_package(CFITSIO REQUIRED)

endif(NOT CCFITS_FOUND)

