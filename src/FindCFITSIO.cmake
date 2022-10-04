# Find the CFITSIO library
#
# This module defines these variables:
#
#   CFITSIO_FOUND
#      True if the CFITSIO library was found.
#   CFITSIO_LIBRARY
#      The location of the CFITSIO library.
#   CFITSIO_INCLUDE_DIR
#      The include path of the CFITSIO library. 

include(FindPackageHandleStandardArgs)

#
# Find the header file
#
FIND_PATH(CFITSIO_INCLUDE_DIR fitsio.h)
#
# Find the library
#
FIND_LIBRARY(CFITSIO_LIBRARY cfitsio)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(CFITSIO DEFAULT_MSG CFITSIO_LIBRARY CFITSIO_INCLUDE_DIR)