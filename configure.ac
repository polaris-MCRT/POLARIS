AC_INIT([polaris], [4.04])
AM_INIT_AUTOMAKE

AC_LANG(C++)
AC_PROG_CXX

AC_CONFIG_FILES([Makefile] [src/Makefile])

CPPFLAGS="$CPPFLAGS -std=c++98"

AC_ARG_ENABLE([fits],
    AS_HELP_STRING([--enable-fits], [Enable fits support]))

AS_IF([test "x$enable_fits" = "xyes"], [
    LDFLAGS="$LDFLAGS -L${ac_pwd}/cfitsio -L${ac_pwd}/CCfits/.libs -lcfitsio -lCCfits"
    CPPFLAGS="$CPPFLAGS -I${ac_pwd}/cfitsio -I${ac_pwd}/CCfits -I${ac_pwd}"

    AC_CHECK_HEADERS([CCfits.h], [],
       [AC_MSG_ERROR([Unable to find the CCfits header files])]
       )

    AC_CHECK_LIB([CCfits], [main], [],
       [AC_MSG_ERROR([Unable to find the CCfits library])]
       )

    AC_CHECK_LIB([cfitsio], [main], [],
       [AC_MSG_ERROR([Unable to find the cfitsio library])]
       )
])

AC_CONFIG_COMMANDS([mkdir], [$MKDIR_P projects])

AC_OUTPUT
