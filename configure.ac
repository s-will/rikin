# -*- Autoconf -*-

# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([RRIKin], [0.3], [])

# special dir for aux config files 
AC_CONFIG_AUX_DIR(config)

AM_INIT_AUTOMAKE


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

AC_LANG_CPLUSPLUS

# Checks for header files.
AC_HEADER_STDC

## for compiling modified ViennaRNA code
if test "$GCC" = yes; then
   AC_DEFINE(UNUSED, __attribute__ ((unused)), avoid warning about unused variables)
else
   AC_DEFINE(UNUSED,, only for gcc)
fi
AC_CHECK_FUNCS(strdup strstr strchr erand48)


## ========================================
## DOXYGEN
#
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN(HybridPF, Doxyfile, doc)


dnl ## ========================================
dnl ## ELL (& BIU)
dnl #

dnl # minimal reqired version of the BIU library
dnl BIU_REQUIRED_VERSION=2.3.5
dnl ELL_REQUIRED_VERSION=3.2.2

dnl # BIU library path support, if not installed in usual directories
dnl AC_ARG_WITH(BIU,
dnl         AC_HELP_STRING(
dnl             [--with-BIU=PREFIX],
dnl             [alternative prefix path to BIU library]
dnl           ),
dnl           BIUPATHSET=1,
dnl           BIUPATHSET=0
dnl )
dnl if test  $BIUPATHSET = 1 ; then
dnl                 CXXFLAGS="-I$with_BIU/include $CXXFLAGS"
dnl                 LDFLAGS="-L$with_BIU/lib $LDFLAGS"
dnl                 export PKG_CONFIG_PATH="$with_BIU/lib/pkgconfig/:$PKG_CONFIG_PATH"
dnl fi

dnl # BIU library always needed
dnl LIBS="-lbiu $LIBS"

dnl # check for biu headers
dnl AC_MSG_CHECKING([for the BIU C++ library])
dnl AC_COMPILE_IFELSE(
dnl      [ #include <biu/VirtualList.hh> ],
dnl      [ 
dnl        AC_MSG_RESULT(yes)
dnl        BIUNOTFOUND=0;
dnl      ],
dnl          [
dnl        AC_MSG_RESULT(no)
dnl        BIUNOTFOUND=1;
dnl      ]
dnl )



dnl # ELL library path support, if not installed in usual directories
dnl AC_ARG_WITH(ELL,
dnl         AC_HELP_STRING(
dnl             [--with-ELL=prefix],
dnl             [try this for the ELL-library prefix install directory. Mandatory for latFold program.]
dnl           ),
dnl           ELLPATHSET=1,
dnl           ELLPATHSET=0
dnl )
dnl if test  $ELLPATHSET = 1 ; then
dnl         CXXFLAGS="$CXXFLAGS -I$with_ELL/include"
dnl         LDFLAGS="-L$with_ELL/lib $LDFLAGS"
dnl         export PKG_CONFIG_PATH="$with_ELL/lib/pkgconfig/:$PKG_CONFIG_PATH"
dnl fi


dnl # ELL library always needed
dnl LIBS="-lell $LIBS"

dnl AC_MSG_CHECKING([for the ELL C++ library])
dnl AC_COMPILE_IFELSE(
dnl      [ #include <ell/State.hh> ],
dnl      [ 
dnl        AC_MSG_RESULT(yes)
dnl        ELLNOTFOUND=0;
dnl      ],
dnl          [
dnl        AC_MSG_RESULT(no)
dnl        ELLNOTFOUND=1;
dnl      ]
dnl )


##################
## Debug option
#
AC_MSG_CHECKING([whether to build with debug information])
debuger=no
AC_ARG_ENABLE([debug],
    AC_HELP_STRING(
        [--enable-debug],
        [enable debug data generation (def=no)]
    ),
    debuger="$enableval"
)
AC_MSG_RESULT($debuger)
if test x"$debuger" = x"yes"; then
    AC_DEFINE([DEBUG], [], [Define to 1 for turning on debugging])
    CPPFLAGS="$CPPFLAGS -g -Wall"
else
    AC_DEFINE([NDEBUG], [], [Define to 1 for turning off debugging])
fi


# configure for use of LocARNA library
PKG_CHECK_MODULES([LOCARNA], [LocARNA-1.1 >= 1.7.5])

CPPFLAGS="$LOCARNA_CFLAGS $CPPFLAGS"
LDFLAGS="$LDFLAGS $LOCARNA_LIBS"
LIBS="-lLocARNA-1.1 $LIBS"


## ========================================
## Vienna RNA
#

# Get Vienna RNA library and include locations
AC_ARG_WITH([vrna],
  [AS_HELP_STRING([--with-vrna],
    [location of the Vienna RNA headers])],
  [VRNA_HOME=$withval;
   VRNA_CFLAGS="-I$withval/include"; 
   VRNA_LDFLAGS="-L$withval/lib"],
  [])

LIBS="-lRNA $LIBS"

#AC_SUBST([VRNA_CFLAGS])
#AC_SUBST([VRNA_LDFLAGS])

LDFLAGS="$VRNA_LDFLAGS $LDFLAGS"
CPPFLAGS="$VRNA_CFLAGS $CPPFLAGS"


AC_CHECK_HEADER("ViennaRNA/fold.h",
	[VRNANOTFOUND=0],
	[echo "ERROR: Cannot find Vienna RNA library (--with-vrna=\"$VRNA_HOME\")",
	VRNANOTFOUND=0],
	[])

# ------------------------------
# handle open mp support
#  this is necessary when linking against the vienna lib
AC_OPENMP
AC_SUBST(OPENMP_CXXFLAGS)

CPPFLAGS="$OPENMP_CXXFLAGS $CPPFLAGS"
LDFLAGS="$OPENMP_CXXFLAGS $LDFLAGS"


####################
# Static linking
#
AC_MSG_CHECKING([whether to link statically])
static_linkage=no
AC_ARG_ENABLE([static],
    AC_HELP_STRING(
        [--enable-static],
        [link statically (def=no)]
    ),
    static_linkage="$enableval"
)
AC_MSG_RESULT($static_linkage)
if test x"$static_linkage" = x"yes"; then
   LDFLAGS="$LDFLAGS -static"
fi


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T

###############################################################################
# BEGIN DO NOT COMPILE CHECKS
###############################################################################

DEPENDENCYNOTFOUND=0;

# error output if ViennaRNA not found
if test "$VRNANOTFOUND" = "1"; then
        AC_MSG_NOTICE()
        AC_MSG_NOTICE(the Vienna RNA C library is required for building the ELL.)
        AC_MSG_NOTICE(it can be obtained from http://www.tbi.univie.ac.at/~ivo/RNA/.)
        AC_MSG_NOTICE()
        DEPENDENCYNOTFOUND=1;
fi

dnl # error output if BIU-headers not found
dnl if test "$BIUNOTFOUND" = "1"; then
dnl         AC_MSG_NOTICE()
dnl         AC_MSG_NOTICE(the BIU C++ library is required for building the ELL.)
dnl         AC_MSG_NOTICE(it can be obtained from http://www.bioinf.uni-freiburg.de/SW/BIU/.)
dnl         AC_MSG_NOTICE()
dnl         if test "$BIUPATHSET" = "1"; then
dnl                 AC_MSG_NOTICE(cant find BIU in given path '$with_BIU'.)
dnl         else
dnl                 AC_MSG_NOTICE(if BIU is installed elsewhere use --with-BIU=PREFIX.)
dnl         fi
dnl         DEPENDENCYNOTFOUND=1;
dnl fi

dnl # error output if ELL-headers not found
dnl if test $ELLNOTFOUND = 1; then
dnl         AC_MSG_NOTICE()
dnl         AC_MSG_NOTICE(the ELL C++ library is required for building the tools.)
dnl         AC_MSG_NOTICE(it can be obtained from http://www.bioinf.uni-freiburg.de/SW/ELL/.)
dnl         AC_MSG_NOTICE()
dnl         if test $ELLPATHSET = 1; then
dnl                 AC_MSG_NOTICE(cant find ELL in given path '$with_ELL'.)
dnl         else
dnl                 AC_MSG_NOTICE(if ELL is installed elsewhere use --with-ELL=prefix.)
dnl         fi
dnl         DEPENDENCYNOTFOUND=1;
dnl fi

# error ABORT if on of the libraries was not found
if test "$DEPENDENCYNOTFOUND" = "1"; then
        AC_MSG_NOTICE()
        AC_MSG_ERROR(library not found! Try './configure --help'.)
fi


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])

AC_OUTPUT([src/Makefile])

