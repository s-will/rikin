# -*- Autoconf -*-

# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([RRIKin], [0.9], [])

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

########################################
## turn on C++11
CXXFLAGS="-std=c++11 $CXXFLAGS"

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
    CXXFLAGS="$CXXFLAGS -g -Wall"
else
    AC_DEFINE([NDEBUG], [], [Define to 1 for turning off debugging])
fi


# configure for use of LocARNA library
PKG_CHECK_MODULES([LOCARNA], [LocARNA-1.4 >= 1.8.7])

CXXFLAGS="$LOCARNA_CFLAGS $CXXFLAGS"
LDFLAGS="$LDFLAGS $LOCARNA_LIBS"
LIBS="-lLocARNA-1.4 $LIBS"

## ========================================
## Vienna RNA
#

AC_ARG_WITH([vrna],
  [AS_HELP_STRING([--with-vrna],
    [location of the Vienna RNA headers])],
  [],
  [])

dnl if with_vrna is set, then use it to set the path for pkg_config
if test "$with_vrna" != "" ; then
    export PKG_CONFIG_PATH="$with_vrna/lib/pkgconfig:$PKG_CONFIG_PATH"
    if test ! -e "$with_vrna/lib/pkgconfig/RNAlib2.pc" ; then
        AC_MSG_NOTICE([Pkg-config of Vienna RNA cannot be found.])
    fi
fi
        
PKG_CHECK_MODULES([VRNA],[RNAlib2 >= 2.1.9],[VRNALIB_OK=true],[VRNALIB_OK=false])
                                  
dnl replace ViennaRNA add end of include path -- use perl, since we require it anyway
dnl This fixes a legacy problem in ViennaRNA <2.2.x pkg-config files
VRNA_CFLAGS=$(echo $VRNA_CFLAGS | perl -ne 's/(-I.*)\/ViennaRNA(\/\s|\s|$)/\1\2/g;print')


LIBS="$VRNA_LIBS $LIBS"
CFLAGS="$VRNA_CFLAGS $CFLAGS"
CXXFLAGS="$VRNA_CFLAGS $CXXFLAGS"

## ========================================
## zlib
AC_CHECK_LIB([z],[main],[HAVE_ZLIB="yes"])
LIBS="-lz $LIBS"

# ------------------------------
# handle open mp support
#  this is necessary when linking against the vienna lib
AC_OPENMP
AC_SUBST(OPENMP_CXXFLAGS)

CXXFLAGS="$OPENMP_CXXFLAGS $CXXFLAGS"
# LDFLAGS="$OPENMP_CXXFLAGS $LDFLAGS"


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
        AC_MSG_NOTICE([The Vienna RNA library is missing.])
        AC_MSG_NOTICE()
        DEPENDENCYNOTFOUND=1;
fi

if test "$HAVE_ZLIB" != "yes" ; then
        DEPENDENCYNOTFOUND=1;   
        AC_MSG_NOTICE()
        AC_MSG_NOTICE([Zlib is missing.])
fi

# error ABORT if on of the libraries was not found
if test "$DEPENDENCYNOTFOUND" = "1"; then
        AC_MSG_NOTICE()
        AC_MSG_ERROR([library not found! Try './configure --help'.])
fi


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile src/Tests/Makefile])

AC_OUTPUT

