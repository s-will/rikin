# -*- Autoconf -*-

# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([RRIKin], [0.2], [])

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


## ========================================
## ELL (& BIU)
#

# minimal reqired version of the BIU library
BIU_REQUIRED_VERSION=2.3.5
ELL_REQUIRED_VERSION=3.2.2

# BIU library path support, if not installed in usual directories
AC_ARG_WITH(BIU,
        AC_HELP_STRING(
            [--with-BIU=PREFIX],
            [alternative prefix path to BIU library]
          ),
          BIUPATHSET=1,
          BIUPATHSET=0
)
if test  $BIUPATHSET = 1 ; then
                CXXFLAGS="-I$with_BIU/include $CXXFLAGS"
                LDFLAGS="-L$with_BIU/lib $LDFLAGS"
                export PKG_CONFIG_PATH="$with_BIU/lib/pkgconfig/:$PKG_CONFIG_PATH"
fi

# BIU library always needed
LIBS="-lbiu $LIBS"

# check for biu headers
AC_MSG_CHECKING([for the BIU C++ library])
AC_COMPILE_IFELSE(
     [ #include <biu/VirtualList.hh> ],
     [ 
       AC_MSG_RESULT(yes)
       BIUNOTFOUND=0;
     ],
         [
       AC_MSG_RESULT(no)
       BIUNOTFOUND=1;
     ]
)



# ELL library path support, if not installed in usual directories
AC_ARG_WITH(ELL,
        AC_HELP_STRING(
            [--with-ELL=prefix],
            [try this for the ELL-library prefix install directory. Mandatory for latFold program.]
          ),
          ELLPATHSET=1,
          ELLPATHSET=0
)
if test  $ELLPATHSET = 1 ; then
        CXXFLAGS="$CXXFLAGS -I$with_ELL/include"
        LDFLAGS="-L$with_ELL/lib $LDFLAGS"
        export PKG_CONFIG_PATH="$with_ELL/lib/pkgconfig/:$PKG_CONFIG_PATH"
fi


# ELL library always needed
LIBS="-lell $LIBS"

AC_MSG_CHECKING([for the ELL C++ library])
AC_COMPILE_IFELSE(
     [ #include <ell/State.hh> ],
     [ 
       AC_MSG_RESULT(yes)
       ELLNOTFOUND=0;
     ],
         [
       AC_MSG_RESULT(no)
       ELLNOTFOUND=1;
     ]
)


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
# PKG_CHECK_MODULES([DEPS], [LocARNA-1.0 >= 1.6.1])

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

LIBS="-lg2 -lRNA $LIBS"

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

# error output if BIU-headers not found
if test "$BIUNOTFOUND" = "1"; then
        AC_MSG_NOTICE()
        AC_MSG_NOTICE(the BIU C++ library is required for building the ELL.)
        AC_MSG_NOTICE(it can be obtained from http://www.bioinf.uni-freiburg.de/SW/BIU/.)
        AC_MSG_NOTICE()
        if test "$BIUPATHSET" = "1"; then
                AC_MSG_NOTICE(cant find BIU in given path '$with_BIU'.)
        else
                AC_MSG_NOTICE(if BIU is installed elsewhere use --with-BIU=PREFIX.)
        fi
        DEPENDENCYNOTFOUND=1;
fi

# error output if ViennaRNA not found
if test "$VRNANOTFOUND" = "1"; then
        AC_MSG_NOTICE()
        AC_MSG_NOTICE(the Vienna RNA C library is required for building the ELL.)
        AC_MSG_NOTICE(it can be obtained from http://www.tbi.univie.ac.at/~ivo/RNA/.)
        AC_MSG_NOTICE()
        DEPENDENCYNOTFOUND=1;
fi

# error output if ELL-headers not found
if test $ELLNOTFOUND = 1; then
        AC_MSG_NOTICE()
        AC_MSG_NOTICE(the ELL C++ library is required for building the tools.)
        AC_MSG_NOTICE(it can be obtained from http://www.bioinf.uni-freiburg.de/SW/ELL/.)
        AC_MSG_NOTICE()
        if test $ELLPATHSET = 1; then
                AC_MSG_NOTICE(cant find ELL in given path '$with_ELL'.)
        else
                AC_MSG_NOTICE(if ELL is installed elsewhere use --with-ELL=prefix.)
        fi
        DEPENDENCYNOTFOUND=1;
fi

# error ABORT if on of the libraries was not found
if test "$DEPENDENCYNOTFOUND" = "1"; then
        AC_MSG_NOTICE()
        AC_MSG_ERROR(library not found! Try './configure --help'.)
fi


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])

AC_OUTPUT([src/Makefile])

