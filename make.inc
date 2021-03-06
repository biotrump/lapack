####################################################################
#  LAPACK make include file.                                       #
#  LAPACK, Version 3.6.0                                           #
#  November 2014                                                   #
####################################################################
#
SHELL = /bin/sh
#
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and
#  desired load options for your machine.
#
#FORTRAN  = gfortran
#-fpie -fPIE for  position independent code can be only linked into executables
#-fpic -fPIc for  position independent code
#OPTS     = -O2 -frecursive -fPIC -fpic
#https://wikis.utexas.edu/display/~cdupree/BLAS,+LAPACK,+ATLAS
OPTS     = -fomit-frame-pointer -mfpmath=sse -Os -fno-schedule-insns2 -fPIC -fPIE
DRVOPTS  = $(OPTS)
#NOOPT    = -O0 -frecursive -fPIC -fpic
NOOPT    = -O0 -fomit-frame-pointer -mfpmath=sse -fno-schedule-insns2 -fPIC -fPIE -g
LOADER   = $(FORTRAN)
LOADOPTS =  -lgfortran -lm -pie
#
# Timer for the SECOND and DSECND routines
#
# Default : SECOND and DSECND will use a call to the EXTERNAL FUNCTION ETIME
# TIMER    = EXT_ETIME
# For RS6K : SECOND and DSECND will use a call to the EXTERNAL FUNCTION ETIME_
# TIMER    = EXT_ETIME_
# For gfortran compiler: SECOND and DSECND will use a call to the INTERNAL FUNCTION ETIME
TIMER    = INT_ETIME
# If your Fortran compiler does not provide etime (like Nag Fortran Compiler, etc...)
# SECOND and DSECND will use a call to the INTERNAL FUNCTION CPU_TIME
# TIMER    = INT_CPU_TIME
# If neither of this works...you can use the NONE value... In that case, SECOND and DSECND will always return 0
# TIMER     = NONE
#
#  Configuration LAPACKE: Native C interface to LAPACK
#  To generate LAPACKE library: type 'make lapackelib'
#  Configuration file: turned off (default)
#  Complex types: C99 (default)
#  Name pattern: mixed case (default)
#  (64-bit) Data model: LP64 (default)
#
# CC is the C compiler, normally invoked with options CFLAGS.
#
#CC = ${CC}
#CFLAGS = -O3 -fPIC
#!!!-DHAVE_LAPACK_CONFIG_H  -DLAPACK_COMPLEX_STRUCTURE : to avoid complex.h missing!!!
CFLAGS += -D_POSIX_C_SOURCE=200112L -std=c99 -O3 -DHAVE_LAPACK_CONFIG_H  -DLAPACK_COMPLEX_STRUCTURE
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
#ARCH     = ar
#ARCHFLAGS= cr
ARCHFLAGS= cru
#RANLIB   = ranlib
#
#  Location of the extended-precision BLAS (XBLAS) Fortran library
#  used for building and testing extended-precision routines.  The
#  relevant routines will be compiled and XBLAS will be linked only if
#  USEXBLAS is defined.
#
# USEXBLAS    = Yes
XBLASLIB     =
# XBLASLIB    = -lxblas
#
#  The location of the libraries to which you will link.  (The
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#lapack stanalone build
#BLASLIB      = ../../librefblas.a
#CBLASLIB     = ../../libcblas.a
#Using ATLAS
BLASLIB		= ../../libf77refblas.a
CBLASLIB     = ../../libcblas.a
LAPACKLIB    = liblapack.a
TMGLIB       = libtmglib.a
LAPACKELIB   = liblapacke.a
