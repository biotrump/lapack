############################################################################
#
#  Program:         XBLAS
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   September 1, 1999, Version 0.0
#
#  Modified:	    June 1, 2000, Version 1.0
############################################################################
#
# Basic include options.
# CC is the C compiler, normally invoked with options OPTS.
# LINKER and LINKOPTS function as CC and OPTS, but for linking stages.
#
VERSION	= 1.0

CC = gcc
OPTS = -O3 -Dx86 -Wall
LINKER = $(CC)
LINKOPTS = 

#
#  The name of the libraries to be created/linked to
#
XBLASLIB = libxblas.a

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib
