############################################################################
#
#  Program:         XBLAS
#
#  Module:          Makefile
#
#  Purpose:         Makefile to create common objects
#
#  Creation date:   June 1, 2000, Version 1.0
#
#  Modified:	    
############################################################################
include ../../../make.inc

LIB_PATH = ../..
HEADER_PATH = ..
COMM_OBJS = BLAS_error.o blas_malloc.o BLAS_fpinfo_x.o BLAS_fpinfo_x-f2c.o

all: $(COMM_OBJS)

lib: $(COMM_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(COMM_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

.c.o:
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ core


