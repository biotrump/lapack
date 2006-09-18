############################################################################
#
#  Program:         XBLAS
#
#  Module:          Makefile
#
#  Purpose:         Makefile to compile source code into the library
#
#  Creation date:   June 1, 2000, Version 1.0
#
#  Modified:	    
############################################################################
include ../../../make.inc

LIB_PATH = ../..
SRC_PATH = .
HEADER_PATH = ..

HEMV_SRCS =\
	$(SRC_PATH)/BLAS_chemv.c\
	$(SRC_PATH)/BLAS_zhemv.c\
	$(SRC_PATH)/BLAS_zhemv_z_c.c\
	$(SRC_PATH)/BLAS_zhemv_c_z.c\
	$(SRC_PATH)/BLAS_zhemv_c_c.c\
	$(SRC_PATH)/BLAS_chemv_c_s.c\
	$(SRC_PATH)/BLAS_zhemv_z_d.c\
	$(SRC_PATH)/BLAS_chemv_x.c\
	$(SRC_PATH)/BLAS_zhemv_x.c\
	$(SRC_PATH)/BLAS_zhemv_z_c_x.c\
	$(SRC_PATH)/BLAS_zhemv_c_z_x.c\
	$(SRC_PATH)/BLAS_zhemv_c_c_x.c\
	$(SRC_PATH)/BLAS_chemv_c_s_x.c\
	$(SRC_PATH)/BLAS_zhemv_z_d_x.c

HEMV_OBJS = $(HEMV_SRCS:.c=.o)

all: $(HEMV_OBJS)

lib: $(HEMV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(HEMV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

