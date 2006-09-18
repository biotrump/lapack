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

SBMV_SRCS =\
	$(SRC_PATH)/BLAS_ssbmv.c\
	$(SRC_PATH)/BLAS_dsbmv.c\
	$(SRC_PATH)/BLAS_csbmv.c\
	$(SRC_PATH)/BLAS_zsbmv.c\
	$(SRC_PATH)/BLAS_dsbmv_d_s.c\
	$(SRC_PATH)/BLAS_dsbmv_s_d.c\
	$(SRC_PATH)/BLAS_dsbmv_s_s.c\
	$(SRC_PATH)/BLAS_zsbmv_z_c.c\
	$(SRC_PATH)/BLAS_zsbmv_c_z.c\
	$(SRC_PATH)/BLAS_zsbmv_c_c.c\
	$(SRC_PATH)/BLAS_csbmv_c_s.c\
	$(SRC_PATH)/BLAS_csbmv_s_c.c\
	$(SRC_PATH)/BLAS_csbmv_s_s.c\
	$(SRC_PATH)/BLAS_zsbmv_z_d.c\
	$(SRC_PATH)/BLAS_zsbmv_d_z.c\
	$(SRC_PATH)/BLAS_zsbmv_d_d.c\
	$(SRC_PATH)/BLAS_ssbmv_x.c\
	$(SRC_PATH)/BLAS_dsbmv_x.c\
	$(SRC_PATH)/BLAS_csbmv_x.c\
	$(SRC_PATH)/BLAS_zsbmv_x.c\
	$(SRC_PATH)/BLAS_dsbmv_d_s_x.c\
	$(SRC_PATH)/BLAS_dsbmv_s_d_x.c\
	$(SRC_PATH)/BLAS_dsbmv_s_s_x.c\
	$(SRC_PATH)/BLAS_zsbmv_z_c_x.c\
	$(SRC_PATH)/BLAS_zsbmv_c_z_x.c\
	$(SRC_PATH)/BLAS_zsbmv_c_c_x.c\
	$(SRC_PATH)/BLAS_csbmv_c_s_x.c\
	$(SRC_PATH)/BLAS_csbmv_s_c_x.c\
	$(SRC_PATH)/BLAS_csbmv_s_s_x.c\
	$(SRC_PATH)/BLAS_zsbmv_z_d_x.c\
	$(SRC_PATH)/BLAS_zsbmv_d_z_x.c\
	$(SRC_PATH)/BLAS_zsbmv_d_d_x.c

SBMV_OBJS = $(SBMV_SRCS:.c=.o)

all: $(SBMV_OBJS)

lib: $(SBMV_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(SBMV_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

