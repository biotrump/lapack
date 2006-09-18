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

DOT_SRCS =\
	$(SRC_PATH)/BLAS_sdot.c\
	$(SRC_PATH)/BLAS_ddot.c\
	$(SRC_PATH)/BLAS_cdot.c\
	$(SRC_PATH)/BLAS_zdot.c\
	$(SRC_PATH)/BLAS_ddot_d_s.c\
	$(SRC_PATH)/BLAS_ddot_s_d.c\
	$(SRC_PATH)/BLAS_ddot_s_s.c\
	$(SRC_PATH)/BLAS_zdot_z_c.c\
	$(SRC_PATH)/BLAS_zdot_c_z.c\
	$(SRC_PATH)/BLAS_zdot_c_c.c\
	$(SRC_PATH)/BLAS_cdot_c_s.c\
	$(SRC_PATH)/BLAS_cdot_s_c.c\
	$(SRC_PATH)/BLAS_cdot_s_s.c\
	$(SRC_PATH)/BLAS_zdot_z_d.c\
	$(SRC_PATH)/BLAS_zdot_d_z.c\
	$(SRC_PATH)/BLAS_zdot_d_d.c\
	$(SRC_PATH)/BLAS_sdot_x.c\
	$(SRC_PATH)/BLAS_ddot_x.c\
	$(SRC_PATH)/BLAS_cdot_x.c\
	$(SRC_PATH)/BLAS_zdot_x.c\
	$(SRC_PATH)/BLAS_ddot_d_s_x.c\
	$(SRC_PATH)/BLAS_ddot_s_d_x.c\
	$(SRC_PATH)/BLAS_ddot_s_s_x.c\
	$(SRC_PATH)/BLAS_zdot_z_c_x.c\
	$(SRC_PATH)/BLAS_zdot_c_z_x.c\
	$(SRC_PATH)/BLAS_zdot_c_c_x.c\
	$(SRC_PATH)/BLAS_cdot_c_s_x.c\
	$(SRC_PATH)/BLAS_cdot_s_c_x.c\
	$(SRC_PATH)/BLAS_cdot_s_s_x.c\
	$(SRC_PATH)/BLAS_zdot_z_d_x.c\
	$(SRC_PATH)/BLAS_zdot_d_z_x.c\
	$(SRC_PATH)/BLAS_zdot_d_d_x.c

DOT_OBJS = $(DOT_SRCS:.c=.o)

all: $(DOT_OBJS)

lib: $(DOT_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(DOT_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

