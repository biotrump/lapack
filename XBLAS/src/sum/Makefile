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
include ../../make.conf
include ../../$(MAKEINC)

LIB_PATH = ../../$(OUTPUT_DIR)
SRC_PATH = .
HEADER_PATH = ..

SUM_SRCS =\
	$(SRC_PATH)/BLAS_ssum.c\
	$(SRC_PATH)/BLAS_dsum.c\
	$(SRC_PATH)/BLAS_csum.c\
	$(SRC_PATH)/BLAS_zsum.c\
	$(SRC_PATH)/BLAS_ssum_x.c\
	$(SRC_PATH)/BLAS_dsum_x.c\
	$(SRC_PATH)/BLAS_csum_x.c\
	$(SRC_PATH)/BLAS_zsum_x.c \
	\
	$(SRC_PATH)/BLAS_ssum-f2c.c\
	$(SRC_PATH)/BLAS_dsum-f2c.c\
	$(SRC_PATH)/BLAS_csum-f2c.c\
	$(SRC_PATH)/BLAS_zsum-f2c.c\
	$(SRC_PATH)/BLAS_ssum_x-f2c.c\
	$(SRC_PATH)/BLAS_dsum_x-f2c.c\
	$(SRC_PATH)/BLAS_csum_x-f2c.c\
	$(SRC_PATH)/BLAS_zsum_x-f2c.c 

SUM_OBJS = $(SUM_SRCS:.c=.o)

all: $(SUM_OBJS)

lib: $(SUM_OBJS)
	$(ARCH) $(ARCHFLAGS) $(LIB_PATH)/$(XBLASLIB) $(SUM_OBJS)
	$(RANLIB) $(LIB_PATH)/$(XBLASLIB)

.c.o: 
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK

