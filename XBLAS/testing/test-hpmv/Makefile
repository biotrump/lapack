############################################################################
#
#  Program:         XBLAS
#
#  Module:          Makefile
#
#  Purpose:         Makefile for HPMV testing
#
#  Creation date:   January 11, 2001, Version 1.0
#
#  Modified:	    
############################################################################
include ../../make.conf
include ../../$(MAKEINC)

DOT_PATH = ../../src/dot
TEST_DOT_PATH = ../test-dot
TEST_HEMM_PATH = ../test-hemm
TEST_HEMV_PATH = ../test-hemv
TEST_SYMV_PATH = ../test-symv
LIB_PATH = ../../$(OUTPUT_DIR)
HEADER_PATH = ../../src

DOT_OBJS = $(DOT_PATH)/BLAS_sdot_x.o $(DOT_PATH)/BLAS_ddot_x.o

HEMM_TEST_OBJS	= $(TEST_HEMM_PATH)/hemm-support.o
HEMV_TEST_OBJS	= $(TEST_HEMV_PATH)/hemv-support.o \
		$(TEST_HEMV_PATH)/BLAS_hemv_testgen.o
SYMV_TEST_OBJS  = $(TEST_SYMV_PATH)/symv-support.o \
		$(TEST_SYMV_PATH)/BLAS_symv_testgen.o

TEST_DOT_OBJS = $(TEST_DOT_PATH)/testgen_BLAS_sdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_ddot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_cdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_zdot.o \
		$(TEST_DOT_PATH)/testgen_aux.o \
		$(TEST_DOT_PATH)/test_dot.o \
		$(TEST_DOT_PATH)/BLAS_dot_testgen.o

TEST_SRCS = do_test_hpmv.c BLAS_hpmv_testgen.c hpmv-support.c 
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(DOT_OBJS) $(HEMM_TEST_OBJS) \
		$(SYMV_TEST_OBJS) $(HEMV_TEST_OBJS) $(COMMON_OBJS)

test: do_test_hpmv
	@echo Testing HPMV
	./do_test_hpmv 7 1 1.0 3 0.01 > hpmv.results

do_test_hpmv: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_hpmv $(EXTRA_LIBS)

test-lib: $(TEST_OBJS)

.c.o:
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK hpmv.results do_test_hpmv core

