############################################################################
#
#  Program:         XBLAS
#
#  Module:          Makefile
#
#  Purpose:         Makefile for GEMV testing
#
#  Creation date:   June 1, 2000, Version 1.0
#
#  Modified:	    
############################################################################
include ../../../make.inc

DOT_PATH = ../../src/dot
TEST_DOT_PATH = ../test-dot
LIB_PATH = ../..
HEADER_PATH = ../../src

DOT_OBJS = $(DOT_PATH)/BLAS_sdot_x.o $(DOT_PATH)/BLAS_ddot_x.o

TEST_DOT_OBJS = $(TEST_DOT_PATH)/testgen_BLAS_sdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_ddot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_cdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_zdot.o \
		$(TEST_DOT_PATH)/testgen_aux.o \
		$(TEST_DOT_PATH)/test_dot.o \
		$(TEST_DOT_PATH)/BLAS_dot_testgen.o

TEST_SRCS = do_test_gemv.c BLAS_gemv_testgen.c gemv-support.c 
TEST_OBJS = $(TEST_SRCS:.c=.o)

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(DOT_OBJS)

test: do_test_gemv
	@echo Testing GEMV - default quick test
	./do_test_gemv 12 1 1.0 0 0.01 > gemv.results

do_test_gemv: $(LIB_PATH)/$(XBLASLIB) $(ALL_OBJS)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_gemv $(EXTRA_LIBS)

.c.o:
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK gemv.results do_test_gemv core

