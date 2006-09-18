############################################################################
#
#  Program:         XBLAS
#
#  Module:          Makefile
#
#  Purpose:         Makefile for SYMM testing
#
#  Creation date:   June 1, 2000, Version 1.0
#
#  Modified:	    
############################################################################
include ../../../make.inc

DOT_PATH = ../../src/dot
TEST_DOT_PATH = ../test-dot
TEST_GEMM_PATH = ../test-gemm
LIB_PATH = ../..
HEADER_PATH = ../../src

DOT_OBJS = $(DOT_PATH)/BLAS_sdot_x.o $(DOT_PATH)/BLAS_ddot_x.o

GEMM_TEST_OBJS	= $(TEST_GEMM_PATH)/gemm-support.o

TEST_DOT_OBJS = $(TEST_DOT_PATH)/testgen_BLAS_sdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_ddot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_cdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_zdot.o \
		$(TEST_DOT_PATH)/testgen_aux.o \
		$(TEST_DOT_PATH)/test_dot.o \
		$(TEST_DOT_PATH)/BLAS_dot_testgen.o

TEST_SRCS = do_test_symm.c BLAS_symm_testgen.c symm-support.c 
TEST_OBJS = $(TEST_SRCS:.c=.o)

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(DOT_OBJS) $(GEMM_TEST_OBJS)

test: do_test_symm
	@echo Testing SYMM - default quick test
	./do_test_symm 7 1 1.0 3 0.001 > symm.results

test-few: do_test_symm
	@echo Testing SYMM - preforming a small set of tests
	./do_test_symm 7 1 1.0 3 $(DO_FEW_TESTS) > symm.results

test-some: do_test_symm
	@echo Testing SYMM - preforming some tests : not all tests 
	./do_test_symm 7 1 1.0 3 $(DO_SOME_TESTS) > symm.results

test-all: do_test_symm
	@echo Testing SYMM Very Thoroughly
	./do_test_symm 7 1 1.0 3 $(DO_ALL_TESTS) > symm.results

test-none: do_test_symm
	@echo NOT Testing SYMM - PREFORMING NO TESTS. NOT TESTING
	./do_test_symm 7 1 1.0 3 $(DO_NO_TESTS) > symm.results

do_test_symm: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LINKOPT) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_symm -lm

test-lib: $(TEST_OBJS)

# $(GEMM_TEST_OBJS):
#	cd $(TEST_GEMM_PATH) && $(MAKE) test-lib

.c.o:
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK symm.results do_test_symm core

