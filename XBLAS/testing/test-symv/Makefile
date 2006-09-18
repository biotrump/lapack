############################################################################
#
#  Program:         XBLAS
#
#  Module:          Makefile
#
#  Purpose:         Makefile for SYMV testing
#
#  Creation date:   Feb 2, 2001, Version 1.0
#
#  Modified:	    
############################################################################
include ../../../make.inc

DOT_PATH = ../../src/dot
TEST_DOT_PATH = ../test-dot
TEST_GEMM_PATH = ../test-gemm
TEST_SYMV_PATH = ../test-symv
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

TEST_SRCS = do_test_symv.c BLAS_symv_testgen.c symv-support.c 
TEST_OBJS = $(TEST_SRCS:.c=.o)

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(DOT_OBJS)

test: do_test_symv
	@echo Testing SYMV - default quick test
	./do_test_symv 7 1 1.0 3 0.01 > symv.results

test-few: do_test_symv
	@echo Testing SYMV - preforming a small set of tests
	./do_test_symv 7 1 1.0 3 $(DO_FEW_TESTS) > symv.results

test-some: do_test_symv
	@echo Testing SYMV - preforming some tests : not all tests 
	./do_test_symv 7 1 1.0 3 $(DO_SOME_TESTS) > symv.results

test-all: do_test_symv
	@echo Testing SYMV Very Thoroughly
	./do_test_symv 7 1 1.0 3 $(DO_ALL_TESTS) > symv.results

test-none: do_test_symv
	@echo NOT Testing SYMV - PREFORMING NO TESTS. NOT TESTING
	./do_test_symv 7 1 1.0 3 $(DO_NO_TESTS) > symv.results

do_test_symv: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LINKOPT) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_symv -lm

test-lib: $(TEST_OBJS)

.c.o:
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK symv.results do_test_symv core

