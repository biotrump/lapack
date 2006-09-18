############################################################################
#
#  Program:         XBLAS
#
#  Module:          Makefile
#
#  Purpose:         Makefile for GE_SUM_MV testing
#
#  Creation date:   Mar 14, 2001, Version 1.0
#
#  Modified:	    
############################################################################
include ../../../make.inc

DOT_PATH = ../../src/dot
TEST_DOT_PATH = ../test-dot
TEST_GEMV_PATH = ../test-gemv
TEST_GEMM_PATH = ../test-gemm
TEST_SYMV_PATH = ../test-symv
TEST_GE_SUM_MV_PATH = ../test-ge_sum_mv
LIB_PATH = ../..
HEADER_PATH = ../../src

DOT_OBJS = $(DOT_PATH)/BLAS_sdot_x.o $(DOT_PATH)/BLAS_ddot_x.o
GEMV_OBJS = $(TEST_GEMV_PATH)/gemv-support.o \
	$(TEST_GEMV_PATH)/BLAS_gemv_testgen.o
GEMM_OBJS = $(TEST_GEMM_PATH)/gemm-support.o \
	$(TEST_GEMM_PATH)/BLAS_gemm_testgen.o
SYMV_OBJS = $(TEST_SYMV_PATH)/symv-support.o

TEST_DOT_OBJS = $(TEST_DOT_PATH)/testgen_BLAS_sdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_ddot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_cdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_zdot.o \
		$(TEST_DOT_PATH)/testgen_aux.o \
		$(TEST_DOT_PATH)/test_dot.o \
		$(TEST_DOT_PATH)/BLAS_dot_testgen.o

TEST_SRCS = do_test_ge_sum_mv.c ge_sum_mv-support.c BLAS_ge_sum_mv_testgen.c 
TEST_OBJS = $(TEST_SRCS:.c=.o)

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(DOT_OBJS) $(GEMV_OBJS) \
	$(SYMV_OBJS) $(GEMM_OBJS)

test: do_test_ge_sum_mv warning
	@echo Testing GE_SUM_MV - default quick test
	./do_test_ge_sum_mv 7 1 1.0 3 0.01 > ge_sum_mv.results

test-few: do_test_ge_sum_mv warning
	@echo Testing GE_SUM_MV - preforming a small set of tests
	./do_test_ge_sum_mv 7 1 1.0 3 $(DO_FEW_TESTS) > ge_sum_mv.results

test-some: do_test_ge_sum_mv warning
	@echo Testing GE_SUM_MV - preforming some tests : not all tests 
	./do_test_ge_sum_mv 7 1 1.0 3 $(DO_SOME_TESTS) > ge_sum_mv.results

test-all: do_test_ge_sum_mv warning
	@echo Testing GE_SUM_MV Very Thoroughly
	./do_test_ge_sum_mv 7 1 1.0 3 $(DO_ALL_TESTS) > ge_sum_mv.results

test-none: do_test_ge_sum_mv
	@echo NOT Testing GE_SUM_MV - PREFORMING NO TESTS. NOT TESTING
	./do_test_ge_sum_mv 7 1 1.0 3 $(DO_NO_TESTS) > ge_sum_mv.results

do_test_ge_sum_mv: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LINKOPT) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_ge_sum_mv -lm

test-lib: $(TEST_OBJS)

warning:
	@echo WARNING -- Not Completely tested.
	@echo Please see warning in /m4/test-ge_sum_mv/WARNING

.c.o:
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK ge_sum_mv.results do_test_ge_sum_mv core

