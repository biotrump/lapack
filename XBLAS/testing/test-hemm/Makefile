############################################################################
#
#  Program:         XBLAS
#
#  Module:          Makefile
#
#  Purpose:         Makefile for HEMM testing
#
#  Creation date:   June 1, 2000, Version 1.0
#
#  Modified:	    
############################################################################
include ../../../make.inc

DOT_PATH = ../../src/dot
TEST_DOT_PATH = ../test-dot
TEST_GEMM_PATH = ../test-gemm
TEST_SYMM_PATH = ../test-symm
LIB_PATH = ../..
HEADER_PATH = ../../src

DOT_OBJS = $(DOT_PATH)/BLAS_sdot_x.o $(DOT_PATH)/BLAS_ddot_x.o

GEMM_TEST_OBJS	= $(TEST_GEMM_PATH)/gemm-support.o
SYMM_TEST_OBJS  = $(TEST_SYMM_PATH)/symm-support.o \
                  $(TEST_SYMM_PATH)/BLAS_symm_testgen.o

TEST_DOT_OBJS = $(TEST_DOT_PATH)/testgen_BLAS_sdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_ddot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_cdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_zdot.o \
		$(TEST_DOT_PATH)/testgen_aux.o \
		$(TEST_DOT_PATH)/test_dot.o \
		$(TEST_DOT_PATH)/BLAS_dot_testgen.o

TEST_SRCS = do_test_hemm.c BLAS_hemm_testgen.c hemm-support.c 
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) \
	   $(DOT_OBJS) $(GEMM_TEST_OBJS) $(SYMM_TEST_OBJS) $(COMMON_OBJS)

test: do_test_hemm
	@echo Testing HEMM - default quick test
	./do_test_hemm 7 1 1.0 0 0.01 > hemm.results

test-few: do_test_hemm
	@echo Testing HEMM - preforming a small set of tests
	./do_test_hemm 7 1 1.0 0 $(DO_FEW_TESTS) > hemm.results

test-some: do_test_hemm
	@echo Testing HEMM - preforming some tests : not all tests 
	./do_test_hemm 7 1 1.0 0 $(DO_SOME_TESTS) > hemm.results

test-all: do_test_hemm
	@echo Testing HEMM Very Thoroughly
	./do_test_hemm 7 1 1.0 0 $(DO_ALL_TESTS) > hemm.results

test-none: do_test_hemm
	@echo NOT Testing HEMM - PREFORMING NO TESTS. NOT TESTING
	./do_test_hemm 7 1 1.0 0 $(DO_NO_TESTS) > hemm.results

do_test_hemm: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_hemm $(EXTRA_LIBS)

# $(GEMM_TEST_OBJS):
#	cd $(TEST_GEMM_PATH) && $(MAKE) test-lib
#
#$(SYMM_TEST_OBJS):
#	cd $(TEST_SYMM_PATH) && $(MAKE) test-lib

.c.o:
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK hemm.results do_test_hemm core

