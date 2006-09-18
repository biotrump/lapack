############################################################################
#
#  Program:         XBLAS
#
#  Module:          Makefile
#
#  Purpose:         Makefile for GEMV2 testing
#
#  Creation date:   August 10, 2005, Version 1.0
#
#  Modified:	    
############################################################################
include ../../../make.inc

TEST_DOT2_PATH = ../test-dot2
LIB_PATH = ../..
HEADER_PATH = ../../src

TEST_DOT2_OBJS = $(TEST_DOT2_PATH)/dot2.o \
		 $(TEST_DOT2_PATH)/testgen_BLAS_sdot2.o \
		 $(TEST_DOT2_PATH)/testgen_BLAS_ddot2.o \
		 $(TEST_DOT2_PATH)/testgen_BLAS_cdot2.o \
		 $(TEST_DOT2_PATH)/testgen_BLAS_zdot2.o \
		 $(TEST_DOT2_PATH)/testgen_aux.o \
		 $(TEST_DOT2_PATH)/test_dot2.o \
		 $(TEST_DOT2_PATH)/r_truth2.o \
		 $(TEST_DOT2_PATH)/BLAS_dot2_testgen.o

TEST_SRCS = do_test_gemv2.c BLAS_gemv2_testgen.c gemv2-support.c 
TEST_OBJS = $(TEST_SRCS:.c=.o)

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT2_OBJS)

test: do_test_gemv2
	@echo Testing GEMV2 - default quick test
	./do_test_gemv2 12 1 1.0 0 0.01 > gemv2.results

do_test_gemv2: $(LIB_PATH)/$(XBLASLIB) $(ALL_OBJS)
	$(LINKER) $(LINKOPT) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_gemv2 -lm

.c.o:
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK gemv2.results do_test_gemv2 core

