############################################################################
#
#  Program:         XBLAS
#
#  Module:          Makefile
#
#  Purpose:         Makefile for WAXPBY testing
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

DOT_OBJS = $(DOT_PATH)/BLAS_sdot_x.o  $(DOT_PATH)/BLAS_ddot_x.o

TEST_SRCS = do_test_waxpby.c 

TEST_DOT_OBJS = $(TEST_DOT_PATH)/test_dot.o \
           $(TEST_DOT_PATH)/testgen_BLAS_sdot.o \
           $(TEST_DOT_PATH)/testgen_BLAS_ddot.o \
           $(TEST_DOT_PATH)/testgen_BLAS_cdot.o \
           $(TEST_DOT_PATH)/testgen_BLAS_zdot.o \
	   $(TEST_DOT_PATH)/testgen_aux.o \
           $(TEST_DOT_PATH)/BLAS_dot_testgen.o

TEST_OBJS = $(TEST_SRCS:.c=.o)

ALL_OBJS = $(TEST_OBJS) $(DOT_OBJS) $(TEST_DOT_OBJS)

test: do_test_waxpby
	@echo Testing WAXPBY
	./do_test_waxpby 20 1 1.0 0 0.1 > waxpby.results

do_test_waxpby: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LINKOPT) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_waxpby -lm

.c.o:
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK waxpby.results do_test_waxpby core

