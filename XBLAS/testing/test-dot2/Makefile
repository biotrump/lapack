############################################################################
#
#  Program:         XBLAS
#
#  Module:          Makefile
#
#  Purpose:         Makefile for DOT2 testing
#
#  Creation date:   July 28, 2005, Version 1.0
#
############################################################################
include ../../../make.inc
LIB_PATH = ../..
HEADER_PATH = ../../src

TEST_SRCS = do_test_dot2.c test_dot2.c BLAS_dot2_testgen.c \
	    dot2.c ../test-dot/testgen_aux.c  r_truth2.c \
	    testgen_BLAS_sdot2.c testgen_BLAS_ddot2.c \
	    testgen_BLAS_zdot2.c testgen_BLAS_cdot2.c

TEST_OBJS = $(TEST_SRCS:.c=.o)

test: do_test_dot2
	@echo Testing DOT2
	./do_test_dot2 20 1 1.0 0 1 > dot2.results

test-lib: $(TEST_OBJS)

do_test_dot2: $(TEST_OBJS)  $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(TEST_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_dot2 $(EXTRA_LIBS)

.c.o:
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK dot2.results do_test_dot2 core

