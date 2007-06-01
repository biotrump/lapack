include ../../make.conf
include ../../$(MAKEINC)

DOT_PATH = ../../src/dot
TEST_DOT_PATH = ../test-dot
TEST_HEMM_PATH = ../test-hemm
TEST_SYMV_PATH = ../test-symv
LIB_PATH = ../../$(OUTPUT_DIR)
HEADER_PATH = ../../src

DOT_OBJS = $(DOT_PATH)/BLAS_sdot_x.o $(DOT_PATH)/BLAS_ddot_x.o

HEMM_TEST_OBJS	= $(TEST_HEMM_PATH)/hemm-support.o
SYMV_TEST_OBJS  = $(TEST_SYMV_PATH)/symv-support.o \
		$(TEST_SYMV_PATH)/BLAS_symv_testgen.o

TEST_DOT_OBJS = $(TEST_DOT_PATH)/testgen_BLAS_sdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_ddot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_cdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_zdot.o \
		$(TEST_DOT_PATH)/testgen_aux.o \
		$(TEST_DOT_PATH)/test_dot.o \
		$(TEST_DOT_PATH)/BLAS_dot_testgen.o

TEST_SRCS = do_test_hemv.c BLAS_hemv_testgen.c hemv-support.c 
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(DOT_OBJS) $(HEMM_TEST_OBJS) \
		$(SYMV_TEST_OBJS) $(COMMON_OBJS)

test: do_test_hemv
	@echo Testing HEMV
	./do_test_hemv 7 1 1.0 3 0.01 > hemv.results

do_test_hemv: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_hemv $(EXTRA_LIBS)

test-lib: $(TEST_OBJS)

.c.o:
	$(CC) $(CFLAGS) -I$(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK hemv.results do_test_hemv core

