include ../../make.conf
include ../../$(MAKEINC)

DOT_PATH = ../../src/dot
TEST_DOT_PATH = ../test-dot
LIB_PATH = ../../$(OUTPUT_DIR)
HEADER_PATH = -I../../src -I..

DOT_OBJS = $(DOT_PATH)/BLAS_sdot_x.o $(DOT_PATH)/BLAS_ddot_x.o

TEST_DOT_OBJS = $(TEST_DOT_PATH)/testgen_BLAS_sdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_ddot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_cdot.o \
		$(TEST_DOT_PATH)/testgen_BLAS_zdot.o \
		$(TEST_DOT_PATH)/testgen_aux.o \
		$(TEST_DOT_PATH)/test_dot.o \
		$(TEST_DOT_PATH)/BLAS_dot_testgen.o \
		$(TEST_DOT_PATH)/print_vector.o \
		$(TEST_DOT_PATH)/copy_vector.o

TEST_SRCS = do_test_trsv.c BLAS_trsv_testgen.c trsv-support.c \
	BLAS_dot_testgen_trsv.c testgen_BLAS_sdot_x.c
TEST_OBJS = $(TEST_SRCS:.c=.o)
COMMON_OBJS = ../common/dummy_main.o

ALL_OBJS = $(TEST_OBJS) $(TEST_DOT_OBJS) $(DOT_OBJS) $(COMMON_OBJS)

test: do_test_trsv
	@echo Testing TRSV - default quick test
	./do_test_trsv 10 1 1.0 3 0.01 > trsv.results

test-few: do_test_trsv
	@echo Testing TRSV - preforming a small set of tests
	./do_test_trsv 10 1 1.0 0 $(DO_FEW_TESTS) > trsv.results

test-some: do_test_trsv
	@echo Testing TRSV - preforming some tests : not all tests 
	./do_test_trsv 10 1 1.0 0 $(DO_SOME_TESTS) > trsv.results

test-all: do_test_trsv
	@echo Testing TRSV Very Thoroughly
	./do_test_trsv 10 1 1.0 0 $(DO_ALL_TESTS) > trsv.results

test-none: do_test_trsv
	@echo NOT Testing TRSV - PREFORMING NO TESTS. NOT TESTING
	./do_test_trsv 10 1 1.0 0 $(DO_NO_TESTS) > trsv.results

do_test_trsv: $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB)
	$(LINKER) $(LDFLAGS) $(ALL_OBJS) $(LIB_PATH)/$(XBLASLIB) \
	-o do_test_trsv $(EXTRA_LIBS)

.c.o:
	$(CC) $(CFLAGS) $(HEADER_PATH) -c -o $@ $<

clean:
	rm -f *.o *~ *.BAK trsv.results do_test_trsv core

