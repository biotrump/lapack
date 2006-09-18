

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cge_sum_mv(enum blas_order_type order, int m, int n,
		     const void *alpha, const void *a, int lda,
		     const void *x, int incx,
		     const void *beta, const void *b, int ldb,
		     void *y, int incy);


extern void FC_FUNC_(blas_cge_sum_mv, BLAS_CGE_SUM_MV)
 
  (int *m, int *n, const void *alpha, const void *a, int *lda, const void *x,
   int *incx, const void *beta, const void *b, int *ldb, void *y, int *incy) {
  BLAS_cge_sum_mv(blas_colmajor, *m, *n, alpha, a, *lda, x, *incx, beta, b,
		  *ldb, y, *incy);
}
