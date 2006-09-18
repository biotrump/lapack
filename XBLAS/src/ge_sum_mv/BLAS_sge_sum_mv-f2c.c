

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_sge_sum_mv(enum blas_order_type order, int m, int n,
		     float alpha, const float *a, int lda,
		     const float *x, int incx,
		     float beta, const float *b, int ldb, float *y, int incy);


extern void FC_FUNC_(blas_sge_sum_mv, BLAS_SGE_SUM_MV)
 
  (int *m, int *n, float *alpha, const float *a, int *lda, const float *x,
   int *incx, float *beta, const float *b, int *ldb, float *y, int *incy) {
  BLAS_sge_sum_mv(blas_colmajor, *m, *n, *alpha, a, *lda, x, *incx, *beta, b,
		  *ldb, y, *incy);
}
