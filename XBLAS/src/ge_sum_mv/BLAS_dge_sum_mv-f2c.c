

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dge_sum_mv(enum blas_order_type order, int m, int n,
		     double alpha, const double *a, int lda,
		     const double *x, int incx,
		     double beta, const double *b, int ldb,
		     double *y, int incy);


extern void FC_FUNC_(blas_dge_sum_mv, BLAS_DGE_SUM_MV)
 
  (int *m, int *n, double *alpha, const double *a, int *lda, const double *x,
   int *incx, double *beta, const double *b, int *ldb, double *y, int *incy) {
  BLAS_dge_sum_mv(blas_colmajor, *m, *n, *alpha, a, *lda, x, *incx, *beta, b,
		  *ldb, y, *incy);
}
