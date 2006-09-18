

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dgbmv(enum blas_order_type order, enum blas_trans_type trans,
		int m, int n, int kl, int ku, double alpha,
		const double *a, int lda, const double *x, int incx,
		double beta, double *y, int incy);


extern void FC_FUNC_(blas_dgbmv, BLAS_DGBMV)
 
  (int *trans, int *m, int *n, int *kl, int *ku, double *alpha,
   const double *a, int *lda, const double *x, int *incx, double *beta,
   double *y, int *incy) {
  BLAS_dgbmv(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, *kl, *ku,
	     *alpha, a, *lda, x, *incx, *beta, y, *incy);
}
