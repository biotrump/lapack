

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_sgbmv(enum blas_order_type order, enum blas_trans_type trans,
		int m, int n, int kl, int ku, float alpha,
		const float *a, int lda, const float *x, int incx,
		float beta, float *y, int incy);


extern void FC_FUNC_(blas_sgbmv, BLAS_SGBMV)
 
  (int *trans, int *m, int *n, int *kl, int *ku, float *alpha, const float *a,
   int *lda, const float *x, int *incx, float *beta, float *y, int *incy) {
  BLAS_sgbmv(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, *kl, *ku,
	     *alpha, a, *lda, x, *incx, *beta, y, *incy);
}
