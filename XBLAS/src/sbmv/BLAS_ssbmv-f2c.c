

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ssbmv(enum blas_order_type order, enum blas_uplo_type uplo,
		int n, int k, float alpha, const float *a, int lda,
		const float *x, int incx, float beta, float *y, int incy);


extern void FC_FUNC_(blas_ssbmv, BLAS_SSBMV)
 
  (int *uplo, int *n, int *k, float *alpha, const float *a, int *lda,
   const float *x, int *incx, float *beta, float *y, int *incy) {
  BLAS_ssbmv(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *k, *alpha, a,
	     *lda, x, *incx, *beta, y, *incy);
}
