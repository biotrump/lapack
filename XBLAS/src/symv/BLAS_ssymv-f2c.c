

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ssymv(enum blas_order_type order, enum blas_uplo_type uplo,
		int n, float alpha, const float *a, int lda,
		const float *x, int incx, float beta, float *y, int incy);


extern void FC_FUNC_(blas_ssymv, BLAS_SSYMV)
 
  (int *uplo, int *n, float *alpha, const float *a, int *lda, const float *x,
   int *incx, float *beta, float *y, int *incy) {
  BLAS_ssymv(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, a, *lda,
	     x, *incx, *beta, y, *incy);
}
