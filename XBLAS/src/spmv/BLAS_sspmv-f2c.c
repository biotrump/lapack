

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_sspmv(enum blas_order_type order, enum blas_uplo_type uplo,
		int n, float alpha, const float *ap,
		const float *x, int incx, float beta, float *y, int incy);


extern void FC_FUNC_(blas_sspmv, BLAS_SSPMV)
 
  (int *uplo, int *n, float *alpha, const float *ap, const float *x,
   int *incx, float *beta, float *y, int *incy) {
  BLAS_sspmv(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, ap, x,
	     *incx, *beta, y, *incy);
}
