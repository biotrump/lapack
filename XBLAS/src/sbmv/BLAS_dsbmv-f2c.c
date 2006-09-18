

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsbmv(enum blas_order_type order, enum blas_uplo_type uplo,
		int n, int k, double alpha, const double *a, int lda,
		const double *x, int incx, double beta, double *y, int incy);


extern void FC_FUNC_(blas_dsbmv, BLAS_DSBMV)
 
  (int *uplo, int *n, int *k, double *alpha, const double *a, int *lda,
   const double *x, int *incx, double *beta, double *y, int *incy) {
  BLAS_dsbmv(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *k, *alpha, a,
	     *lda, x, *incx, *beta, y, *incy);
}
