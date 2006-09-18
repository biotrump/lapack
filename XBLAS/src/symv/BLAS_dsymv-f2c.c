

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsymv(enum blas_order_type order, enum blas_uplo_type uplo,
		int n, double alpha, const double *a, int lda,
		const double *x, int incx, double beta, double *y, int incy);


extern void FC_FUNC_(blas_dsymv, BLAS_DSYMV)
 
  (int *uplo, int *n, double *alpha, const double *a, int *lda,
   const double *x, int *incx, double *beta, double *y, int *incy) {
  BLAS_dsymv(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, a, *lda,
	     x, *incx, *beta, y, *incy);
}
