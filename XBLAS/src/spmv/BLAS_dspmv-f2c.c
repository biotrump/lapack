

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dspmv(enum blas_order_type order, enum blas_uplo_type uplo,
		int n, double alpha, const double *ap,
		const double *x, int incx, double beta, double *y, int incy);


extern void FC_FUNC_(blas_dspmv, BLAS_DSPMV)
 
  (int *uplo, int *n, double *alpha, const double *ap, const double *x,
   int *incx, double *beta, double *y, int *incy) {
  BLAS_dspmv(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, ap, x,
	     *incx, *beta, y, *incy);
}
