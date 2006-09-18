

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zsymv_d_z_x(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, const void *alpha, const double *a, int lda,
		      const void *x, int incx, const void *beta,
		      void *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_zsymv_d_z_x, BLAS_ZSYMV_D_Z_X)
 
  (int *uplo, int *n, const void *alpha, const double *a, int *lda,
   const void *x, int *incx, const void *beta, void *y, int *incy,
   int *prec) {
  BLAS_zsymv_d_z_x(blas_colmajor, (enum blas_uplo_type) *uplo, *n, alpha, a,
		   *lda, x, *incx, beta, y, *incy,
		   (enum blas_prec_type) *prec);
}
