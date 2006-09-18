

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dtrsv(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, enum blas_diag_type diag,
		int n, double alpha, const double *T, int ldt,
		double *x, int incx);


extern void FC_FUNC_(blas_dtrsv, BLAS_DTRSV)
 
  (int *uplo, int *trans, int *diag, int *n, double *alpha, const double *T,
   int *ldt, double *x, int *incx) {
  BLAS_dtrsv(blas_colmajor, (enum blas_uplo_type) *uplo,
	     (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	     *alpha, T, *ldt, x, *incx);
}
