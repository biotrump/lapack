

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dtbsv(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, enum blas_diag_type diag,
		int n, int k, double alpha, const double *t, int ldt,
		double *x, int incx);


extern void FC_FUNC_(blas_dtbsv, BLAS_DTBSV)
 
  (int *uplo, int *trans, int *diag, int *n, int *k, double *alpha,
   const double *t, int *ldt, double *x, int *incx) {
  BLAS_dtbsv(blas_colmajor, (enum blas_uplo_type) *uplo,
	     (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	     *k, *alpha, t, *ldt, x, *incx);
}
