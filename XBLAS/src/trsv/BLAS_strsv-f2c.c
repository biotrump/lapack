

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_strsv(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, enum blas_diag_type diag,
		int n, float alpha, const float *T, int ldt,
		float *x, int incx);


extern void FC_FUNC_(blas_strsv, BLAS_STRSV)
 
  (int *uplo, int *trans, int *diag, int *n, float *alpha, const float *T,
   int *ldt, float *x, int *incx) {
  BLAS_strsv(blas_colmajor, (enum blas_uplo_type) *uplo,
	     (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	     *alpha, T, *ldt, x, *incx);
}
