

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_stpmv(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, enum blas_diag_type diag,
		int n, float alpha, const float *tp, float *x, int incx);


extern void FC_FUNC_(blas_stpmv, BLAS_STPMV)
 
  (int *uplo, int *trans, int *diag, int *n, float *alpha, const float *tp,
   float *x, int *incx) {
  BLAS_stpmv(blas_colmajor, (enum blas_uplo_type) *uplo,
	     (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	     *alpha, tp, x, *incx);
}
