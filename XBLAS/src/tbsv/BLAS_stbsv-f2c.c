

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_stbsv(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, enum blas_diag_type diag,
		int n, int k, float alpha, const float *t, int ldt,
		float *x, int incx);


extern void FC_FUNC_(blas_stbsv, BLAS_STBSV)
 
  (int *uplo, int *trans, int *diag, int *n, int *k, float *alpha,
   const float *t, int *ldt, float *x, int *incx) {
  BLAS_stbsv(blas_colmajor, (enum blas_uplo_type) *uplo,
	     (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	     *k, *alpha, t, *ldt, x, *incx);
}
