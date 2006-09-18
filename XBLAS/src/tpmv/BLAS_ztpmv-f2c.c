

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ztpmv(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, enum blas_diag_type diag,
		int n, const void *alpha, const void *tp, void *x, int incx);


extern void FC_FUNC_(blas_ztpmv, BLAS_ZTPMV)
 
  (int *uplo, int *trans, int *diag, int *n, const void *alpha,
   const void *tp, void *x, int *incx) {
  BLAS_ztpmv(blas_colmajor, (enum blas_uplo_type) *uplo,
	     (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	     alpha, tp, x, *incx);
}
