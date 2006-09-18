

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ztrmv(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, enum blas_diag_type diag, int n,
		const void *alpha, const void *T, int ldt, void *x, int incx);


extern void FC_FUNC_(blas_ztrmv, BLAS_ZTRMV)
 
  (int *uplo, int *trans, int *diag, int *n, const void *alpha, const void *T,
   int *ldt, void *x, int *incx) {
  BLAS_ztrmv(blas_colmajor, (enum blas_uplo_type) *uplo,
	     (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	     alpha, T, *ldt, x, *incx);
}
