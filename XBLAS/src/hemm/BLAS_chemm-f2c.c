

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_chemm(enum blas_order_type order, enum blas_side_type side,
		enum blas_uplo_type uplo, int m, int n,
		const void *alpha, const void *a, int lda,
		const void *b, int ldb, const void *beta, void *c, int ldc);


extern void FC_FUNC_(blas_chemm, BLAS_CHEMM)
 
  (int *side, int *uplo, int *m, int *n, const void *alpha, const void *a,
   int *lda, const void *b, int *ldb, const void *beta, void *c, int *ldc) {
  BLAS_chemm(blas_colmajor, (enum blas_side_type) *side,
	     (enum blas_uplo_type) *uplo, *m, *n, alpha, a, *lda, b, *ldb,
	     beta, c, *ldc);
}
