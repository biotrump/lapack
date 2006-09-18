

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ssymm(enum blas_order_type order, enum blas_side_type side,
		enum blas_uplo_type uplo, int m, int n,
		float alpha, const float *a, int lda,
		const float *b, int ldb, float beta, float *c, int ldc);


extern void FC_FUNC_(blas_ssymm, BLAS_SSYMM)
 
  (int *side, int *uplo, int *m, int *n, float *alpha, const float *a,
   int *lda, const float *b, int *ldb, float *beta, float *c, int *ldc) {
  BLAS_ssymm(blas_colmajor, (enum blas_side_type) *side,
	     (enum blas_uplo_type) *uplo, *m, *n, *alpha, a, *lda, b, *ldb,
	     *beta, c, *ldc);
}
