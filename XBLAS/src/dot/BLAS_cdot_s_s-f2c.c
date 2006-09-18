

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_cdot_s_s(enum blas_conj_type conj, int n, const void *alpha,
		   const float *x, int incx, const void *beta,
		   const float *y, int incy, void *r);


extern void FC_FUNC_(blas_cdot_s_s, BLAS_CDOT_S_S)
 
  (int *conj, int *n, const void *alpha, const float *x, int *incx,
   const void *beta, const float *y, int *incy, void *r) {
  BLAS_cdot_s_s((enum blas_conj_type) *conj, *n, alpha, x, *incx, beta, y,
		*incy, r);
}
