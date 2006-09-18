

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ssum(int n, const float *x, int incx, float *sum);


extern void FC_FUNC_(blas_ssum, BLAS_SSUM)
  (int *n, const float *x, int *incx, float *sum) {
  BLAS_ssum(*n, x, *incx, sum);
}
