

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_csum(int n, const void *x, int incx, void *sum);


extern void FC_FUNC_(blas_csum, BLAS_CSUM)
  (int *n, const void *x, int *incx, void *sum) {
  BLAS_csum(*n, x, *incx, sum);
}
