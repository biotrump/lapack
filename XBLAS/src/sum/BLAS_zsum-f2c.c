

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zsum(int n, const void *x, int incx, void *sum);


extern void FC_FUNC_(blas_zsum, BLAS_ZSUM)
  (int *n, const void *x, int *incx, void *sum) {
  BLAS_zsum(*n, x, *incx, sum);
}
