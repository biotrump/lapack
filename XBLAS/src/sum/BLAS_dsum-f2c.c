

#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_dsum(int n, const double *x, int incx, double *sum);


extern void FC_FUNC_(blas_dsum, BLAS_DSUM)
  (int *n, const double *x, int *incx, double *sum) {
  BLAS_dsum(*n, x, *incx, sum);
}
