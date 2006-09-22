#include "blas_extended.h"
#include "blas_fpu.h"
void BLAS_ssum(int n, const float *x, int incx, float *sum)

/*
 * Purpose
 * =======
 * 
 * This routine computes the summation:
 * 
 *     sum <- SUM_{i=0, n-1} x[i].
 * 
 * Arguments
 * =========
 *
 * n      (input) int
 *        The length of vector x.
 * 
 * x      (input) const float*
 *        Array of length n.
 * 
 * incx   (input) int
 *        The stride used to access components x[i].
 *
 * sum    (output) float*
 * 
 */
{
  static const char routine_name[] = "BLAS_ssum";

  int i, xi;
  float *sum_i = sum;
  const float *x_i = x;
  float x_elem;
  float tmp;


  /* Test the input parameters. */
  if (n < 0)
    BLAS_error(routine_name, -1, n, NULL);
  if (incx == 0)
    BLAS_error(routine_name, -3, incx, NULL);

  /* Immediate return. */
  if (n <= 0) {
    *sum_i = 0.0;
    return;
  }



  tmp = 0.0;


  if (incx < 0)
    xi = -(n - 1) * incx;
  else
    xi = 0;

  for (i = 0; i < n; i++, xi += incx) {
    x_elem = x_i[xi];
    tmp = tmp + x_elem;
  }
  *sum = tmp;



}
