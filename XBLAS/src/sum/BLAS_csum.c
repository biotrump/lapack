#include "blas_extended.h"
#include "blas_fpu.h"
void BLAS_csum(int n, const void *x, int incx, void *sum)

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
 * x      (input) const void*
 *        Array of length n.
 * 
 * incx   (input) int
 *        The stride used to access components x[i].
 *
 * sum    (output) void*
 * 
 */
{
  const char routine_name[] = "BLAS_csum";

  int i, xi;
  float *sum_i = (float *) sum;
  const float *x_i = (float *) x;
  float x_elem[2];
  float tmp[2];


  /* Test the input parameters. */
  if (n < 0)
    BLAS_error(routine_name, -1, n, NULL);
  if (incx == 0)
    BLAS_error(routine_name, -3, incx, NULL);

  /* Immediate return. */
  if (n <= 0) {
    sum_i[0] = sum_i[1] = 0.0;
    return;
  }



  tmp[0] = tmp[1] = 0.0;

  incx *= 2;
  if (incx < 0)
    xi = -(n - 1) * incx;
  else
    xi = 0;

  for (i = 0; i < n; i++, xi += incx) {
    x_elem[0] = x_i[xi];
    x_elem[1] = x_i[xi + 1];
    tmp[0] = tmp[0] + x_elem[0];
    tmp[1] = tmp[1] + x_elem[1];
  }
  ((float *) sum)[0] = tmp[0];
  ((float *) sum)[1] = tmp[1];



}
