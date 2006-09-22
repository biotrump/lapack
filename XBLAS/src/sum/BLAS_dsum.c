#include "blas_extended.h"
#include "blas_fpu.h"
void BLAS_dsum(int n, const double *x, int incx, double *sum)

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
 * x      (input) const double*
 *        Array of length n.
 * 
 * incx   (input) int
 *        The stride used to access components x[i].
 *
 * sum    (output) double*
 * 
 */
{
  static const char routine_name[] = "BLAS_dsum";

  int i, xi;
  double *sum_i = sum;
  const double *x_i = x;
  double x_elem;
  double tmp;


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
