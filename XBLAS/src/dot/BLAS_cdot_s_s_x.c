#include "blas_extended.h"
#include "blas_fpu.h"
void BLAS_cdot_s_s_x(enum blas_conj_type conj, int n, const void *alpha,
		     const float *x, int incx, const void *beta,
		     const float *y, int incy,
		     void *r, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 * 
 * This routine computes the inner product:
 * 
 *     r <- beta * r + alpha * SUM_{i=0, n-1} x[i] * y[i].
 * 
 * Arguments
 * =========
 *  
 * conj   (input) enum blas_conj_type
 *        When x and y are complex vectors, specifies whether vector
 *        components x[i] are used unconjugated or conjugated. 
 * 
 * n      (input) int
 *        The length of vectors x and y.
 * 
 * alpha  (input) const void*
 * 
 * x      (input) const float*
 *        Array of length n.
 * 
 * incx   (input) int
 *        The stride used to access components x[i].
 *
 * beta   (input) const void*
 *
 * y      (input) const float*
 *        Array of length n.
 *      
 * incy   (input) int
 *        The stride used to access components y[i].
 *
 * r      (input/output) void*
 * 
 * prec   (input) enum blas_prec_type
 *        Specifies the internal precision to be used.
 *        = blas_prec_single: single precision.
 *        = blas_prec_double: double precision.
 *        = blas_prec_extra : anything at least 1.5 times as accurate
 *                            than double, and wider than 80-bits.
 *                            We use double-double in our implementation.
 *
 */
{
  static const char routine_name[] = "BLAS_cdot_s_s_x";

  switch (prec) {
  case blas_prec_single:
    {
      int i, ix = 0, iy = 0;
      float *r_i = (float *) r;
      const float *x_i = x;
      const float *y_i = y;
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;
      float x_ii;
      float y_ii;
      float r_v[2];
      float prod;
      float sum;
      float tmp1[2];
      float tmp2[2];


      /* Test the input parameters. */
      if (n < 0)
	BLAS_error(routine_name, -2, n, NULL);
      else if (incx == 0)
	BLAS_error(routine_name, -5, incx, NULL);
      else if (incy == 0)
	BLAS_error(routine_name, -8, incy, NULL);

      /* Immediate return. */
      if (((beta_i[0] == 1.0 && beta_i[1] == 0.0))
	  && (n == 0 || (alpha_i[0] == 0.0 && alpha_i[1] == 0.0)))
	return;



      r_v[0] = r_i[0];
      r_v[1] = r_i[0 + 1];
      sum = 0.0;


      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii = y_i[iy];

	prod = x_ii * y_ii;	/* prod = x[i]*y[i] */
	sum = sum + prod;	/* sum = sum+prod */
	ix += incx;
	iy += incy;
      }				/* endfor */


      {
	tmp1[0] = alpha_i[0] * sum;
	tmp1[1] = alpha_i[1] * sum;
      }				/* tmp1 = sum*alpha */
      {
	tmp2[0] = r_v[0] * beta_i[0] - r_v[1] * beta_i[1];
	tmp2[1] = r_v[0] * beta_i[1] + r_v[1] * beta_i[0];
      }
      /* tmp2 = r*beta */
      tmp1[0] = tmp1[0] + tmp2[0];
      tmp1[1] = tmp1[1] + tmp2[1];	/* tmp1 = tmp1+tmp2 */
      ((float *) r)[0] = tmp1[0];
      ((float *) r)[1] = tmp1[1];	/* r = tmp1 */


    }
    break;
  case blas_prec_double:
  case blas_prec_indigenous:
    {
      int i, ix = 0, iy = 0;
      float *r_i = (float *) r;
      const float *x_i = x;
      const float *y_i = y;
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;
      float x_ii;
      float y_ii;
      float r_v[2];
      double prod;
      double sum;
      double tmp1[2];
      double tmp2[2];


      /* Test the input parameters. */
      if (n < 0)
	BLAS_error(routine_name, -2, n, NULL);
      else if (incx == 0)
	BLAS_error(routine_name, -5, incx, NULL);
      else if (incy == 0)
	BLAS_error(routine_name, -8, incy, NULL);

      /* Immediate return. */
      if (((beta_i[0] == 1.0 && beta_i[1] == 0.0))
	  && (n == 0 || (alpha_i[0] == 0.0 && alpha_i[1] == 0.0)))
	return;



      r_v[0] = r_i[0];
      r_v[1] = r_i[0 + 1];
      sum = 0.0;


      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii = y_i[iy];

	prod = (double) x_ii *y_ii;	/* prod = x[i]*y[i] */
	sum = sum + prod;	/* sum = sum+prod */
	ix += incx;
	iy += incy;
      }				/* endfor */


      {
	tmp1[0] = alpha_i[0] * sum;
	tmp1[1] = alpha_i[1] * sum;
      }				/* tmp1 = sum*alpha */
      {
	tmp2[0] = (double) r_v[0] * beta_i[0] - (double) r_v[1] * beta_i[1];
	tmp2[1] = (double) r_v[0] * beta_i[1] + (double) r_v[1] * beta_i[0];
      }				/* tmp2 = r*beta */
      tmp1[0] = tmp1[0] + tmp2[0];
      tmp1[1] = tmp1[1] + tmp2[1];	/* tmp1 = tmp1+tmp2 */
      ((float *) r)[0] = tmp1[0];
      ((float *) r)[1] = tmp1[1];	/* r = tmp1 */


    }
    break;
  case blas_prec_extra:
    {
      int i, ix = 0, iy = 0;
      float *r_i = (float *) r;
      const float *x_i = x;
      const float *y_i = y;
      float *alpha_i = (float *) alpha;
      float *beta_i = (float *) beta;
      float x_ii;
      float y_ii;
      float r_v[2];
      double head_prod, tail_prod;
      double head_sum, tail_sum;
      double head_tmp1[2], tail_tmp1[2];
      double head_tmp2[2], tail_tmp2[2];
      FPU_FIX_DECL;

      /* Test the input parameters. */
      if (n < 0)
	BLAS_error(routine_name, -2, n, NULL);
      else if (incx == 0)
	BLAS_error(routine_name, -5, incx, NULL);
      else if (incy == 0)
	BLAS_error(routine_name, -8, incy, NULL);

      /* Immediate return. */
      if (((beta_i[0] == 1.0 && beta_i[1] == 0.0))
	  && (n == 0 || (alpha_i[0] == 0.0 && alpha_i[1] == 0.0)))
	return;

      FPU_FIX_START;

      r_v[0] = r_i[0];
      r_v[1] = r_i[0 + 1];
      head_sum = tail_sum = 0.0;


      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii = y_i[iy];

	head_prod = (double) x_ii *y_ii;
	tail_prod = 0.0;	/* prod = x[i]*y[i] */
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_sum + head_prod;
	  bv = s1 - head_sum;
	  s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_sum + tail_prod;
	  bv = t1 - tail_sum;
	  t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_sum = t1 + t2;
	  tail_sum = t2 - (head_sum - t1);
	}			/* sum = sum+prod */
	ix += incx;
	iy += incy;
      }				/* endfor */


      {
	double head_e1, tail_e1;
	double dt;
	dt = (double) alpha_i[0];
	{
	  /* Compute double-double = double-double * double. */
	  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	  con = head_sum * split;
	  a11 = con - head_sum;
	  a11 = con - a11;
	  a21 = head_sum - a11;
	  con = dt * split;
	  b1 = con - dt;
	  b1 = con - b1;
	  b2 = dt - b1;

	  c11 = head_sum * dt;
	  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	  c2 = tail_sum * dt;
	  t1 = c11 + c2;
	  t2 = (c2 - (t1 - c11)) + c21;

	  head_e1 = t1 + t2;
	  tail_e1 = t2 - (head_e1 - t1);
	}
	head_tmp1[0] = head_e1;
	tail_tmp1[0] = tail_e1;
	dt = (double) alpha_i[1];
	{
	  /* Compute double-double = double-double * double. */
	  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	  con = head_sum * split;
	  a11 = con - head_sum;
	  a11 = con - a11;
	  a21 = head_sum - a11;
	  con = dt * split;
	  b1 = con - dt;
	  b1 = con - b1;
	  b2 = dt - b1;

	  c11 = head_sum * dt;
	  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	  c2 = tail_sum * dt;
	  t1 = c11 + c2;
	  t2 = (c2 - (t1 - c11)) + c21;

	  head_e1 = t1 + t2;
	  tail_e1 = t2 - (head_e1 - t1);
	}
	head_tmp1[1] = head_e1;
	tail_tmp1[1] = tail_e1;
      }				/* tmp1 = sum*alpha */
      {
	double head_e1, tail_e1;
	double d1;
	double d2;
	/* Real part */
	d1 = (double) r_v[0] * beta_i[0];
	d2 = (double) -r_v[1] * beta_i[1];
	{
	  /* Compute double-double = double + double. */
	  double e, t1, t2;

	  /* Knuth trick. */
	  t1 = d1 + d2;
	  e = t1 - d1;
	  t2 = ((d2 - e) + (d1 - (t1 - e)));

	  /* The result is t1 + t2, after normalization. */
	  head_e1 = t1 + t2;
	  tail_e1 = t2 - (head_e1 - t1);
	}
	head_tmp2[0] = head_e1;
	tail_tmp2[0] = tail_e1;
	/* imaginary part */
	d1 = (double) r_v[0] * beta_i[1];
	d2 = (double) r_v[1] * beta_i[0];
	{
	  /* Compute double-double = double + double. */
	  double e, t1, t2;

	  /* Knuth trick. */
	  t1 = d1 + d2;
	  e = t1 - d1;
	  t2 = ((d2 - e) + (d1 - (t1 - e)));

	  /* The result is t1 + t2, after normalization. */
	  head_e1 = t1 + t2;
	  tail_e1 = t2 - (head_e1 - t1);
	}
	head_tmp2[1] = head_e1;
	tail_tmp2[1] = tail_e1;
      }				/* tmp2 = r*beta */
      {
	double head_t, tail_t;
	double head_a, tail_a;
	double head_b, tail_b;
	/* Real part */
	head_a = head_tmp1[0];
	tail_a = tail_tmp1[0];
	head_b = head_tmp2[0];
	tail_b = tail_tmp2[0];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_tmp1[0] = head_t;
	tail_tmp1[0] = tail_t;
	/* Imaginary part */
	head_a = head_tmp1[1];
	tail_a = tail_tmp1[1];
	head_b = head_tmp2[1];
	tail_b = tail_tmp2[1];
	{
	  /* Compute double-double = double-double + double-double. */
	  double bv;
	  double s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_a + head_b;
	  bv = s1 - head_a;
	  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_a + tail_b;
	  bv = t1 - tail_a;
	  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_t = t1 + t2;
	  tail_t = t2 - (head_t - t1);
	}
	head_tmp1[1] = head_t;
	tail_tmp1[1] = tail_t;
      }				/* tmp1 = tmp1+tmp2 */
      ((float *) r)[0] = head_tmp1[0];
      ((float *) r)[1] = head_tmp1[1];	/* r = tmp1 */

      FPU_FIX_STOP;
    }
    break;
  }
}
