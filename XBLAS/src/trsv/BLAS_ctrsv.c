#include "blas_extended.h"
#include "blas_fpu.h"

void BLAS_ctrsv(enum blas_order_type order, enum blas_uplo_type uplo,
		enum blas_trans_type trans, enum blas_diag_type diag,
		int n, const void *alpha, const void *T, int ldt,
		void *x, int incx)

/*
 * Purpose
 * =======
 * 
 * This routine solve :
 * 
 *     x <- alpha * inverse(T) * x
 * 
 * Arguments
 * =========
 * 
 * order  (input) enum blas_order_type
 *        column major, row major
 *
 * uplo   (input) enum blas_uplo_type
 *        upper, lower
 *
 * trans  (input) enum blas_trans_type
 *        no trans, trans, conj trans
 * 
 * diag   (input) enum blas_diag_type
 *        unit, non unit
 *
 * n      (input) int
 *        the dimension of T
 * 
 * alpha  (input) const void*
 * 
 * T      (input) void*
 *        Triangular matrix
 *
 * x      (input) const void*
 *           Array of length n.
 * 
 * incx   (input) int
 *           The stride used to access components x[i].
 *
 */
{
  char *routine_name = "BLAS_ctrsv";

  int i, j;			/* used to idx matrix */
  int ix, jx;			/* used to idx vector x */
  int start_x;			/* used as the starting idx to vector x */
  const float *T_i = (float *) T;	/* internal matrix T */
  float *x_i = (float *) x;	/* internal x */
  float *alpha_i = (float *) alpha;	/* internal alpha */
  float T_element[2];		/* temporary variable for an element of matrix A */
  int incT = 1;			/* internal ldt */

  if ((order != blas_rowmajor && order != blas_colmajor) ||
      (uplo != blas_upper && uplo != blas_lower) ||
      (trans != blas_trans && trans !=
       blas_no_trans && trans != blas_conj_trans) ||
      (diag != blas_non_unit_diag && diag != blas_unit_diag) ||
      (ldt < n) || (incx == 0)) {
    BLAS_error(routine_name, 0, 0, NULL);
  }

  if (n <= 0)
    return;

  incT *= 2;
  incx *= 2;
  /* configuring the vector starting idx */
  if (incx <= 0) {
    start_x = -(n - 1) * incx;
  } else {
    start_x = 0;
  }

  /* if alpha is zero, then return x as a zero vector */
  if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
    ix = start_x;
    for (i = 0; i < n; i++) {
      x_i[ix] = 0.0;
      x_i[ix + 1] = 0.0;
      ix += incx;
    }
    return;
  }

  {
    float temp1[2];		/* temporary variable for calculations */
    float temp2[2];		/* temporary variable for calculations */
    float temp3[2];		/* temporary variable for calculations */

    if ((order == blas_rowmajor &&
	 trans == blas_no_trans && uplo == blas_upper) ||
	(order == blas_colmajor &&
	 trans != blas_no_trans && uplo == blas_lower)) {
      if (trans == blas_conj_trans) {

	jx = start_x + (n - 1) * incx;
	for (j = n - 1; j >= 0; j--) {

	  /* compute Xj = alpha*Xj - SUM Tij(or Tji) * Xi
	     i=j+1 to n-1           */
	  temp3[0] = x_i[jx];
	  temp3[1] = x_i[1 + jx];
	  {
	    temp1[0] = temp3[0] * alpha_i[0] - temp3[1] * alpha_i[1];
	    temp1[1] = temp3[0] * alpha_i[1] + temp3[1] * alpha_i[0];
	  }


	  ix = start_x + (n - 1) * incx;
	  for (i = n - 1; i >= j + 1; i--) {
	    T_element[0] = T_i[i * incT + j * ldt * incT];
	    T_element[1] = T_i[i * incT + j * ldt * incT + 1];
	    T_element[1] = -T_element[1];
	    temp3[0] = x_i[ix];
	    temp3[1] = x_i[1 + ix];
	    {
	      temp2[0] = temp3[0] * T_element[0] - temp3[1] * T_element[1];
	      temp2[1] = temp3[0] * T_element[1] + temp3[1] * T_element[0];
	    }

	    temp1[0] = temp1[0] - temp2[0];
	    temp1[1] = temp1[1] - temp2[1];
	    ix -= incx;
	  }			/* for j<n */

	  /* if the diagonal entry is not equal to one, then divide Xj by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element[0] = T_i[j * incT + j * ldt * incT];
	    T_element[1] = T_i[j * incT + j * ldt * incT + 1];
	    T_element[1] = -T_element[1];

	    {
	      double S = 1.0, eps, ov, un;
	      double abs_a, abs_b, abs_c, abs_d, ab, cd;
	      double r;
	      double t;
	      double q[2];

	      eps = pow(2.0, -24.0);
	      un = pow(2.0, -126.0);
	      ov = pow(2.0, 128.0) * (1 - eps);
	      abs_a = fabs((double) temp1[0]);
	      abs_b = fabs((double) temp1[1]);
	      abs_c = fabs((double) T_element[0]);
	      abs_d = fabs((double) T_element[1]);
	      ab = MAX(abs_a, abs_b);
	      cd = MAX(abs_c, abs_d);

	      /* Scaling */
	      if (ab > ov / 16) {	/* scale down a, b */
		temp1[0] /= 16;
		temp1[1] /= 16;
		S = S * 16;
	      }
	      if (cd > ov / 16) {	/* scale down c, d */
		T_element[0] /= 16;
		T_element[1] /= 16;
		S = S / 16;
	      }
	      if (ab < un / eps * 2) {	/* scale up a, b */
		t = 2.0 / (eps * eps);
		temp1[0] *= t;
		temp1[1] *= t;
		S = S / t;
	      }
	      if (cd < un / eps * 2) {	/* scale up c, d */
		t = 2.0 / (eps * eps);
		T_element[0] *= t;
		T_element[1] *= t;
		S = S * t;
	      }

	      /* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
	      if (abs_c > abs_d) {
		r = T_element[1] / T_element[0];
		t = 1 / (T_element[0] + T_element[1] * r);
		q[0] = (temp1[0] + temp1[1] * r) * t;
		q[1] = (temp1[1] - temp1[0] * r) * t;
	      } else {
		r = T_element[0] / T_element[1];
		t = 1 / (T_element[1] + T_element[0] * r);
		q[0] = (temp1[1] + temp1[0] * r) * t;
		q[1] = (-temp1[0] + temp1[1] * r) * t;
	      }
	      /* Scale back */
	      temp1[0] = q[0] * S;
	      temp1[1] = q[1] * S;
	    }

	  }
	  /* if (diag == blas_non_unit_diag) */
	  x_i[jx] = temp1[0];
	  x_i[jx + 1] = temp1[1];

	  jx -= incx;
	}			/* for j>=0 */
      } else {

	jx = start_x + (n - 1) * incx;
	for (j = n - 1; j >= 0; j--) {

	  /* compute Xj = alpha*Xj - SUM Tij(or Tji) * Xi
	     i=j+1 to n-1           */
	  temp3[0] = x_i[jx];
	  temp3[1] = x_i[1 + jx];
	  {
	    temp1[0] = temp3[0] * alpha_i[0] - temp3[1] * alpha_i[1];
	    temp1[1] = temp3[0] * alpha_i[1] + temp3[1] * alpha_i[0];
	  }


	  ix = start_x + (n - 1) * incx;
	  for (i = n - 1; i >= j + 1; i--) {
	    T_element[0] = T_i[i * incT + j * ldt * incT];
	    T_element[1] = T_i[i * incT + j * ldt * incT + 1];

	    temp3[0] = x_i[ix];
	    temp3[1] = x_i[1 + ix];
	    {
	      temp2[0] = temp3[0] * T_element[0] - temp3[1] * T_element[1];
	      temp2[1] = temp3[0] * T_element[1] + temp3[1] * T_element[0];
	    }

	    temp1[0] = temp1[0] - temp2[0];
	    temp1[1] = temp1[1] - temp2[1];
	    ix -= incx;
	  }			/* for j<n */

	  /* if the diagonal entry is not equal to one, then divide Xj by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element[0] = T_i[j * incT + j * ldt * incT];
	    T_element[1] = T_i[j * incT + j * ldt * incT + 1];


	    {
	      double S = 1.0, eps, ov, un;
	      double abs_a, abs_b, abs_c, abs_d, ab, cd;
	      double r;
	      double t;
	      double q[2];

	      eps = pow(2.0, -24.0);
	      un = pow(2.0, -126.0);
	      ov = pow(2.0, 128.0) * (1 - eps);
	      abs_a = fabs((double) temp1[0]);
	      abs_b = fabs((double) temp1[1]);
	      abs_c = fabs((double) T_element[0]);
	      abs_d = fabs((double) T_element[1]);
	      ab = MAX(abs_a, abs_b);
	      cd = MAX(abs_c, abs_d);

	      /* Scaling */
	      if (ab > ov / 16) {	/* scale down a, b */
		temp1[0] /= 16;
		temp1[1] /= 16;
		S = S * 16;
	      }
	      if (cd > ov / 16) {	/* scale down c, d */
		T_element[0] /= 16;
		T_element[1] /= 16;
		S = S / 16;
	      }
	      if (ab < un / eps * 2) {	/* scale up a, b */
		t = 2.0 / (eps * eps);
		temp1[0] *= t;
		temp1[1] *= t;
		S = S / t;
	      }
	      if (cd < un / eps * 2) {	/* scale up c, d */
		t = 2.0 / (eps * eps);
		T_element[0] *= t;
		T_element[1] *= t;
		S = S * t;
	      }

	      /* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
	      if (abs_c > abs_d) {
		r = T_element[1] / T_element[0];
		t = 1 / (T_element[0] + T_element[1] * r);
		q[0] = (temp1[0] + temp1[1] * r) * t;
		q[1] = (temp1[1] - temp1[0] * r) * t;
	      } else {
		r = T_element[0] / T_element[1];
		t = 1 / (T_element[1] + T_element[0] * r);
		q[0] = (temp1[1] + temp1[0] * r) * t;
		q[1] = (-temp1[0] + temp1[1] * r) * t;
	      }
	      /* Scale back */
	      temp1[0] = q[0] * S;
	      temp1[1] = q[1] * S;
	    }

	  }
	  /* if (diag == blas_non_unit_diag) */
	  x_i[jx] = temp1[0];
	  x_i[jx + 1] = temp1[1];

	  jx -= incx;
	}			/* for j>=0 */
      }
    } else if ((order == blas_rowmajor &&
		trans == blas_no_trans && uplo == blas_lower) ||
	       (order == blas_colmajor &&
		trans != blas_no_trans && uplo == blas_upper)) {
      if (trans == blas_conj_trans) {

	jx = start_x;
	for (j = 0; j < n; j++) {

	  /* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
	     i=j+1 to n-1           */
	  temp3[0] = x_i[jx];
	  temp3[1] = x_i[1 + jx];
	  /* multiply by alpha */
	  {
	    temp1[0] = temp3[0] * alpha_i[0] - temp3[1] * alpha_i[1];
	    temp1[1] = temp3[0] * alpha_i[1] + temp3[1] * alpha_i[0];
	  }


	  ix = start_x;
	  for (i = 0; i < j; i++) {
	    T_element[0] = T_i[i * incT + j * ldt * incT];
	    T_element[1] = T_i[i * incT + j * ldt * incT + 1];
	    T_element[1] = -T_element[1];
	    temp3[0] = x_i[ix];
	    temp3[1] = x_i[1 + ix];
	    {
	      temp2[0] = temp3[0] * T_element[0] - temp3[1] * T_element[1];
	      temp2[1] = temp3[0] * T_element[1] + temp3[1] * T_element[0];
	    }

	    temp1[0] = temp1[0] - temp2[0];
	    temp1[1] = temp1[1] - temp2[1];
	    ix += incx;
	  }			/* for i<j */

	  /* if the diagonal entry is not equal to one, then divide Xj by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element[0] = T_i[j * incT + j * ldt * incT];
	    T_element[1] = T_i[j * incT + j * ldt * incT + 1];
	    T_element[1] = -T_element[1];

	    {
	      double S = 1.0, eps, ov, un;
	      double abs_a, abs_b, abs_c, abs_d, ab, cd;
	      double r;
	      double t;
	      double q[2];

	      eps = pow(2.0, -24.0);
	      un = pow(2.0, -126.0);
	      ov = pow(2.0, 128.0) * (1 - eps);
	      abs_a = fabs((double) temp1[0]);
	      abs_b = fabs((double) temp1[1]);
	      abs_c = fabs((double) T_element[0]);
	      abs_d = fabs((double) T_element[1]);
	      ab = MAX(abs_a, abs_b);
	      cd = MAX(abs_c, abs_d);

	      /* Scaling */
	      if (ab > ov / 16) {	/* scale down a, b */
		temp1[0] /= 16;
		temp1[1] /= 16;
		S = S * 16;
	      }
	      if (cd > ov / 16) {	/* scale down c, d */
		T_element[0] /= 16;
		T_element[1] /= 16;
		S = S / 16;
	      }
	      if (ab < un / eps * 2) {	/* scale up a, b */
		t = 2.0 / (eps * eps);
		temp1[0] *= t;
		temp1[1] *= t;
		S = S / t;
	      }
	      if (cd < un / eps * 2) {	/* scale up c, d */
		t = 2.0 / (eps * eps);
		T_element[0] *= t;
		T_element[1] *= t;
		S = S * t;
	      }

	      /* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
	      if (abs_c > abs_d) {
		r = T_element[1] / T_element[0];
		t = 1 / (T_element[0] + T_element[1] * r);
		q[0] = (temp1[0] + temp1[1] * r) * t;
		q[1] = (temp1[1] - temp1[0] * r) * t;
	      } else {
		r = T_element[0] / T_element[1];
		t = 1 / (T_element[1] + T_element[0] * r);
		q[0] = (temp1[1] + temp1[0] * r) * t;
		q[1] = (-temp1[0] + temp1[1] * r) * t;
	      }
	      /* Scale back */
	      temp1[0] = q[0] * S;
	      temp1[1] = q[1] * S;
	    }

	  }
	  /* if (diag == blas_non_unit_diag) */
	  x_i[jx] = temp1[0];
	  x_i[jx + 1] = temp1[1];
	  jx += incx;
	}			/* for j<n */
      } else {

	jx = start_x;
	for (j = 0; j < n; j++) {

	  /* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
	     i=j+1 to n-1           */
	  temp3[0] = x_i[jx];
	  temp3[1] = x_i[1 + jx];
	  /* multiply by alpha */
	  {
	    temp1[0] = temp3[0] * alpha_i[0] - temp3[1] * alpha_i[1];
	    temp1[1] = temp3[0] * alpha_i[1] + temp3[1] * alpha_i[0];
	  }


	  ix = start_x;
	  for (i = 0; i < j; i++) {
	    T_element[0] = T_i[i * incT + j * ldt * incT];
	    T_element[1] = T_i[i * incT + j * ldt * incT + 1];

	    temp3[0] = x_i[ix];
	    temp3[1] = x_i[1 + ix];
	    {
	      temp2[0] = temp3[0] * T_element[0] - temp3[1] * T_element[1];
	      temp2[1] = temp3[0] * T_element[1] + temp3[1] * T_element[0];
	    }

	    temp1[0] = temp1[0] - temp2[0];
	    temp1[1] = temp1[1] - temp2[1];
	    ix += incx;
	  }			/* for i<j */

	  /* if the diagonal entry is not equal to one, then divide Xj by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element[0] = T_i[j * incT + j * ldt * incT];
	    T_element[1] = T_i[j * incT + j * ldt * incT + 1];


	    {
	      double S = 1.0, eps, ov, un;
	      double abs_a, abs_b, abs_c, abs_d, ab, cd;
	      double r;
	      double t;
	      double q[2];

	      eps = pow(2.0, -24.0);
	      un = pow(2.0, -126.0);
	      ov = pow(2.0, 128.0) * (1 - eps);
	      abs_a = fabs((double) temp1[0]);
	      abs_b = fabs((double) temp1[1]);
	      abs_c = fabs((double) T_element[0]);
	      abs_d = fabs((double) T_element[1]);
	      ab = MAX(abs_a, abs_b);
	      cd = MAX(abs_c, abs_d);

	      /* Scaling */
	      if (ab > ov / 16) {	/* scale down a, b */
		temp1[0] /= 16;
		temp1[1] /= 16;
		S = S * 16;
	      }
	      if (cd > ov / 16) {	/* scale down c, d */
		T_element[0] /= 16;
		T_element[1] /= 16;
		S = S / 16;
	      }
	      if (ab < un / eps * 2) {	/* scale up a, b */
		t = 2.0 / (eps * eps);
		temp1[0] *= t;
		temp1[1] *= t;
		S = S / t;
	      }
	      if (cd < un / eps * 2) {	/* scale up c, d */
		t = 2.0 / (eps * eps);
		T_element[0] *= t;
		T_element[1] *= t;
		S = S * t;
	      }

	      /* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
	      if (abs_c > abs_d) {
		r = T_element[1] / T_element[0];
		t = 1 / (T_element[0] + T_element[1] * r);
		q[0] = (temp1[0] + temp1[1] * r) * t;
		q[1] = (temp1[1] - temp1[0] * r) * t;
	      } else {
		r = T_element[0] / T_element[1];
		t = 1 / (T_element[1] + T_element[0] * r);
		q[0] = (temp1[1] + temp1[0] * r) * t;
		q[1] = (-temp1[0] + temp1[1] * r) * t;
	      }
	      /* Scale back */
	      temp1[0] = q[0] * S;
	      temp1[1] = q[1] * S;
	    }

	  }
	  /* if (diag == blas_non_unit_diag) */
	  x_i[jx] = temp1[0];
	  x_i[jx + 1] = temp1[1];
	  jx += incx;
	}			/* for j<n */
      }
    } else if ((order == blas_rowmajor &&
		trans != blas_no_trans && uplo == blas_lower) ||
	       (order == blas_colmajor &&
		trans == blas_no_trans && uplo == blas_upper)) {
      if (trans == blas_conj_trans) {

	jx = start_x + (n - 1) * incx;
	for (j = n - 1; j >= 0; j--) {

	  /* compute Xj = alpha*Xj - SUM Tij(or Tji) * Xi
	     i=j+1 to n-1           */
	  temp3[0] = x_i[jx];
	  temp3[1] = x_i[1 + jx];
	  {
	    temp1[0] = temp3[0] * alpha_i[0] - temp3[1] * alpha_i[1];
	    temp1[1] = temp3[0] * alpha_i[1] + temp3[1] * alpha_i[0];
	  }


	  ix = start_x + (n - 1) * incx;
	  for (i = n - 1; i >= j + 1; i--) {
	    T_element[0] = T_i[j * incT + i * ldt * incT];
	    T_element[1] = T_i[j * incT + i * ldt * incT + 1];
	    T_element[1] = -T_element[1];
	    temp3[0] = x_i[ix];
	    temp3[1] = x_i[1 + ix];
	    {
	      temp2[0] = temp3[0] * T_element[0] - temp3[1] * T_element[1];
	      temp2[1] = temp3[0] * T_element[1] + temp3[1] * T_element[0];
	    }

	    temp1[0] = temp1[0] - temp2[0];
	    temp1[1] = temp1[1] - temp2[1];
	    ix -= incx;
	  }			/* for j<n */

	  /* if the diagonal entry is not equal to one, then divide Xj by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element[0] = T_i[j * incT + j * ldt * incT];
	    T_element[1] = T_i[j * incT + j * ldt * incT + 1];
	    T_element[1] = -T_element[1];

	    {
	      double S = 1.0, eps, ov, un;
	      double abs_a, abs_b, abs_c, abs_d, ab, cd;
	      double r;
	      double t;
	      double q[2];

	      eps = pow(2.0, -24.0);
	      un = pow(2.0, -126.0);
	      ov = pow(2.0, 128.0) * (1 - eps);
	      abs_a = fabs((double) temp1[0]);
	      abs_b = fabs((double) temp1[1]);
	      abs_c = fabs((double) T_element[0]);
	      abs_d = fabs((double) T_element[1]);
	      ab = MAX(abs_a, abs_b);
	      cd = MAX(abs_c, abs_d);

	      /* Scaling */
	      if (ab > ov / 16) {	/* scale down a, b */
		temp1[0] /= 16;
		temp1[1] /= 16;
		S = S * 16;
	      }
	      if (cd > ov / 16) {	/* scale down c, d */
		T_element[0] /= 16;
		T_element[1] /= 16;
		S = S / 16;
	      }
	      if (ab < un / eps * 2) {	/* scale up a, b */
		t = 2.0 / (eps * eps);
		temp1[0] *= t;
		temp1[1] *= t;
		S = S / t;
	      }
	      if (cd < un / eps * 2) {	/* scale up c, d */
		t = 2.0 / (eps * eps);
		T_element[0] *= t;
		T_element[1] *= t;
		S = S * t;
	      }

	      /* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
	      if (abs_c > abs_d) {
		r = T_element[1] / T_element[0];
		t = 1 / (T_element[0] + T_element[1] * r);
		q[0] = (temp1[0] + temp1[1] * r) * t;
		q[1] = (temp1[1] - temp1[0] * r) * t;
	      } else {
		r = T_element[0] / T_element[1];
		t = 1 / (T_element[1] + T_element[0] * r);
		q[0] = (temp1[1] + temp1[0] * r) * t;
		q[1] = (-temp1[0] + temp1[1] * r) * t;
	      }
	      /* Scale back */
	      temp1[0] = q[0] * S;
	      temp1[1] = q[1] * S;
	    }

	  }
	  /* if (diag == blas_non_unit_diag) */
	  x_i[jx] = temp1[0];
	  x_i[jx + 1] = temp1[1];

	  jx -= incx;
	}			/* for j>=0 */
      } else {

	jx = start_x + (n - 1) * incx;
	for (j = n - 1; j >= 0; j--) {

	  /* compute Xj = alpha*Xj - SUM Tij(or Tji) * Xi
	     i=j+1 to n-1           */
	  temp3[0] = x_i[jx];
	  temp3[1] = x_i[1 + jx];
	  {
	    temp1[0] = temp3[0] * alpha_i[0] - temp3[1] * alpha_i[1];
	    temp1[1] = temp3[0] * alpha_i[1] + temp3[1] * alpha_i[0];
	  }


	  ix = start_x + (n - 1) * incx;
	  for (i = n - 1; i >= j + 1; i--) {
	    T_element[0] = T_i[j * incT + i * ldt * incT];
	    T_element[1] = T_i[j * incT + i * ldt * incT + 1];

	    temp3[0] = x_i[ix];
	    temp3[1] = x_i[1 + ix];
	    {
	      temp2[0] = temp3[0] * T_element[0] - temp3[1] * T_element[1];
	      temp2[1] = temp3[0] * T_element[1] + temp3[1] * T_element[0];
	    }

	    temp1[0] = temp1[0] - temp2[0];
	    temp1[1] = temp1[1] - temp2[1];
	    ix -= incx;
	  }			/* for j<n */

	  /* if the diagonal entry is not equal to one, then divide Xj by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element[0] = T_i[j * incT + j * ldt * incT];
	    T_element[1] = T_i[j * incT + j * ldt * incT + 1];


	    {
	      double S = 1.0, eps, ov, un;
	      double abs_a, abs_b, abs_c, abs_d, ab, cd;
	      double r;
	      double t;
	      double q[2];

	      eps = pow(2.0, -24.0);
	      un = pow(2.0, -126.0);
	      ov = pow(2.0, 128.0) * (1 - eps);
	      abs_a = fabs((double) temp1[0]);
	      abs_b = fabs((double) temp1[1]);
	      abs_c = fabs((double) T_element[0]);
	      abs_d = fabs((double) T_element[1]);
	      ab = MAX(abs_a, abs_b);
	      cd = MAX(abs_c, abs_d);

	      /* Scaling */
	      if (ab > ov / 16) {	/* scale down a, b */
		temp1[0] /= 16;
		temp1[1] /= 16;
		S = S * 16;
	      }
	      if (cd > ov / 16) {	/* scale down c, d */
		T_element[0] /= 16;
		T_element[1] /= 16;
		S = S / 16;
	      }
	      if (ab < un / eps * 2) {	/* scale up a, b */
		t = 2.0 / (eps * eps);
		temp1[0] *= t;
		temp1[1] *= t;
		S = S / t;
	      }
	      if (cd < un / eps * 2) {	/* scale up c, d */
		t = 2.0 / (eps * eps);
		T_element[0] *= t;
		T_element[1] *= t;
		S = S * t;
	      }

	      /* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
	      if (abs_c > abs_d) {
		r = T_element[1] / T_element[0];
		t = 1 / (T_element[0] + T_element[1] * r);
		q[0] = (temp1[0] + temp1[1] * r) * t;
		q[1] = (temp1[1] - temp1[0] * r) * t;
	      } else {
		r = T_element[0] / T_element[1];
		t = 1 / (T_element[1] + T_element[0] * r);
		q[0] = (temp1[1] + temp1[0] * r) * t;
		q[1] = (-temp1[0] + temp1[1] * r) * t;
	      }
	      /* Scale back */
	      temp1[0] = q[0] * S;
	      temp1[1] = q[1] * S;
	    }

	  }
	  /* if (diag == blas_non_unit_diag) */
	  x_i[jx] = temp1[0];
	  x_i[jx + 1] = temp1[1];

	  jx -= incx;
	}			/* for j>=0 */
      }
    } else if ((order == blas_rowmajor &&
		trans != blas_no_trans && uplo == blas_upper) ||
	       (order == blas_colmajor &&
		trans == blas_no_trans && uplo == blas_lower)) {
      if (trans == blas_conj_trans) {

	jx = start_x;
	for (j = 0; j < n; j++) {

	  /* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
	     i=j+1 to n-1           */
	  temp3[0] = x_i[jx];
	  temp3[1] = x_i[1 + jx];
	  /* multiply by alpha */
	  {
	    temp1[0] = temp3[0] * alpha_i[0] - temp3[1] * alpha_i[1];
	    temp1[1] = temp3[0] * alpha_i[1] + temp3[1] * alpha_i[0];
	  }


	  ix = start_x;
	  for (i = 0; i < j; i++) {
	    T_element[0] = T_i[j * incT + i * ldt * incT];
	    T_element[1] = T_i[j * incT + i * ldt * incT + 1];
	    T_element[1] = -T_element[1];
	    temp3[0] = x_i[ix];
	    temp3[1] = x_i[1 + ix];
	    {
	      temp2[0] = temp3[0] * T_element[0] - temp3[1] * T_element[1];
	      temp2[1] = temp3[0] * T_element[1] + temp3[1] * T_element[0];
	    }

	    temp1[0] = temp1[0] - temp2[0];
	    temp1[1] = temp1[1] - temp2[1];
	    ix += incx;
	  }			/* for i<j */

	  /* if the diagonal entry is not equal to one, then divide Xj by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element[0] = T_i[j * incT + j * ldt * incT];
	    T_element[1] = T_i[j * incT + j * ldt * incT + 1];
	    T_element[1] = -T_element[1];

	    {
	      double S = 1.0, eps, ov, un;
	      double abs_a, abs_b, abs_c, abs_d, ab, cd;
	      double r;
	      double t;
	      double q[2];

	      eps = pow(2.0, -24.0);
	      un = pow(2.0, -126.0);
	      ov = pow(2.0, 128.0) * (1 - eps);
	      abs_a = fabs((double) temp1[0]);
	      abs_b = fabs((double) temp1[1]);
	      abs_c = fabs((double) T_element[0]);
	      abs_d = fabs((double) T_element[1]);
	      ab = MAX(abs_a, abs_b);
	      cd = MAX(abs_c, abs_d);

	      /* Scaling */
	      if (ab > ov / 16) {	/* scale down a, b */
		temp1[0] /= 16;
		temp1[1] /= 16;
		S = S * 16;
	      }
	      if (cd > ov / 16) {	/* scale down c, d */
		T_element[0] /= 16;
		T_element[1] /= 16;
		S = S / 16;
	      }
	      if (ab < un / eps * 2) {	/* scale up a, b */
		t = 2.0 / (eps * eps);
		temp1[0] *= t;
		temp1[1] *= t;
		S = S / t;
	      }
	      if (cd < un / eps * 2) {	/* scale up c, d */
		t = 2.0 / (eps * eps);
		T_element[0] *= t;
		T_element[1] *= t;
		S = S * t;
	      }

	      /* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
	      if (abs_c > abs_d) {
		r = T_element[1] / T_element[0];
		t = 1 / (T_element[0] + T_element[1] * r);
		q[0] = (temp1[0] + temp1[1] * r) * t;
		q[1] = (temp1[1] - temp1[0] * r) * t;
	      } else {
		r = T_element[0] / T_element[1];
		t = 1 / (T_element[1] + T_element[0] * r);
		q[0] = (temp1[1] + temp1[0] * r) * t;
		q[1] = (-temp1[0] + temp1[1] * r) * t;
	      }
	      /* Scale back */
	      temp1[0] = q[0] * S;
	      temp1[1] = q[1] * S;
	    }

	  }
	  /* if (diag == blas_non_unit_diag) */
	  x_i[jx] = temp1[0];
	  x_i[jx + 1] = temp1[1];
	  jx += incx;
	}			/* for j<n */
      } else {

	jx = start_x;
	for (j = 0; j < n; j++) {

	  /* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
	     i=j+1 to n-1           */
	  temp3[0] = x_i[jx];
	  temp3[1] = x_i[1 + jx];
	  /* multiply by alpha */
	  {
	    temp1[0] = temp3[0] * alpha_i[0] - temp3[1] * alpha_i[1];
	    temp1[1] = temp3[0] * alpha_i[1] + temp3[1] * alpha_i[0];
	  }


	  ix = start_x;
	  for (i = 0; i < j; i++) {
	    T_element[0] = T_i[j * incT + i * ldt * incT];
	    T_element[1] = T_i[j * incT + i * ldt * incT + 1];

	    temp3[0] = x_i[ix];
	    temp3[1] = x_i[1 + ix];
	    {
	      temp2[0] = temp3[0] * T_element[0] - temp3[1] * T_element[1];
	      temp2[1] = temp3[0] * T_element[1] + temp3[1] * T_element[0];
	    }

	    temp1[0] = temp1[0] - temp2[0];
	    temp1[1] = temp1[1] - temp2[1];
	    ix += incx;
	  }			/* for i<j */

	  /* if the diagonal entry is not equal to one, then divide Xj by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element[0] = T_i[j * incT + j * ldt * incT];
	    T_element[1] = T_i[j * incT + j * ldt * incT + 1];


	    {
	      double S = 1.0, eps, ov, un;
	      double abs_a, abs_b, abs_c, abs_d, ab, cd;
	      double r;
	      double t;
	      double q[2];

	      eps = pow(2.0, -24.0);
	      un = pow(2.0, -126.0);
	      ov = pow(2.0, 128.0) * (1 - eps);
	      abs_a = fabs((double) temp1[0]);
	      abs_b = fabs((double) temp1[1]);
	      abs_c = fabs((double) T_element[0]);
	      abs_d = fabs((double) T_element[1]);
	      ab = MAX(abs_a, abs_b);
	      cd = MAX(abs_c, abs_d);

	      /* Scaling */
	      if (ab > ov / 16) {	/* scale down a, b */
		temp1[0] /= 16;
		temp1[1] /= 16;
		S = S * 16;
	      }
	      if (cd > ov / 16) {	/* scale down c, d */
		T_element[0] /= 16;
		T_element[1] /= 16;
		S = S / 16;
	      }
	      if (ab < un / eps * 2) {	/* scale up a, b */
		t = 2.0 / (eps * eps);
		temp1[0] *= t;
		temp1[1] *= t;
		S = S / t;
	      }
	      if (cd < un / eps * 2) {	/* scale up c, d */
		t = 2.0 / (eps * eps);
		T_element[0] *= t;
		T_element[1] *= t;
		S = S * t;
	      }

	      /* Now un/eps*2 <= (a, b, c, d) >= ov/16 */
	      if (abs_c > abs_d) {
		r = T_element[1] / T_element[0];
		t = 1 / (T_element[0] + T_element[1] * r);
		q[0] = (temp1[0] + temp1[1] * r) * t;
		q[1] = (temp1[1] - temp1[0] * r) * t;
	      } else {
		r = T_element[0] / T_element[1];
		t = 1 / (T_element[1] + T_element[0] * r);
		q[0] = (temp1[1] + temp1[0] * r) * t;
		q[1] = (-temp1[0] + temp1[1] * r) * t;
	      }
	      /* Scale back */
	      temp1[0] = q[0] * S;
	      temp1[1] = q[1] * S;
	    }

	  }
	  /* if (diag == blas_non_unit_diag) */
	  x_i[jx] = temp1[0];
	  x_i[jx + 1] = temp1[1];
	  jx += incx;
	}			/* for j<n */
      }
    }
  }
}
