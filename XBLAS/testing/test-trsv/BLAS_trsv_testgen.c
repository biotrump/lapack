#include "blas_extended.h"
#include "cblas_test.h"


extern double xrand(int *);








void BLAS_strsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, float *alpha,
			int alpha_flag, float *T, int lda, float *x,
			int *seed, double *r_true_l, double *r_true_t,
			int row, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) float*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) float*
 *
 * x            (input/output) float*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  int start;
  int length;
  int i, j;
  float alpha_i;
  float minus_one;
  float Tii;
  float *temp;
  float *xtemp2;

  temp = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }

  xtemp2 = NULL;
  if (prec != blas_prec_extra) {
    xtemp2 = (float *) blas_malloc(n * sizeof(float));
    if (n > 0 && xtemp2 == NULL) {
      printf("malloc failed\n");
      exit(-1);
    }
  }

  minus_one = -1.0;

  /* if alpha_flag=0, gives a random value to alpha */
  if (alpha_flag == 0) {
    alpha_i = xrand(seed);
    *alpha = alpha_i;
    alpha_flag = 1;
  }

  for (i = 0; i < 4 * n * n; i++) {
    T[i] = 0.0;
  }

  for (i = 0; i < n; i++) {

    if (i != row) {
      if (diag == blas_non_unit_diag) {
	Tii = xrand(seed);
	T[i * lda + i] = Tii;
      } else {
	Tii = 1.0;
	T[i * lda + i] = 1.0;
      }

      x[i] = xrand(seed);

      switch (prec) {
      case blas_prec_single:
	{
	  float multemp;
	  float divtemp;

	  multemp = x[i] * *alpha;
	  divtemp = multemp / Tii;
	  r_true_l[i] = divtemp;
	  r_true_t[i] = 0.0;
	  xtemp2[i] = divtemp;
	  break;
	}
      case blas_prec_double:
      case blas_prec_indigenous:
	{
	  double multemp;
	  double divtemp;

	  multemp = (double) x[i] * *alpha;
	  divtemp = (double) multemp / Tii;
	  r_true_l[i] = divtemp;
	  r_true_t[i] = 0.0;
	  break;
	}
      case blas_prec_extra:
	{
	  double head_multemp, tail_multemp;
	  double head_divtemp, tail_divtemp;

	  head_multemp = (double) x[i] * *alpha;
	  tail_multemp = 0.0;
	  {
	    double dt = (double) Tii;
	    {
	      /* Compute double-double = double-double / double,
	         using a Newton iteration scheme. */
	      double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

	      /* Compute a DP approximation to the quotient. */
	      t1 = head_multemp / dt;

	      /* Split t1 and b into two parts with at most 26 bits each,
	         using the Dekker-Veltkamp method. */
	      con = t1 * split;
	      t11 = con - (con - t1);
	      t21 = t1 - t11;
	      con = dt * split;
	      b1 = con - (con - dt);
	      b2 = dt - b1;

	      /* Compute t1 * b using Dekker method. */
	      t12 = t1 * dt;
	      t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

	      /* Compute dda - (t12, t22) using Knuth trick. */
	      t11 = head_multemp - t12;
	      e = t11 - head_multemp;
	      t21 =
		((-t12 - e) + (head_multemp - (t11 - e))) + tail_multemp -
		t22;

	      /* Compute high-order word of (t11, t21) and divide by b. */
	      t2 = (t11 + t21) / dt;

	      /* The result is t1 + t2, after normalization. */
	      head_divtemp = t1 + t2;
	      tail_divtemp = t2 - (head_divtemp - t1);
	    }
	  }
	  r_true_l[i] = head_divtemp;
	  r_true_t[i] = tail_divtemp;
	  break;
	}
      }				/* case */
    }				/* if */
  }				/* for */

  for (j = 0; j < n; j++) {
    temp[j] = 0.0;
  }

  T[row * lda + row] = 1.0;

  if ((uplo == blas_lower && trans == blas_no_trans) ||
      (uplo == blas_upper && trans != blas_no_trans)) {
    length = row;
    start = 0;
  } else {
    length = n - row - 1;
    start = row + 1;
  }

  if (length != 0) {

    switch (prec) {
    case blas_prec_single:
      BLAS_sdot_testgen(length, 0, length, norm,
			blas_no_conj, &minus_one, 1, alpha, 1,
			&xtemp2[start], temp, seed, &x[row],
			&r_true_l[row], &r_true_t[row]);
      break;
    case blas_prec_double:
    case blas_prec_indigenous:
    case blas_prec_extra:
      BLAS_sdot_x_testgen(length, 0, length, norm,
			  blas_no_conj, &minus_one, 1, alpha, 1,
			  &r_true_l[start], &r_true_t[start], temp,
			  seed, &x[row], &r_true_l[row], &r_true_t[row]);
      break;
    }
    strsv_commit(order, uplo, trans, length, T, lda, temp, row);
  } else {
    x[row] = xrand(seed);

    switch (prec) {
    case blas_prec_single:
      {
	float multemp;

	multemp = x[row] * *alpha;
	r_true_l[row] = multemp;
	r_true_t[row] = 0.0;
	break;
      }
    case blas_prec_indigenous:
    case blas_prec_double:
      {
	double multemp;

	multemp = (double) x[row] * *alpha;
	r_true_l[row] = multemp;
	r_true_t[row] = 0.0;
	break;
      }
    case blas_prec_extra:
      {
	double head_multemp, tail_multemp;

	head_multemp = (double) x[row] * *alpha;
	tail_multemp = 0.0;
	r_true_l[row] = head_multemp;
	r_true_t[row] = tail_multemp;
	break;
      }
    }
  }

  blas_free(temp);

  if (prec != blas_prec_extra)
    blas_free(xtemp2);
}				/* end of BLAS_strsv_testgen */

void BLAS_dtrsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, double *alpha,
			int alpha_flag, double *T, int lda, double *x,
			int *seed, double *r_true_l, double *r_true_t,
			int row, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) double*
 *
 * x            (input/output) double*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  int start;
  int length;
  int i, j;
  float alpha_i;
  double minus_one;
  double Tii;
  double *temp;
  double *xtemp2;

  temp = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }

  xtemp2 = NULL;
  if (prec != blas_prec_extra) {
    xtemp2 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && xtemp2 == NULL) {
      printf("malloc failed\n");
      exit(-1);
    }
  }

  minus_one = -1.0;

  /* if alpha_flag=0, gives a random value to alpha */
  if (alpha_flag == 0) {
    alpha_i = xrand(seed);
    *alpha = alpha_i;
    alpha_flag = 1;
  }

  for (i = 0; i < 4 * n * n; i++) {
    T[i] = 0.0;
  }

  for (i = 0; i < n; i++) {

    if (i != row) {
      if (diag == blas_non_unit_diag) {
	Tii = xrand(seed);
	T[i * lda + i] = Tii;
      } else {
	Tii = 1.0;
	T[i * lda + i] = 1.0;
      }

      x[i] = xrand(seed);

      switch (prec) {
      case blas_prec_single:
	{
	  double multemp;
	  double divtemp;

	  multemp = x[i] * *alpha;
	  divtemp = multemp / Tii;
	  r_true_l[i] = divtemp;
	  r_true_t[i] = 0.0;
	  xtemp2[i] = divtemp;
	  break;
	}
      case blas_prec_double:
      case blas_prec_indigenous:
	{
	  double multemp;
	  double divtemp;

	  multemp = x[i] * *alpha;
	  divtemp = multemp / Tii;
	  r_true_l[i] = divtemp;
	  r_true_t[i] = 0.0;
	  xtemp2[i] = divtemp;
	  break;
	}
      case blas_prec_extra:
	{
	  double head_multemp, tail_multemp;
	  double head_divtemp, tail_divtemp;

	  {
	    /* Compute double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = x[i] * split;
	    a1 = con - x[i];
	    a1 = con - a1;
	    a2 = x[i] - a1;
	    con = *alpha * split;
	    b1 = con - *alpha;
	    b1 = con - b1;
	    b2 = *alpha - b1;

	    head_multemp = x[i] * *alpha;
	    tail_multemp =
	      (((a1 * b1 - head_multemp) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute double-double = double-double / double,
	       using a Newton iteration scheme. */
	    double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

	    /* Compute a DP approximation to the quotient. */
	    t1 = head_multemp / Tii;

	    /* Split t1 and b into two parts with at most 26 bits each,
	       using the Dekker-Veltkamp method. */
	    con = t1 * split;
	    t11 = con - (con - t1);
	    t21 = t1 - t11;
	    con = Tii * split;
	    b1 = con - (con - Tii);
	    b2 = Tii - b1;

	    /* Compute t1 * b using Dekker method. */
	    t12 = t1 * Tii;
	    t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

	    /* Compute dda - (t12, t22) using Knuth trick. */
	    t11 = head_multemp - t12;
	    e = t11 - head_multemp;
	    t21 =
	      ((-t12 - e) + (head_multemp - (t11 - e))) + tail_multemp - t22;

	    /* Compute high-order word of (t11, t21) and divide by b. */
	    t2 = (t11 + t21) / Tii;

	    /* The result is t1 + t2, after normalization. */
	    head_divtemp = t1 + t2;
	    tail_divtemp = t2 - (head_divtemp - t1);
	  }
	  r_true_l[i] = head_divtemp;
	  r_true_t[i] = tail_divtemp;
	  break;
	}
      }				/* case */
    }				/* if */
  }				/* for */

  for (j = 0; j < n; j++) {
    temp[j] = 0.0;
  }

  T[row * lda + row] = 1.0;

  if ((uplo == blas_lower && trans == blas_no_trans) ||
      (uplo == blas_upper && trans != blas_no_trans)) {
    length = row;
    start = 0;
  } else {
    length = n - row - 1;
    start = row + 1;
  }

  if (length != 0) {

    switch (prec) {
    case blas_prec_single:
    case blas_prec_double:
    case blas_prec_indigenous:
      BLAS_ddot_testgen(length, 0, length, norm,
			blas_no_conj, &minus_one, 1, alpha, 1,
			&xtemp2[start], temp, seed, &x[row],
			&r_true_l[row], &r_true_t[row]);
      break;
    case blas_prec_extra:
      BLAS_ddot_x_testgen(length, 0, length, norm,
			  blas_no_conj, &minus_one, 1, alpha, 1,
			  &r_true_l[start], &r_true_t[start], temp,
			  seed, &x[row], &r_true_l[row], &r_true_t[row]);
      break;

    }
    dtrsv_commit(order, uplo, trans, length, T, lda, temp, row);
  } else {
    x[row] = xrand(seed);

    switch (prec) {
    case blas_prec_single:
      {
	double multemp;

	multemp = x[row] * *alpha;
	r_true_l[row] = multemp;
	r_true_t[row] = 0.0;
	break;
      }
    case blas_prec_indigenous:
    case blas_prec_double:
      {
	double multemp;

	multemp = x[row] * *alpha;
	r_true_l[row] = multemp;
	r_true_t[row] = 0.0;
	break;
      }
    case blas_prec_extra:
      {
	double head_multemp, tail_multemp;

	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x[row] * split;
	  a1 = con - x[row];
	  a1 = con - a1;
	  a2 = x[row] - a1;
	  con = *alpha * split;
	  b1 = con - *alpha;
	  b1 = con - b1;
	  b2 = *alpha - b1;

	  head_multemp = x[row] * *alpha;
	  tail_multemp =
	    (((a1 * b1 - head_multemp) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	r_true_l[row] = head_multemp;
	r_true_t[row] = tail_multemp;
	break;
      }
    }
  }

  blas_free(temp);

  if (prec != blas_prec_extra)
    blas_free(xtemp2);
}				/* end of BLAS_dtrsv_testgen */

void BLAS_dtrsv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, double *alpha,
			  int alpha_flag, float *T, int lda, double *x,
			  int *seed, double *r_true_l, double *r_true_t,
			  int row, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) float*
 *
 * x            (input/output) double*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  int start;
  int length;
  int i, j;
  float alpha_i;
  double minus_one;
  float Tii;
  float *temp;
  double *xtemp2;

  temp = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }

  xtemp2 = NULL;
  if (prec != blas_prec_extra) {
    xtemp2 = (double *) blas_malloc(n * sizeof(double));
    if (n > 0 && xtemp2 == NULL) {
      printf("malloc failed\n");
      exit(-1);
    }
  }

  minus_one = -1.0;

  /* if alpha_flag=0, gives a random value to alpha */
  if (alpha_flag == 0) {
    alpha_i = xrand(seed);
    *alpha = alpha_i;
    alpha_flag = 1;
  }

  for (i = 0; i < 4 * n * n; i++) {
    T[i] = 0.0;
  }

  for (i = 0; i < n; i++) {

    if (i != row) {
      if (diag == blas_non_unit_diag) {
	Tii = xrand(seed);
	T[i * lda + i] = Tii;
      } else {
	Tii = 1.0;
	T[i * lda + i] = 1.0;
      }

      x[i] = xrand(seed);

      switch (prec) {
      case blas_prec_single:
	{
	  double multemp;
	  double divtemp;

	  multemp = x[i] * *alpha;
	  divtemp = multemp / Tii;
	  r_true_l[i] = divtemp;
	  r_true_t[i] = 0.0;
	  xtemp2[i] = divtemp;
	  break;
	}
      case blas_prec_double:
      case blas_prec_indigenous:
	{
	  double multemp;
	  double divtemp;

	  multemp = x[i] * *alpha;
	  divtemp = multemp / Tii;
	  r_true_l[i] = divtemp;
	  r_true_t[i] = 0.0;
	  xtemp2[i] = divtemp;
	  break;
	}
      case blas_prec_extra:
	{
	  double head_multemp, tail_multemp;
	  double head_divtemp, tail_divtemp;

	  {
	    /* Compute double_double = double * double. */
	    double a1, a2, b1, b2, con;

	    con = x[i] * split;
	    a1 = con - x[i];
	    a1 = con - a1;
	    a2 = x[i] - a1;
	    con = *alpha * split;
	    b1 = con - *alpha;
	    b1 = con - b1;
	    b2 = *alpha - b1;

	    head_multemp = x[i] * *alpha;
	    tail_multemp =
	      (((a1 * b1 - head_multemp) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    double dt = (double) Tii;
	    {
	      /* Compute double-double = double-double / double,
	         using a Newton iteration scheme. */
	      double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

	      /* Compute a DP approximation to the quotient. */
	      t1 = head_multemp / dt;

	      /* Split t1 and b into two parts with at most 26 bits each,
	         using the Dekker-Veltkamp method. */
	      con = t1 * split;
	      t11 = con - (con - t1);
	      t21 = t1 - t11;
	      con = dt * split;
	      b1 = con - (con - dt);
	      b2 = dt - b1;

	      /* Compute t1 * b using Dekker method. */
	      t12 = t1 * dt;
	      t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

	      /* Compute dda - (t12, t22) using Knuth trick. */
	      t11 = head_multemp - t12;
	      e = t11 - head_multemp;
	      t21 =
		((-t12 - e) + (head_multemp - (t11 - e))) + tail_multemp -
		t22;

	      /* Compute high-order word of (t11, t21) and divide by b. */
	      t2 = (t11 + t21) / dt;

	      /* The result is t1 + t2, after normalization. */
	      head_divtemp = t1 + t2;
	      tail_divtemp = t2 - (head_divtemp - t1);
	    }
	  }
	  r_true_l[i] = head_divtemp;
	  r_true_t[i] = tail_divtemp;
	  break;
	}
      }				/* case */
    }				/* if */
  }				/* for */

  for (j = 0; j < n; j++) {
    temp[j] = 0.0;
  }

  T[row * lda + row] = 1.0;

  if ((uplo == blas_lower && trans == blas_no_trans) ||
      (uplo == blas_upper && trans != blas_no_trans)) {
    length = row;
    start = 0;
  } else {
    length = n - row - 1;
    start = row + 1;
  }

  if (length != 0) {

    switch (prec) {
    case blas_prec_single:
    case blas_prec_double:
    case blas_prec_indigenous:
      /*BLAS_ddot_s_x_testgen(length, 0, length, norm, 
         blas_no_conj, &minus_one, 1, alpha, 1, 
         &r_true_l[start], &r_true_t[start], temp, 
         seed, &x[row], &r_true_l[row], &r_true_t[row]);
         break; */
    case blas_prec_extra:
      BLAS_ddot_s_x_testgen(length, 0, length, norm,
			    blas_no_conj, &minus_one, 1, alpha, 1,
			    &r_true_l[start], &r_true_t[start], temp,
			    seed, &x[row], &r_true_l[row], &r_true_t[row]);
      break;
    }
    strsv_commit(order, uplo, trans, length, T, lda, temp, row);
  } else {
    x[row] = xrand(seed);

    switch (prec) {
    case blas_prec_single:
      {
	double multemp;

	multemp = x[row] * *alpha;
	r_true_l[row] = multemp;
	r_true_t[row] = 0.0;
	break;
      }
    case blas_prec_indigenous:
    case blas_prec_double:
      {
	double multemp;

	multemp = x[row] * *alpha;
	r_true_l[row] = multemp;
	r_true_t[row] = 0.0;
	break;
      }
    case blas_prec_extra:
      {
	double head_multemp, tail_multemp;

	{
	  /* Compute double_double = double * double. */
	  double a1, a2, b1, b2, con;

	  con = x[row] * split;
	  a1 = con - x[row];
	  a1 = con - a1;
	  a2 = x[row] - a1;
	  con = *alpha * split;
	  b1 = con - *alpha;
	  b1 = con - b1;
	  b2 = *alpha - b1;

	  head_multemp = x[row] * *alpha;
	  tail_multemp =
	    (((a1 * b1 - head_multemp) + a1 * b2) + a2 * b1) + a2 * b2;
	}
	r_true_l[row] = head_multemp;
	r_true_t[row] = tail_multemp;
	break;
      }
    }
  }

  blas_free(temp);

  if (prec != blas_prec_extra)
    blas_free(xtemp2);
}				/* end of BLAS_dtrsv_s_testgen */

void BLAS_ctrsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, void *alpha,
			int alpha_flag, void *T, int lda, void *x,
			int *seed, double *r_true_l, double *r_true_t,
			int row, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) void*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  float *x_i = (float *) x;
  float *alpha_i = (float *) alpha;
  float *T_i = (float *) T;
  float alpha_r;
  float *T_r;
  float *x_r;
  double *r_true_l_r, *r_true_t_r;
  int i, inc = 2, length;

  T_r = (float *) blas_malloc(4 * n * n * sizeof(float));
  if (4 * n * n > 0 && T_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  x_r = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && x_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  r_true_l_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && r_true_l_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  r_true_t_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && r_true_t_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }

  if (alpha_flag == 1) {
    alpha_r = alpha_i[0];
  }

  if ((uplo == blas_lower && trans == blas_no_trans) ||
      (uplo == blas_upper && trans != blas_no_trans)) {
    length = row;
  } else {
    length = n - row - 1;
  }

  BLAS_strsv_testgen(norm, order, uplo, trans, diag, n, &alpha_r,
		     alpha_flag, T_r, lda, x_r, seed, r_true_l_r, r_true_t_r,
		     row, prec);

  alpha_i[0] = alpha_r;
  alpha_i[1] = alpha_r;

  if (diag == blas_non_unit_diag) {
    for (i = 0; i < n; i++) {
      x_i[i * inc] = 0.0;
      x_i[i * inc + 1] = x_r[i];

      if (i != row) {
	r_true_l[i * inc] = 0.0;
	r_true_l[i * inc + 1] = r_true_l_r[i];
	r_true_t[i * inc] = 0.0;
	r_true_t[i * inc + 1] = r_true_t_r[i];
      } else {
	r_true_l[i * inc] = -r_true_l_r[i];
	r_true_l[i * inc + 1] = r_true_l_r[i];
	r_true_t[i * inc] = -r_true_t_r[i];
	r_true_t[i * inc + 1] = r_true_t_r[i];
      }
    }

    for (i = 0; i < 4 * n * n; i++) {
      T_i[i * inc] = T_r[i];

      if (trans != blas_conj_trans)
	T_i[i * inc + 1] = T_r[i];
      else
	T_i[i * inc + 1] = -T_r[i];
    }

    T_i[(row + lda * row) * inc + 1] = 0.0;
  } else {
    for (i = 0; i < n; i++) {
      x_i[i * inc] = 0.0;
      x_i[i * inc + 1] = x_r[i];

      if (i != row || length == 0) {
	r_true_l[i * inc] = -r_true_l_r[i];
	r_true_l[i * inc + 1] = r_true_l_r[i];
	r_true_t[i * inc] = -r_true_t_r[i];
	r_true_t[i * inc + 1] = r_true_t_r[i];
      } else {
	x_i[i * inc] = x_r[i];
	x_i[i * inc + 1] = x_r[i];

	r_true_l[i * inc] = 0.0;
	r_true_l[i * inc + 1] = 2 * r_true_l_r[i];
	r_true_t[i * inc] = 0.0;
	r_true_t[i * inc + 1] = 2 * r_true_t_r[i];
      }
    }

    for (i = 0; i < 4 * n * n; i++) {
      T_i[i * inc] = T_r[i];

      if (trans != blas_conj_trans)
	T_i[i * inc + 1] = -T_r[i];
      else
	T_i[i * inc + 1] = T_r[i];
    }

    for (i = 0; i < n; i++) {
      T_i[(i + lda * i) * inc + 1] = 0.0;
    }
  }

  blas_free(T_r);
  blas_free(x_r);
  blas_free(r_true_l_r);
  blas_free(r_true_t_r);
}				/* end of BLAS_ctrsv_testgen */

void BLAS_ztrsv_c_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, void *T, int lda, void *x,
			  int *seed, double *r_true_l, double *r_true_t,
			  int row, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) void*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  double *x_i = (double *) x;
  double *alpha_i = (double *) alpha;
  float *T_i = (float *) T;
  double alpha_r;
  float *T_r;
  double *x_r;
  double *r_true_l_r, *r_true_t_r;
  int i, inc = 2, length;

  T_r = (float *) blas_malloc(4 * n * n * sizeof(float));
  if (4 * n * n > 0 && T_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  x_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && x_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  r_true_l_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && r_true_l_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  r_true_t_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && r_true_t_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }

  if (alpha_flag == 1) {
    alpha_r = alpha_i[0];
  }

  if ((uplo == blas_lower && trans == blas_no_trans) ||
      (uplo == blas_upper && trans != blas_no_trans)) {
    length = row;
  } else {
    length = n - row - 1;
  }

  BLAS_dtrsv_s_testgen(norm, order, uplo, trans, diag, n, &alpha_r,
		       alpha_flag, T_r, lda, x_r, seed, r_true_l_r,
		       r_true_t_r, row, prec);

  alpha_i[0] = alpha_r;
  alpha_i[1] = alpha_r;

  if (diag == blas_non_unit_diag) {
    for (i = 0; i < n; i++) {
      x_i[i * inc] = 0.0;
      x_i[i * inc + 1] = x_r[i];

      if (i != row) {
	r_true_l[i * inc] = 0.0;
	r_true_l[i * inc + 1] = r_true_l_r[i];
	r_true_t[i * inc] = 0.0;
	r_true_t[i * inc + 1] = r_true_t_r[i];
      } else {
	r_true_l[i * inc] = -r_true_l_r[i];
	r_true_l[i * inc + 1] = r_true_l_r[i];
	r_true_t[i * inc] = -r_true_t_r[i];
	r_true_t[i * inc + 1] = r_true_t_r[i];
      }
    }

    for (i = 0; i < 4 * n * n; i++) {
      T_i[i * inc] = T_r[i];

      if (trans != blas_conj_trans)
	T_i[i * inc + 1] = T_r[i];
      else
	T_i[i * inc + 1] = -T_r[i];
    }

    T_i[(row + lda * row) * inc + 1] = 0.0;
  } else {
    for (i = 0; i < n; i++) {
      x_i[i * inc] = 0.0;
      x_i[i * inc + 1] = x_r[i];

      if (i != row || length == 0) {
	r_true_l[i * inc] = -r_true_l_r[i];
	r_true_l[i * inc + 1] = r_true_l_r[i];
	r_true_t[i * inc] = -r_true_t_r[i];
	r_true_t[i * inc + 1] = r_true_t_r[i];
      } else {
	x_i[i * inc] = x_r[i];
	x_i[i * inc + 1] = x_r[i];

	r_true_l[i * inc] = 0.0;
	r_true_l[i * inc + 1] = 2 * r_true_l_r[i];
	r_true_t[i * inc] = 0.0;
	r_true_t[i * inc + 1] = 2 * r_true_t_r[i];
      }
    }

    for (i = 0; i < 4 * n * n; i++) {
      T_i[i * inc] = T_r[i];

      if (trans != blas_conj_trans)
	T_i[i * inc + 1] = -T_r[i];
      else
	T_i[i * inc + 1] = T_r[i];
    }

    for (i = 0; i < n; i++) {
      T_i[(i + lda * i) * inc + 1] = 0.0;
    }
  }

  blas_free(T_r);
  blas_free(x_r);
  blas_free(r_true_l_r);
  blas_free(r_true_t_r);
}				/* end of BLAS_ztrsv_c_testgen */

void BLAS_ztrsv_testgen(int norm, enum blas_order_type order,
			enum blas_uplo_type uplo, enum blas_trans_type trans,
			enum blas_diag_type diag, int n, void *alpha,
			int alpha_flag, void *T, int lda, void *x,
			int *seed, double *r_true_l, double *r_true_t,
			int row, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) void*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  double *x_i = (double *) x;
  double *alpha_i = (double *) alpha;
  double *T_i = (double *) T;
  double alpha_r;
  double *T_r;
  double *x_r;
  double *r_true_l_r, *r_true_t_r;
  int i, inc = 2, length;

  T_r = (double *) blas_malloc(4 * n * n * sizeof(double));
  if (4 * n * n > 0 && T_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  x_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && x_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  r_true_l_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && r_true_l_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  r_true_t_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && r_true_t_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }

  if (alpha_flag == 1) {
    alpha_r = alpha_i[0];
  }

  if ((uplo == blas_lower && trans == blas_no_trans) ||
      (uplo == blas_upper && trans != blas_no_trans)) {
    length = row;
  } else {
    length = n - row - 1;
  }

  BLAS_dtrsv_testgen(norm, order, uplo, trans, diag, n, &alpha_r,
		     alpha_flag, T_r, lda, x_r, seed, r_true_l_r, r_true_t_r,
		     row, prec);

  alpha_i[0] = alpha_r;
  alpha_i[1] = alpha_r;

  if (diag == blas_non_unit_diag) {
    for (i = 0; i < n; i++) {
      x_i[i * inc] = 0.0;
      x_i[i * inc + 1] = x_r[i];

      if (i != row) {
	r_true_l[i * inc] = 0.0;
	r_true_l[i * inc + 1] = r_true_l_r[i];
	r_true_t[i * inc] = 0.0;
	r_true_t[i * inc + 1] = r_true_t_r[i];
      } else {
	r_true_l[i * inc] = -r_true_l_r[i];
	r_true_l[i * inc + 1] = r_true_l_r[i];
	r_true_t[i * inc] = -r_true_t_r[i];
	r_true_t[i * inc + 1] = r_true_t_r[i];
      }
    }

    for (i = 0; i < 4 * n * n; i++) {
      T_i[i * inc] = T_r[i];

      if (trans != blas_conj_trans)
	T_i[i * inc + 1] = T_r[i];
      else
	T_i[i * inc + 1] = -T_r[i];
    }

    T_i[(row + lda * row) * inc + 1] = 0.0;
  } else {
    for (i = 0; i < n; i++) {
      x_i[i * inc] = 0.0;
      x_i[i * inc + 1] = x_r[i];

      if (i != row || length == 0) {
	r_true_l[i * inc] = -r_true_l_r[i];
	r_true_l[i * inc + 1] = r_true_l_r[i];
	r_true_t[i * inc] = -r_true_t_r[i];
	r_true_t[i * inc + 1] = r_true_t_r[i];
      } else {
	x_i[i * inc] = x_r[i];
	x_i[i * inc + 1] = x_r[i];

	r_true_l[i * inc] = 0.0;
	r_true_l[i * inc + 1] = 2 * r_true_l_r[i];
	r_true_t[i * inc] = 0.0;
	r_true_t[i * inc + 1] = 2 * r_true_t_r[i];
      }
    }

    for (i = 0; i < 4 * n * n; i++) {
      T_i[i * inc] = T_r[i];

      if (trans != blas_conj_trans)
	T_i[i * inc + 1] = -T_r[i];
      else
	T_i[i * inc + 1] = T_r[i];
    }

    for (i = 0; i < n; i++) {
      T_i[(i + lda * i) * inc + 1] = 0.0;
    }
  }

  blas_free(T_r);
  blas_free(x_r);
  blas_free(r_true_l_r);
  blas_free(r_true_t_r);
}				/* end of BLAS_ztrsv_testgen */

void BLAS_ctrsv_s_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, float *T, int lda, void *x,
			  int *seed, double *r_true_l, double *r_true_t,
			  int row, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) float*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  float *x_i = (float *) x;
  float *alpha_i = (float *) alpha;
  float *T_i = T;
  float alpha_r;
  float *x_r;
  double *r_true_l_r, *r_true_t_r;
  int i, inc = 2;

  x_r = (float *) blas_malloc(n * sizeof(float));
  if (n > 0 && x_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  r_true_l_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && r_true_l_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  r_true_t_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && r_true_t_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }

  if (alpha_flag == 1) {
    alpha_r = alpha_i[0];
  }

  BLAS_strsv_testgen(norm, order, uplo, trans, diag, n, &alpha_r,
		     alpha_flag, T_i, lda, x_r, seed, r_true_l_r, r_true_t_r,
		     row, prec);

  alpha_i[0] = alpha_r;
  alpha_i[1] = alpha_r;

  for (i = 0; i < n; i++) {
    x_i[i * inc] = 0.0;
    x_i[i * inc + 1] = x_r[i];

    r_true_l[i * inc] = -r_true_l_r[i];
    r_true_l[i * inc + 1] = r_true_l_r[i];
    r_true_t[i * inc] = -r_true_t_r[i];
    r_true_t[i * inc + 1] = r_true_t_r[i];
  }

  blas_free(x_r);
  blas_free(r_true_l_r);
  blas_free(r_true_t_r);
}				/* end of BLAS_ctrsv_s_testgen */

void BLAS_ztrsv_d_testgen(int norm, enum blas_order_type order,
			  enum blas_uplo_type uplo,
			  enum blas_trans_type trans,
			  enum blas_diag_type diag, int n, void *alpha,
			  int alpha_flag, double *T, int lda, void *x,
			  int *seed, double *r_true_l, double *r_true_t,
			  int row, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Generates alpha, x and T, where T is a triangular matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              No trans, trans, conj trans
 *
 * diag         (input) blas_diag_type
 *              non unit, unit
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * T            (output) double*
 *
 * x            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 * row          (input) int
 *              The true row being generated
 *
 * prec         (input) blas_prec_type
 *              single, double, or extra precision   
 *
 */
{
  double *x_i = (double *) x;
  double *alpha_i = (double *) alpha;
  double *T_i = T;
  double alpha_r;
  double *x_r;
  double *r_true_l_r, *r_true_t_r;
  int i, inc = 2;

  x_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && x_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  r_true_l_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && r_true_l_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  r_true_t_r = (double *) blas_malloc(n * sizeof(double));
  if (n > 0 && r_true_t_r == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }

  if (alpha_flag == 1) {
    alpha_r = alpha_i[0];
  }

  BLAS_dtrsv_testgen(norm, order, uplo, trans, diag, n, &alpha_r,
		     alpha_flag, T_i, lda, x_r, seed, r_true_l_r, r_true_t_r,
		     row, prec);

  alpha_i[0] = alpha_r;
  alpha_i[1] = alpha_r;

  for (i = 0; i < n; i++) {
    x_i[i * inc] = 0.0;
    x_i[i * inc + 1] = x_r[i];

    r_true_l[i * inc] = -r_true_l_r[i];
    r_true_l[i * inc + 1] = r_true_l_r[i];
    r_true_t[i * inc] = -r_true_t_r[i];
    r_true_t[i * inc + 1] = r_true_t_r[i];
  }

  blas_free(x_r);
  blas_free(r_true_l_r);
  blas_free(r_true_t_r);
}				/* end of BLAS_ztrsv_d_testgen */