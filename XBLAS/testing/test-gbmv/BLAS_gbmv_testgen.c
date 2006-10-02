#include "blas_extended.h"
#include "cblas_test.h"





void BLAS_sgbmv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n, int kl,
			int ku, float *alpha, int alpha_flag, float *AB,
			int lda, float *x, float *beta, int beta_flag,
			float *y, int *seed, double *r_true_l,
			double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) float*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) float*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) float*
 *
 * beta         (input/output) float*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) float*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  float *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;



  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    AB[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    sgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_sdot_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
		      alpha_flag, beta, beta_flag, x, temp, seed,
		      &y[i * incy], &((double *) r_true_l)[i * incy],
		      &((double *) r_true_t)[i * incy]);



    /* copy temp to AB */
    sgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_sgbmv_testgen */

void BLAS_dgbmv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n, int kl,
			int ku, double *alpha, int alpha_flag, double *AB,
			int lda, double *x, double *beta, int beta_flag,
			double *y, int *seed, double *r_true_l,
			double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) double*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) double*
 *
 * beta         (input/output) double*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) double*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  double *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;



  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    AB[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    dgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
		      alpha_flag, beta, beta_flag, x, temp, seed,
		      &y[i * incy], &((double *) r_true_l)[i * incy],
		      &((double *) r_true_t)[i * incy]);



    /* copy temp to AB */
    dgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_dgbmv_testgen */

void BLAS_dgbmv_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, double *alpha, int alpha_flag, double *AB,
			    int lda, float *x, double *beta, int beta_flag,
			    double *y, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) double*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) float*
 *
 * beta         (input/output) double*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) double*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  double *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;



  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    AB[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    dgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot_s_d_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed,
			  &y[i * incy], &((double *) r_true_l)[i * incy],
			  &((double *) r_true_t)[i * incy]);



    /* copy temp to AB */
    dgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_dgbmv_d_s_testgen */

void BLAS_dgbmv_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, double *alpha, int alpha_flag, float *AB,
			    int lda, double *x, double *beta, int beta_flag,
			    double *y, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) float*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) double*
 *
 * beta         (input/output) double*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) double*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  float *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;



  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    AB[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    sgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot_d_s_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed,
			  &y[i * incy], &((double *) r_true_l)[i * incy],
			  &((double *) r_true_t)[i * incy]);



    /* copy temp to AB */
    sgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_dgbmv_s_d_testgen */

void BLAS_dgbmv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, double *alpha, int alpha_flag, float *AB,
			    int lda, float *x, double *beta, int beta_flag,
			    double *y, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) float*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) float*
 *
 * beta         (input/output) double*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) double*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  float *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;



  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    AB[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    sgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot_s_s_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed,
			  &y[i * incy], &((double *) r_true_l)[i * incy],
			  &((double *) r_true_t)[i * incy]);



    /* copy temp to AB */
    sgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_dgbmv_s_s_testgen */


void BLAS_cgbmv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n, int kl,
			int ku, void *alpha, int alpha_flag, void *AB,
			int lda, void *x, void *beta, int beta_flag, void *y,
			int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) void*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) void*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  int j;
  float *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;
  incy *= 2;
  incAB *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    ((float *) AB)[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    cgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
		      alpha_flag, beta, beta_flag, x, temp, seed,
		      &((float *) y)[i * incy],
		      &((double *) r_true_l)[i * incy],
		      &((double *) r_true_t)[i * incy]);

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incAB; j += 2) {
	temp[j + 1] = -temp[j + 1];
      }
    }

    /* copy temp to AB */
    cgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_cgbmv_testgen */

void BLAS_cgbmv_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, float *AB,
			    int lda, float *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) float*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) float*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  float *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;
  incy *= 2;


  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    AB[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    sgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot_s_s_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed,
			  &((float *) y)[i * incy],
			  &((double *) r_true_l)[i * incy],
			  &((double *) r_true_t)[i * incy]);



    /* copy temp to AB */
    sgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_cgbmv_s_s_testgen */

void BLAS_cgbmv_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, void *AB,
			    int lda, float *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) void*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) float*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  int j;
  float *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;
  incy *= 2;
  incAB *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    ((float *) AB)[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    cgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot_s_c_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed,
			  &((float *) y)[i * incy],
			  &((double *) r_true_l)[i * incy],
			  &((double *) r_true_t)[i * incy]);

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incAB; j += 2) {
	temp[j + 1] = -temp[j + 1];
      }
    }

    /* copy temp to AB */
    cgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_cgbmv_c_s_testgen */

void BLAS_cgbmv_s_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, float *AB,
			    int lda, void *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) float*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) void*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  float *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;
  incy *= 2;


  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    AB[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    sgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot_c_s_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed,
			  &((float *) y)[i * incy],
			  &((double *) r_true_l)[i * incy],
			  &((double *) r_true_t)[i * incy]);



    /* copy temp to AB */
    sgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_cgbmv_s_c_testgen */


void BLAS_zgbmv_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type trans, int m, int n, int kl,
			int ku, void *alpha, int alpha_flag, void *AB,
			int lda, void *x, void *beta, int beta_flag, void *y,
			int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) void*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) void*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  int j;
  double *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;
  incy *= 2;
  incAB *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    ((double *) AB)[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    zgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
		      alpha_flag, beta, beta_flag, x, temp, seed,
		      &((double *) y)[i * incy],
		      &((double *) r_true_l)[i * incy],
		      &((double *) r_true_t)[i * incy]);

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incAB; j += 2) {
	temp[j + 1] = -temp[j + 1];
      }
    }

    /* copy temp to AB */
    zgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_zgbmv_testgen */

void BLAS_zgbmv_d_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, double *AB,
			    int lda, double *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) double*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) double*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  double *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;
  incy *= 2;


  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    AB[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    dgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot_d_d_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed,
			  &((double *) y)[i * incy],
			  &((double *) r_true_l)[i * incy],
			  &((double *) r_true_t)[i * incy]);



    /* copy temp to AB */
    dgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_zgbmv_d_d_testgen */

void BLAS_zgbmv_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, void *AB,
			    int lda, double *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) void*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) double*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  int j;
  double *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;
  incy *= 2;
  incAB *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    ((double *) AB)[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    zgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot_d_z_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed,
			  &((double *) y)[i * incy],
			  &((double *) r_true_l)[i * incy],
			  &((double *) r_true_t)[i * incy]);

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incAB; j += 2) {
	temp[j + 1] = -temp[j + 1];
      }
    }

    /* copy temp to AB */
    zgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_zgbmv_z_d_testgen */

void BLAS_zgbmv_d_z_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, double *AB,
			    int lda, void *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) double*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) void*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;

  double *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;
  incy *= 2;


  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * sizeof(double));
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    AB[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    dgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot_z_d_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed,
			  &((double *) y)[i * incy],
			  &((double *) r_true_l)[i * incy],
			  &((double *) r_true_t)[i * incy]);



    /* copy temp to AB */
    dgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_zgbmv_d_z_testgen */


void BLAS_zgbmv_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, void *AB,
			    int lda, void *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) void*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) void*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  int j;
  float *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;
  incy *= 2;
  incAB *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    ((float *) AB)[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    cgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot_c_c_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed,
			  &((double *) y)[i * incy],
			  &((double *) r_true_l)[i * incy],
			  &((double *) r_true_t)[i * incy]);

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incAB; j += 2) {
	temp[j + 1] = -temp[j + 1];
      }
    }

    /* copy temp to AB */
    cgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_zgbmv_c_c_testgen */

void BLAS_zgbmv_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, void *AB,
			    int lda, void *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) void*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) void*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  int j;
  double *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;
  incy *= 2;
  incAB *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * sizeof(double) * 2);
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    ((double *) AB)[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    zgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot_c_z_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed,
			  &((double *) y)[i * incy],
			  &((double *) r_true_l)[i * incy],
			  &((double *) r_true_t)[i * incy]);

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incAB; j += 2) {
	temp[j + 1] = -temp[j + 1];
      }
    }

    /* copy temp to AB */
    zgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_zgbmv_z_c_testgen */

void BLAS_zgbmv_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type trans, int m, int n, int kl,
			    int ku, void *alpha, int alpha_flag, void *AB,
			    int lda, void *x, void *beta, int beta_flag,
			    void *y, int *seed, double *r_true_l,
			    double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, AB, x, beta, and y, where AB is a banded
 * matrix; and computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of AB; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether AB is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 * 
 * kl           (input) int
 *              The number of subdiagonals
 *
 * ku           (input) int
 *              The number of superdiagonals
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * AB           (output) void*
 *              Matrix A in the banded storage.
 *
 *
 * lda          (input) int
 *              The first dimension of AB
 *
 * x            (input/output) void*
 *
 * beta         (input/output) void*
 *              If beta_flag = 1, beta is input.
 *              If beta_flag = 0, beta is output.
 *
 * beta_flag    (input) int
 *              = 0 : beta is free, and is output.
 *              = 1 : beta is fixed on input.
 *
 * y            (input/output) void*
 *
 * seed         (input/output) int
 *
 * r_true_l     (output) double*
 *              The leading part of the truth in double-double.
 *
 * r_true_t     (output) double*
 *              The trailing part of the truth in double-double.
 *
 */
{
  int n_fix2;
  int n_mix;
  int ysize;
  int i;
  int j;
  float *temp;
  int m_i, n_i;
  int max_mn;
  int incy, incAB;

  max_mn = MAX(m, n);
  incy = incAB = 1;
  incy *= 2;
  incAB *= 2;

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * sizeof(float) * 2);
  if (max_mn > 0 && temp == NULL) {
    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
  }
  for (i = 0; i < max_mn * incAB; i += incAB) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }

  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incAB; i++) {
    ((float *) AB)[i] = 0.0;
  }

  /* calling dot_testgen n time. in each iteration, one row of AB, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  for (i = 0; i < m_i; i++) {
    /* copy AB to temp */
    cgbmv_prepare(order, trans, m, n, kl, ku, AB, lda, temp, i,
		  &n_fix2, &n_mix, &ysize);

    if (i == 1) {
      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot_z_c_testgen(ysize, n_fix2, n_mix, norm, blas_no_conj, alpha,
			  alpha_flag, beta, beta_flag, x, temp, seed,
			  &((double *) y)[i * incy],
			  &((double *) r_true_l)[i * incy],
			  &((double *) r_true_t)[i * incy]);

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incAB; j += 2) {
	temp[j + 1] = -temp[j + 1];
      }
    }

    /* copy temp to AB */
    cgbmv_commit(order, trans, m, n, kl, ku, AB, lda, temp, i);
  }

  blas_free(temp);
}				/* end of BLAS_zgbmv_c_z_testgen */
