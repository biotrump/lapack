#include "blas_extended.h"
#include "cblas_test.h"




void BLAS_sgemv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n,
			 float *alpha, int alpha_flag, float *A, int lda,
			 float *head_x, float *tail_x, float *beta,
			 int beta_flag, float *y, int *seed,
			 double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) float*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) float*
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
  int i;
  float *temp;
  int m_i, n_i;
  int max_mn;

  int incy, incA;

  incy = incA = 1;



  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * incA * sizeof(float));
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    A[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_sdot2_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
		       alpha, alpha_flag, beta, beta_flag,
		       head_x, tail_x, temp, seed,
		       &y[i * incy],
		       &((double *) r_true_l)[i * incy],
		       &((double *) r_true_t)[i * incy]);



    /* copy temp to A */
    sgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_sgemv2_testgen */

void BLAS_dgemv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n,
			 double *alpha, int alpha_flag, double *A, int lda,
			 double *head_x, double *tail_x, double *beta,
			 int beta_flag, double *y, int *seed,
			 double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) double*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) double*
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
  int i;
  double *temp;
  int m_i, n_i;
  int max_mn;

  int incy, incA;

  incy = incA = 1;



  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * incA * sizeof(double));
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    A[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot2_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
		       alpha, alpha_flag, beta, beta_flag,
		       head_x, tail_x, temp, seed,
		       &y[i * incy],
		       &((double *) r_true_l)[i * incy],
		       &((double *) r_true_t)[i * incy]);



    /* copy temp to A */
    dgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_dgemv2_testgen */

void BLAS_cgemv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n,
			 void *alpha, int alpha_flag, void *A, int lda,
			 void *head_x, void *tail_x, void *beta,
			 int beta_flag, void *y, int *seed,
			 double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) void*
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
  int i;
  float *temp;
  int m_i, n_i;
  int max_mn;
  int j;
  int incy, incA;

  incy = incA = 1;
  incy *= 2;
  incA *= 2;

  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * incA * sizeof(float) * 2);
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    ((float *) A)[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot2_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
		       alpha, alpha_flag, beta, beta_flag,
		       head_x, tail_x, temp, seed,
		       &((float *) y)[i * incy],
		       &((double *) r_true_l)[i * incy],
		       &((double *) r_true_t)[i * incy]);

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incA; j += 2) {
	temp[j + 1] = -temp[j + 1];
      }
    }

    /* copy temp to A */
    cgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_cgemv2_testgen */

void BLAS_zgemv2_testgen(int norm, enum blas_order_type order,
			 enum blas_trans_type trans, int m, int n,
			 void *alpha, int alpha_flag, void *A, int lda,
			 void *head_x, void *tail_x, void *beta,
			 int beta_flag, void *y, int *seed,
			 double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) void*
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
  int i;
  double *temp;
  int m_i, n_i;
  int max_mn;
  int j;
  int incy, incA;

  incy = incA = 1;
  incy *= 2;
  incA *= 2;

  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * incA * sizeof(double) * 2);
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    ((double *) A)[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot2_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
		       alpha, alpha_flag, beta, beta_flag,
		       head_x, tail_x, temp, seed,
		       &((double *) y)[i * incy],
		       &((double *) r_true_l)[i * incy],
		       &((double *) r_true_t)[i * incy]);

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incA; j += 2) {
	temp[j + 1] = -temp[j + 1];
      }
    }

    /* copy temp to A */
    zgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_zgemv2_testgen */

void BLAS_dgemv2_s_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     double *alpha, int alpha_flag, float *A, int lda,
			     float *head_x, float *tail_x, double *beta,
			     int beta_flag, double *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) float*
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
  int i;
  float *temp;
  int m_i, n_i;
  int max_mn;

  int incy, incA;

  incy = incA = 1;



  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * incA * sizeof(float));
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    A[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot2_s_s_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
			   alpha, alpha_flag, beta, beta_flag,
			   head_x, tail_x, temp, seed,
			   &y[i * incy],
			   &((double *) r_true_l)[i * incy],
			   &((double *) r_true_t)[i * incy]);



    /* copy temp to A */
    sgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_dgemv2_s_s_testgen */

void BLAS_dgemv2_s_d_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     double *alpha, int alpha_flag, float *A, int lda,
			     double *head_x, double *tail_x, double *beta,
			     int beta_flag, double *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) double*
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
  int i;
  float *temp;
  int m_i, n_i;
  int max_mn;

  int incy, incA;

  incy = incA = 1;



  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * incA * sizeof(float));
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    A[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot2_d_s_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
			   alpha, alpha_flag, beta, beta_flag,
			   head_x, tail_x, temp, seed,
			   &y[i * incy],
			   &((double *) r_true_l)[i * incy],
			   &((double *) r_true_t)[i * incy]);



    /* copy temp to A */
    sgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_dgemv2_s_d_testgen */

void BLAS_dgemv2_d_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     double *alpha, int alpha_flag, double *A,
			     int lda, float *head_x, float *tail_x,
			     double *beta, int beta_flag, double *y,
			     int *seed, double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) double*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) double*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) float*
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
  int i;
  double *temp;
  int m_i, n_i;
  int max_mn;

  int incy, incA;

  incy = incA = 1;



  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * incA * sizeof(double));
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    A[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_ddot2_s_d_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
			   alpha, alpha_flag, beta, beta_flag,
			   head_x, tail_x, temp, seed,
			   &y[i * incy],
			   &((double *) r_true_l)[i * incy],
			   &((double *) r_true_t)[i * incy]);



    /* copy temp to A */
    dgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_dgemv2_d_s_testgen */

void BLAS_zgemv2_c_c_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, void *A, int lda,
			     void *head_x, void *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) void*
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
  int i;
  float *temp;
  int m_i, n_i;
  int max_mn;
  int j;
  int incy, incA;

  incy = incA = 1;
  incy *= 2;
  incA *= 2;

  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * incA * sizeof(float) * 2);
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    ((float *) A)[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot2_c_c_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
			   alpha, alpha_flag, beta, beta_flag,
			   head_x, tail_x, temp, seed,
			   &((double *) y)[i * incy],
			   &((double *) r_true_l)[i * incy],
			   &((double *) r_true_t)[i * incy]);

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incA; j += 2) {
	temp[j + 1] = -temp[j + 1];
      }
    }

    /* copy temp to A */
    cgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_zgemv2_c_c_testgen */

void BLAS_zgemv2_c_z_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, void *A, int lda,
			     void *head_x, void *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) void*
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
  int i;
  float *temp;
  int m_i, n_i;
  int max_mn;
  int j;
  int incy, incA;

  incy = incA = 1;
  incy *= 2;
  incA *= 2;

  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * incA * sizeof(float) * 2);
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    ((float *) A)[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot2_z_c_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
			   alpha, alpha_flag, beta, beta_flag,
			   head_x, tail_x, temp, seed,
			   &((double *) y)[i * incy],
			   &((double *) r_true_l)[i * incy],
			   &((double *) r_true_t)[i * incy]);

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incA; j += 2) {
	temp[j + 1] = -temp[j + 1];
      }
    }

    /* copy temp to A */
    cgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_zgemv2_c_z_testgen */

void BLAS_zgemv2_z_c_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, void *A, int lda,
			     void *head_x, void *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) void*
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
  int i;
  double *temp;
  int m_i, n_i;
  int max_mn;
  int j;
  int incy, incA;

  incy = incA = 1;
  incy *= 2;
  incA *= 2;

  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * incA * sizeof(double) * 2);
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    ((double *) A)[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot2_c_z_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
			   alpha, alpha_flag, beta, beta_flag,
			   head_x, tail_x, temp, seed,
			   &((double *) y)[i * incy],
			   &((double *) r_true_l)[i * incy],
			   &((double *) r_true_t)[i * incy]);

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incA; j += 2) {
	temp[j + 1] = -temp[j + 1];
      }
    }

    /* copy temp to A */
    zgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_zgemv2_z_c_testgen */

void BLAS_cgemv2_s_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, float *A, int lda,
			     float *head_x, float *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) float*
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
  int i;
  float *temp;
  int m_i, n_i;
  int max_mn;

  int incy, incA;

  incy = incA = 1;
  incy *= 2;


  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * incA * sizeof(float));
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    A[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot2_s_s_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
			   alpha, alpha_flag, beta, beta_flag,
			   head_x, tail_x, temp, seed,
			   &((float *) y)[i * incy],
			   &((double *) r_true_l)[i * incy],
			   &((double *) r_true_t)[i * incy]);



    /* copy temp to A */
    sgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_cgemv2_s_s_testgen */

void BLAS_cgemv2_s_c_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, float *A, int lda,
			     void *head_x, void *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) float*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) void*
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
  int i;
  float *temp;
  int m_i, n_i;
  int max_mn;

  int incy, incA;

  incy = incA = 1;
  incy *= 2;


  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * incA * sizeof(float));
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    A[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot2_c_s_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
			   alpha, alpha_flag, beta, beta_flag,
			   head_x, tail_x, temp, seed,
			   &((float *) y)[i * incy],
			   &((double *) r_true_l)[i * incy],
			   &((double *) r_true_t)[i * incy]);



    /* copy temp to A */
    sgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_cgemv2_s_c_testgen */

void BLAS_cgemv2_c_s_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, void *A, int lda,
			     float *head_x, float *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) float*
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
  int i;
  float *temp;
  int m_i, n_i;
  int max_mn;
  int j;
  int incy, incA;

  incy = incA = 1;
  incy *= 2;
  incA *= 2;

  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (float *) blas_malloc(max_mn * incA * sizeof(float) * 2);
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    ((float *) A)[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_cdot2_s_c_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
			   alpha, alpha_flag, beta, beta_flag,
			   head_x, tail_x, temp, seed,
			   &((float *) y)[i * incy],
			   &((double *) r_true_l)[i * incy],
			   &((double *) r_true_t)[i * incy]);

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incA; j += 2) {
	temp[j + 1] = -temp[j + 1];
      }
    }

    /* copy temp to A */
    cgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_cgemv2_c_s_testgen */

void BLAS_zgemv2_d_d_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, double *A, int lda,
			     double *head_x, double *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) double*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) double*
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
  int i;
  double *temp;
  int m_i, n_i;
  int max_mn;

  int incy, incA;

  incy = incA = 1;
  incy *= 2;


  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * incA * sizeof(double));
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    A[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot2_d_d_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
			   alpha, alpha_flag, beta, beta_flag,
			   head_x, tail_x, temp, seed,
			   &((double *) y)[i * incy],
			   &((double *) r_true_l)[i * incy],
			   &((double *) r_true_t)[i * incy]);



    /* copy temp to A */
    dgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_zgemv2_d_d_testgen */

void BLAS_zgemv2_d_z_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, double *A, int lda,
			     void *head_x, void *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) double*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) void*
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
  int i;
  double *temp;
  int m_i, n_i;
  int max_mn;

  int incy, incA;

  incy = incA = 1;
  incy *= 2;


  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * incA * sizeof(double));
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    A[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot2_z_d_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
			   alpha, alpha_flag, beta, beta_flag,
			   head_x, tail_x, temp, seed,
			   &((double *) y)[i * incy],
			   &((double *) r_true_l)[i * incy],
			   &((double *) r_true_t)[i * incy]);



    /* copy temp to A */
    dgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_zgemv2_d_z_testgen */

void BLAS_zgemv2_z_d_testgen(int norm, enum blas_order_type order,
			     enum blas_trans_type trans, int m, int n,
			     void *alpha, int alpha_flag, void *A, int lda,
			     double *head_x, double *tail_x, void *beta,
			     int beta_flag, void *y, int *seed,
			     double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 *
 * Generates alpha, A, x, beta, and y, where A is a general
 * matrix, and x has two parts: (head_x, tail_x); Computes r_true.
 *
 * Arguments
 * =========
 *
 * norm         (input) blas_norm_type 
 *
 * order        (input) blas_order_type
 *              Order of A; row or column major
 *
 * trans         (input) blas_trans_type
 *              Whether A is no trans, trans, or conj trans
 *
 * m            (input) int
 *              The number of rows 
 *
 * n            (input) int
 *              The number of columns
 *
 * alpha        (input/output) void*
 *              If alpha_flag = 1, alpha is input.
 *              If alpha_flag = 0, alpha is output.
 *
 * alpha_flag   (input) int
 *              = 0 : alpha is free, and is output.
 *              = 1 : alpha is fixed on input.
 *              
 * A            (output) void*
 *              Matrix A 
 *
 * lda          (input) int
 *              The first dimension of A
 *
 * head_x
 * tail_x       (input/output) double*
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
  int i;
  double *temp;
  int m_i, n_i;
  int max_mn;
  int j;
  int incy, incA;

  incy = incA = 1;
  incy *= 2;
  incA *= 2;

  max_mn = MAX(m, n);

  if (trans == blas_no_trans) {
    m_i = m;
    n_i = n;
  } else {
    m_i = n;
    n_i = m;
  }

  temp = (double *) blas_malloc(max_mn * incA * sizeof(double) * 2);
  if (max_mn * incA > 0 && temp == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < max_mn * incA; i += incA) {
    temp[i] = 0.0;
    temp[i + 1] = 0.0;
  }
  for (i = 0; i < (m - 1 + n - 1 + 1) * max_mn * 2 * incA; i++) {
    ((double *) A)[i] = 0.0;
  }

  /* calling dot2_testgen n times. in each iteration, one row of A, and one 
     element of y are produced. the vector x is produced at the first 
     iteration only */
  n_fix2 = n_mix = 0;
  for (i = 0; i < m_i; i++) {
    if (i == 0) {
      n_fix2 = 0;
      n_mix = 0;
    } else if (i == 1) {
      /* from now on, x is fixed */
      n_mix = n_i;

      /* from now on, fix alpha and beta */
      alpha_flag = 1;
      beta_flag = 1;
    }

    BLAS_zdot2_d_z_testgen(n_i, n_fix2, n_mix, norm, blas_no_conj,
			   alpha, alpha_flag, beta, beta_flag,
			   head_x, tail_x, temp, seed,
			   &((double *) y)[i * incy],
			   &((double *) r_true_l)[i * incy],
			   &((double *) r_true_t)[i * incy]);

    if (trans == blas_conj_trans) {
      for (j = 0; j < n_i * incA; j += 2) {
	temp[j + 1] = -temp[j + 1];
      }
    }

    /* copy temp to A */
    zgemv2_commit(order, trans, m, n, A, lda, temp, i);
  }


  blas_free(temp);
}				/* end of BLAS_zgemv2_z_d_testgen */