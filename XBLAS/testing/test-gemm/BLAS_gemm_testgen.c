#include <stdio.h>
#include <stdlib.h>
#include "blas_extended.h"
#include "cblas_test.h"














void BLAS_sgemm_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type transa,
			enum blas_trans_type transb, int m, int n, int k,
			int randomize, float *alpha, int alpha_flag, float *a,
			int lda, float *beta, int beta_flag, float *b,
			int ldb, float *c, int ldc, int *seed,
			double *r_true_l, double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  float c_elem;
  float a_elem;
  float b_elem;

  double r_elem;
  double r_true_t_elem;
  double r_true_l_elem;

  float *c_i = c;
  float *b_i = b;
  float *a_i = a;
  float *alpha_i = alpha;
  float *beta_i = beta;

  /* Temporary storage space for vectors of length k */
  float *a_vec;
  float *b_vec;

  inca = incb = 1;



  a_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }




  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_sdot_testgen(k, 0, 0, norm, blas_no_conj,
		      alpha, alpha_flag, beta, beta_flag,
		      b_vec, a_vec, seed, &c_elem,
		      &r_true_l_elem, &r_true_t_elem);
    c_i[cij] = c_elem;
    r_true_l[cij] = r_true_l_elem;
    r_true_t[cij] = r_true_t_elem;

    /* Copy a_vec to the first row of A */
    sgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      sgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_sdot_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			beta, 1, b_vec, a_vec, seed,
			&c_elem, &r_true_l_elem, &r_true_t_elem);

      c_i[cij] = c_elem;
      r_true_l[cij] = r_true_l_elem;
      r_true_t[cij] = r_true_t_elem;
      sgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem = c_i[ci];
	c_i[cij] = c_elem;
	r_elem = r_true_l[ci];
	r_true_l[cij] = r_elem;
	r_elem = r_true_t[ci];
	r_true_t[cij] = r_elem;
      }
    }
  } else {








    if (alpha_flag == 0) {
      c_elem = (float) drand48();
      alpha_i[0] = c_elem;
    }
    if (beta_flag == 0) {
      c_elem = (float) drand48();
      beta_i[0] = c_elem;
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }






    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	b_elem = (float) drand48();
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	a_elem = (float) drand48();
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      sgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);



      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	sgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);



	BLAS_sdot_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1, b_vec, a_vec, seed,
			  &c_elem, &r_true_l_elem, &r_true_t_elem);

	c_i[cij] = c_elem;
	r_true_l[cij] = r_true_l_elem;
	r_true_t[cij] = r_true_t_elem;

      }
    }




  }

  free(a_vec);
  free(b_vec);

}
void BLAS_dgemm_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type transa,
			enum blas_trans_type transb, int m, int n, int k,
			int randomize, double *alpha, int alpha_flag,
			double *a, int lda, double *beta, int beta_flag,
			double *b, int ldb, double *c, int ldc, int *seed,
			double *r_true_l, double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  double c_elem;
  double a_elem;
  double b_elem;

  double r_elem;
  double r_true_t_elem;
  double r_true_l_elem;

  double *c_i = c;
  double *b_i = b;
  double *a_i = a;
  double *alpha_i = alpha;
  double *beta_i = beta;

  /* Temporary storage space for vectors of length k */
  double *a_vec;
  double *b_vec;

  inca = incb = 1;



  a_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }




  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_ddot_testgen(k, 0, 0, norm, blas_no_conj,
		      alpha, alpha_flag, beta, beta_flag,
		      b_vec, a_vec, seed, &c_elem,
		      &r_true_l_elem, &r_true_t_elem);
    c_i[cij] = c_elem;
    r_true_l[cij] = r_true_l_elem;
    r_true_t[cij] = r_true_t_elem;

    /* Copy a_vec to the first row of A */
    dgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      dgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_ddot_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			beta, 1, b_vec, a_vec, seed,
			&c_elem, &r_true_l_elem, &r_true_t_elem);

      c_i[cij] = c_elem;
      r_true_l[cij] = r_true_l_elem;
      r_true_t[cij] = r_true_t_elem;
      dgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem = c_i[ci];
	c_i[cij] = c_elem;
	r_elem = r_true_l[ci];
	r_true_l[cij] = r_elem;
	r_elem = r_true_t[ci];
	r_true_t[cij] = r_elem;
      }
    }
  } else {








    if (alpha_flag == 0) {
      c_elem = (float) drand48();
      alpha_i[0] = c_elem;
    }
    if (beta_flag == 0) {
      c_elem = (float) drand48();
      beta_i[0] = c_elem;
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }






    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	b_elem = (float) drand48();
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	a_elem = (float) drand48();
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      dgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);



      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	dgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);



	BLAS_ddot_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1, b_vec, a_vec, seed,
			  &c_elem, &r_true_l_elem, &r_true_t_elem);

	c_i[cij] = c_elem;
	r_true_l[cij] = r_true_l_elem;
	r_true_t[cij] = r_true_t_elem;

      }
    }




  }

  free(a_vec);
  free(b_vec);

}
void BLAS_dgemm_d_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, double *alpha, int alpha_flag,
			    double *a, int lda, double *beta, int beta_flag,
			    float *b, int ldb, double *c, int ldc, int *seed,
			    double *r_true_l, double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  double c_elem;
  double a_elem;
  float b_elem;

  double r_elem;
  double r_true_t_elem;
  double r_true_l_elem;

  double *c_i = c;
  float *b_i = b;
  double *a_i = a;
  double *alpha_i = alpha;
  double *beta_i = beta;

  /* Temporary storage space for vectors of length k */
  double *a_vec;
  float *b_vec;

  inca = incb = 1;



  a_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }




  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_ddot_s_d_testgen(k, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, &c_elem,
			  &r_true_l_elem, &r_true_t_elem);
    c_i[cij] = c_elem;
    r_true_l[cij] = r_true_l_elem;
    r_true_t[cij] = r_true_t_elem;

    /* Copy a_vec to the first row of A */
    dgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      sgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_ddot_s_d_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    &c_elem, &r_true_l_elem, &r_true_t_elem);

      c_i[cij] = c_elem;
      r_true_l[cij] = r_true_l_elem;
      r_true_t[cij] = r_true_t_elem;
      dgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem = c_i[ci];
	c_i[cij] = c_elem;
	r_elem = r_true_l[ci];
	r_true_l[cij] = r_elem;
	r_elem = r_true_t[ci];
	r_true_t[cij] = r_elem;
      }
    }
  } else {








    if (alpha_flag == 0) {
      c_elem = (float) drand48();
      alpha_i[0] = c_elem;
    }
    if (beta_flag == 0) {
      c_elem = (float) drand48();
      beta_i[0] = c_elem;
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }






    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	b_elem = (float) drand48();
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	a_elem = (float) drand48();
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      dgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);



      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	sgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);



	BLAS_ddot_s_d_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			      beta, 1, b_vec, a_vec, seed,
			      &c_elem, &r_true_l_elem, &r_true_t_elem);

	c_i[cij] = c_elem;
	r_true_l[cij] = r_true_l_elem;
	r_true_t[cij] = r_true_t_elem;

      }
    }




  }

  free(a_vec);
  free(b_vec);

}
void BLAS_dgemm_s_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, double *alpha, int alpha_flag,
			    float *a, int lda, double *beta, int beta_flag,
			    double *b, int ldb, double *c, int ldc, int *seed,
			    double *r_true_l, double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  double c_elem;
  float a_elem;
  double b_elem;

  double r_elem;
  double r_true_t_elem;
  double r_true_l_elem;

  double *c_i = c;
  double *b_i = b;
  float *a_i = a;
  double *alpha_i = alpha;
  double *beta_i = beta;

  /* Temporary storage space for vectors of length k */
  float *a_vec;
  double *b_vec;

  inca = incb = 1;



  a_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }




  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_ddot_d_s_testgen(k, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, &c_elem,
			  &r_true_l_elem, &r_true_t_elem);
    c_i[cij] = c_elem;
    r_true_l[cij] = r_true_l_elem;
    r_true_t[cij] = r_true_t_elem;

    /* Copy a_vec to the first row of A */
    sgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      dgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_ddot_d_s_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    &c_elem, &r_true_l_elem, &r_true_t_elem);

      c_i[cij] = c_elem;
      r_true_l[cij] = r_true_l_elem;
      r_true_t[cij] = r_true_t_elem;
      sgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem = c_i[ci];
	c_i[cij] = c_elem;
	r_elem = r_true_l[ci];
	r_true_l[cij] = r_elem;
	r_elem = r_true_t[ci];
	r_true_t[cij] = r_elem;
      }
    }
  } else {








    if (alpha_flag == 0) {
      c_elem = (float) drand48();
      alpha_i[0] = c_elem;
    }
    if (beta_flag == 0) {
      c_elem = (float) drand48();
      beta_i[0] = c_elem;
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }






    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	b_elem = (float) drand48();
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	a_elem = (float) drand48();
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      sgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);



      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	dgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);



	BLAS_ddot_d_s_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			      beta, 1, b_vec, a_vec, seed,
			      &c_elem, &r_true_l_elem, &r_true_t_elem);

	c_i[cij] = c_elem;
	r_true_l[cij] = r_true_l_elem;
	r_true_t[cij] = r_true_t_elem;

      }
    }




  }

  free(a_vec);
  free(b_vec);

}
void BLAS_dgemm_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, double *alpha, int alpha_flag,
			    float *a, int lda, double *beta, int beta_flag,
			    float *b, int ldb, double *c, int ldc, int *seed,
			    double *r_true_l, double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  double c_elem;
  float a_elem;
  float b_elem;

  double r_elem;
  double r_true_t_elem;
  double r_true_l_elem;

  double *c_i = c;
  float *b_i = b;
  float *a_i = a;
  double *alpha_i = alpha;
  double *beta_i = beta;

  /* Temporary storage space for vectors of length k */
  float *a_vec;
  float *b_vec;

  inca = incb = 1;



  a_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }




  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_ddot_s_s_testgen(k, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, &c_elem,
			  &r_true_l_elem, &r_true_t_elem);
    c_i[cij] = c_elem;
    r_true_l[cij] = r_true_l_elem;
    r_true_t[cij] = r_true_t_elem;

    /* Copy a_vec to the first row of A */
    sgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      sgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_ddot_s_s_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    &c_elem, &r_true_l_elem, &r_true_t_elem);

      c_i[cij] = c_elem;
      r_true_l[cij] = r_true_l_elem;
      r_true_t[cij] = r_true_t_elem;
      sgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem = c_i[ci];
	c_i[cij] = c_elem;
	r_elem = r_true_l[ci];
	r_true_l[cij] = r_elem;
	r_elem = r_true_t[ci];
	r_true_t[cij] = r_elem;
      }
    }
  } else {








    if (alpha_flag == 0) {
      c_elem = (float) drand48();
      alpha_i[0] = c_elem;
    }
    if (beta_flag == 0) {
      c_elem = (float) drand48();
      beta_i[0] = c_elem;
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }






    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	b_elem = (float) drand48();
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	a_elem = (float) drand48();
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      sgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);



      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	sgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);



	BLAS_ddot_s_s_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			      beta, 1, b_vec, a_vec, seed,
			      &c_elem, &r_true_l_elem, &r_true_t_elem);

	c_i[cij] = c_elem;
	r_true_l[cij] = r_true_l_elem;
	r_true_t[cij] = r_true_t_elem;

      }
    }




  }

  free(a_vec);
  free(b_vec);

}

void BLAS_cgemm_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type transa,
			enum blas_trans_type transb, int m, int n, int k,
			int randomize, void *alpha, int alpha_flag, void *a,
			int lda, void *beta, int beta_flag, void *b, int ldb,
			void *c, int ldc, int *seed, double *r_true_l,
			double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  float c_elem[2];
  float a_elem[2];
  float b_elem[2];

  double r_elem[2];
  double r_true_t_elem[2];
  double r_true_l_elem[2];

  float *c_i = (float *) c;
  float *b_i = (float *) b;
  float *a_i = (float *) a;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;

  /* Temporary storage space for vectors of length k */
  float *a_vec;
  float *b_vec;

  inca = incb = 1;
  inca *= 2;
  incb *= 2;

  a_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }
  b_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
    b_vec[i + 1] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;

  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_cdot_testgen(k, 0, 0, norm, blas_no_conj,
		      alpha, alpha_flag, beta, beta_flag,
		      b_vec, a_vec, seed, c_elem,
		      r_true_l_elem, r_true_t_elem);
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    r_true_l[cij] = r_true_l_elem[0];
    r_true_l[cij + 1] = r_true_l_elem[1];
    r_true_t[cij] = r_true_t_elem[0];
    r_true_t[cij + 1] = r_true_t_elem[1];

    /* Copy a_vec to the first row of A */
    cgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      cgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_cdot_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			beta, 1, b_vec, a_vec, seed,
			c_elem, r_true_l_elem, r_true_t_elem);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      r_true_l[cij] = r_true_l_elem[0];
      r_true_l[cij + 1] = r_true_l_elem[1];
      r_true_t[cij] = r_true_t_elem[0];
      r_true_t[cij + 1] = r_true_t_elem[1];
      cgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_elem[0] = r_true_l[ci];
	r_elem[1] = r_true_l[ci + 1];
	r_true_l[cij] = r_elem[0];
	r_true_l[cij + 1] = r_elem[1];
	r_elem[0] = r_true_t[ci];
	r_elem[1] = r_true_t[ci + 1];
	r_true_t[cij] = r_elem[0];
	r_true_t[cij + 1] = r_elem[1];
      }
    }
  } else {








    if (alpha_flag == 0) {
      ((float *) c_elem)[0] = (float) drand48();
      ((float *) c_elem)[1] = (float) drand48();
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      ((float *) c_elem)[0] = (float) drand48();
      ((float *) c_elem)[1] = (float) drand48();
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incbi *= 2;
    incbij *= 2;
    incai *= 2;
    incaij *= 2;

    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	((float *) b_elem)[0] = (float) drand48();
	((float *) b_elem)[1] = (float) drand48();
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	((float *) a_elem)[0] = (float) drand48();
	((float *) a_elem)[1] = (float) drand48();
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      cgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);



      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	cgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);



	BLAS_cdot_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1, b_vec, a_vec, seed,
			  c_elem, r_true_l_elem, r_true_t_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_true_l[cij] = r_true_l_elem[0];
	r_true_l[cij + 1] = r_true_l_elem[1];
	r_true_t[cij] = r_true_t_elem[0];
	r_true_t[cij + 1] = r_true_t_elem[1];

      }
    }




  }

  free(a_vec);
  free(b_vec);

}
void BLAS_zgemm_testgen(int norm, enum blas_order_type order,
			enum blas_trans_type transa,
			enum blas_trans_type transb, int m, int n, int k,
			int randomize, void *alpha, int alpha_flag, void *a,
			int lda, void *beta, int beta_flag, void *b, int ldb,
			void *c, int ldc, int *seed, double *r_true_l,
			double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  double c_elem[2];
  double a_elem[2];
  double b_elem[2];

  double r_elem[2];
  double r_true_t_elem[2];
  double r_true_l_elem[2];

  double *c_i = (double *) c;
  double *b_i = (double *) b;
  double *a_i = (double *) a;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;

  /* Temporary storage space for vectors of length k */
  double *a_vec;
  double *b_vec;

  inca = incb = 1;
  inca *= 2;
  incb *= 2;

  a_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }
  b_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
    b_vec[i + 1] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;

  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_zdot_testgen(k, 0, 0, norm, blas_no_conj,
		      alpha, alpha_flag, beta, beta_flag,
		      b_vec, a_vec, seed, c_elem,
		      r_true_l_elem, r_true_t_elem);
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    r_true_l[cij] = r_true_l_elem[0];
    r_true_l[cij + 1] = r_true_l_elem[1];
    r_true_t[cij] = r_true_t_elem[0];
    r_true_t[cij + 1] = r_true_t_elem[1];

    /* Copy a_vec to the first row of A */
    zgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      zgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_zdot_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			beta, 1, b_vec, a_vec, seed,
			c_elem, r_true_l_elem, r_true_t_elem);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      r_true_l[cij] = r_true_l_elem[0];
      r_true_l[cij + 1] = r_true_l_elem[1];
      r_true_t[cij] = r_true_t_elem[0];
      r_true_t[cij + 1] = r_true_t_elem[1];
      zgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_elem[0] = r_true_l[ci];
	r_elem[1] = r_true_l[ci + 1];
	r_true_l[cij] = r_elem[0];
	r_true_l[cij + 1] = r_elem[1];
	r_elem[0] = r_true_t[ci];
	r_elem[1] = r_true_t[ci + 1];
	r_true_t[cij] = r_elem[0];
	r_true_t[cij + 1] = r_elem[1];
      }
    }
  } else {








    if (alpha_flag == 0) {
      ((double *) c_elem)[0] = (float) drand48();
      ((double *) c_elem)[1] = (float) drand48();
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      ((double *) c_elem)[0] = (float) drand48();
      ((double *) c_elem)[1] = (float) drand48();
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incbi *= 2;
    incbij *= 2;
    incai *= 2;
    incaij *= 2;

    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	((double *) b_elem)[0] = (float) drand48();
	((double *) b_elem)[1] = (float) drand48();
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	((double *) a_elem)[0] = (float) drand48();
	((double *) a_elem)[1] = (float) drand48();
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      zgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);



      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	zgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);



	BLAS_zdot_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1, b_vec, a_vec, seed,
			  c_elem, r_true_l_elem, r_true_t_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_true_l[cij] = r_true_l_elem[0];
	r_true_l[cij + 1] = r_true_l_elem[1];
	r_true_t[cij] = r_true_t_elem[0];
	r_true_t[cij + 1] = r_true_t_elem[1];

      }
    }




  }

  free(a_vec);
  free(b_vec);

}
void BLAS_zgemm_z_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    void *a, int lda, void *beta, int beta_flag,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *r_true_l, double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  double c_elem[2];
  double a_elem[2];
  float b_elem[2];

  double r_elem[2];
  double r_true_t_elem[2];
  double r_true_l_elem[2];

  double *c_i = (double *) c;
  float *b_i = (float *) b;
  double *a_i = (double *) a;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;

  /* Temporary storage space for vectors of length k */
  double *a_vec;
  float *b_vec;

  inca = incb = 1;
  inca *= 2;
  incb *= 2;

  a_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }
  b_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
    b_vec[i + 1] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;

  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_zdot_c_z_testgen(k, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  r_true_l_elem, r_true_t_elem);
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    r_true_l[cij] = r_true_l_elem[0];
    r_true_l[cij + 1] = r_true_l_elem[1];
    r_true_t[cij] = r_true_t_elem[0];
    r_true_t[cij + 1] = r_true_t_elem[1];

    /* Copy a_vec to the first row of A */
    zgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      cgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_zdot_c_z_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, r_true_l_elem, r_true_t_elem);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      r_true_l[cij] = r_true_l_elem[0];
      r_true_l[cij + 1] = r_true_l_elem[1];
      r_true_t[cij] = r_true_t_elem[0];
      r_true_t[cij + 1] = r_true_t_elem[1];
      zgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_elem[0] = r_true_l[ci];
	r_elem[1] = r_true_l[ci + 1];
	r_true_l[cij] = r_elem[0];
	r_true_l[cij + 1] = r_elem[1];
	r_elem[0] = r_true_t[ci];
	r_elem[1] = r_true_t[ci + 1];
	r_true_t[cij] = r_elem[0];
	r_true_t[cij + 1] = r_elem[1];
      }
    }
  } else {








    if (alpha_flag == 0) {
      ((double *) c_elem)[0] = (float) drand48();
      ((double *) c_elem)[1] = (float) drand48();
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      ((double *) c_elem)[0] = (float) drand48();
      ((double *) c_elem)[1] = (float) drand48();
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incbi *= 2;
    incbij *= 2;
    incai *= 2;
    incaij *= 2;

    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	((float *) b_elem)[0] = (float) drand48();
	((float *) b_elem)[1] = (float) drand48();
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	((double *) a_elem)[0] = (float) drand48();
	((double *) a_elem)[1] = (float) drand48();
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      zgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);



      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	cgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);



	BLAS_zdot_c_z_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			      beta, 1, b_vec, a_vec, seed,
			      c_elem, r_true_l_elem, r_true_t_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_true_l[cij] = r_true_l_elem[0];
	r_true_l[cij + 1] = r_true_l_elem[1];
	r_true_t[cij] = r_true_t_elem[0];
	r_true_t[cij + 1] = r_true_t_elem[1];

      }
    }




  }

  free(a_vec);
  free(b_vec);

}
void BLAS_zgemm_c_z_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    void *a, int lda, void *beta, int beta_flag,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *r_true_l, double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  double c_elem[2];
  float a_elem[2];
  double b_elem[2];

  double r_elem[2];
  double r_true_t_elem[2];
  double r_true_l_elem[2];

  double *c_i = (double *) c;
  double *b_i = (double *) b;
  float *a_i = (float *) a;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;

  /* Temporary storage space for vectors of length k */
  float *a_vec;
  double *b_vec;

  inca = incb = 1;
  inca *= 2;
  incb *= 2;

  a_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }
  b_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
    b_vec[i + 1] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;

  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_zdot_z_c_testgen(k, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  r_true_l_elem, r_true_t_elem);
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    r_true_l[cij] = r_true_l_elem[0];
    r_true_l[cij + 1] = r_true_l_elem[1];
    r_true_t[cij] = r_true_t_elem[0];
    r_true_t[cij + 1] = r_true_t_elem[1];

    /* Copy a_vec to the first row of A */
    cgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      zgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_zdot_z_c_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, r_true_l_elem, r_true_t_elem);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      r_true_l[cij] = r_true_l_elem[0];
      r_true_l[cij + 1] = r_true_l_elem[1];
      r_true_t[cij] = r_true_t_elem[0];
      r_true_t[cij + 1] = r_true_t_elem[1];
      cgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_elem[0] = r_true_l[ci];
	r_elem[1] = r_true_l[ci + 1];
	r_true_l[cij] = r_elem[0];
	r_true_l[cij + 1] = r_elem[1];
	r_elem[0] = r_true_t[ci];
	r_elem[1] = r_true_t[ci + 1];
	r_true_t[cij] = r_elem[0];
	r_true_t[cij + 1] = r_elem[1];
      }
    }
  } else {








    if (alpha_flag == 0) {
      ((double *) c_elem)[0] = (float) drand48();
      ((double *) c_elem)[1] = (float) drand48();
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      ((double *) c_elem)[0] = (float) drand48();
      ((double *) c_elem)[1] = (float) drand48();
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incbi *= 2;
    incbij *= 2;
    incai *= 2;
    incaij *= 2;

    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	((double *) b_elem)[0] = (float) drand48();
	((double *) b_elem)[1] = (float) drand48();
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	((float *) a_elem)[0] = (float) drand48();
	((float *) a_elem)[1] = (float) drand48();
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      cgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);



      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	zgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);



	BLAS_zdot_z_c_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			      beta, 1, b_vec, a_vec, seed,
			      c_elem, r_true_l_elem, r_true_t_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_true_l[cij] = r_true_l_elem[0];
	r_true_l[cij + 1] = r_true_l_elem[1];
	r_true_t[cij] = r_true_t_elem[0];
	r_true_t[cij + 1] = r_true_t_elem[1];

      }
    }




  }

  free(a_vec);
  free(b_vec);

}
void BLAS_zgemm_c_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    void *a, int lda, void *beta, int beta_flag,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *r_true_l, double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  double c_elem[2];
  float a_elem[2];
  float b_elem[2];

  double r_elem[2];
  double r_true_t_elem[2];
  double r_true_l_elem[2];

  double *c_i = (double *) c;
  float *b_i = (float *) b;
  float *a_i = (float *) a;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;

  /* Temporary storage space for vectors of length k */
  float *a_vec;
  float *b_vec;

  inca = incb = 1;
  inca *= 2;
  incb *= 2;

  a_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }
  b_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
    b_vec[i + 1] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;

  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_zdot_c_c_testgen(k, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  r_true_l_elem, r_true_t_elem);
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    r_true_l[cij] = r_true_l_elem[0];
    r_true_l[cij + 1] = r_true_l_elem[1];
    r_true_t[cij] = r_true_t_elem[0];
    r_true_t[cij + 1] = r_true_t_elem[1];

    /* Copy a_vec to the first row of A */
    cgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      cgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_zdot_c_c_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, r_true_l_elem, r_true_t_elem);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      r_true_l[cij] = r_true_l_elem[0];
      r_true_l[cij + 1] = r_true_l_elem[1];
      r_true_t[cij] = r_true_t_elem[0];
      r_true_t[cij + 1] = r_true_t_elem[1];
      cgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_elem[0] = r_true_l[ci];
	r_elem[1] = r_true_l[ci + 1];
	r_true_l[cij] = r_elem[0];
	r_true_l[cij + 1] = r_elem[1];
	r_elem[0] = r_true_t[ci];
	r_elem[1] = r_true_t[ci + 1];
	r_true_t[cij] = r_elem[0];
	r_true_t[cij + 1] = r_elem[1];
      }
    }
  } else {








    if (alpha_flag == 0) {
      ((double *) c_elem)[0] = (float) drand48();
      ((double *) c_elem)[1] = (float) drand48();
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      ((double *) c_elem)[0] = (float) drand48();
      ((double *) c_elem)[1] = (float) drand48();
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incbi *= 2;
    incbij *= 2;
    incai *= 2;
    incaij *= 2;

    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	((float *) b_elem)[0] = (float) drand48();
	((float *) b_elem)[1] = (float) drand48();
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	((float *) a_elem)[0] = (float) drand48();
	((float *) a_elem)[1] = (float) drand48();
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      cgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);



      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	cgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);



	BLAS_zdot_c_c_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			      beta, 1, b_vec, a_vec, seed,
			      c_elem, r_true_l_elem, r_true_t_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_true_l[cij] = r_true_l_elem[0];
	r_true_l[cij + 1] = r_true_l_elem[1];
	r_true_t[cij] = r_true_t_elem[0];
	r_true_t[cij + 1] = r_true_t_elem[1];

      }
    }




  }

  free(a_vec);
  free(b_vec);

}
void BLAS_cgemm_c_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    void *a, int lda, void *beta, int beta_flag,
			    float *b, int ldb, void *c, int ldc, int *seed,
			    double *r_true_l, double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  float c_elem[2];
  float a_elem[2];
  float b_elem;

  double r_elem[2];
  double r_true_t_elem[2];
  double r_true_l_elem[2];

  float *c_i = (float *) c;
  float *b_i = b;
  float *a_i = (float *) a;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;

  /* Temporary storage space for vectors of length k */
  float *a_vec;
  float *b_vec;

  inca = incb = 1;
  inca *= 2;


  a_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }
  b_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;

  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_cdot_s_c_testgen(k, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  r_true_l_elem, r_true_t_elem);
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    r_true_l[cij] = r_true_l_elem[0];
    r_true_l[cij + 1] = r_true_l_elem[1];
    r_true_t[cij] = r_true_t_elem[0];
    r_true_t[cij + 1] = r_true_t_elem[1];

    /* Copy a_vec to the first row of A */
    cgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      sgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_cdot_s_c_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, r_true_l_elem, r_true_t_elem);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      r_true_l[cij] = r_true_l_elem[0];
      r_true_l[cij + 1] = r_true_l_elem[1];
      r_true_t[cij] = r_true_t_elem[0];
      r_true_t[cij + 1] = r_true_t_elem[1];
      cgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_elem[0] = r_true_l[ci];
	r_elem[1] = r_true_l[ci + 1];
	r_true_l[cij] = r_elem[0];
	r_true_l[cij + 1] = r_elem[1];
	r_elem[0] = r_true_t[ci];
	r_elem[1] = r_true_t[ci + 1];
	r_true_t[cij] = r_elem[0];
	r_true_t[cij + 1] = r_elem[1];
      }
    }
  } else {


    float *bb_vec;


    bb_vec = (float *) blas_malloc(k * sizeof(float) * 2);
    if (k > 0 && bb_vec == NULL) {
      printf("malloc failed\n");
      exit(-1);
    }


    if (alpha_flag == 0) {
      ((float *) c_elem)[0] = (float) drand48();
      ((float *) c_elem)[1] = (float) drand48();
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      ((float *) c_elem)[0] = (float) drand48();
      ((float *) c_elem)[1] = (float) drand48();
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }



    incai *= 2;
    incaij *= 2;

    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	b_elem = (float) drand48();
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	((float *) a_elem)[0] = (float) drand48();
	((float *) a_elem)[1] = (float) drand48();
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      cgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);



      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	sgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);

	{
	  int r;
	  for (r = 0; r < k; r++) {
	    bb_vec[2 * r] = b_vec[r];
	    bb_vec[2 * r + 1] = 0.0;
	  }
	}

	BLAS_cdot_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1,
			  bb_vec,
			  a_vec, seed, c_elem, r_true_l_elem, r_true_t_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_true_l[cij] = r_true_l_elem[0];
	r_true_l[cij + 1] = r_true_l_elem[1];
	r_true_t[cij] = r_true_t_elem[0];
	r_true_t[cij + 1] = r_true_t_elem[1];

      }
    }


    free(bb_vec);

  }

  free(a_vec);
  free(b_vec);

}
void BLAS_cgemm_s_c_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    float *a, int lda, void *beta, int beta_flag,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *r_true_l, double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  float c_elem[2];
  float a_elem;
  float b_elem[2];

  double r_elem[2];
  double r_true_t_elem[2];
  double r_true_l_elem[2];

  float *c_i = (float *) c;
  float *b_i = (float *) b;
  float *a_i = a;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;

  /* Temporary storage space for vectors of length k */
  float *a_vec;
  float *b_vec;

  inca = incb = 1;

  incb *= 2;

  a_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (float *) blas_malloc(k * sizeof(float) * 2);
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
    b_vec[i + 1] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;

  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_cdot_c_s_testgen(k, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  r_true_l_elem, r_true_t_elem);
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    r_true_l[cij] = r_true_l_elem[0];
    r_true_l[cij + 1] = r_true_l_elem[1];
    r_true_t[cij] = r_true_t_elem[0];
    r_true_t[cij + 1] = r_true_t_elem[1];

    /* Copy a_vec to the first row of A */
    sgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      cgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_cdot_c_s_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, r_true_l_elem, r_true_t_elem);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      r_true_l[cij] = r_true_l_elem[0];
      r_true_l[cij + 1] = r_true_l_elem[1];
      r_true_t[cij] = r_true_t_elem[0];
      r_true_t[cij + 1] = r_true_t_elem[1];
      sgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_elem[0] = r_true_l[ci];
	r_elem[1] = r_true_l[ci + 1];
	r_true_l[cij] = r_elem[0];
	r_true_l[cij + 1] = r_elem[1];
	r_elem[0] = r_true_t[ci];
	r_elem[1] = r_true_t[ci + 1];
	r_true_t[cij] = r_elem[0];
	r_true_t[cij + 1] = r_elem[1];
      }
    }
  } else {

    float *aa_vec;


    aa_vec = (float *) blas_malloc(k * sizeof(float) * 2);
    if (k > 0 && aa_vec == NULL) {
      printf("malloc failed\n");
      exit(-1);
    }



    if (alpha_flag == 0) {
      ((float *) c_elem)[0] = (float) drand48();
      ((float *) c_elem)[1] = (float) drand48();
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      ((float *) c_elem)[0] = (float) drand48();
      ((float *) c_elem)[1] = (float) drand48();
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incbi *= 2;
    incbij *= 2;



    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	((float *) b_elem)[0] = (float) drand48();
	((float *) b_elem)[1] = (float) drand48();
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	a_elem = (float) drand48();
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      sgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);

      {
	int r;
	for (r = 0; r < k; r++) {
	  aa_vec[2 * r] = a_vec[r];
	  aa_vec[2 * r + 1] = 0.0;
	}
      }

      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	cgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);



	BLAS_cdot_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1,
			  b_vec,
			  aa_vec, seed, c_elem, r_true_l_elem, r_true_t_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_true_l[cij] = r_true_l_elem[0];
	r_true_l[cij + 1] = r_true_l_elem[1];
	r_true_t[cij] = r_true_t_elem[0];
	r_true_t[cij + 1] = r_true_t_elem[1];

      }
    }

    free(aa_vec);


  }

  free(a_vec);
  free(b_vec);

}
void BLAS_cgemm_s_s_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    float *a, int lda, void *beta, int beta_flag,
			    float *b, int ldb, void *c, int ldc, int *seed,
			    double *r_true_l, double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  float c_elem[2];
  float a_elem;
  float b_elem;

  double r_elem[2];
  double r_true_t_elem[2];
  double r_true_l_elem[2];

  float *c_i = (float *) c;
  float *b_i = b;
  float *a_i = a;
  float *alpha_i = (float *) alpha;
  float *beta_i = (float *) beta;

  /* Temporary storage space for vectors of length k */
  float *a_vec;
  float *b_vec;

  inca = incb = 1;



  a_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (float *) blas_malloc(k * sizeof(float));
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;

  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_cdot_s_s_testgen(k, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  r_true_l_elem, r_true_t_elem);
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    r_true_l[cij] = r_true_l_elem[0];
    r_true_l[cij + 1] = r_true_l_elem[1];
    r_true_t[cij] = r_true_t_elem[0];
    r_true_t[cij + 1] = r_true_t_elem[1];

    /* Copy a_vec to the first row of A */
    sgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      sgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_cdot_s_s_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, r_true_l_elem, r_true_t_elem);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      r_true_l[cij] = r_true_l_elem[0];
      r_true_l[cij + 1] = r_true_l_elem[1];
      r_true_t[cij] = r_true_t_elem[0];
      r_true_t[cij + 1] = r_true_t_elem[1];
      sgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_elem[0] = r_true_l[ci];
	r_elem[1] = r_true_l[ci + 1];
	r_true_l[cij] = r_elem[0];
	r_true_l[cij + 1] = r_elem[1];
	r_elem[0] = r_true_t[ci];
	r_elem[1] = r_true_t[ci + 1];
	r_true_t[cij] = r_elem[0];
	r_true_t[cij + 1] = r_elem[1];
      }
    }
  } else {

    float *aa_vec;
    float *bb_vec;

    aa_vec = (float *) blas_malloc(k * sizeof(float) * 2);
    if (k > 0 && aa_vec == NULL) {
      printf("malloc failed\n");
      exit(-1);
    }
    bb_vec = (float *) blas_malloc(k * sizeof(float) * 2);
    if (k > 0 && bb_vec == NULL) {
      printf("malloc failed\n");
      exit(-1);
    }


    if (alpha_flag == 0) {
      ((float *) c_elem)[0] = (float) drand48();
      ((float *) c_elem)[1] = (float) drand48();
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      ((float *) c_elem)[0] = (float) drand48();
      ((float *) c_elem)[1] = (float) drand48();
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }






    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	b_elem = (float) drand48();
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	a_elem = (float) drand48();
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      sgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);

      {
	int r;
	for (r = 0; r < k; r++) {
	  aa_vec[2 * r] = a_vec[r];
	  aa_vec[2 * r + 1] = 0.0;
	}
      }

      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	sgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);

	{
	  int r;
	  for (r = 0; r < k; r++) {
	    bb_vec[2 * r] = b_vec[r];
	    bb_vec[2 * r + 1] = 0.0;
	  }
	}

	BLAS_cdot_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1,
			  bb_vec,
			  aa_vec, seed, c_elem, r_true_l_elem, r_true_t_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_true_l[cij] = r_true_l_elem[0];
	r_true_l[cij + 1] = r_true_l_elem[1];
	r_true_t[cij] = r_true_t_elem[0];
	r_true_t[cij + 1] = r_true_t_elem[1];

      }
    }

    free(aa_vec);
    free(bb_vec);

  }

  free(a_vec);
  free(b_vec);

}
void BLAS_zgemm_z_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    void *a, int lda, void *beta, int beta_flag,
			    double *b, int ldb, void *c, int ldc, int *seed,
			    double *r_true_l, double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  double c_elem[2];
  double a_elem[2];
  double b_elem;

  double r_elem[2];
  double r_true_t_elem[2];
  double r_true_l_elem[2];

  double *c_i = (double *) c;
  double *b_i = b;
  double *a_i = (double *) a;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;

  /* Temporary storage space for vectors of length k */
  double *a_vec;
  double *b_vec;

  inca = incb = 1;
  inca *= 2;


  a_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
    a_vec[i + 1] = 0.0;
  }
  b_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;

  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_zdot_d_z_testgen(k, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  r_true_l_elem, r_true_t_elem);
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    r_true_l[cij] = r_true_l_elem[0];
    r_true_l[cij + 1] = r_true_l_elem[1];
    r_true_t[cij] = r_true_t_elem[0];
    r_true_t[cij + 1] = r_true_t_elem[1];

    /* Copy a_vec to the first row of A */
    zgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      dgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_zdot_d_z_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, r_true_l_elem, r_true_t_elem);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      r_true_l[cij] = r_true_l_elem[0];
      r_true_l[cij + 1] = r_true_l_elem[1];
      r_true_t[cij] = r_true_t_elem[0];
      r_true_t[cij + 1] = r_true_t_elem[1];
      zgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_elem[0] = r_true_l[ci];
	r_elem[1] = r_true_l[ci + 1];
	r_true_l[cij] = r_elem[0];
	r_true_l[cij + 1] = r_elem[1];
	r_elem[0] = r_true_t[ci];
	r_elem[1] = r_true_t[ci + 1];
	r_true_t[cij] = r_elem[0];
	r_true_t[cij + 1] = r_elem[1];
      }
    }
  } else {


    double *bb_vec;


    bb_vec = (double *) blas_malloc(k * sizeof(double) * 2);
    if (k > 0 && bb_vec == NULL) {
      printf("malloc failed\n");
      exit(-1);
    }


    if (alpha_flag == 0) {
      ((double *) c_elem)[0] = (float) drand48();
      ((double *) c_elem)[1] = (float) drand48();
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      ((double *) c_elem)[0] = (float) drand48();
      ((double *) c_elem)[1] = (float) drand48();
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }



    incai *= 2;
    incaij *= 2;

    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	b_elem = (float) drand48();
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	((double *) a_elem)[0] = (float) drand48();
	((double *) a_elem)[1] = (float) drand48();
	a_i[aij] = a_elem[0];
	a_i[aij + 1] = a_elem[1];
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      zgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);



      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	dgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);

	{
	  int r;
	  for (r = 0; r < k; r++) {
	    bb_vec[2 * r] = b_vec[r];
	    bb_vec[2 * r + 1] = 0.0;
	  }
	}

	BLAS_zdot_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1,
			  bb_vec,
			  a_vec, seed, c_elem, r_true_l_elem, r_true_t_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_true_l[cij] = r_true_l_elem[0];
	r_true_l[cij + 1] = r_true_l_elem[1];
	r_true_t[cij] = r_true_t_elem[0];
	r_true_t[cij + 1] = r_true_t_elem[1];

      }
    }


    free(bb_vec);

  }

  free(a_vec);
  free(b_vec);

}
void BLAS_zgemm_d_z_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    double *a, int lda, void *beta, int beta_flag,
			    void *b, int ldb, void *c, int ldc, int *seed,
			    double *r_true_l, double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  double c_elem[2];
  double a_elem;
  double b_elem[2];

  double r_elem[2];
  double r_true_t_elem[2];
  double r_true_l_elem[2];

  double *c_i = (double *) c;
  double *b_i = (double *) b;
  double *a_i = a;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;

  /* Temporary storage space for vectors of length k */
  double *a_vec;
  double *b_vec;

  inca = incb = 1;

  incb *= 2;

  a_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (double *) blas_malloc(k * sizeof(double) * 2);
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
    b_vec[i + 1] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;

  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_zdot_z_d_testgen(k, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  r_true_l_elem, r_true_t_elem);
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    r_true_l[cij] = r_true_l_elem[0];
    r_true_l[cij + 1] = r_true_l_elem[1];
    r_true_t[cij] = r_true_t_elem[0];
    r_true_t[cij + 1] = r_true_t_elem[1];

    /* Copy a_vec to the first row of A */
    dgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      zgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_zdot_z_d_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, r_true_l_elem, r_true_t_elem);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      r_true_l[cij] = r_true_l_elem[0];
      r_true_l[cij + 1] = r_true_l_elem[1];
      r_true_t[cij] = r_true_t_elem[0];
      r_true_t[cij + 1] = r_true_t_elem[1];
      dgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_elem[0] = r_true_l[ci];
	r_elem[1] = r_true_l[ci + 1];
	r_true_l[cij] = r_elem[0];
	r_true_l[cij + 1] = r_elem[1];
	r_elem[0] = r_true_t[ci];
	r_elem[1] = r_true_t[ci + 1];
	r_true_t[cij] = r_elem[0];
	r_true_t[cij + 1] = r_elem[1];
      }
    }
  } else {

    double *aa_vec;


    aa_vec = (double *) blas_malloc(k * sizeof(double) * 2);
    if (k > 0 && aa_vec == NULL) {
      printf("malloc failed\n");
      exit(-1);
    }



    if (alpha_flag == 0) {
      ((double *) c_elem)[0] = (float) drand48();
      ((double *) c_elem)[1] = (float) drand48();
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      ((double *) c_elem)[0] = (float) drand48();
      ((double *) c_elem)[1] = (float) drand48();
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }

    incbi *= 2;
    incbij *= 2;



    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	((double *) b_elem)[0] = (float) drand48();
	((double *) b_elem)[1] = (float) drand48();
	b_i[bij] = b_elem[0];
	b_i[bij + 1] = b_elem[1];
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	a_elem = (float) drand48();
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      dgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);

      {
	int r;
	for (r = 0; r < k; r++) {
	  aa_vec[2 * r] = a_vec[r];
	  aa_vec[2 * r + 1] = 0.0;
	}
      }

      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	zgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);



	BLAS_zdot_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1,
			  b_vec,
			  aa_vec, seed, c_elem, r_true_l_elem, r_true_t_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_true_l[cij] = r_true_l_elem[0];
	r_true_l[cij + 1] = r_true_l_elem[1];
	r_true_t[cij] = r_true_t_elem[0];
	r_true_t[cij + 1] = r_true_t_elem[1];

      }
    }

    free(aa_vec);


  }

  free(a_vec);
  free(b_vec);

}
void BLAS_zgemm_d_d_testgen(int norm, enum blas_order_type order,
			    enum blas_trans_type transa,
			    enum blas_trans_type transb, int m, int n, int k,
			    int randomize, void *alpha, int alpha_flag,
			    double *a, int lda, void *beta, int beta_flag,
			    double *b, int ldb, void *c, int ldc, int *seed,
			    double *r_true_l, double *r_true_t)
{

  int i, j;
  int cij, ci;
  int inccij, incci;
  int bij, bi;
  int incbij, incbi;
  int aij, ai;
  int incaij, incai;
  int inca, incb;

  double c_elem[2];
  double a_elem;
  double b_elem;

  double r_elem[2];
  double r_true_t_elem[2];
  double r_true_l_elem[2];

  double *c_i = (double *) c;
  double *b_i = b;
  double *a_i = a;
  double *alpha_i = (double *) alpha;
  double *beta_i = (double *) beta;

  /* Temporary storage space for vectors of length k */
  double *a_vec;
  double *b_vec;

  inca = incb = 1;



  a_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && a_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * inca; i += inca) {
    a_vec[i] = 0.0;
  }
  b_vec = (double *) blas_malloc(k * sizeof(double));
  if (k > 0 && b_vec == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  for (i = 0; i < k * incb; i += incb) {
    b_vec[i] = 0.0;
  }

  srand48(*seed);

  if (order == blas_colmajor) {
    incci = 1;
    inccij = ldc;
  } else {
    incci = ldc;
    inccij = 1;
  }

  incci *= 2;
  inccij *= 2;

  if (randomize == 0) {

    /* First fill in the first row of A and the first column of B */
    cij = 0;
    BLAS_zdot_d_d_testgen(k, 0, 0, norm, blas_no_conj,
			  alpha, alpha_flag, beta, beta_flag,
			  b_vec, a_vec, seed, c_elem,
			  r_true_l_elem, r_true_t_elem);
    c_i[cij] = c_elem[0];
    c_i[cij + 1] = c_elem[1];
    r_true_l[cij] = r_true_l_elem[0];
    r_true_l[cij + 1] = r_true_l_elem[1];
    r_true_t[cij] = r_true_t_elem[0];
    r_true_t[cij + 1] = r_true_t_elem[1];

    /* Copy a_vec to the first row of A */
    dgemm_commit_row(order, transa, m, k, a, lda, a_vec, 0);

    /* set every column of B to be b_vec */
    for (j = 0; j < n; j++)
      dgemm_commit_col(order, transb, k, n, b, ldb, b_vec, j);

    /* Next fill in the rest of matrix A.
       Note that alpha and beta are now fixed.   */
    cij = incci;
    for (i = 1; i < m; i++, cij += incci) {
      BLAS_zdot_d_d_testgen(k, 0, k, norm, blas_no_conj, alpha, 1,
			    beta, 1, b_vec, a_vec, seed,
			    c_elem, r_true_l_elem, r_true_t_elem);

      c_i[cij] = c_elem[0];
      c_i[cij + 1] = c_elem[1];
      r_true_l[cij] = r_true_l_elem[0];
      r_true_l[cij + 1] = r_true_l_elem[1];
      r_true_t[cij] = r_true_t_elem[0];
      r_true_t[cij + 1] = r_true_t_elem[1];
      dgemm_commit_row(order, transa, m, k, a, lda, a_vec, i);
    }

    /* Now fill in c and r_true */
    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      for (j = 1, cij = ci + inccij; j < n; j++, cij += inccij) {
	c_elem[0] = c_i[ci];
	c_elem[1] = c_i[ci + 1];
	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_elem[0] = r_true_l[ci];
	r_elem[1] = r_true_l[ci + 1];
	r_true_l[cij] = r_elem[0];
	r_true_l[cij + 1] = r_elem[1];
	r_elem[0] = r_true_t[ci];
	r_elem[1] = r_true_t[ci + 1];
	r_true_t[cij] = r_elem[0];
	r_true_t[cij + 1] = r_elem[1];
      }
    }
  } else {

    double *aa_vec;
    double *bb_vec;

    aa_vec = (double *) blas_malloc(k * sizeof(double) * 2);
    if (k > 0 && aa_vec == NULL) {
      printf("malloc failed\n");
      exit(-1);
    }
    bb_vec = (double *) blas_malloc(k * sizeof(double) * 2);
    if (k > 0 && bb_vec == NULL) {
      printf("malloc failed\n");
      exit(-1);
    }


    if (alpha_flag == 0) {
      ((double *) c_elem)[0] = (float) drand48();
      ((double *) c_elem)[1] = (float) drand48();
      alpha_i[0] = c_elem[0];
      alpha_i[0 + 1] = c_elem[1];
    }
    if (beta_flag == 0) {
      ((double *) c_elem)[0] = (float) drand48();
      ((double *) c_elem)[1] = (float) drand48();
      beta_i[0] = c_elem[0];
      beta_i[0 + 1] = c_elem[1];
    }

    if ((order == blas_colmajor && transb == blas_no_trans) ||
	(order == blas_rowmajor && transb != blas_no_trans)) {
      incbi = 1;
      incbij = ldb;
    } else {
      incbi = ldb;
      incbij = 1;
    }

    if ((order == blas_colmajor && transa == blas_no_trans) ||
	(order == blas_rowmajor && transa != blas_no_trans)) {
      incai = 1;
      incaij = lda;
    } else {
      incai = lda;
      incaij = 1;
    }






    for (i = 0, bi = 0; i < k; i++, bi += incbi) {
      for (j = 0, bij = bi; j < n; j++, bij += incbij) {
	b_elem = (float) drand48();
	b_i[bij] = b_elem;
      }
    }

    for (i = 0, ai = 0; i < m; i++, ai += incai) {
      for (j = 0, aij = ai; j < k; j++, aij += incaij) {
	a_elem = (float) drand48();
	a_i[aij] = a_elem;
      }
    }

    for (i = 0, ci = 0; i < m; i++, ci += incci) {
      dgemm_copy_row(order, transa, m, k, a, lda, a_vec, i);

      {
	int r;
	for (r = 0; r < k; r++) {
	  aa_vec[2 * r] = a_vec[r];
	  aa_vec[2 * r + 1] = 0.0;
	}
      }

      for (j = 0, cij = ci; j < n; j++, cij += inccij) {
	dgemm_copy_col(order, transb, k, n, b, ldb, b_vec, j);

	{
	  int r;
	  for (r = 0; r < k; r++) {
	    bb_vec[2 * r] = b_vec[r];
	    bb_vec[2 * r + 1] = 0.0;
	  }
	}

	BLAS_zdot_testgen(k, k, 0, norm, blas_no_conj, alpha, 1,
			  beta, 1,
			  bb_vec,
			  aa_vec, seed, c_elem, r_true_l_elem, r_true_t_elem);

	c_i[cij] = c_elem[0];
	c_i[cij + 1] = c_elem[1];
	r_true_l[cij] = r_true_l_elem[0];
	r_true_l[cij + 1] = r_true_l_elem[1];
	r_true_t[cij] = r_true_t_elem[0];
	r_true_t[cij + 1] = r_true_t_elem[1];

      }
    }

    free(aa_vec);
    free(bb_vec);

  }

  free(a_vec);
  free(b_vec);

}
