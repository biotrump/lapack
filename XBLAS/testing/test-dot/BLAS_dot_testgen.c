#include <stdio.h>
#include "blas_extended.h"
#include "cblas_test.h"

void BLAS_sdot_testgen(int n, int n_fix2, int n_mix, int norm,
		       enum blas_conj_type conj,
		       float *alpha, int alpha_flag,
		       float *beta, int beta_flag,
		       float *x, float *y, int *seed, float *r,
		       double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_sdot{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) float*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) float*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) float*
 *
 * y       (input/output) float*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) float*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  float alpha_i;
  float beta_i;
  float r_i;
  float *x_i;
  float *y_i;

  alpha_i = *alpha;
  beta_i = *beta;

  x_i = (float *) blas_malloc(2 * n * sizeof(float));
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    x_i[i] = x[i];
    y_i[i] = y[i];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    x_i[i] = x[i];
  }

  /* Call generator now. */
  testgen_BLAS_sdot(n, n_fix2, n_mix, norm, conj,
		    &alpha_i, alpha_flag,
		    &beta_i, beta_flag,
		    x_i, y_i, seed, &r_i, r_true_l, r_true_t);

  *alpha = alpha_i;
  *beta = beta_i;
  *r = r_i;
  for (i = 0; i < inc * n; i += inc) {
    x[i] = x_i[i];
    y[i] = y_i[i];
  }

  blas_free(x_i);		/* also y_i */
}				/* end BLAS_sdot_testgen */

void BLAS_ddot_testgen(int n, int n_fix2, int n_mix, int norm,
		       enum blas_conj_type conj,
		       double *alpha, int alpha_flag,
		       double *beta, int beta_flag,
		       double *x, double *y, int *seed, double *r,
		       double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_ddot{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) double*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) double*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) double*
 *
 * y       (input/output) double*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) double*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  double alpha_i;
  double beta_i;
  double r_i;
  double *x_i;
  double *y_i;

  alpha_i = *alpha;
  beta_i = *beta;

  x_i = (double *) blas_malloc(2 * n * sizeof(double));
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    x_i[i] = x[i];
    y_i[i] = y[i];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    x_i[i] = x[i];
  }

  /* Call generator now. */
  testgen_BLAS_ddot(n, n_fix2, n_mix, norm, conj,
		    &alpha_i, alpha_flag,
		    &beta_i, beta_flag,
		    x_i, y_i, seed, &r_i, r_true_l, r_true_t);

  *alpha = alpha_i;
  *beta = beta_i;
  *r = r_i;
  for (i = 0; i < inc * n; i += inc) {
    x[i] = x_i[i];
    y[i] = y_i[i];
  }

  blas_free(x_i);		/* also y_i */
}				/* end BLAS_ddot_testgen */

void BLAS_ddot_s_s_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj,
			   double *alpha, int alpha_flag,
			   double *beta, int beta_flag,
			   float *x, float *y, int *seed, double *r,
			   double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_ddot_s_s{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) double*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) double*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) float*
 *
 * y       (input/output) float*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) double*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  float alpha_i;
  float beta_i;
  float r_i;
  float *x_i;
  float *y_i;

  alpha_i = *alpha;
  beta_i = *beta;

  x_i = (float *) blas_malloc(2 * n * sizeof(float));
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    x_i[i] = x[i];
    y_i[i] = y[i];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    x_i[i] = x[i];
  }

  /* Call generator now. */
  testgen_BLAS_sdot(n, n_fix2, n_mix, norm, conj,
		    &alpha_i, alpha_flag,
		    &beta_i, beta_flag,
		    x_i, y_i, seed, &r_i, r_true_l, r_true_t);

  *alpha = alpha_i;
  *beta = beta_i;
  *r = r_i;
  for (i = 0; i < inc * n; i += inc) {
    x[i] = x_i[i];
    y[i] = y_i[i];
  }

  blas_free(x_i);		/* also y_i */
}				/* end BLAS_ddot_s_s_testgen */

void BLAS_ddot_s_d_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj,
			   double *alpha, int alpha_flag,
			   double *beta, int beta_flag,
			   float *x, double *y, int *seed, double *r,
			   double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_ddot_s_d{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) double*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) double*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) float*
 *
 * y       (input/output) double*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) double*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  float alpha_i;
  float beta_i;
  float r_i;
  float *x_i;
  float *y_i;

  alpha_i = *alpha;
  beta_i = *beta;

  x_i = (float *) blas_malloc(2 * n * sizeof(float));
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    x_i[i] = x[i];
    y_i[i] = y[i];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    x_i[i] = x[i];
  }

  /* Call generator now. */
  testgen_BLAS_sdot(n, n_fix2, n_mix, norm, conj,
		    &alpha_i, alpha_flag,
		    &beta_i, beta_flag,
		    x_i, y_i, seed, &r_i, r_true_l, r_true_t);

  *alpha = alpha_i;
  *beta = beta_i;
  *r = r_i;
  for (i = 0; i < inc * n; i += inc) {
    x[i] = x_i[i];
    y[i] = y_i[i];
  }

  blas_free(x_i);		/* also y_i */
}				/* end BLAS_ddot_s_d_testgen */

void BLAS_ddot_d_s_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj,
			   double *alpha, int alpha_flag,
			   double *beta, int beta_flag,
			   double *x, float *y, int *seed, double *r,
			   double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_ddot_d_s{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) double*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) double*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) double*
 *
 * y       (input/output) float*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) double*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  float alpha_i;
  float beta_i;
  float r_i;
  float *x_i;
  float *y_i;

  alpha_i = *alpha;
  beta_i = *beta;

  x_i = (float *) blas_malloc(2 * n * sizeof(float));
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    x_i[i] = x[i];
    y_i[i] = y[i];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    x_i[i] = x[i];
  }

  /* Call generator now. */
  testgen_BLAS_sdot(n, n_fix2, n_mix, norm, conj,
		    &alpha_i, alpha_flag,
		    &beta_i, beta_flag,
		    x_i, y_i, seed, &r_i, r_true_l, r_true_t);

  *alpha = alpha_i;
  *beta = beta_i;
  *r = r_i;
  for (i = 0; i < inc * n; i += inc) {
    x[i] = x_i[i];
    y[i] = y_i[i];
  }

  blas_free(x_i);		/* also y_i */
}				/* end BLAS_ddot_d_s_testgen */

void BLAS_cdot_testgen(int n, int n_fix2, int n_mix, int norm,
		       enum blas_conj_type conj,
		       void *alpha, int alpha_flag,
		       void *beta, int beta_flag,
		       void *x, void *y, int *seed, void *r,
		       double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_cdot{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  float alpha_i[2];
  float beta_i[2];
  float r_i[2];
  float *x_i;
  float *y_i;

  alpha_i[0] = ((float *) alpha)[0];
  alpha_i[1] = ((float *) alpha)[1];
  beta_i[0] = ((float *) beta)[0];
  beta_i[1] = ((float *) beta)[1];
  inc *= 2;
  x_i = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    ((float *) x_i)[i] = ((float *) x)[i];
    ((float *) x_i)[i + 1] = ((float *) x)[i + 1];
    ((float *) y_i)[i] = ((float *) y)[i];
    ((float *) y_i)[i + 1] = ((float *) y)[i + 1];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    ((float *) x_i)[i] = ((float *) x)[i];
    ((float *) x_i)[i + 1] = ((float *) x)[i + 1];
  }

  /* Call generator now. */
  testgen_BLAS_cdot(n, n_fix2, n_mix, norm, conj,
		    alpha_i, alpha_flag,
		    beta_i, beta_flag,
		    x_i, y_i, seed, r_i, r_true_l, r_true_t);

  ((float *) alpha)[0] = alpha_i[0];
  ((float *) alpha)[1] = alpha_i[1];
  ((float *) beta)[0] = beta_i[0];
  ((float *) beta)[1] = beta_i[1];
  ((float *) r)[0] = r_i[0];
  ((float *) r)[1] = r_i[1];
  for (i = 0; i < inc * n; i += inc) {
    ((float *) x)[i] = ((float *) x_i)[i];
    ((float *) x)[i + 1] = ((float *) x_i)[i + 1];
    ((float *) y)[i] = ((float *) y_i)[i];
    ((float *) y)[i + 1] = ((float *) y_i)[i + 1];
  }

  blas_free(x_i);		/* also y_i */
}				/* end BLAS_cdot_testgen */

void BLAS_zdot_testgen(int n, int n_fix2, int n_mix, int norm,
		       enum blas_conj_type conj,
		       void *alpha, int alpha_flag,
		       void *beta, int beta_flag,
		       void *x, void *y, int *seed, void *r,
		       double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zdot{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  double alpha_i[2];
  double beta_i[2];
  double r_i[2];
  double *x_i;
  double *y_i;

  alpha_i[0] = ((double *) alpha)[0];
  alpha_i[1] = ((double *) alpha)[1];
  beta_i[0] = ((double *) beta)[0];
  beta_i[1] = ((double *) beta)[1];
  inc *= 2;
  x_i = (double *) blas_malloc(2 * n * sizeof(double) * 2);
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    ((double *) x_i)[i] = ((double *) x)[i];
    ((double *) x_i)[i + 1] = ((double *) x)[i + 1];
    ((double *) y_i)[i] = ((double *) y)[i];
    ((double *) y_i)[i + 1] = ((double *) y)[i + 1];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    ((double *) x_i)[i] = ((double *) x)[i];
    ((double *) x_i)[i + 1] = ((double *) x)[i + 1];
  }

  /* Call generator now. */
  testgen_BLAS_zdot(n, n_fix2, n_mix, norm, conj,
		    alpha_i, alpha_flag,
		    beta_i, beta_flag,
		    x_i, y_i, seed, r_i, r_true_l, r_true_t);

  ((double *) alpha)[0] = alpha_i[0];
  ((double *) alpha)[1] = alpha_i[1];
  ((double *) beta)[0] = beta_i[0];
  ((double *) beta)[1] = beta_i[1];
  ((double *) r)[0] = r_i[0];
  ((double *) r)[1] = r_i[1];
  for (i = 0; i < inc * n; i += inc) {
    ((double *) x)[i] = ((double *) x_i)[i];
    ((double *) x)[i + 1] = ((double *) x_i)[i + 1];
    ((double *) y)[i] = ((double *) y_i)[i];
    ((double *) y)[i + 1] = ((double *) y_i)[i + 1];
  }

  blas_free(x_i);		/* also y_i */
}				/* end BLAS_zdot_testgen */

void BLAS_zdot_c_c_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj,
			   void *alpha, int alpha_flag,
			   void *beta, int beta_flag,
			   void *x, void *y, int *seed, void *r,
			   double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zdot_c_c{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  float alpha_i[2];
  float beta_i[2];
  float r_i[2];
  float *x_i;
  float *y_i;

  alpha_i[0] = ((double *) alpha)[0];
  alpha_i[1] = ((double *) alpha)[1];
  beta_i[0] = ((double *) beta)[0];
  beta_i[1] = ((double *) beta)[1];
  inc *= 2;
  x_i = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    ((float *) x_i)[i] = ((float *) x)[i];
    ((float *) x_i)[i + 1] = ((float *) x)[i + 1];
    ((float *) y_i)[i] = ((float *) y)[i];
    ((float *) y_i)[i + 1] = ((float *) y)[i + 1];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    ((float *) x_i)[i] = ((float *) x)[i];
    ((float *) x_i)[i + 1] = ((float *) x)[i + 1];
  }

  /* Call generator now. */
  testgen_BLAS_cdot(n, n_fix2, n_mix, norm, conj,
		    alpha_i, alpha_flag,
		    beta_i, beta_flag,
		    x_i, y_i, seed, r_i, r_true_l, r_true_t);

  ((double *) alpha)[0] = alpha_i[0];
  ((double *) alpha)[1] = alpha_i[1];
  ((double *) beta)[0] = beta_i[0];
  ((double *) beta)[1] = beta_i[1];
  ((double *) r)[0] = r_i[0];
  ((double *) r)[1] = r_i[1];
  for (i = 0; i < inc * n; i += inc) {
    ((float *) x)[i] = ((float *) x_i)[i];
    ((float *) x)[i + 1] = ((float *) x_i)[i + 1];
    ((float *) y)[i] = ((float *) y_i)[i];
    ((float *) y)[i + 1] = ((float *) y_i)[i + 1];
  }

  blas_free(x_i);		/* also y_i */
}				/* end BLAS_zdot_c_c_testgen */

void BLAS_zdot_c_z_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj,
			   void *alpha, int alpha_flag,
			   void *beta, int beta_flag,
			   void *x, void *y, int *seed, void *r,
			   double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zdot_c_z{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  float alpha_i[2];
  float beta_i[2];
  float r_i[2];
  float *x_i;
  float *y_i;

  alpha_i[0] = ((double *) alpha)[0];
  alpha_i[1] = ((double *) alpha)[1];
  beta_i[0] = ((double *) beta)[0];
  beta_i[1] = ((double *) beta)[1];
  inc *= 2;
  x_i = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    ((float *) x_i)[i] = ((float *) x)[i];
    ((float *) x_i)[i + 1] = ((float *) x)[i + 1];
    ((float *) y_i)[i] = ((double *) y)[i];
    ((float *) y_i)[i + 1] = ((double *) y)[i + 1];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    ((float *) x_i)[i] = ((float *) x)[i];
    ((float *) x_i)[i + 1] = ((float *) x)[i + 1];
  }

  /* Call generator now. */
  testgen_BLAS_cdot(n, n_fix2, n_mix, norm, conj,
		    alpha_i, alpha_flag,
		    beta_i, beta_flag,
		    x_i, y_i, seed, r_i, r_true_l, r_true_t);

  ((double *) alpha)[0] = alpha_i[0];
  ((double *) alpha)[1] = alpha_i[1];
  ((double *) beta)[0] = beta_i[0];
  ((double *) beta)[1] = beta_i[1];
  ((double *) r)[0] = r_i[0];
  ((double *) r)[1] = r_i[1];
  for (i = 0; i < inc * n; i += inc) {
    ((float *) x)[i] = ((float *) x_i)[i];
    ((float *) x)[i + 1] = ((float *) x_i)[i + 1];
    ((double *) y)[i] = ((float *) y_i)[i];
    ((double *) y)[i + 1] = ((float *) y_i)[i + 1];
  }

  blas_free(x_i);		/* also y_i */
}				/* end BLAS_zdot_c_z_testgen */

void BLAS_zdot_z_c_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj,
			   void *alpha, int alpha_flag,
			   void *beta, int beta_flag,
			   void *x, void *y, int *seed, void *r,
			   double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zdot_z_c{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i, inc = 1;
  float alpha_i[2];
  float beta_i[2];
  float r_i[2];
  float *x_i;
  float *y_i;

  alpha_i[0] = ((double *) alpha)[0];
  alpha_i[1] = ((double *) alpha)[1];
  beta_i[0] = ((double *) beta)[0];
  beta_i[1] = ((double *) beta)[1];
  inc *= 2;
  x_i = (float *) blas_malloc(2 * n * sizeof(float) * 2);
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + inc * n;
  for (i = 0; i < inc * n_fix2; i += inc) {
    ((float *) x_i)[i] = ((double *) x)[i];
    ((float *) x_i)[i + 1] = ((double *) x)[i + 1];
    ((float *) y_i)[i] = ((float *) y)[i];
    ((float *) y_i)[i + 1] = ((float *) y)[i + 1];
  }
  for (; i < inc * (n_fix2 + n_mix); i += inc) {
    ((float *) x_i)[i] = ((double *) x)[i];
    ((float *) x_i)[i + 1] = ((double *) x)[i + 1];
  }

  /* Call generator now. */
  testgen_BLAS_cdot(n, n_fix2, n_mix, norm, conj,
		    alpha_i, alpha_flag,
		    beta_i, beta_flag,
		    x_i, y_i, seed, r_i, r_true_l, r_true_t);

  ((double *) alpha)[0] = alpha_i[0];
  ((double *) alpha)[1] = alpha_i[1];
  ((double *) beta)[0] = beta_i[0];
  ((double *) beta)[1] = beta_i[1];
  ((double *) r)[0] = r_i[0];
  ((double *) r)[1] = r_i[1];
  for (i = 0; i < inc * n; i += inc) {
    ((double *) x)[i] = ((float *) x_i)[i];
    ((double *) x)[i + 1] = ((float *) x_i)[i + 1];
    ((float *) y)[i] = ((float *) y_i)[i];
    ((float *) y)[i + 1] = ((float *) y_i)[i + 1];
  }

  blas_free(x_i);		/* also y_i */
}				/* end BLAS_zdot_z_c_testgen */

void BLAS_cdot_s_s_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj,
			   void *alpha, int alpha_flag,
			   void *beta, int beta_flag,
			   float *x, float *y, int *seed, void *r,
			   double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_cdot_s_s{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) float*
 *
 * y       (input/output) float*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i;
  float alpha_i_r;
  float alpha_i_i;
  float beta_i_r;
  float beta_i_i;
  float r_i;
  float *x_i;
  float *y_i;

  alpha_i_r = ((float *) alpha)[0];
  alpha_i_i = ((float *) alpha)[1];
  beta_i_r = ((float *) beta)[0];
  beta_i_i = ((float *) beta)[1];
  x_i = (float *) blas_malloc(2 * n * sizeof(float));
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + n;
  for (i = 0; i < n_fix2; i++) {
    x_i[i] = x[i];
    y_i[i] = y[i];
  }
  for (; i < n_fix2 + n_mix; i++) {
    x_i[i] = x[i];
  }

  /* Call generator now. */
  testgen_BLAS_sdot(n, n_fix2, n_mix, norm, conj,
		    &alpha_i_r, alpha_flag,
		    &beta_i_r, beta_flag,
		    x_i, y_i, seed, &r_i, &r_true_l[0], &r_true_t[0]);

  if (alpha_flag == 1) {	/* alpha is fixed */
    if (alpha_i_r == 1.0 && alpha_i_i == 0.) {	/* alpha == 1.0 */
      if (beta_flag == 1 && ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.))) {	/* beta == 0 or 1 */
	((float *) r)[0] = r_i;
	((float *) r)[1] = 0.0;
      } else {			/* beta *= (1-i), r *= (1+i)/2 --> prod = 1 */
	((float *) beta)[0] = beta_i_r;
	((float *) beta)[1] = -beta_i_r;
	((float *) r)[0] = r_i / 2.;
	((float *) r)[1] = r_i / 2.;
      }
      r_true_l[1] = r_true_t[1] = 0.0;
    } else if (alpha_i_r == 0. && alpha_i_i == 0.) {	/* alpha == 0.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta == 0 or 1 --> r *= (1+i) */
	((float *) r)[0] = r_i;
	((float *) r)[1] = r_i;
      } else {			/* beta *= (1-i), r *= i --> prod = 1+i */
	((float *) beta)[0] = beta_i_r;
	((float *) beta)[1] = -beta_i_r;
	((float *) r)[0] = 0.0;
	((float *) r)[1] = r_i;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else {			/* alpha is a fixed multiple of (1+i) */
      ((float *) alpha)[0] = alpha_i_r;
      ((float *) alpha)[1] = alpha_i_r;
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta == 0 or 1 --> r *= (1+i) */
	((float *) r)[0] = r_i;
	((float *) r)[1] = r_i;
      } else {			/* beta *= (1-i), r *= i --> prod = 1+i */
	((float *) beta)[0] = beta_i_r;
	((float *) beta)[1] = -beta_i_r;
	((float *) r)[0] = 0.0;
	((float *) r)[1] = r_i;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    }
  } else if (beta_flag == 1) {	/* alpha is free, beta is fixed */
    /* alpha *= (1+i) */
    ((float *) alpha)[0] = alpha_i_r;
    ((float *) alpha)[1] = alpha_i_r;
    if ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.)) {	/* r *= (1+i) */
      ((float *) r)[0] = r_i;
      ((float *) r)[1] = r_i;
    } else {			/* beta *= (1-i), r *= i */
      ((float *) beta)[0] = beta_i_r;
      ((float *) beta)[1] = -beta_i_r;
      ((float *) r)[0] = 0.;
      ((float *) r)[1] = r_i;
    }
    r_true_l[1] = r_true_l[0];
    r_true_t[1] = r_true_t[0];
  } else {			/* both alpha and beta are free */
    assert(alpha_flag == 0 && beta_flag == 0);
    ((float *) alpha)[0] = alpha_i_r;
    ((float *) alpha)[1] = alpha_i_r;
    ((float *) beta)[0] = beta_i_r;
    ((float *) beta)[1] = -beta_i_r;
    ((float *) r)[0] = 0;
    ((float *) r)[1] = r_i;
    /* imaginary part of r_true */
    r_true_l[1] = r_true_l[0];
    r_true_t[1] = r_true_t[0];
  }

  for (i = 0; i < n; ++i) {
    x[i] = x_i[i];
    y[i] = y_i[i];
  }

  blas_free(x_i);		/* also y_i */
}				/* end BLAS_cdot_s_s_testgen */

void BLAS_cdot_s_c_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj,
			   void *alpha, int alpha_flag,
			   void *beta, int beta_flag,
			   float *x, void *y, int *seed, void *r,
			   double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_cdot_s_c{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) float*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i;
  float alpha_i_r;
  float alpha_i_i;
  float beta_i_r;
  float beta_i_i;
  float r_i;
  float *x_i;
  float *y_i;

  alpha_i_r = ((float *) alpha)[0];
  alpha_i_i = ((float *) alpha)[1];
  beta_i_r = ((float *) beta)[0];
  beta_i_i = ((float *) beta)[1];
  x_i = (float *) blas_malloc(2 * n * sizeof(float));
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + n;
  for (i = 0; i < n_fix2; i++) {
    x_i[i] = x[i];
    ((float *) y_i)[i] = ((float *) y)[2 * i];
  }
  for (; i < n_fix2 + n_mix; i++) {
    x_i[i] = x[i];
  }

  /* Call generator now. */
  testgen_BLAS_sdot(n, n_fix2, n_mix, norm, conj,
		    &alpha_i_r, alpha_flag,
		    &beta_i_r, beta_flag,
		    x_i, y_i, seed, &r_i, &r_true_l[0], &r_true_t[0]);

  if (alpha_flag == 1) {	/* alpha is fixed */
    if (alpha_i_r == 1.0 && alpha_i_i == 0.) {	/* alpha == 1.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta == 0 or 1 --> r *= (1+i) */
	((float *) r)[0] = r_i;
	((float *) r)[1] = r_i;
      } else {			/* beta *= (1-i), r *= i --> prod = 1+i */
	((float *) beta)[0] = beta_i_r;
	((float *) beta)[1] = -beta_i_r;
	((float *) r)[0] = 0.0;
	((float *) r)[1] = r_i;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else if (alpha_i_r == 0. && alpha_i_i == 0.) {	/* alpha == 0.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta == 0 or 1 --> r *= (1+i) */
	((float *) r)[0] = r_i;
	((float *) r)[1] = r_i;
      } else {			/* beta *= (1-i), r *= i --> prod = 1+i */
	((float *) beta)[0] = beta_i_r;
	((float *) beta)[1] = -beta_i_r;
	((float *) r)[0] = 0.0;
	((float *) r)[1] = r_i;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else {			/* alpha is a fixed multiple of (1+i) */
      ((float *) alpha)[0] = alpha_i_r;
      ((float *) alpha)[1] = alpha_i_r;
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta is 0 or 1 --> r *= 2i --> prod = 2i */
	((float *) r)[0] = 0.0;
	((float *) r)[1] = 2.0 * r_i;
      } else {			/* beta *= (1+i), r *= (1+i) --> prod = 2i */
	((float *) beta)[0] = beta_i_r;
	((float *) beta)[1] = beta_i_r;
	((float *) r)[0] = r_i;
	((float *) r)[1] = r_i;
      }
      r_true_l[1] = 2.0 * r_true_l[0];
      r_true_t[1] = 2.0 * r_true_t[0];
      r_true_l[0] = r_true_t[0] = 0.0;
    }
  } else if (beta_flag == 1) {	/* alpha is free, beta is fixed */
    /* alpha *= (1+i) */
    ((float *) alpha)[0] = alpha_i_r;
    ((float *) alpha)[1] = alpha_i_r;
    if ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.)) {	/* r*=2i --> prod = 2i */
      ((float *) r)[0] = 0.0;
      ((float *) r)[1] = 2.0 * r_i;
    } else {			/* beta *= (1+i), r *= (1+i) */
      ((float *) beta)[0] = beta_i_r;
      ((float *) beta)[1] = beta_i_r;
      ((float *) r)[0] = r_i;
      ((float *) r)[1] = r_i;
    }
    r_true_l[1] = 2.0 * r_true_l[0];
    r_true_t[1] = 2.0 * r_true_t[0];
    r_true_l[0] = r_true_t[0] = 0.0;
  } else {			/* both alpha and beta are free */
    assert(alpha_flag == 0 && beta_flag == 0);
    ((float *) alpha)[0] = alpha_i_r;
    ((float *) alpha)[1] = alpha_i_r;
    ((float *) beta)[0] = beta_i_r;
    ((float *) beta)[1] = beta_i_r;
    ((float *) r)[0] = r_i;
    ((float *) r)[1] = r_i;
    /* imaginary part of r_true */
    ddmuld(r_true_l[0], r_true_t[0], 2.0, &r_true_l[1], &r_true_t[1]);
    /* real part of r_true */
    r_true_l[0] = 0.;
    r_true_t[0] = 0.;
  }

  for (i = 0; i < n; ++i) {
    x[i] = x_i[i];
    ((float *) y)[2 * i] = ((float *) y_i)[i];
    ((float *) y)[2 * i + 1] = ((float *) y_i)[i];
  }

  blas_free(x_i);		/* also y_i */
}				/* end BLAS_cdot_s_c_testgen */

void BLAS_cdot_c_s_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj,
			   void *alpha, int alpha_flag,
			   void *beta, int beta_flag,
			   void *x, float *y, int *seed, void *r,
			   double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_cdot_c_s{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) float*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i;
  float alpha_i_r;
  float alpha_i_i;
  float beta_i_r;
  float beta_i_i;
  float r_i;
  float *x_i;
  float *y_i;

  alpha_i_r = ((float *) alpha)[0];
  alpha_i_i = ((float *) alpha)[1];
  beta_i_r = ((float *) beta)[0];
  beta_i_i = ((float *) beta)[1];
  x_i = (float *) blas_malloc(2 * n * sizeof(float));
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + n;
  for (i = 0; i < n_fix2; i++) {
    ((float *) x_i)[i] = ((float *) x)[2 * i];
    y_i[i] = y[i];
  }
  for (; i < n_fix2 + n_mix; i++) {
    ((float *) x_i)[i] = ((float *) x)[2 * i];
  }

  /* Call generator now. */
  testgen_BLAS_sdot(n, n_fix2, n_mix, norm, conj,
		    &alpha_i_r, alpha_flag,
		    &beta_i_r, beta_flag,
		    x_i, y_i, seed, &r_i, &r_true_l[0], &r_true_t[0]);

  if (alpha_flag == 1) {	/* alpha is fixed */
    if (alpha_i_r == 1.0 && alpha_i_i == 0.) {	/* alpha == 1.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta == 0 or 1 --> r *= (1+i) */
	((float *) r)[0] = r_i;
	((float *) r)[1] = r_i;
      } else {			/* beta *= (1-i), r *= i --> prod = 1+i */
	((float *) beta)[0] = beta_i_r;
	((float *) beta)[1] = -beta_i_r;
	((float *) r)[0] = 0.0;
	((float *) r)[1] = r_i;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else if (alpha_i_r == 0. && alpha_i_i == 0.) {	/* alpha == 0.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta == 0 or 1 --> r *= (1+i) */
	((float *) r)[0] = r_i;
	((float *) r)[1] = r_i;
      } else {			/* beta *= (1-i), r *= i --> prod = 1+i */
	((float *) beta)[0] = beta_i_r;
	((float *) beta)[1] = -beta_i_r;
	((float *) r)[0] = 0.0;
	((float *) r)[1] = r_i;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else {			/* alpha is a fixed multiple of (1+i) */
      ((float *) alpha)[0] = alpha_i_r;
      ((float *) alpha)[1] = alpha_i_r;
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta is 0 or 1 --> r *= 2i --> prod = 2i */
	((float *) r)[0] = 0.0;
	((float *) r)[1] = 2.0 * r_i;
      } else {			/* beta *= (1+i), r *= (1+i) --> prod = 2i */
	((float *) beta)[0] = beta_i_r;
	((float *) beta)[1] = beta_i_r;
	((float *) r)[0] = r_i;
	((float *) r)[1] = r_i;
      }
      r_true_l[1] = 2.0 * r_true_l[0];
      r_true_t[1] = 2.0 * r_true_t[0];
      r_true_l[0] = r_true_t[0] = 0.0;
    }
  } else if (beta_flag == 1) {	/* alpha is free, beta is fixed */
    /* alpha *= (1+i) */
    ((float *) alpha)[0] = alpha_i_r;
    ((float *) alpha)[1] = alpha_i_r;
    if ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.)) {	/* r*=2i --> prod = 2i */
      ((float *) r)[0] = 0.0;
      ((float *) r)[1] = 2.0 * r_i;
    } else {			/* beta *= (1+i), r *= (1+i) */
      ((float *) beta)[0] = beta_i_r;
      ((float *) beta)[1] = beta_i_r;
      ((float *) r)[0] = r_i;
      ((float *) r)[1] = r_i;
    }
    r_true_l[1] = 2.0 * r_true_l[0];
    r_true_t[1] = 2.0 * r_true_t[0];
    r_true_l[0] = r_true_t[0] = 0.0;
  } else {			/* both alpha and beta are free */
    assert(alpha_flag == 0 && beta_flag == 0);
    ((float *) alpha)[0] = alpha_i_r;
    ((float *) alpha)[1] = alpha_i_r;
    ((float *) beta)[0] = beta_i_r;
    ((float *) beta)[1] = beta_i_r;
    ((float *) r)[0] = r_i;
    ((float *) r)[1] = r_i;
    /* imaginary part of r_true */
    ddmuld(r_true_l[0], r_true_t[0], 2.0, &r_true_l[1], &r_true_t[1]);
    /* real part of r_true */
    r_true_l[0] = 0.;
    r_true_t[0] = 0.;
  }

  for (i = 0; i < n; ++i) {
    ((float *) x)[2 * i] = ((float *) x_i)[i];
    ((float *) x)[2 * i + 1] = ((float *) x_i)[i];
    y[i] = y_i[i];
  }
  if (conj == blas_conj) {
    for (i = 0; i < n; ++i)
      ((float *) x)[2 * i + 1] = -((float *) x)[2 * i + 1];
  }
  blas_free(x_i);		/* also y_i */
}				/* end BLAS_cdot_c_s_testgen */

void BLAS_zdot_d_d_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj,
			   void *alpha, int alpha_flag,
			   void *beta, int beta_flag,
			   double *x, double *y, int *seed, void *r,
			   double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zdot_d_d{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) double*
 *
 * y       (input/output) double*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i;
  double alpha_i_r;
  double alpha_i_i;
  double beta_i_r;
  double beta_i_i;
  double r_i;
  double *x_i;
  double *y_i;

  alpha_i_r = ((double *) alpha)[0];
  alpha_i_i = ((double *) alpha)[1];
  beta_i_r = ((double *) beta)[0];
  beta_i_i = ((double *) beta)[1];
  x_i = (double *) blas_malloc(2 * n * sizeof(double));
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + n;
  for (i = 0; i < n_fix2; i++) {
    x_i[i] = x[i];
    y_i[i] = y[i];
  }
  for (; i < n_fix2 + n_mix; i++) {
    x_i[i] = x[i];
  }

  /* Call generator now. */
  testgen_BLAS_ddot(n, n_fix2, n_mix, norm, conj,
		    &alpha_i_r, alpha_flag,
		    &beta_i_r, beta_flag,
		    x_i, y_i, seed, &r_i, &r_true_l[0], &r_true_t[0]);

  if (alpha_flag == 1) {	/* alpha is fixed */
    if (alpha_i_r == 1.0 && alpha_i_i == 0.) {	/* alpha == 1.0 */
      if (beta_flag == 1 && ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.))) {	/* beta == 0 or 1 */
	((double *) r)[0] = r_i;
	((double *) r)[1] = 0.0;
      } else {			/* beta *= (1-i), r *= (1+i)/2 --> prod = 1 */
	((double *) beta)[0] = beta_i_r;
	((double *) beta)[1] = -beta_i_r;
	((double *) r)[0] = r_i / 2.;
	((double *) r)[1] = r_i / 2.;
      }
      r_true_l[1] = r_true_t[1] = 0.0;
    } else if (alpha_i_r == 0. && alpha_i_i == 0.) {	/* alpha == 0.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta == 0 or 1 --> r *= (1+i) */
	((double *) r)[0] = r_i;
	((double *) r)[1] = r_i;
      } else {			/* beta *= (1-i), r *= i --> prod = 1+i */
	((double *) beta)[0] = beta_i_r;
	((double *) beta)[1] = -beta_i_r;
	((double *) r)[0] = 0.0;
	((double *) r)[1] = r_i;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else {			/* alpha is a fixed multiple of (1+i) */
      ((double *) alpha)[0] = alpha_i_r;
      ((double *) alpha)[1] = alpha_i_r;
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta == 0 or 1 --> r *= (1+i) */
	((double *) r)[0] = r_i;
	((double *) r)[1] = r_i;
      } else {			/* beta *= (1-i), r *= i --> prod = 1+i */
	((double *) beta)[0] = beta_i_r;
	((double *) beta)[1] = -beta_i_r;
	((double *) r)[0] = 0.0;
	((double *) r)[1] = r_i;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    }
  } else if (beta_flag == 1) {	/* alpha is free, beta is fixed */
    /* alpha *= (1+i) */
    ((double *) alpha)[0] = alpha_i_r;
    ((double *) alpha)[1] = alpha_i_r;
    if ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.)) {	/* r *= (1+i) */
      ((double *) r)[0] = r_i;
      ((double *) r)[1] = r_i;
    } else {			/* beta *= (1-i), r *= i */
      ((double *) beta)[0] = beta_i_r;
      ((double *) beta)[1] = -beta_i_r;
      ((double *) r)[0] = 0.;
      ((double *) r)[1] = r_i;
    }
    r_true_l[1] = r_true_l[0];
    r_true_t[1] = r_true_t[0];
  } else {			/* both alpha and beta are free */
    assert(alpha_flag == 0 && beta_flag == 0);
    ((double *) alpha)[0] = alpha_i_r;
    ((double *) alpha)[1] = alpha_i_r;
    ((double *) beta)[0] = beta_i_r;
    ((double *) beta)[1] = -beta_i_r;
    ((double *) r)[0] = 0;
    ((double *) r)[1] = r_i;
    /* imaginary part of r_true */
    r_true_l[1] = r_true_l[0];
    r_true_t[1] = r_true_t[0];
  }

  for (i = 0; i < n; ++i) {
    x[i] = x_i[i];
    y[i] = y_i[i];
  }

  blas_free(x_i);		/* also y_i */
}

	/* end BLAS_zdot_d_d_testgen */

void BLAS_zdot_z_d_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj,
			   void *alpha, int alpha_flag,
			   void *beta, int beta_flag,
			   void *x, double *y, int *seed, void *r,
			   double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zdot_z_d{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) void*
 *
 * y       (input/output) double*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i;
  double alpha_i_r;
  double alpha_i_i;
  double beta_i_r;
  double beta_i_i;
  double r_i;
  double *x_i;
  double *y_i;

  alpha_i_r = ((double *) alpha)[0];
  alpha_i_i = ((double *) alpha)[1];
  beta_i_r = ((double *) beta)[0];
  beta_i_i = ((double *) beta)[1];
  x_i = (double *) blas_malloc(2 * n * sizeof(double));
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + n;
  for (i = 0; i < n_fix2; i++) {
    ((double *) x_i)[i] = ((double *) x)[2 * i];
    y_i[i] = y[i];
  }
  for (; i < n_fix2 + n_mix; i++) {
    ((double *) x_i)[i] = ((double *) x)[2 * i];
  }

  /* Call generator now. */
  testgen_BLAS_ddot(n, n_fix2, n_mix, norm, conj,
		    &alpha_i_r, alpha_flag,
		    &beta_i_r, beta_flag,
		    x_i, y_i, seed, &r_i, &r_true_l[0], &r_true_t[0]);

  if (alpha_flag == 1) {	/* alpha is fixed */
    if (alpha_i_r == 1.0 && alpha_i_i == 0.) {	/* alpha == 1.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta == 0 or 1 --> r *= (1+i) */
	((double *) r)[0] = r_i;
	((double *) r)[1] = r_i;
      } else {			/* beta *= (1-i), r *= i --> prod = 1+i */
	((double *) beta)[0] = beta_i_r;
	((double *) beta)[1] = -beta_i_r;
	((double *) r)[0] = 0.0;
	((double *) r)[1] = r_i;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else if (alpha_i_r == 0. && alpha_i_i == 0.) {	/* alpha == 0.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta == 0 or 1 --> r *= (1+i) */
	((double *) r)[0] = r_i;
	((double *) r)[1] = r_i;
      } else {			/* beta *= (1-i), r *= i --> prod = 1+i */
	((double *) beta)[0] = beta_i_r;
	((double *) beta)[1] = -beta_i_r;
	((double *) r)[0] = 0.0;
	((double *) r)[1] = r_i;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else {			/* alpha is a fixed multiple of (1+i) */
      ((double *) alpha)[0] = alpha_i_r;
      ((double *) alpha)[1] = alpha_i_r;
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta is 0 or 1 --> r *= 2i --> prod = 2i */
	((double *) r)[0] = 0.0;
	((double *) r)[1] = 2.0 * r_i;
      } else {			/* beta *= (1+i), r *= (1+i) --> prod = 2i */
	((double *) beta)[0] = beta_i_r;
	((double *) beta)[1] = beta_i_r;
	((double *) r)[0] = r_i;
	((double *) r)[1] = r_i;
      }
      r_true_l[1] = 2.0 * r_true_l[0];
      r_true_t[1] = 2.0 * r_true_t[0];
      r_true_l[0] = r_true_t[0] = 0.0;
    }
  } else if (beta_flag == 1) {	/* alpha is free, beta is fixed */
    /* alpha *= (1+i) */
    ((double *) alpha)[0] = alpha_i_r;
    ((double *) alpha)[1] = alpha_i_r;
    if ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.)) {	/* r*=2i --> prod = 2i */
      ((double *) r)[0] = 0.0;
      ((double *) r)[1] = 2.0 * r_i;
    } else {			/* beta *= (1+i), r *= (1+i) */
      ((double *) beta)[0] = beta_i_r;
      ((double *) beta)[1] = beta_i_r;
      ((double *) r)[0] = r_i;
      ((double *) r)[1] = r_i;
    }
    r_true_l[1] = 2.0 * r_true_l[0];
    r_true_t[1] = 2.0 * r_true_t[0];
    r_true_l[0] = r_true_t[0] = 0.0;
  } else {			/* both alpha and beta are free */
    assert(alpha_flag == 0 && beta_flag == 0);
    ((double *) alpha)[0] = alpha_i_r;
    ((double *) alpha)[1] = alpha_i_r;
    ((double *) beta)[0] = beta_i_r;
    ((double *) beta)[1] = beta_i_r;
    ((double *) r)[0] = r_i;
    ((double *) r)[1] = r_i;
    /* imaginary part of r_true */
    ddmuld(r_true_l[0], r_true_t[0], 2.0, &r_true_l[1], &r_true_t[1]);
    /* real part of r_true */
    r_true_l[0] = 0.;
    r_true_t[0] = 0.;
  }

  for (i = 0; i < n; ++i) {
    ((double *) x)[2 * i] = ((double *) x_i)[i];
    ((double *) x)[2 * i + 1] = ((double *) x_i)[i];
    y[i] = y_i[i];
  }
  if (conj == blas_conj) {
    for (i = 0; i < n; ++i)
      ((double *) x)[2 * i + 1] = -((double *) x)[2 * i + 1];
  }
  blas_free(x_i);		/* also y_i */
}

	/* end BLAS_zdot_z_d_testgen */

void BLAS_zdot_d_z_testgen(int n, int n_fix2, int n_mix, int norm,
			   enum blas_conj_type conj,
			   void *alpha, int alpha_flag,
			   void *beta, int beta_flag,
			   double *x, void *y, int *seed, void *r,
			   double *r_true_l, double *r_true_t)

/*
 * Purpose
 * =======
 * 
 * This routine generates the test inputs to BLAS_zdot_d_z{_x}.
 * 
 * Arguments
 * =========
 *  
 * n       (input) int
 *         The length of the vectors X and Y.
 * 
 * n_fix2  (input) int
 *         Number of pairs in the vectors X and Y that are fixed in value,
 *
 * n_mix   (input) int
 *         Number of pairs in the vectors X and Y with X(i) fixed
 *         and Y(i) free in value.
 *
 * norm    (input) int
 *         = -1 : the vectors are scaled with norms near underflow.
 *         = 0  : the vectors have norms of order 1.
 *         = 1  : the vectors are scaled with norms near overflow.
 *
 * conj    (input) enum blas_conj_type
 *
 * alpha   (input/output) void*
 *         If alpha_flag = 1, alpha is input.
 *         If alpha_flag = 0, alpha is output.
 *
 * alpha_flag (input) int
 *         = 0 : alpha is free, and is output.
 *         = 1 : alpha is fixed on input.
 *
 * beta    (input) void*
 *         If beta_flag = 1, beta is input.
 *         If beta_flag = 0, beta is output.
 *
 * beta_flag (input) int
 *         = 0 : beta is free, and is output.
 *         = 1 : beta is fixed on input.
 *
 * x       (input/output) double*
 *
 * y       (input/output) void*
 *
 * seed    (input/output) int*
 *         The seed for the random number generator.
 * 
 * r       (output) void*
 *         The generated scalar r that will be used as an input to DOT.
 *
 * r_true_l (output) double*
 *         The leading part of the truth in double-double.
 *
 * r_true_t (output) double*
 *         The trailing part of the truth in double-double.
 *
 */
{
  int i;
  double alpha_i_r;
  double alpha_i_i;
  double beta_i_r;
  double beta_i_i;
  double r_i;
  double *x_i;
  double *y_i;

  alpha_i_r = ((double *) alpha)[0];
  alpha_i_i = ((double *) alpha)[1];
  beta_i_r = ((double *) beta)[0];
  beta_i_i = ((double *) beta)[1];
  x_i = (double *) blas_malloc(2 * n * sizeof(double));
  if (2 * n > 0 && x_i == NULL) {
    printf("malloc failed\n");
    exit(-1);
  }
  y_i = x_i + n;
  for (i = 0; i < n_fix2; i++) {
    x_i[i] = x[i];
    ((double *) y_i)[i] = ((double *) y)[2 * i];
  }
  for (; i < n_fix2 + n_mix; i++) {
    x_i[i] = x[i];
  }

  /* Call generator now. */
  testgen_BLAS_ddot(n, n_fix2, n_mix, norm, conj,
		    &alpha_i_r, alpha_flag,
		    &beta_i_r, beta_flag,
		    x_i, y_i, seed, &r_i, &r_true_l[0], &r_true_t[0]);

  if (alpha_flag == 1) {	/* alpha is fixed */
    if (alpha_i_r == 1.0 && alpha_i_i == 0.) {	/* alpha == 1.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta == 0 or 1 --> r *= (1+i) */
	((double *) r)[0] = r_i;
	((double *) r)[1] = r_i;
      } else {			/* beta *= (1-i), r *= i --> prod = 1+i */
	((double *) beta)[0] = beta_i_r;
	((double *) beta)[1] = -beta_i_r;
	((double *) r)[0] = 0.0;
	((double *) r)[1] = r_i;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else if (alpha_i_r == 0. && alpha_i_i == 0.) {	/* alpha == 0.0 */
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta == 0 or 1 --> r *= (1+i) */
	((double *) r)[0] = r_i;
	((double *) r)[1] = r_i;
      } else {			/* beta *= (1-i), r *= i --> prod = 1+i */
	((double *) beta)[0] = beta_i_r;
	((double *) beta)[1] = -beta_i_r;
	((double *) r)[0] = 0.0;
	((double *) r)[1] = r_i;
      }
      r_true_l[1] = r_true_l[0];
      r_true_t[1] = r_true_t[0];
    } else {			/* alpha is a fixed multiple of (1+i) */
      ((double *) alpha)[0] = alpha_i_r;
      ((double *) alpha)[1] = alpha_i_r;
      if (beta_flag == 1 &&
	  ((beta_i_r == 0. && beta_i_i == 0.) ||
	   (beta_i_r == 1. && beta_i_i == 0.))) {
	/* beta is 0 or 1 --> r *= 2i --> prod = 2i */
	((double *) r)[0] = 0.0;
	((double *) r)[1] = 2.0 * r_i;
      } else {			/* beta *= (1+i), r *= (1+i) --> prod = 2i */
	((double *) beta)[0] = beta_i_r;
	((double *) beta)[1] = beta_i_r;
	((double *) r)[0] = r_i;
	((double *) r)[1] = r_i;
      }
      r_true_l[1] = 2.0 * r_true_l[0];
      r_true_t[1] = 2.0 * r_true_t[0];
      r_true_l[0] = r_true_t[0] = 0.0;
    }
  } else if (beta_flag == 1) {	/* alpha is free, beta is fixed */
    /* alpha *= (1+i) */
    ((double *) alpha)[0] = alpha_i_r;
    ((double *) alpha)[1] = alpha_i_r;
    if ((beta_i_r == 0. && beta_i_i == 0.) || (beta_i_r == 1. && beta_i_i == 0.)) {	/* r*=2i --> prod = 2i */
      ((double *) r)[0] = 0.0;
      ((double *) r)[1] = 2.0 * r_i;
    } else {			/* beta *= (1+i), r *= (1+i) */
      ((double *) beta)[0] = beta_i_r;
      ((double *) beta)[1] = beta_i_r;
      ((double *) r)[0] = r_i;
      ((double *) r)[1] = r_i;
    }
    r_true_l[1] = 2.0 * r_true_l[0];
    r_true_t[1] = 2.0 * r_true_t[0];
    r_true_l[0] = r_true_t[0] = 0.0;
  } else {			/* both alpha and beta are free */
    assert(alpha_flag == 0 && beta_flag == 0);
    ((double *) alpha)[0] = alpha_i_r;
    ((double *) alpha)[1] = alpha_i_r;
    ((double *) beta)[0] = beta_i_r;
    ((double *) beta)[1] = beta_i_r;
    ((double *) r)[0] = r_i;
    ((double *) r)[1] = r_i;
    /* imaginary part of r_true */
    ddmuld(r_true_l[0], r_true_t[0], 2.0, &r_true_l[1], &r_true_t[1]);
    /* real part of r_true */
    r_true_l[0] = 0.;
    r_true_t[0] = 0.;
  }

  for (i = 0; i < n; ++i) {
    x[i] = x_i[i];
    ((double *) y)[2 * i] = ((double *) y_i)[i];
    ((double *) y)[2 * i + 1] = ((double *) y_i)[i];
  }

  blas_free(x_i);		/* also y_i */
}

	/* end BLAS_zdot_d_z_testgen */
