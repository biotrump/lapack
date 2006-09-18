#include "blas_extended.h"



void strsv_copy(enum blas_order_type order,
		enum blas_uplo_type uplo,
		enum blas_trans_type trans,
		int n, const float *T, int lda, float *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from T to y
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              Whether T is no trans, trans, or conj trans 
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * T            (input) float*
 *              The triangular matrix T
 *
 * lda          (input) int
 *              Leading dimension
 *
 * y            (output) float*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  int i, inc = 1;
  float tmp;
  const float *T_i = T;
  float *y_i = y;


  if (((order == blas_rowmajor) && (trans == blas_no_trans)) ||
      ((order == blas_colmajor) && (trans != blas_no_trans))) {
    for (i = 0; i < n; i++) {
      tmp = T_i[(row * lda + i) * inc];
      y_i[i * inc] = tmp;
    }
  } else {
    for (i = 0; i < n; i++) {
      tmp = T_i[(row + i * lda) * inc];
      y_i[i * inc] = tmp;
    }
  }

}

void strsv_commit(enum blas_order_type order,
		  enum blas_uplo_type uplo,
		  enum blas_trans_type trans,
		  int length, float *T, int lda, const float *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy y to T
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              Whether T is no trans, trans, or conj trans 
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * T            (input) float*
 *              The triangular matrix T
 *
 * lda          (input) int
 *              Leading dimension
 *
 * y            (output) float*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  int i;
  float tmp;

  if (order == blas_rowmajor && uplo == blas_upper && trans == blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[(row * lda) + (row + i + 1)] = tmp;
    }
  } else
    if (order == blas_rowmajor && uplo == blas_lower
	&& trans == blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[row * lda + i] = tmp;
    }
  } else
    if (order == blas_rowmajor && uplo == blas_lower
	&& trans != blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[(row + i + 1) * lda + row] = tmp;
    }
  } else
    if (order == blas_rowmajor && uplo == blas_upper
	&& trans != blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[i * lda + row] = tmp;
    }
  } else
    if (order == blas_colmajor && uplo == blas_upper
	&& trans == blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[(row + i + 1) * lda + row] = tmp;
    }
  } else
    if (order == blas_colmajor && uplo == blas_lower
	&& trans == blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[i * lda + row] = tmp;
    }
  } else
    if (order == blas_colmajor && uplo == blas_upper
	&& trans != blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[row * lda + i] = tmp;
    }
  } else
    if (order == blas_colmajor && uplo == blas_lower
	&& trans != blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[row * lda + (row + i + 1)] = tmp;
    }
  }
}

void dtrsv_copy(enum blas_order_type order,
		enum blas_uplo_type uplo,
		enum blas_trans_type trans,
		int n, const double *T, int lda, double *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from T to y
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              Whether T is no trans, trans, or conj trans 
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * T            (input) double*
 *              The triangular matrix T
 *
 * lda          (input) int
 *              Leading dimension
 *
 * y            (output) double*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  int i, inc = 1;
  double tmp;
  const double *T_i = T;
  double *y_i = y;


  if (((order == blas_rowmajor) && (trans == blas_no_trans)) ||
      ((order == blas_colmajor) && (trans != blas_no_trans))) {
    for (i = 0; i < n; i++) {
      tmp = T_i[(row * lda + i) * inc];
      y_i[i * inc] = tmp;
    }
  } else {
    for (i = 0; i < n; i++) {
      tmp = T_i[(row + i * lda) * inc];
      y_i[i * inc] = tmp;
    }
  }

}

void dtrsv_commit(enum blas_order_type order,
		  enum blas_uplo_type uplo,
		  enum blas_trans_type trans,
		  int length, double *T, int lda, const double *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy y to T
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              Whether T is no trans, trans, or conj trans 
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * T            (input) double*
 *              The triangular matrix T
 *
 * lda          (input) int
 *              Leading dimension
 *
 * y            (output) double*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  int i;
  double tmp;

  if (order == blas_rowmajor && uplo == blas_upper && trans == blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[(row * lda) + (row + i + 1)] = tmp;
    }
  } else
    if (order == blas_rowmajor && uplo == blas_lower
	&& trans == blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[row * lda + i] = tmp;
    }
  } else
    if (order == blas_rowmajor && uplo == blas_lower
	&& trans != blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[(row + i + 1) * lda + row] = tmp;
    }
  } else
    if (order == blas_rowmajor && uplo == blas_upper
	&& trans != blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[i * lda + row] = tmp;
    }
  } else
    if (order == blas_colmajor && uplo == blas_upper
	&& trans == blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[(row + i + 1) * lda + row] = tmp;
    }
  } else
    if (order == blas_colmajor && uplo == blas_lower
	&& trans == blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[i * lda + row] = tmp;
    }
  } else
    if (order == blas_colmajor && uplo == blas_upper
	&& trans != blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[row * lda + i] = tmp;
    }
  } else
    if (order == blas_colmajor && uplo == blas_lower
	&& trans != blas_no_trans) {
    for (i = 0; i < length; i++) {
      tmp = y[i];
      T[row * lda + (row + i + 1)] = tmp;
    }
  }
}

void ctrsv_copy(enum blas_order_type order,
		enum blas_uplo_type uplo,
		enum blas_trans_type trans,
		int n, const void *T, int lda, void *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from T to y
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              Whether T is no trans, trans, or conj trans 
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * T            (input) void*
 *              The triangular matrix T
 *
 * lda          (input) int
 *              Leading dimension
 *
 * y            (output) void*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  int i, inc = 1;
  float tmp[2];
  const float *T_i = (float *) T;
  float *y_i = (float *) y;
  inc *= 2;

  if (((order == blas_rowmajor) && (trans == blas_no_trans)) ||
      ((order == blas_colmajor) && (trans != blas_no_trans))) {
    for (i = 0; i < n; i++) {
      tmp[0] = T_i[(row * lda + i) * inc];
      tmp[1] = T_i[(row * lda + i) * inc + 1];
      y_i[i * inc] = tmp[0];
      y_i[i * inc + 1] = tmp[1];
    }
  } else {
    for (i = 0; i < n; i++) {
      tmp[0] = T_i[(row + i * lda) * inc];
      tmp[1] = T_i[(row + i * lda) * inc + 1];
      y_i[i * inc] = tmp[0];
      y_i[i * inc + 1] = tmp[1];
    }
  }

}

void ztrsv_copy(enum blas_order_type order,
		enum blas_uplo_type uplo,
		enum blas_trans_type trans,
		int n, const void *T, int lda, void *y, int row)
/*
 * Purpose
 * =======
 *
 * Copy a row from T to y
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of T; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether T is upper or lower
 *
 * trans        (input) blas_trans_type
 *              Whether T is no trans, trans, or conj trans 
 *
 * n            (input) int
 *              Dimension of AP and the length of vector x
 *
 * T            (input) void*
 *              The triangular matrix T
 *
 * lda          (input) int
 *              Leading dimension
 *
 * y            (output) void*
 *              The vector y
 *
 * row          (input) int
 *              The row to be copyied to y
 */
{
  int i, inc = 1;
  double tmp[2];
  const double *T_i = (double *) T;
  double *y_i = (double *) y;
  inc *= 2;

  if (((order == blas_rowmajor) && (trans == blas_no_trans)) ||
      ((order == blas_colmajor) && (trans != blas_no_trans))) {
    for (i = 0; i < n; i++) {
      tmp[0] = T_i[(row * lda + i) * inc];
      tmp[1] = T_i[(row * lda + i) * inc + 1];
      y_i[i * inc] = tmp[0];
      y_i[i * inc + 1] = tmp[1];
    }
  } else {
    for (i = 0; i < n; i++) {
      tmp[0] = T_i[(row + i * lda) * inc];
      tmp[1] = T_i[(row + i * lda) * inc + 1];
      y_i[i * inc] = tmp[0];
      y_i[i * inc + 1] = tmp[1];
    }
  }

}
