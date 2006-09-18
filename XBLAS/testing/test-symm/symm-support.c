#include <stdio.h>
#include "blas_extended.h"




















void ssymm_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_side_type side, int n, float *a, int lda,
		      float *a_vec, int row)

/*
 *  Copies the given vector into the given row of matrix a.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;

  float a_elem;
  float *a_i = a;
  const float *a_vec_i = a_vec;


  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }






  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    a_elem = a_vec_i[vi];
    a_i[ai] = a_elem;
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    a_elem = a_vec_i[vi];
    a_i[ai] = a_elem;
  }
}
void dsymm_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_side_type side, int n, double *a, int lda,
		      double *a_vec, int row)

/*
 *  Copies the given vector into the given row of matrix a.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;

  double a_elem;
  double *a_i = a;
  const double *a_vec_i = a_vec;


  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }






  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    a_elem = a_vec_i[vi];
    a_i[ai] = a_elem;
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    a_elem = a_vec_i[vi];
    a_i[ai] = a_elem;
  }
}
void csymm_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_side_type side, int n, void *a, int lda,
		      void *a_vec, int row)

/*
 *  Copies the given vector into the given row of matrix a.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;

  float a_elem[2];
  float *a_i = (float *) a;
  const float *a_vec_i = (float *) a_vec;


  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }

  incai1 *= 2;
  incai2 *= 2;
  ai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    a_elem[0] = a_vec_i[vi];
    a_elem[1] = a_vec_i[vi + 1];
    a_i[ai] = a_elem[0];
    a_i[ai + 1] = a_elem[1];
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    a_elem[0] = a_vec_i[vi];
    a_elem[1] = a_vec_i[vi + 1];
    a_i[ai] = a_elem[0];
    a_i[ai + 1] = a_elem[1];
  }
}
void zsymm_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_side_type side, int n, void *a, int lda,
		      void *a_vec, int row)

/*
 *  Copies the given vector into the given row of matrix a.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;

  double a_elem[2];
  double *a_i = (double *) a;
  const double *a_vec_i = (double *) a_vec;


  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }

  incai1 *= 2;
  incai2 *= 2;
  ai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    a_elem[0] = a_vec_i[vi];
    a_elem[1] = a_vec_i[vi + 1];
    a_i[ai] = a_elem[0];
    a_i[ai + 1] = a_elem[1];
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    a_elem[0] = a_vec_i[vi];
    a_elem[1] = a_vec_i[vi + 1];
    a_i[ai] = a_elem[0];
    a_i[ai + 1] = a_elem[1];
  }
}

void ssymm_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_side_type side, int n, float *a, int lda,
		    float *a_vec, int row)

/*
 *  Copies the given row of matrix a into the supplied vector.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;

  float a_elem;
  const float *a_i = a;
  float *a_vec_i = a_vec;


  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }






  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    a_elem = a_i[ai];
    a_vec_i[vi] = a_elem;
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    a_elem = a_i[ai];
    a_vec_i[vi] = a_elem;
  }
}
void dsymm_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_side_type side, int n, double *a, int lda,
		    double *a_vec, int row)

/*
 *  Copies the given row of matrix a into the supplied vector.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;

  double a_elem;
  const double *a_i = a;
  double *a_vec_i = a_vec;


  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }






  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    a_elem = a_i[ai];
    a_vec_i[vi] = a_elem;
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    a_elem = a_i[ai];
    a_vec_i[vi] = a_elem;
  }
}
void csymm_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_side_type side, int n, void *a, int lda,
		    void *a_vec, int row)

/*
 *  Copies the given row of matrix a into the supplied vector.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;

  float a_elem[2];
  const float *a_i = (float *) a;
  float *a_vec_i = (float *) a_vec;


  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }

  incai1 *= 2;
  incai2 *= 2;
  ai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }
}
void zsymm_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_side_type side, int n, void *a, int lda,
		    void *a_vec, int row)

/*
 *  Copies the given row of matrix a into the supplied vector.
 */
{

  int ai;
  int i;
  int incai1;
  int incai2;
  int vi;
  int incvi = 1;

  double a_elem[2];
  const double *a_i = (double *) a;
  double *a_vec_i = (double *) a_vec;


  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }

  incai1 *= 2;
  incai2 *= 2;
  ai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }
}

void sprint_symm_matrix(float *a, int n, int lda,
			enum blas_order_type order, enum blas_uplo_type uplo)
{
  int ai, aij;
  int incai, incaij1, incaij2;
  int i, j;

  float a_elem;
  const float *a_i = a;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai = lda;
    incaij1 = 1;
    incaij2 = lda;
  } else {
    incai = 1;
    incaij1 = lda;
    incaij2 = 1;
  }





  for (i = 0, ai = 0; i < n; i++, ai += incai) {

    for (j = 0, aij = ai; j < i; j++, aij += incaij1) {
      a_elem = a_i[aij];
      printf(" %.12g ", a_elem);

    }
    for (; j < n; j++, aij += incaij2) {
      a_elem = a_i[aij];
      printf(" %.12g ", a_elem);

    }
    printf("\n");
  }

}
void dprint_symm_matrix(double *a, int n, int lda,
			enum blas_order_type order, enum blas_uplo_type uplo)
{
  int ai, aij;
  int incai, incaij1, incaij2;
  int i, j;

  double a_elem;
  const double *a_i = a;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai = lda;
    incaij1 = 1;
    incaij2 = lda;
  } else {
    incai = 1;
    incaij1 = lda;
    incaij2 = 1;
  }





  for (i = 0, ai = 0; i < n; i++, ai += incai) {

    for (j = 0, aij = ai; j < i; j++, aij += incaij1) {
      a_elem = a_i[aij];
      printf(" %.24g ", a_elem);

    }
    for (; j < n; j++, aij += incaij2) {
      a_elem = a_i[aij];
      printf(" %.24g ", a_elem);

    }
    printf("\n");
  }

}
void cprint_symm_matrix(void *a, int n, int lda,
			enum blas_order_type order, enum blas_uplo_type uplo)
{
  int ai, aij;
  int incai, incaij1, incaij2;
  int i, j;

  float a_elem[2];
  const float *a_i = (float *) a;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai = lda;
    incaij1 = 1;
    incaij2 = lda;
  } else {
    incai = 1;
    incaij1 = lda;
    incaij2 = 1;
  }

  incai *= 2;
  incaij1 *= 2;
  incaij2 *= 2;

  for (i = 0, ai = 0; i < n; i++, ai += incai) {

    for (j = 0, aij = ai; j < i; j++, aij += incaij1) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      printf(" %.12g + %.12gi", ((float *) a_elem)[0], ((float *) a_elem)[1]);

    }
    for (; j < n; j++, aij += incaij2) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      printf(" %.12g + %.12gi", ((float *) a_elem)[0], ((float *) a_elem)[1]);

    }
    printf("\n");
  }

}
void zprint_symm_matrix(void *a, int n, int lda,
			enum blas_order_type order, enum blas_uplo_type uplo)
{
  int ai, aij;
  int incai, incaij1, incaij2;
  int i, j;

  double a_elem[2];
  const double *a_i = (double *) a;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai = lda;
    incaij1 = 1;
    incaij2 = lda;
  } else {
    incai = 1;
    incaij1 = lda;
    incaij2 = 1;
  }

  incai *= 2;
  incaij1 *= 2;
  incaij2 *= 2;

  for (i = 0, ai = 0; i < n; i++, ai += incai) {

    for (j = 0, aij = ai; j < i; j++, aij += incaij1) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      printf(" %.24g + %.24gi ", ((double *) a_elem)[0],
	     ((double *) a_elem)[1]);

    }
    for (; j < n; j++, aij += incaij2) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      printf(" %.24g + %.24gi ", ((double *) a_elem)[0],
	     ((double *) a_elem)[1]);

    }
    printf("\n");
  }

}