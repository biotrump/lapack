#include <stdio.h>
#include "blas_extended.h"


































void ssymv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, float *a, int lda, float *a_vec, int row)

/*
 *  Copies the given vector into the given row of symmetric matrix a.
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
void dsymv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, double *a, int lda, double *a_vec, int row)

/*
 *  Copies the given vector into the given row of symmetric matrix a.
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
void csymv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, int lda, void *a_vec, int row)

/*
 *  Copies the given vector into the given row of symmetric matrix a.
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
void zsymv_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      int n, void *a, int lda, void *a_vec, int row)

/*
 *  Copies the given vector into the given row of symmetric matrix a.
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

void ssymv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, float *a, int lda, float *a_vec, int row)

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
void dsymv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, double *a, int lda, double *a_vec, int row)

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
void csymv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, void *a, int lda, void *a_vec, int row)

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
void zsymv_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, void *a, int lda, void *a_vec, int row)

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

void ssymv_copy_vector(int n, float *x_to, int incx_to,
		       float *x_from, int incx_from)

/*
 *  Copies the given vector x_from into the supplied vector x_to.
 */
{

  int x_toi, x_fromi;
  int incx_toi, incx_fromi;
  int i, n_i;
  int x_starti_to, x_starti_from;

  float x_elem;
  const float *x_from_i = x_from;
  float *x_to_i = x_to;

  incx_toi = incx_to;
  incx_fromi = incx_from;



  if (incx_to < 0) {
    x_starti_to = (-n + 1) * incx_toi;
  } else {
    x_starti_to = 0;
  }
  if (incx_from < 0) {
    x_starti_from = (-n + 1) * incx_fromi;
  } else {
    x_starti_from = 0;
  }


  n_i = n;

  for (i = 0, x_toi = x_starti_to, x_fromi = x_starti_from;
       i < n_i; i++, x_toi += incx_toi, x_fromi += incx_fromi) {
    x_elem = x_from_i[x_fromi];
    x_to_i[x_toi] = x_elem;
  }

}
void dsymv_copy_vector(int n, double *x_to, int incx_to,
		       double *x_from, int incx_from)

/*
 *  Copies the given vector x_from into the supplied vector x_to.
 */
{

  int x_toi, x_fromi;
  int incx_toi, incx_fromi;
  int i, n_i;
  int x_starti_to, x_starti_from;

  double x_elem;
  const double *x_from_i = x_from;
  double *x_to_i = x_to;

  incx_toi = incx_to;
  incx_fromi = incx_from;



  if (incx_to < 0) {
    x_starti_to = (-n + 1) * incx_toi;
  } else {
    x_starti_to = 0;
  }
  if (incx_from < 0) {
    x_starti_from = (-n + 1) * incx_fromi;
  } else {
    x_starti_from = 0;
  }


  n_i = n;

  for (i = 0, x_toi = x_starti_to, x_fromi = x_starti_from;
       i < n_i; i++, x_toi += incx_toi, x_fromi += incx_fromi) {
    x_elem = x_from_i[x_fromi];
    x_to_i[x_toi] = x_elem;
  }

}
void csymv_copy_vector(int n, void *x_to, int incx_to,
		       void *x_from, int incx_from)

/*
 *  Copies the given vector x_from into the supplied vector x_to.
 */
{

  int x_toi, x_fromi;
  int incx_toi, incx_fromi;
  int i, n_i;
  int x_starti_to, x_starti_from;

  float x_elem[2];
  const float *x_from_i = (float *) x_from;
  float *x_to_i = (float *) x_to;

  incx_toi = incx_to;
  incx_fromi = incx_from;

  incx_toi *= 2;
  incx_fromi *= 2;
  if (incx_to < 0) {
    x_starti_to = (-n + 1) * incx_toi;
  } else {
    x_starti_to = 0;
  }
  if (incx_from < 0) {
    x_starti_from = (-n + 1) * incx_fromi;
  } else {
    x_starti_from = 0;
  }


  n_i = n;

  for (i = 0, x_toi = x_starti_to, x_fromi = x_starti_from;
       i < n_i; i++, x_toi += incx_toi, x_fromi += incx_fromi) {
    x_elem[0] = x_from_i[x_fromi];
    x_elem[1] = x_from_i[x_fromi + 1];
    x_to_i[x_toi] = x_elem[0];
    x_to_i[x_toi + 1] = x_elem[1];
  }

}
void zsymv_copy_vector(int n, void *x_to, int incx_to,
		       void *x_from, int incx_from)

/*
 *  Copies the given vector x_from into the supplied vector x_to.
 */
{

  int x_toi, x_fromi;
  int incx_toi, incx_fromi;
  int i, n_i;
  int x_starti_to, x_starti_from;

  double x_elem[2];
  const double *x_from_i = (double *) x_from;
  double *x_to_i = (double *) x_to;

  incx_toi = incx_to;
  incx_fromi = incx_from;

  incx_toi *= 2;
  incx_fromi *= 2;
  if (incx_to < 0) {
    x_starti_to = (-n + 1) * incx_toi;
  } else {
    x_starti_to = 0;
  }
  if (incx_from < 0) {
    x_starti_from = (-n + 1) * incx_fromi;
  } else {
    x_starti_from = 0;
  }


  n_i = n;

  for (i = 0, x_toi = x_starti_to, x_fromi = x_starti_from;
       i < n_i; i++, x_toi += incx_toi, x_fromi += incx_fromi) {
    x_elem[0] = x_from_i[x_fromi];
    x_elem[1] = x_from_i[x_fromi + 1];
    x_to_i[x_toi] = x_elem[0];
    x_to_i[x_toi + 1] = x_elem[1];
  }

}

void sprint_symv_matrix(float *a, int n, int lda,
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
void dprint_symv_matrix(double *a, int n, int lda,
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
void cprint_symv_matrix(void *a, int n, int lda,
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
void zprint_symv_matrix(void *a, int n, int lda,
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

void sprint_vector(float *x, int n, int incx)
{
  int xi;
  int incxi;
  int x_starti;
  int i;

  float x_elem;
  const float *x_i = x;

  incxi = incx;

  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  printf("\n");

  for (i = 0, xi = x_starti; i < n; i++, xi += incxi) {
    x_elem = x_i[xi];
    printf(" %.12g ", x_elem);

  }
  printf("\n");


}
void dprint_vector(double *x, int n, int incx)
{
  int xi;
  int incxi;
  int x_starti;
  int i;

  double x_elem;
  const double *x_i = x;

  incxi = incx;

  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  printf("\n");

  for (i = 0, xi = x_starti; i < n; i++, xi += incxi) {
    x_elem = x_i[xi];
    printf(" %.24g ", x_elem);

  }
  printf("\n");


}
void cprint_vector(void *x, int n, int incx)
{
  int xi;
  int incxi;
  int x_starti;
  int i;

  float x_elem[2];
  const float *x_i = (float *) x;

  incxi = incx;
  incxi *= 2;
  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  printf("\n");

  for (i = 0, xi = x_starti; i < n; i++, xi += incxi) {
    x_elem[0] = x_i[xi];
    x_elem[1] = x_i[xi + 1];
    printf(" %.12g + %.12gi", ((float *) x_elem)[0], ((float *) x_elem)[1]);

  }
  printf("\n");


}
void zprint_vector(void *x, int n, int incx)
{
  int xi;
  int incxi;
  int x_starti;
  int i;

  double x_elem[2];
  const double *x_i = (double *) x;

  incxi = incx;
  incxi *= 2;
  if (incxi < 0) {
    x_starti = (-n + 1) * incxi;
  } else {
    x_starti = 0;
  }

  printf("\n");

  for (i = 0, xi = x_starti; i < n; i++, xi += incxi) {
    x_elem[0] = x_i[xi];
    x_elem[1] = x_i[xi + 1];
    printf(" %.24g + %.24gi ", ((double *) x_elem)[0],
	   ((double *) x_elem)[1]);

  }
  printf("\n");


}
