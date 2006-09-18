#include <stdio.h>
#include "blas_extended.h"




































void sgemm_commit_row(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, float *a, int lda, float *a_vec, int row)
{
  int ai;
  int i;
  int incai;
  int incvi = 1;
  int vi;

  float a_elem;
  float *a_i = a;
  const float *a_vec_i = a_vec;

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      incai = lda;
      ai = row;
    } else {
      incai = 1;
      ai = row * lda;
    }
  } else {
    if (trans == blas_no_trans) {
      incai = 1;
      ai = row * lda;
    } else {
      incai = lda;
      ai = row;
    }
  }





  for (i = 0, vi = 0; i < n; i++, vi += incvi, ai += incai) {
    a_elem = a_vec_i[vi];
    if (trans == blas_conj_trans) {

    }
    a_i[ai] = a_elem;
  }

}
void dgemm_commit_row(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, double *a, int lda, double *a_vec,
		      int row)
{
  int ai;
  int i;
  int incai;
  int incvi = 1;
  int vi;

  double a_elem;
  double *a_i = a;
  const double *a_vec_i = a_vec;

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      incai = lda;
      ai = row;
    } else {
      incai = 1;
      ai = row * lda;
    }
  } else {
    if (trans == blas_no_trans) {
      incai = 1;
      ai = row * lda;
    } else {
      incai = lda;
      ai = row;
    }
  }





  for (i = 0, vi = 0; i < n; i++, vi += incvi, ai += incai) {
    a_elem = a_vec_i[vi];
    if (trans == blas_conj_trans) {

    }
    a_i[ai] = a_elem;
  }

}
void cgemm_commit_row(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, void *a, int lda, void *a_vec, int row)
{
  int ai;
  int i;
  int incai;
  int incvi = 1;
  int vi;

  float a_elem[2];
  float *a_i = (float *) a;
  const float *a_vec_i = (float *) a_vec;

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      incai = lda;
      ai = row;
    } else {
      incai = 1;
      ai = row * lda;
    }
  } else {
    if (trans == blas_no_trans) {
      incai = 1;
      ai = row * lda;
    } else {
      incai = lda;
      ai = row;
    }
  }

  ai *= 2;
  incai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < n; i++, vi += incvi, ai += incai) {
    a_elem[0] = a_vec_i[vi];
    a_elem[1] = a_vec_i[vi + 1];
    if (trans == blas_conj_trans) {
      a_elem[1] = -a_elem[1];
    }
    a_i[ai] = a_elem[0];
    a_i[ai + 1] = a_elem[1];
  }

}
void zgemm_commit_row(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, void *a, int lda, void *a_vec, int row)
{
  int ai;
  int i;
  int incai;
  int incvi = 1;
  int vi;

  double a_elem[2];
  double *a_i = (double *) a;
  const double *a_vec_i = (double *) a_vec;

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      incai = lda;
      ai = row;
    } else {
      incai = 1;
      ai = row * lda;
    }
  } else {
    if (trans == blas_no_trans) {
      incai = 1;
      ai = row * lda;
    } else {
      incai = lda;
      ai = row;
    }
  }

  ai *= 2;
  incai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < n; i++, vi += incvi, ai += incai) {
    a_elem[0] = a_vec_i[vi];
    a_elem[1] = a_vec_i[vi + 1];
    if (trans == blas_conj_trans) {
      a_elem[1] = -a_elem[1];
    }
    a_i[ai] = a_elem[0];
    a_i[ai + 1] = a_elem[1];
  }

}

void sgemm_commit_col(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, float *a, int lda, float *a_vec, int col)
{
  int ai;
  int i;
  int incai;
  int incvi = 1;
  int vi;

  float a_elem;
  float *a_i = a;
  const float *a_vec_i = a_vec;

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      incai = 1;
      ai = col * lda;
    } else {
      incai = lda;
      ai = col;
    }
  } else {
    if (trans == blas_no_trans) {
      incai = lda;
      ai = col;
    } else {
      incai = 1;
      ai = col * lda;
    }
  }





  for (i = 0, vi = 0; i < m; i++, vi += incvi, ai += incai) {
    a_elem = a_vec_i[vi];
    if (trans == blas_conj_trans) {

    }
    a_i[ai] = a_elem;
  }

}
void dgemm_commit_col(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, double *a, int lda, double *a_vec,
		      int col)
{
  int ai;
  int i;
  int incai;
  int incvi = 1;
  int vi;

  double a_elem;
  double *a_i = a;
  const double *a_vec_i = a_vec;

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      incai = 1;
      ai = col * lda;
    } else {
      incai = lda;
      ai = col;
    }
  } else {
    if (trans == blas_no_trans) {
      incai = lda;
      ai = col;
    } else {
      incai = 1;
      ai = col * lda;
    }
  }





  for (i = 0, vi = 0; i < m; i++, vi += incvi, ai += incai) {
    a_elem = a_vec_i[vi];
    if (trans == blas_conj_trans) {

    }
    a_i[ai] = a_elem;
  }

}
void cgemm_commit_col(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, void *a, int lda, void *a_vec, int col)
{
  int ai;
  int i;
  int incai;
  int incvi = 1;
  int vi;

  float a_elem[2];
  float *a_i = (float *) a;
  const float *a_vec_i = (float *) a_vec;

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      incai = 1;
      ai = col * lda;
    } else {
      incai = lda;
      ai = col;
    }
  } else {
    if (trans == blas_no_trans) {
      incai = lda;
      ai = col;
    } else {
      incai = 1;
      ai = col * lda;
    }
  }

  ai *= 2;
  incai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < m; i++, vi += incvi, ai += incai) {
    a_elem[0] = a_vec_i[vi];
    a_elem[1] = a_vec_i[vi + 1];
    if (trans == blas_conj_trans) {
      a_elem[1] = -a_elem[1];
    }
    a_i[ai] = a_elem[0];
    a_i[ai + 1] = a_elem[1];
  }

}
void zgemm_commit_col(enum blas_order_type order, enum blas_trans_type trans,
		      int m, int n, void *a, int lda, void *a_vec, int col)
{
  int ai;
  int i;
  int incai;
  int incvi = 1;
  int vi;

  double a_elem[2];
  double *a_i = (double *) a;
  const double *a_vec_i = (double *) a_vec;

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      incai = 1;
      ai = col * lda;
    } else {
      incai = lda;
      ai = col;
    }
  } else {
    if (trans == blas_no_trans) {
      incai = lda;
      ai = col;
    } else {
      incai = 1;
      ai = col * lda;
    }
  }

  ai *= 2;
  incai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < m; i++, vi += incvi, ai += incai) {
    a_elem[0] = a_vec_i[vi];
    a_elem[1] = a_vec_i[vi + 1];
    if (trans == blas_conj_trans) {
      a_elem[1] = -a_elem[1];
    }
    a_i[ai] = a_elem[0];
    a_i[ai + 1] = a_elem[1];
  }

}

void sgemm_copy_row
  (enum blas_order_type order, enum blas_trans_type trans,
   int m, int n, float *a, int lda, float *a_vec, int row) {

  int ai;
  int incai;
  int i;
  int incvi = 1;
  int vi;

  float a_elem;
  const float *a_i = a;
  float *a_vec_i = a_vec;

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      ai = row;
      incai = lda;
    } else {
      ai = row * lda;
      incai = 1;
    }
  } else {
    if (trans == blas_no_trans) {
      ai = row * lda;
      incai = 1;
    } else {
      ai = row;
      incai = lda;
    }

  }





  for (i = 0, vi = 0; i < n; i++, vi += incvi, ai += incai) {
    a_elem = a_i[ai];
    if (trans == blas_conj_trans) {

    }
    a_vec_i[vi] = a_elem;
  }


}
void dgemm_copy_row
  (enum blas_order_type order, enum blas_trans_type trans,
   int m, int n, double *a, int lda, double *a_vec, int row) {

  int ai;
  int incai;
  int i;
  int incvi = 1;
  int vi;

  double a_elem;
  const double *a_i = a;
  double *a_vec_i = a_vec;

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      ai = row;
      incai = lda;
    } else {
      ai = row * lda;
      incai = 1;
    }
  } else {
    if (trans == blas_no_trans) {
      ai = row * lda;
      incai = 1;
    } else {
      ai = row;
      incai = lda;
    }

  }





  for (i = 0, vi = 0; i < n; i++, vi += incvi, ai += incai) {
    a_elem = a_i[ai];
    if (trans == blas_conj_trans) {

    }
    a_vec_i[vi] = a_elem;
  }


}
void cgemm_copy_row
  (enum blas_order_type order, enum blas_trans_type trans,
   int m, int n, void *a, int lda, void *a_vec, int row) {

  int ai;
  int incai;
  int i;
  int incvi = 1;
  int vi;

  float a_elem[2];
  const float *a_i = (float *) a;
  float *a_vec_i = (float *) a_vec;

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      ai = row;
      incai = lda;
    } else {
      ai = row * lda;
      incai = 1;
    }
  } else {
    if (trans == blas_no_trans) {
      ai = row * lda;
      incai = 1;
    } else {
      ai = row;
      incai = lda;
    }

  }

  ai *= 2;
  incai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < n; i++, vi += incvi, ai += incai) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    if (trans == blas_conj_trans) {
      a_elem[1] = -a_elem[1];
    }
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }


}
void zgemm_copy_row
  (enum blas_order_type order, enum blas_trans_type trans,
   int m, int n, void *a, int lda, void *a_vec, int row) {

  int ai;
  int incai;
  int i;
  int incvi = 1;
  int vi;

  double a_elem[2];
  const double *a_i = (double *) a;
  double *a_vec_i = (double *) a_vec;

  if (order == blas_colmajor) {
    if (trans == blas_no_trans) {
      ai = row;
      incai = lda;
    } else {
      ai = row * lda;
      incai = 1;
    }
  } else {
    if (trans == blas_no_trans) {
      ai = row * lda;
      incai = 1;
    } else {
      ai = row;
      incai = lda;
    }

  }

  ai *= 2;
  incai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < n; i++, vi += incvi, ai += incai) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    if (trans == blas_conj_trans) {
      a_elem[1] = -a_elem[1];
    }
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }


}

void sgemm_copy_col
  (enum blas_order_type order, enum blas_trans_type trans,
   int m, int n, float *a, int lda, float *a_vec, int col) {

  int ai;
  int incai;
  int i;
  int incvi = 1;
  int vi;

  float a_elem;
  const float *a_i = a;
  float *a_vec_i = a_vec;

  if (order == blas_rowmajor) {
    if (trans == blas_no_trans) {
      ai = col;
      incai = lda;
    } else {
      ai = col * lda;
      incai = 1;
    }
  } else {
    if (trans == blas_no_trans) {
      ai = col * lda;
      incai = 1;
    } else {
      ai = col;
      incai = lda;
    }
  }





  for (i = 0, vi = 0; i < m; i++, vi += incvi, ai += incai) {
    a_elem = a_i[ai];
    if (trans == blas_conj_trans) {

    }
    a_vec_i[vi] = a_elem;
  }


}
void dgemm_copy_col
  (enum blas_order_type order, enum blas_trans_type trans,
   int m, int n, double *a, int lda, double *a_vec, int col) {

  int ai;
  int incai;
  int i;
  int incvi = 1;
  int vi;

  double a_elem;
  const double *a_i = a;
  double *a_vec_i = a_vec;

  if (order == blas_rowmajor) {
    if (trans == blas_no_trans) {
      ai = col;
      incai = lda;
    } else {
      ai = col * lda;
      incai = 1;
    }
  } else {
    if (trans == blas_no_trans) {
      ai = col * lda;
      incai = 1;
    } else {
      ai = col;
      incai = lda;
    }
  }





  for (i = 0, vi = 0; i < m; i++, vi += incvi, ai += incai) {
    a_elem = a_i[ai];
    if (trans == blas_conj_trans) {

    }
    a_vec_i[vi] = a_elem;
  }


}
void cgemm_copy_col
  (enum blas_order_type order, enum blas_trans_type trans,
   int m, int n, void *a, int lda, void *a_vec, int col) {

  int ai;
  int incai;
  int i;
  int incvi = 1;
  int vi;

  float a_elem[2];
  const float *a_i = (float *) a;
  float *a_vec_i = (float *) a_vec;

  if (order == blas_rowmajor) {
    if (trans == blas_no_trans) {
      ai = col;
      incai = lda;
    } else {
      ai = col * lda;
      incai = 1;
    }
  } else {
    if (trans == blas_no_trans) {
      ai = col * lda;
      incai = 1;
    } else {
      ai = col;
      incai = lda;
    }
  }

  ai *= 2;
  incai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < m; i++, vi += incvi, ai += incai) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    if (trans == blas_conj_trans) {
      a_elem[1] = -a_elem[1];
    }
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }


}
void zgemm_copy_col
  (enum blas_order_type order, enum blas_trans_type trans,
   int m, int n, void *a, int lda, void *a_vec, int col) {

  int ai;
  int incai;
  int i;
  int incvi = 1;
  int vi;

  double a_elem[2];
  const double *a_i = (double *) a;
  double *a_vec_i = (double *) a_vec;

  if (order == blas_rowmajor) {
    if (trans == blas_no_trans) {
      ai = col;
      incai = lda;
    } else {
      ai = col * lda;
      incai = 1;
    }
  } else {
    if (trans == blas_no_trans) {
      ai = col * lda;
      incai = 1;
    } else {
      ai = col;
      incai = lda;
    }
  }

  ai *= 2;
  incai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < m; i++, vi += incvi, ai += incai) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    if (trans == blas_conj_trans) {
      a_elem[1] = -a_elem[1];
    }
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }


}

void sgemm_copy(enum blas_order_type order, int m, int n,
		float *a, int lda, float *b, int ldb)
{

  int ai, aij;
  int bi, bij;
  int incai, incaij;
  int incbi, incbij;

  int i, j;

  float elem;
  float *a_i = a;
  const float *b_i = b;

  if (order == blas_colmajor) {
    incai = 1;
    incbi = 1;
    incaij = lda;
    incbij = ldb;
  } else {
    incai = lda;
    incbi = ldb;
    incaij = 1;
    incbij = 1;
  }






  for (i = 0, ai = 0, bi = 0; i < m; i++, ai += incai, bi += incbi) {
    for (j = 0, aij = ai, bij = bi; j < n; j++, aij += incaij, bij += incbij) {
      elem = b_i[bij];
      a_i[aij] = elem;
    }
  }


}
void dgemm_copy(enum blas_order_type order, int m, int n,
		double *a, int lda, double *b, int ldb)
{

  int ai, aij;
  int bi, bij;
  int incai, incaij;
  int incbi, incbij;

  int i, j;

  double elem;
  double *a_i = a;
  const double *b_i = b;

  if (order == blas_colmajor) {
    incai = 1;
    incbi = 1;
    incaij = lda;
    incbij = ldb;
  } else {
    incai = lda;
    incbi = ldb;
    incaij = 1;
    incbij = 1;
  }






  for (i = 0, ai = 0, bi = 0; i < m; i++, ai += incai, bi += incbi) {
    for (j = 0, aij = ai, bij = bi; j < n; j++, aij += incaij, bij += incbij) {
      elem = b_i[bij];
      a_i[aij] = elem;
    }
  }


}
void cgemm_copy(enum blas_order_type order, int m, int n,
		void *a, int lda, void *b, int ldb)
{

  int ai, aij;
  int bi, bij;
  int incai, incaij;
  int incbi, incbij;

  int i, j;

  float elem[2];
  float *a_i = (float *) a;
  const float *b_i = (float *) b;

  if (order == blas_colmajor) {
    incai = 1;
    incbi = 1;
    incaij = lda;
    incbij = ldb;
  } else {
    incai = lda;
    incbi = ldb;
    incaij = 1;
    incbij = 1;
  }

  incai *= 2;
  incaij *= 2;
  incbi *= 2;
  incbij *= 2;

  for (i = 0, ai = 0, bi = 0; i < m; i++, ai += incai, bi += incbi) {
    for (j = 0, aij = ai, bij = bi; j < n; j++, aij += incaij, bij += incbij) {
      elem[0] = b_i[bij];
      elem[1] = b_i[bij + 1];
      a_i[aij] = elem[0];
      a_i[aij + 1] = elem[1];
    }
  }


}
void zgemm_copy(enum blas_order_type order, int m, int n,
		void *a, int lda, void *b, int ldb)
{

  int ai, aij;
  int bi, bij;
  int incai, incaij;
  int incbi, incbij;

  int i, j;

  double elem[2];
  double *a_i = (double *) a;
  const double *b_i = (double *) b;

  if (order == blas_colmajor) {
    incai = 1;
    incbi = 1;
    incaij = lda;
    incbij = ldb;
  } else {
    incai = lda;
    incbi = ldb;
    incaij = 1;
    incbij = 1;
  }

  incai *= 2;
  incaij *= 2;
  incbi *= 2;
  incbij *= 2;

  for (i = 0, ai = 0, bi = 0; i < m; i++, ai += incai, bi += incbi) {
    for (j = 0, aij = ai, bij = bi; j < n; j++, aij += incaij, bij += incbij) {
      elem[0] = b_i[bij];
      elem[1] = b_i[bij + 1];
      a_i[aij] = elem[0];
      a_i[aij + 1] = elem[1];
    }
  }


}

void sgemm_zero(enum blas_order_type order, int m, int n, float *a, int lda)
{
  int i, j;
  int ai, aij;
  int incai, incaij;

  float a_elem;
  float *a_i = a;

  a_elem = 0.0;

  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }




  for (i = 0, ai = 0; i < m; i++, ai += incai) {
    for (j = 0; aij = ai, j < n; j++, aij += incaij) {
      a_i[aij] = a_elem;
    }
  }
}
void dgemm_zero(enum blas_order_type order, int m, int n, double *a, int lda)
{
  int i, j;
  int ai, aij;
  int incai, incaij;

  double a_elem;
  double *a_i = a;

  a_elem = 0.0;

  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }




  for (i = 0, ai = 0; i < m; i++, ai += incai) {
    for (j = 0; aij = ai, j < n; j++, aij += incaij) {
      a_i[aij] = a_elem;
    }
  }
}
void cgemm_zero(enum blas_order_type order, int m, int n, void *a, int lda)
{
  int i, j;
  int ai, aij;
  int incai, incaij;

  float a_elem[2];
  float *a_i = (float *) a;

  a_elem[0] = a_elem[1] = 0.0;

  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }

  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m; i++, ai += incai) {
    for (j = 0; aij = ai, j < n; j++, aij += incaij) {
      a_i[aij] = a_elem[0];
      a_i[aij + 1] = a_elem[1];
    }
  }
}
void zgemm_zero(enum blas_order_type order, int m, int n, void *a, int lda)
{
  int i, j;
  int ai, aij;
  int incai, incaij;

  double a_elem[2];
  double *a_i = (double *) a;

  a_elem[0] = a_elem[1] = 0.0;

  if (order == blas_colmajor) {
    incai = 1;
    incaij = lda;
  } else {
    incai = lda;
    incaij = 1;
  }

  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m; i++, ai += incai) {
    for (j = 0; aij = ai, j < n; j++, aij += incaij) {
      a_i[aij] = a_elem[0];
      a_i[aij + 1] = a_elem[1];
    }
  }
}

void sprint_matrix(float *a, int m, int n, int lda,
		   enum blas_order_type order)
{
  int ai, aij;
  int incai, incaij;
  int i, j;

  float a_elem;
  const float *a_i = a;

  if (order == blas_rowmajor) {
    incai = lda;
    incaij = 1;
  } else {
    incai = 1;
    incaij = lda;
  }




  for (i = 0, ai = 0; i < m; i++, ai += incai) {
    for (j = 0, aij = ai; j < n; j++, aij += incaij) {
      a_elem = a_i[aij];
      printf(" %.12g ", a_elem);

    }
    printf("\n");
  }

}
void dprint_matrix(double *a, int m, int n, int lda,
		   enum blas_order_type order)
{
  int ai, aij;
  int incai, incaij;
  int i, j;

  double a_elem;
  const double *a_i = a;

  if (order == blas_rowmajor) {
    incai = lda;
    incaij = 1;
  } else {
    incai = 1;
    incaij = lda;
  }




  for (i = 0, ai = 0; i < m; i++, ai += incai) {
    for (j = 0, aij = ai; j < n; j++, aij += incaij) {
      a_elem = a_i[aij];
      printf(" %.24g ", a_elem);

    }
    printf("\n");
  }

}
void cprint_matrix(void *a, int m, int n, int lda, enum blas_order_type order)
{
  int ai, aij;
  int incai, incaij;
  int i, j;

  float a_elem[2];
  const float *a_i = (float *) a;

  if (order == blas_rowmajor) {
    incai = lda;
    incaij = 1;
  } else {
    incai = 1;
    incaij = lda;
  }

  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m; i++, ai += incai) {
    for (j = 0, aij = ai; j < n; j++, aij += incaij) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      printf(" %.12g + %.12gi", ((float *) a_elem)[0], ((float *) a_elem)[1]);

    }
    printf("\n");
  }

}
void zprint_matrix(void *a, int m, int n, int lda, enum blas_order_type order)
{
  int ai, aij;
  int incai, incaij;
  int i, j;

  double a_elem[2];
  const double *a_i = (double *) a;

  if (order == blas_rowmajor) {
    incai = lda;
    incaij = 1;
  } else {
    incai = 1;
    incaij = lda;
  }

  incai *= 2;
  incaij *= 2;

  for (i = 0, ai = 0; i < m; i++, ai += incai) {
    for (j = 0, aij = ai; j < n; j++, aij += incaij) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      printf(" %.24g + %.24gi ", ((double *) a_elem)[0],
	     ((double *) a_elem)[1]);

    }
    printf("\n");
  }

}
