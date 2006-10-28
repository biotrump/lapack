      SUBROUTINE CSTEGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
     $           ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
     $           LIWORK, INFO )

      IMPLICIT NONE
*
*
*  -- LAPACK computational routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.
*     October 7, 2006
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE
      INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, M, N
      REAL             ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            ISUPPZ( * ), IWORK( * )
      REAL               D( * ), E( * ), W( * ), WORK( * )
      COMPLEX            Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  CSTEGR computes selected eigenvalues and, optionally, eigenvectors
*  of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
*  a well defined set of pairwise different real eigenvalues, the corresponding
*  real eigenvectors are pairwise orthogonal.
*
*  The spectrum may be computed either completely or partially by specifying
*  either an interval (VL,VU] or a range of indices IL:IU for the desired
*  eigenvalues.
*
*  CSTEGR is a compatability wrapper around the improved CSTEMR routine.
*  See SSTEMR for further details.
*
*  One important change is that the ABSTOL parameter no longer provides any
*  benefit and hence is no longer used.
*
*  Note : CSTEGR and CSTEMR work only on machines which follow
*  IEEE-754 floating-point standard in their handling of infinities and
*  NaNs.  Normal execution may create these exceptiona values and hence
*  may abort due to a floating point exception in environments which
*  do not conform to the IEEE-754 standard.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the half-open interval (VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) REAL array, dimension (N)
*          On entry, the N diagonal elements of the tridiagonal matrix
*          T. On exit, D is overwritten.
*
*  E       (input/output) REAL array, dimension (N)
*          On entry, the (N-1) subdiagonal elements of the tridiagonal
*          matrix T in elements 1 to N-1 of E. E(N) need not be set on
*          input, but is used internally as workspace.
*          On exit, E is overwritten.
*
*  VL      (input) REAL
*  VU      (input) REAL
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues. VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  ABSTOL  (input) REAL
*          Unused.  Was the absolute error tolerance for the
*          eigenvalues/eigenvectors in previous versions.
*
*  M       (output) INTEGER
*          The total number of eigenvalues found.  0 <= M <= N.
*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*
*  W       (output) REAL array, dimension (N)
*          The first M elements contain the selected eigenvalues in
*          ascending order.
*
*  Z       (output) COMPLEX array, dimension (LDZ, max(1,M) )
*          If JOBZ = 'V', and if INFO = 0, then the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix T
*          corresponding to the selected eigenvalues, with the i-th
*          column of Z holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*          Note: the user must ensure that at least max(1,M) columns are
*          supplied in the array Z; if RANGE = 'V', the exact value of M
*          is not known in advance and an upper bound must be used.
*          Supplying N columns is always safe.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', then LDZ >= max(1,N).
*
*  ISUPPZ  (output) INTEGER ARRAY, dimension ( 2*max(1,M) )
*          The support of the eigenvectors in Z, i.e., the indices
*          indicating the nonzero elements in Z. The i-th computed eigenvector
*          is nonzero only in elements ISUPPZ( 2*i-1 ) through
*          ISUPPZ( 2*i ). This is relevant in the case when the matrix
*          is split. ISUPPZ is only set if N>2.
*
*  WORK    (workspace/output) REAL array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal
*          (and minimal) LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,18*N)
*          if JOBZ = 'V', and LWORK >= max(1,12*N) if JOBZ = 'N'.
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.  LIWORK >= max(1,10*N)
*          if the eigenvectors are desired, and LIWORK >= max(1,8*N)
*          if only the eigenvalues are to be computed.
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (input/output) INTEGER
*          On entry, INFO
*          <>0: indicates that the code should check whether
*               the tridiagonal matrix defines its eigenvalues to high
*               relative accuracy.
*               If this is the case, the code uses relative-accuracy
*               preserving algorithms that might be (a bit) slower depending
*               on the matrix.
*               If the eigenvalues are not defined to high relative accuracy,
*               can use possibly faster algorithms.
*          = 0: the code is not required to guarantee relatively accurate
*               eigenvalues and can use the fastest possible techniques.
*          On exit, INFO
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = 1X, internal error in SLARRE,
*                if INFO = 2X, internal error in CLARRV.
*                Here, the digit X = ABS( IINFO ) < 10, where IINFO is
*                the nonzero error code returned by SLARRE or
*                CLARRV, respectively.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Inderjit Dhillon, IBM Almaden, USA
*     Osni Marques, LBNL/NERSC, USA
*     Christof Voemel, LBNL/NERSC, USA
*
*  =====================================================================
*
*     .. External Subroutines ..
      EXTERNAL CSTEMR
*     ..

      INFO = 0

      CALL CSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
     $                   M, W, Z, LDZ, N, ISUPPZ, WORK, LWORK, IWORK,
     $                   LIWORK, INFO )
*
*     End of CSTEGR
*
      END
