      SUBROUTINE SGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK,
     $                   INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LWORK, M, N, P
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( * ), D( * ),
     $                   WORK( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  SGGLSE solves the linear equality-constrained least squares (LSE)
*  problem:
*
*          minimize || c - A*x ||_2   subject to   B*x = d
*
*  where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
*  M-vector, and d is a given P-vector. It is assumed that
*  P <= N <= M+P, and
*
*           rank(B) = P and  rank( (A) ) = N.
*                                ( (B) )
*
*  These conditions ensure that the LSE problem has a unique solution,
*  which is obtained using a generalized RQ factorization of the
*  matrices (B, A) given by
*
*     B = (0 R)*Q,   A = Z*T*Q.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrices A and B. N >= 0.
*
*  P       (input) INTEGER
*          The number of rows of the matrix B. 0 <= P <= N <= M+P.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(M,N)-by-N upper trapezoidal matrix T.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  B       (input/output) REAL array, dimension (LDB,N)
*          On entry, the P-by-N matrix B.
*          On exit, the upper triangle of the subarray B(1:P,N-P+1:N)
*          contains the P-by-P upper triangular matrix R.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,P).
*
*  C       (input/output) REAL array, dimension (M)
*          On entry, C contains the right hand side vector for the
*          least squares part of the LSE problem.
*          On exit, the residual sum of squares for the solution
*          is given by the sum of squares of elements N-P+1 to M of
*          vector C.
*
*  D       (input/output) REAL array, dimension (P)
*          On entry, D contains the right hand side vector for the
*          constrained equation.
*          On exit, D is destroyed.
*
*  X       (output) REAL array, dimension (N)
*          On exit, X is the solution of the LSE problem.
*
*  WORK    (workspace/output) REAL array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,M+N+P).
*          For optimum performance LWORK >= P+min(M,N)+max(M,N)*NB,
*          where NB is an upper bound for the optimal blocksizes for
*          SGEQRF, SGERQF, SORMQR and SORMRQ.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          = 1:  the upper triangular factor R associated with B in the
*                generalized RQ factorization of the pair (B, A) is
*                singular, so that rank(B) < P; the least squares
*                solution could not be computed.
*          = 2:  the (N-P) by (N-P) part of the upper trapezoidal factor
*                T associated with A in the generalized RQ factorization
*                of the pair (B, A) is singular, so that
*                rank( (A) ) < N; the least squares solution could not
*                    ( (B) )
*                be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            LOPT, LWKMIN, LWKOPT, MN, NB, NB1, NB2, NB3,
     $                   NB4, NR
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SCOPY, SGEMV, SGGRQF, SORMQR, SORMRQ,
     $                   STRMV, STRTRS, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV 
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 .OR. P.GT.N .OR. P.LT.N-M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -7
      END IF
*
*     Calculate workspace
*
      IF( INFO.EQ.0) THEN
         IF( N.EQ.0 ) THEN
            LWKMIN = 1
            LWKOPT = 1
         ELSE
            NB1 = ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 )
            NB2 = ILAENV( 1, 'SGERQF', ' ', M, N, -1, -1 )
            NB3 = ILAENV( 1, 'SORMQR', ' ', M, N, P, -1 )
            NB4 = ILAENV( 1, 'SORMRQ', ' ', M, N, P, -1 )
            NB = MAX( NB1, NB2, NB3, NB4 )
            LWKMIN = M + N + P
            LWKOPT = P + MN + MAX( M, N )*NB
         END IF
         WORK( 1 ) = LWKOPT
*
         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGGLSE', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Compute the GRQ factorization of matrices B and A:
*
*            B*Q' = (  0  T12 ) P   Z'*A*Q' = ( R11 R12 ) N-P
*                     N-P  P                  (  0  R22 ) M+P-N
*                                               N-P  P
*
*     where T12 and R11 are upper triangular, and Q and Z are
*     orthogonal.
*
      CALL SGGRQF( P, M, N, B, LDB, WORK, A, LDA, WORK( P+1 ),
     $             WORK( P+MN+1 ), LWORK-P-MN, INFO )
      LOPT = WORK( P+MN+1 )
*
*     Update c = Z'*c = ( c1 ) N-P
*                       ( c2 ) M+P-N
*
      CALL SORMQR( 'Left', 'Transpose', M, 1, MN, A, LDA, WORK( P+1 ),
     $             C, MAX( 1, M ), WORK( P+MN+1 ), LWORK-P-MN, INFO )
      LOPT = MAX( LOPT, INT( WORK( P+MN+1 ) ) )
*
*     Solve T12*x2 = d for x2
*
      CALL STRTRS( 'Upper', 'No transpose', 'Non-unit', P, 1,
     $             B( 1, N-P+1 ), LDB, D, P, INFO )
*
      IF( INFO.GT.0 ) THEN
         INFO = 1
         RETURN
      END IF
*
*     Update c1
*
      CALL SGEMV( 'No transpose', N-P, P, -ONE, A( 1, N-P+1 ), LDA, D,
     $            1, ONE, C, 1 )
*
*     Sovle R11*x1 = c1 for x1
*
      CALL STRTRS( 'Upper', 'No transpose', 'Non-unit', N-P, 1, A, LDA,
     $             C, N-P, INFO )
*
      IF( INFO.GT.0 ) THEN
         INFO = 2
         RETURN
      END IF
*
*     Put the solutions in X
*
      CALL SCOPY( N-P, C, 1, X, 1 )
      CALL SCOPY( P, D, 1, X( N-P+1 ), 1 )
*
*     Compute the residual vector:
*
      IF( M.LT.N ) THEN
         NR = M + P - N
         IF( NR.GT.0 )
     $      CALL SGEMV( 'No transpose', NR, N-M, -ONE, A( N-P+1, M+1 ),
     $                  LDA, D( NR+1 ), 1, ONE, C( N-P+1 ), 1 )
      ELSE
         NR = P
      END IF
      IF( NR.GT.0 ) THEN
         CALL STRMV( 'Upper', 'No transpose', 'Non unit', NR,
     $               A( N-P+1, N-P+1 ), LDA, D, 1 )
         CALL SAXPY( NR, -ONE, D, 1, C( N-P+1 ), 1 )
      END IF
*
*     Backward transformation x = Q'*x
*
      CALL SORMRQ( 'Left', 'Transpose', N, 1, P, B, LDB, WORK( 1 ), X,
     $             N, WORK( P+MN+1 ), LWORK-P-MN, INFO )
      WORK( 1 ) = P + MN + MAX( LOPT, INT( WORK( P+MN+1 ) ) )
*
      RETURN
*
*     End of SGGLSE
*
      END
