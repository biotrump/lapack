      REAL             FUNCTION SRZT02( M, N, AF, LDA, TAU, WORK,
     $                 LWORK )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL               AF( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  SRZT02 returns
*       || I - Q'*Q || / ( M * eps)
*  where the matrix Q is defined by the Householder transformations
*  generated by STZRZF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix AF.
*
*  N       (input) INTEGER
*          The number of columns of the matrix AF.
*
*  AF      (input) REAL array, dimension (LDA,N)
*          The output of STZRZF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array AF.
*
*  TAU     (input) REAL array, dimension (M)
*          Details of the Householder transformations as returned by
*          STZRZF.
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          length of WORK array. LWORK >= N*N+N*NB.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO
*     ..
*     .. Local Arrays ..
      REAL               RWORK( 1 )
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE
      EXTERNAL           SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLASET, SORMRZ, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, REAL
*     ..
*     .. Executable Statements ..
*
      SRZT02 = ZERO
*
      IF( LWORK.LT.N*N+N ) THEN
         CALL XERBLA( 'SRZT02', 7 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
*     Q := I
*
      CALL SLASET( 'Full', N, N, ZERO, ONE, WORK, N )
*
*     Q := P(1) * ... * P(m) * Q
*
      CALL SORMRZ( 'Left', 'No transpose', N, N, M, N-M, AF, LDA, TAU,
     $             WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO )
*
*     Q := P(m) * ... * P(1) * Q
*
      CALL SORMRZ( 'Left', 'Transpose', N, N, M, N-M, AF, LDA, TAU,
     $             WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO )
*
*     Q := Q - I
*
      DO 10 I = 1, N
         WORK( ( I-1 )*N+I ) = WORK( ( I-1 )*N+I ) - ONE
   10 CONTINUE
*
      SRZT02 = SLANGE( 'One-norm', N, N, WORK, N, RWORK ) /
     $         ( SLAMCH( 'Epsilon' )*REAL( MAX( M, N ) ) )
      RETURN
*
*     End of SRZT02
*
      END
