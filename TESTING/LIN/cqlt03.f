      SUBROUTINE CQLT03( M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK,
     $                   RWORK, RESULT )
*
*  -- LAPACK test routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL               RESULT( * ), RWORK( * )
      COMPLEX            AF( LDA, * ), C( LDA, * ), CC( LDA, * ),
     $                   Q( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  CQLT03 tests CUNMQL, which computes Q*C, Q'*C, C*Q or C*Q'.
*
*  CQLT03 compares the results of a call to CUNMQL with the results of
*  forming Q explicitly by a call to CUNGQL and then performing matrix
*  multiplication by a call to CGEMM.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The order of the orthogonal matrix Q.  M >= 0.
*
*  N       (input) INTEGER
*          The number of rows or columns of the matrix C; C is m-by-n if
*          Q is applied from the left, or n-by-m if Q is applied from
*          the right.  N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          orthogonal matrix Q.  M >= K >= 0.
*
*  AF      (input) COMPLEX array, dimension (LDA,N)
*          Details of the QL factorization of an m-by-n matrix, as
*          returned by CGEQLF. See CGEQLF for further details.
*
*  C       (workspace) COMPLEX array, dimension (LDA,N)
*
*  CC      (workspace) COMPLEX array, dimension (LDA,N)
*
*  Q       (workspace) COMPLEX array, dimension (LDA,M)
*
*  LDA     (input) INTEGER
*          The leading dimension of the arrays AF, C, CC, and Q.
*
*  TAU     (input) COMPLEX array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors corresponding
*          to the QL factorization in AF.
*
*  WORK    (workspace) COMPLEX array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The length of WORK.  LWORK must be at least M, and should be
*          M*NB, where NB is the blocksize for this environment.
*
*  RWORK   (workspace) REAL array, dimension (M)
*
*  RESULT  (output) REAL array, dimension (4)
*          The test ratios compare two techniques for multiplying a
*          random matrix C by an m-by-m orthogonal matrix Q.
*          RESULT(1) = norm( Q*C - Q*C )  / ( M * norm(C) * EPS )
*          RESULT(2) = norm( C*Q - C*Q )  / ( M * norm(C) * EPS )
*          RESULT(3) = norm( Q'*C - Q'*C )/ ( M * norm(C) * EPS )
*          RESULT(4) = norm( C*Q' - C*Q' )/ ( M * norm(C) * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            ROGUE
      PARAMETER          ( ROGUE = ( -1.0E+10, -1.0E+10 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, ISIDE, ITRANS, J, MC, MINMN, NC
      REAL               CNORM, EPS, RESID
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANGE, SLAMCH
      EXTERNAL           LSAME, CLANGE, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CLACPY, CLARNV, CLASET, CUNGQL, CUNMQL
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, MAX, MIN, REAL
*     ..
*     .. Scalars in Common ..
      CHARACTER(32)      SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               ISEED / 1988, 1989, 1990, 1991 /
*     ..
*     .. Executable Statements ..
*
      EPS = SLAMCH( 'Epsilon' )
      MINMN = MIN( M, N )
*
*     Quick return if possible
*
      IF( MINMN.EQ.0 ) THEN
         RESULT( 1 ) = ZERO
         RESULT( 2 ) = ZERO
         RESULT( 3 ) = ZERO
         RESULT( 4 ) = ZERO
         RETURN
      ENDIF
*
*     Copy the last k columns of the factorization to the array Q
*
      CALL CLASET( 'Full', M, M, ROGUE, ROGUE, Q, LDA )
      IF( K.GT.0 .AND. M.GT.K )
     $   CALL CLACPY( 'Full', M-K, K, AF( 1, N-K+1 ), LDA,
     $                Q( 1, M-K+1 ), LDA )
      IF( K.GT.1 )
     $   CALL CLACPY( 'Upper', K-1, K-1, AF( M-K+1, N-K+2 ), LDA,
     $                Q( M-K+1, M-K+2 ), LDA )
*
*     Generate the m-by-m matrix Q
*
      SRNAMT = 'CUNGQL'
      CALL CUNGQL( M, M, K, Q, LDA, TAU( MINMN-K+1 ), WORK, LWORK,
     $             INFO )
*
      DO 30 ISIDE = 1, 2
         IF( ISIDE.EQ.1 ) THEN
            SIDE = 'L'
            MC = M
            NC = N
         ELSE
            SIDE = 'R'
            MC = N
            NC = M
         END IF
*
*        Generate MC by NC matrix C
*
         DO 10 J = 1, NC
            CALL CLARNV( 2, ISEED, MC, C( 1, J ) )
   10    CONTINUE
         CNORM = CLANGE( '1', MC, NC, C, LDA, RWORK )
         IF( CNORM.EQ.ZERO )
     $      CNORM = ONE
*
         DO 20 ITRANS = 1, 2
            IF( ITRANS.EQ.1 ) THEN
               TRANS = 'N'
            ELSE
               TRANS = 'C'
            END IF
*
*           Copy C
*
            CALL CLACPY( 'Full', MC, NC, C, LDA, CC, LDA )
*
*           Apply Q or Q' to C
*
            SRNAMT = 'CUNMQL'
            IF( K.GT.0 )
     $         CALL CUNMQL( SIDE, TRANS, MC, NC, K, AF( 1, N-K+1 ),
     $                      LDA, TAU( MINMN-K+1 ), CC, LDA, WORK,
     $                      LWORK, INFO )
*
*           Form explicit product and subtract
*
            IF( LSAME( SIDE, 'L' ) ) THEN
               CALL CGEMM( TRANS, 'No transpose', MC, NC, MC,
     $                     CMPLX( -ONE ), Q, LDA, C, LDA, CMPLX( ONE ),
     $                     CC, LDA )
            ELSE
               CALL CGEMM( 'No transpose', TRANS, MC, NC, NC,
     $                     CMPLX( -ONE ), C, LDA, Q, LDA, CMPLX( ONE ),
     $                     CC, LDA )
            END IF
*
*           Compute error in the difference
*
            RESID = CLANGE( '1', MC, NC, CC, LDA, RWORK )
            RESULT( ( ISIDE-1 )*2+ITRANS ) = RESID /
     $         ( REAL( MAX( 1, M ) )*CNORM*EPS )
*
   20    CONTINUE
   30 CONTINUE
*
      RETURN
*
*     End of CQLT03
*
      END
