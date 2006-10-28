      FUNCTION SLANEG( N, D, LLD, SIGMA, PIVMIN, R )
      INTEGER SLANEG
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.
*     October 7, 2006
*
*     .. Scalar Arguments ..
      INTEGER            N, R
      REAL               PIVMIN, SIGMA
*     ..
*     .. Array Arguments ..
      REAL               D( * ), LLD( * )
*     ..
*
*  Purpose
*  =======
*
*  SLANEG computes the Sturm count, the number of negative pivots
*  encountered while factoring tridiagonal T - sigma I = L D L^T.
*  This implementation works directly on the factors without forming
*  the tridiagonal matrix T.  The Sturm count is also the number of
*  eigenvalues of T less than sigma.
*
*  This routine is called from SLARRB.
*
*  Note : This routine requires IEEE-754 propagation of Infinities and
*  NaNs as well as default handling of 0/0 => NaN.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  D       (input) REAL             array, dimension (N)
*          The N diagonal elements of the diagonal matrix D.
*
*  LLD     (input) REAL             array, dimension (N-1)
*          The (N-1) elements L(i)*L(i)*D(i).
*
*  SIGMA   (input) REAL            
*          Shift amount in T - sigma I = L D L^T.
*
*  PIVMIN  (input) REAL            
*          The minimum pivot in the Sturm sequence.  Used when zero
*          pivots are encountered.
*
*  R       (input) INTEGER
*          The twist index for the twisted factorization that is used
*          for the negcount.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER        ( ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, NEG1, NEG2, NEGCNT
      REAL               DMINUS, DPLUS, GAMMA, P, S, T, TMP
      LOGICAL SAWNAN
*     ..
*     .. External Functions ..
      LOGICAL SISNAN
      EXTERNAL SISNAN
*     .. Executable Statements ..

      NEGCNT = 0

*     I) upper part: L D L^T - SIGMA I = L+ D+ L+^T
      NEG1 = 0
      S = ZERO
      DO 21 J = 1, R - 1
         T = S - SIGMA
         DPLUS = D( J ) + T
         S = T*LLD( J ) / DPLUS
         IF( DPLUS.LT.ZERO ) NEG1 = NEG1 + 1
 21   CONTINUE
      SAWNAN = SISNAN( S )
*     Run a slower version of the above loop if a NaN is detected
      IF( SAWNAN ) THEN
         NEG1 = 0
         S = ZERO
         T = -SIGMA
         DO 22 J = 1, R - 1
            DPLUS = D( J ) + T
            IF(ABS(DPLUS).LT.PIVMIN) DPLUS = -PIVMIN
            TMP = LLD( J ) / DPLUS
            IF( DPLUS.LT.ZERO ) NEG1 = NEG1 + 1
            S = T*TMP
            IF( TMP.EQ.ZERO ) S = LLD( J )
            T = S - SIGMA
 22      CONTINUE
      END IF
      NEGCNT = NEGCNT + NEG1
*
*     II) lower part: L D L^T - SIGMA I = U- D- U-^T
      NEG2 = 0
      P = D( N ) - SIGMA
      DO 23 J = N - 1, R, -1
         DMINUS = LLD( J ) + P
         P = P*D( J )/DMINUS - SIGMA
         IF( DMINUS.LT.ZERO ) NEG2 = NEG2 + 1
 23   CONTINUE
      SAWNAN = SISNAN( P )
      IF( SAWNAN ) THEN
         NEG2 = 0
         P = D( N ) - SIGMA
         DO 24 J = N - 1, R, -1
            DMINUS = LLD( J ) + P
            IF(ABS(DMINUS).LT.PIVMIN) DMINUS = -PIVMIN
            TMP = D( J ) / DMINUS
            IF( DMINUS.LT.ZERO ) NEG2 = NEG2 + 1
            P = P*TMP - SIGMA
            IF( TMP.EQ.ZERO ) P = D( J ) - SIGMA
 24      CONTINUE
      END IF
      NEGCNT = NEGCNT + NEG2
*
*     III) Twist index
      GAMMA = S + P
      IF( GAMMA.LT.ZERO ) NEGCNT = NEGCNT+1

      SLANEG = NEGCNT
      END
