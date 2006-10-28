      SUBROUTINE SLARRR( N, D, E, INFO )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.
*     October 7, 2006
*
*     .. Scalar Arguments ..
      INTEGER            N, INFO
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * )
*     ..
*
*
*  Purpose
*  =======
*
*  Perform tests to decide whether the symmetric tridiagonal matrix T
*  warrants expensive computations which guarantee high relative accuracy
*  in the eigenvalues.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix. N > 0.
*
*  D       (input) REAL             array, dimension (N)
*          The N diagonal elements of the tridiagonal matrix T.
*
*  E       (input/output) REAL             array, dimension (N)
*          On entry, the first (N-1) entries contain the subdiagonal
*          elements of the tridiagonal matrix T; E(N) is set to ZERO.
*
*  INFO    (output) INTEGER
*          INFO = 0(default) : the matrix warrants computations preserving
*                              relative accuracy.
*          INFO = 1          : the matrix warrants computations guaranteeing
*                              only absolute accuracy.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Beresford Parlett, University of California, Berkeley, USA
*     Jim Demmel, University of California, Berkeley, USA
*     Inderjit Dhillon, University of Texas, Austin, USA
*     Osni Marques, LBNL/NERSC, USA
*     Christof Voemel, University of California, Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, RELCOND
      PARAMETER          ( ZERO = 0.0E0,
     $                     RELCOND = 0.999E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      LOGICAL            YESREL
      REAL               EPS, SAFMIN, SMLNUM, RMIN, TMP, TMP2,
     $          OFFDIG, OFFDIG2

*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     As a default, do NOT go for relative-accuracy preserving computations.
      INFO = 1

      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      RMIN = SQRT( SMLNUM )

*     Tests for relative accuracy
*
*     Test for scaled diagonal dominance
*     Scale the diagonal entries to one and check whether the sum of the
*     off-diagonals is less than one
*
      YESREL = .TRUE.
      OFFDIG = ZERO
      TMP = SQRT(ABS(D(1)))
      IF (TMP.LT.RMIN) YESREL = .FALSE.
      IF(.NOT.YESREL) GOTO 11
      DO 10 I = 2, N
         TMP2 = SQRT(ABS(D(I)))
         IF (TMP2.LT.RMIN) YESREL = .FALSE.
         IF(.NOT.YESREL) GOTO 11
         OFFDIG2 = ABS(E(I-1))/(TMP*TMP2)
         IF(OFFDIG+OFFDIG2.GE.RELCOND) YESREL = .FALSE.
         IF(.NOT.YESREL) GOTO 11
         TMP = TMP2
         OFFDIG = OFFDIG2
 10   CONTINUE
 11   CONTINUE

      IF( YESREL ) THEN
         INFO = 0
         RETURN
      ELSE
      ENDIF
*

*
*     *** MORE TO BE IMPLEMENTED ***
*

*
*     Test if the lower bidiagonal matrix L from T = L D L^T
*     (zero shift facto) is well conditioned
*

*
*     Test if the upper bidiagonal matrix U from T = U D U^T
*     (zero shift facto) is well conditioned.
*     In this case, the matrix needs to be flipped and, at the end
*     of the eigenvector computation, the flip needs to be applied
*     to the computed eigenvectors (and the support)
*

*
      RETURN
*
*     END OF SLARRR
*
      END
