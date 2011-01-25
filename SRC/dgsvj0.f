      SUBROUTINE DGSVJ0( JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS,
     +                   SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.3.1)                                  --
*
*  -- Contributed by Zlatko Drmac of the University of Zagreb and     --
*  -- Kresimir Veselic of the Fernuniversitaet Hagen                  --
*     January 2011
*
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
* This routine is also part of SIGMA (version 1.23, October 23. 2008.)
* SIGMA is a library of algorithms for highly accurate algorithms for
* computation of SVD, PSVD, QSVD, (H,K)-SVD, and for solution of the
* eigenvalue problems Hx = lambda M x, H M x = lambda x with H, M > 0.
*
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDV, LWORK, M, MV, N, NSWEEP
      DOUBLE PRECISION   EPS, SFMIN, TOL
      CHARACTER*1        JOBV
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), SVA( N ), D( N ), V( LDV, * ),
     +                   WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  DGSVJ0 is called from DGESVJ as a pre-processor and that is its main
*  purpose. It applies Jacobi rotations in the same way as DGESVJ does, but
*  it does not check convergence (stopping criterion). Few tuning
*  parameters (marked by [TP]) are available for the implementer.
*
*  Further Details
*  ~~~~~~~~~~~~~~~
*  DGSVJ0 is used just to enable SGESVJ to call a simplified version of
*  itself to work on a submatrix of the original matrix.
*
*  Contributors
*  ~~~~~~~~~~~~
*  Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany)
*
*  Bugs, Examples and Comments
*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~
*  Please report all bugs and send interesting test examples and comments to
*  drmac@math.hr. Thank you.
*
*  Arguments
*  =========
*
*  JOBV    (input) CHARACTER*1
*          Specifies whether the output from this procedure is used
*          to compute the matrix V:
*          = 'V': the product of the Jacobi rotations is accumulated
*                 by postmulyiplying the N-by-N array V.
*                (See the description of V.)
*          = 'A': the product of the Jacobi rotations is accumulated
*                 by postmulyiplying the MV-by-N array V.
*                (See the descriptions of MV and V.)
*          = 'N': the Jacobi rotations are not accumulated.
*
*  M       (input) INTEGER
*          The number of rows of the input matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the input matrix A.
*          M >= N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, M-by-N matrix A, such that A*diag(D) represents
*          the input matrix.
*          On exit,
*          A_onexit * D_onexit represents the input matrix A*diag(D)
*          post-multiplied by a sequence of Jacobi rotations, where the
*          rotation threshold and the total number of sweeps are given in
*          TOL and NSWEEP, respectively.
*          (See the descriptions of D, TOL and NSWEEP.)
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  D       (input/workspace/output) DOUBLE PRECISION array, dimension (N)
*          The array D accumulates the scaling factors from the fast scaled
*          Jacobi rotations.
*          On entry, A*diag(D) represents the input matrix.
*          On exit, A_onexit*diag(D_onexit) represents the input matrix
*          post-multiplied by a sequence of Jacobi rotations, where the
*          rotation threshold and the total number of sweeps are given in
*          TOL and NSWEEP, respectively.
*          (See the descriptions of A, TOL and NSWEEP.)
*
*  SVA     (input/workspace/output) DOUBLE PRECISION array, dimension (N)
*          On entry, SVA contains the Euclidean norms of the columns of
*          the matrix A*diag(D).
*          On exit, SVA contains the Euclidean norms of the columns of
*          the matrix onexit*diag(D_onexit).
*
*  MV      (input) INTEGER
*          If JOBV .EQ. 'A', then MV rows of V are post-multipled by a
*                           sequence of Jacobi rotations.
*          If JOBV = 'N',   then MV is not referenced.
*
*  V       (input/output) DOUBLE PRECISION array, dimension (LDV,N)
*          If JOBV .EQ. 'V' then N rows of V are post-multipled by a
*                           sequence of Jacobi rotations.
*          If JOBV .EQ. 'A' then MV rows of V are post-multipled by a
*                           sequence of Jacobi rotations.
*          If JOBV = 'N',   then V is not referenced.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V,  LDV >= 1.
*          If JOBV = 'V', LDV .GE. N.
*          If JOBV = 'A', LDV .GE. MV.
*
*  EPS     (input) DOUBLE PRECISION
*          EPS = DLAMCH('Epsilon')
*
*  SFMIN   (input) DOUBLE PRECISION
*          SFMIN = DLAMCH('Safe Minimum')
*
*  TOL     (input) DOUBLE PRECISION
*          TOL is the threshold for Jacobi rotations. For a pair
*          A(:,p), A(:,q) of pivot columns, the Jacobi rotation is
*          applied only if DABS(COS(angle(A(:,p),A(:,q)))) .GT. TOL.
*
*  NSWEEP  (input) INTEGER
*          NSWEEP is the number of sweeps of Jacobi rotations to be
*          performed.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          LWORK is the dimension of WORK. LWORK .GE. M.
*
*  INFO    (output) INTEGER
*          = 0 : successful exit.
*          < 0 : if INFO = -i, then the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Local Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0,
     +                   TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AAPP, AAPP0, AAPQ, AAQQ, APOAQ, AQOAP, BIG,
     +                   BIGTHETA, CS, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS,
     +                   ROOTSFMIN, ROOTTOL, SMALL, SN, T, TEMP1, THETA,
     +                   THSIGN
      INTEGER            BLSKIP, EMPTSW, i, ibr, IERR, igl, IJBLSK, ir1,
     +                   ISWROT, jbc, jgl, KBL, LKAHEAD, MVL, NBL,
     +                   NOTROT, p, PSKIPPED, q, ROWSKIP, SWBAND
      LOGICAL            APPLV, ROTOK, RSVEC
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   FASTR( 5 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DABS, DMAX1, DBLE, MIN0, DSIGN, DSQRT
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DDOT, DNRM2
      INTEGER            IDAMAX
      LOGICAL            LSAME
      EXTERNAL           IDAMAX, LSAME, DDOT, DNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DLASCL, DLASSQ, DROTM, DSWAP
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      APPLV = LSAME( JOBV, 'A' )
      RSVEC = LSAME( JOBV, 'V' )
      IF( .NOT.( RSVEC .OR. APPLV .OR. LSAME( JOBV, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( N.LT.0 ) .OR. ( N.GT.M ) ) THEN
         INFO = -3
      ELSE IF( LDA.LT.M ) THEN
         INFO = -5
      ELSE IF( ( RSVEC.OR.APPLV ) .AND. ( MV.LT.0 ) ) THEN
         INFO = -8
      ELSE IF( ( RSVEC.AND.( LDV.LT.N ) ).OR. 
     &         ( APPLV.AND.( LDV.LT.MV ) ) ) THEN
         INFO = -10
      ELSE IF( TOL.LE.EPS ) THEN
         INFO = -13
      ELSE IF( NSWEEP.LT.0 ) THEN
         INFO = -14
      ELSE IF( LWORK.LT.M ) THEN
         INFO = -16
      ELSE
         INFO = 0
      END IF
*
*     #:(
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGSVJ0', -INFO )
         RETURN
      END IF
*
      IF( RSVEC ) THEN
         MVL = N
      ELSE IF( APPLV ) THEN
         MVL = MV
      END IF
      RSVEC = RSVEC .OR. APPLV

      ROOTEPS = DSQRT( EPS )
      ROOTSFMIN = DSQRT( SFMIN )
      SMALL = SFMIN / EPS
      BIG = ONE / SFMIN
      ROOTBIG = ONE / ROOTSFMIN
      BIGTHETA = ONE / ROOTEPS
      ROOTTOL = DSQRT( TOL )
*
*     -#- Row-cyclic Jacobi SVD algorithm with column pivoting -#-
*
      EMPTSW = ( N*( N-1 ) ) / 2
      NOTROT = 0
      FASTR( 1 ) = ZERO
*
*     -#- Row-cyclic pivot strategy with de Rijk's pivoting -#-
*

      SWBAND = 0
*[TP] SWBAND is a tuning parameter. It is meaningful and effective
*     if SGESVJ is used as a computational routine in the preconditioned
*     Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure
*     ......

      KBL = MIN0( 8, N )
*[TP] KBL is a tuning parameter that defines the tile size in the
*     tiling of the p-q loops of pivot pairs. In general, an optimal
*     value of KBL depends on the matrix dimensions and on the
*     parameters of the computer's memory.
*
      NBL = N / KBL
      IF( ( NBL*KBL ).NE.N )NBL = NBL + 1

      BLSKIP = ( KBL**2 ) + 1
*[TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.

      ROWSKIP = MIN0( 5, KBL )
*[TP] ROWSKIP is a tuning parameter.

      LKAHEAD = 1
*[TP] LKAHEAD is a tuning parameter.
      SWBAND = 0
      PSKIPPED = 0
*
      DO 1993 i = 1, NSWEEP
*     .. go go go ...
*
         MXAAPQ = ZERO
         MXSINJ = ZERO
         ISWROT = 0
*
         NOTROT = 0
         PSKIPPED = 0
*
         DO 2000 ibr = 1, NBL

            igl = ( ibr-1 )*KBL + 1
*
            DO 1002 ir1 = 0, MIN0( LKAHEAD, NBL-ibr )
*
               igl = igl + ir1*KBL
*
               DO 2001 p = igl, MIN0( igl+KBL-1, N-1 )

*     .. de Rijk's pivoting
                  q = IDAMAX( N-p+1, SVA( p ), 1 ) + p - 1
                  IF( p.NE.q ) THEN
                     CALL DSWAP( M, A( 1, p ), 1, A( 1, q ), 1 )
                     IF( RSVEC )CALL DSWAP( MVL, V( 1, p ), 1,
     +                                      V( 1, q ), 1 )
                     TEMP1 = SVA( p )
                     SVA( p ) = SVA( q )
                     SVA( q ) = TEMP1
                     TEMP1 = D( p )
                     D( p ) = D( q )
                     D( q ) = TEMP1
                  END IF
*
                  IF( ir1.EQ.0 ) THEN
*
*        Column norms are periodically updated by explicit
*        norm computation.
*        Caveat:
*        Some BLAS implementations compute DNRM2(M,A(1,p),1)
*        as DSQRT(DDOT(M,A(1,p),1,A(1,p),1)), which may result in
*        overflow for ||A(:,p)||_2 > DSQRT(overflow_threshold), and
*        undeflow for ||A(:,p)||_2 < DSQRT(underflow_threshold).
*        Hence, DNRM2 cannot be trusted, not even in the case when
*        the true norm is far from the under(over)flow boundaries.
*        If properly implemented DNRM2 is available, the IF-THEN-ELSE
*        below should read "AAPP = DNRM2( M, A(1,p), 1 ) * D(p)".
*
                     IF( ( SVA( p ).LT.ROOTBIG ) .AND.
     +                   ( SVA( p ).GT.ROOTSFMIN ) ) THEN
                        SVA( p ) = DNRM2( M, A( 1, p ), 1 )*D( p )
                     ELSE
                        TEMP1 = ZERO
                        AAPP = ONE
                        CALL DLASSQ( M, A( 1, p ), 1, TEMP1, AAPP )
                        SVA( p ) = TEMP1*DSQRT( AAPP )*D( p )
                     END IF
                     AAPP = SVA( p )
                  ELSE
                     AAPP = SVA( p )
                  END IF

*
                  IF( AAPP.GT.ZERO ) THEN
*
                     PSKIPPED = 0
*
                     DO 2002 q = p + 1, MIN0( igl+KBL-1, N )
*
                        AAQQ = SVA( q )

                        IF( AAQQ.GT.ZERO ) THEN
*
                           AAPP0 = AAPP
                           IF( AAQQ.GE.ONE ) THEN
                              ROTOK = ( SMALL*AAPP ).LE.AAQQ
                              IF( AAPP.LT.( BIG / AAQQ ) ) THEN
                                 AAPQ = ( DDOT( M, A( 1, p ), 1, A( 1,
     +                                  q ), 1 )*D( p )*D( q ) / AAQQ )
     +                                  / AAPP
                              ELSE
                                 CALL DCOPY( M, A( 1, p ), 1, WORK, 1 )
                                 CALL DLASCL( 'G', 0, 0, AAPP, D( p ),
     +                                        M, 1, WORK, LDA, IERR )
                                 AAPQ = DDOT( M, WORK, 1, A( 1, q ),
     +                                  1 )*D( q ) / AAQQ
                              END IF
                           ELSE
                              ROTOK = AAPP.LE.( AAQQ / SMALL )
                              IF( AAPP.GT.( SMALL / AAQQ ) ) THEN
                                 AAPQ = ( DDOT( M, A( 1, p ), 1, A( 1,
     +                                  q ), 1 )*D( p )*D( q ) / AAQQ )
     +                                  / AAPP
                              ELSE
                                 CALL DCOPY( M, A( 1, q ), 1, WORK, 1 )
                                 CALL DLASCL( 'G', 0, 0, AAQQ, D( q ),
     +                                        M, 1, WORK, LDA, IERR )
                                 AAPQ = DDOT( M, WORK, 1, A( 1, p ),
     +                                  1 )*D( p ) / AAPP
                              END IF
                           END IF
*
                           MXAAPQ = DMAX1( MXAAPQ, DABS( AAPQ ) )
*
*        TO rotate or NOT to rotate, THAT is the question ...
*
                           IF( DABS( AAPQ ).GT.TOL ) THEN
*
*           .. rotate
*           ROTATED = ROTATED + ONE
*
                              IF( ir1.EQ.0 ) THEN
                                 NOTROT = 0
                                 PSKIPPED = 0
                                 ISWROT = ISWROT + 1
                              END IF
*
                              IF( ROTOK ) THEN
*
                                 AQOAP = AAQQ / AAPP
                                 APOAQ = AAPP / AAQQ
                                 THETA = -HALF*DABS( AQOAP-APOAQ )/AAPQ
*
                                 IF( DABS( THETA ).GT.BIGTHETA ) THEN
*
                                    T = HALF / THETA
                                    FASTR( 3 ) = T*D( p ) / D( q )
                                    FASTR( 4 ) = -T*D( q ) / D( p )
                                    CALL DROTM( M, A( 1, p ), 1,
     +                                          A( 1, q ), 1, FASTR )
                                    IF( RSVEC )CALL DROTM( MVL,
     +                                              V( 1, p ), 1,
     +                                              V( 1, q ), 1,
     +                                              FASTR )
                                    SVA( q ) = AAQQ*DSQRT( DMAX1( ZERO,
     +                                         ONE+T*APOAQ*AAPQ ) )
                                    AAPP = AAPP*DSQRT( DMAX1( ZERO, 
     +                                     ONE-T*AQOAP*AAPQ ) )
                                    MXSINJ = DMAX1( MXSINJ, DABS( T ) )
*
                                 ELSE
*
*                 .. choose correct signum for THETA and rotate
*
                                    THSIGN = -DSIGN( ONE, AAPQ )
                                    T = ONE / ( THETA+THSIGN*
     +                                  DSQRT( ONE+THETA*THETA ) )
                                    CS = DSQRT( ONE / ( ONE+T*T ) )
                                    SN = T*CS
*
                                    MXSINJ = DMAX1( MXSINJ, DABS( SN ) )
                                    SVA( q ) = AAQQ*DSQRT( DMAX1( ZERO,
     +                                         ONE+T*APOAQ*AAPQ ) )
                                    AAPP = AAPP*DSQRT( DMAX1( ZERO,
     +                                     ONE-T*AQOAP*AAPQ ) )
*
                                    APOAQ = D( p ) / D( q )
                                    AQOAP = D( q ) / D( p )
                                    IF( D( p ).GE.ONE ) THEN
                                       IF( D( q ).GE.ONE ) THEN
                                          FASTR( 3 ) = T*APOAQ
                                          FASTR( 4 ) = -T*AQOAP
                                          D( p ) = D( p )*CS
                                          D( q ) = D( q )*CS
                                          CALL DROTM( M, A( 1, p ), 1,
     +                                                A( 1, q ), 1,
     +                                                FASTR )
                                          IF( RSVEC )CALL DROTM( MVL,
     +                                        V( 1, p ), 1, V( 1, q ),
     +                                        1, FASTR )
                                       ELSE
                                          CALL DAXPY( M, -T*AQOAP,
     +                                                A( 1, q ), 1,
     +                                                A( 1, p ), 1 )
                                          CALL DAXPY( M, CS*SN*APOAQ,
     +                                                A( 1, p ), 1,
     +                                                A( 1, q ), 1 )
                                          D( p ) = D( p )*CS
                                          D( q ) = D( q ) / CS
                                          IF( RSVEC ) THEN
                                             CALL DAXPY( MVL, -T*AQOAP,
     +                                                   V( 1, q ), 1,
     +                                                   V( 1, p ), 1 )
                                             CALL DAXPY( MVL,
     +                                                   CS*SN*APOAQ,
     +                                                   V( 1, p ), 1,
     +                                                   V( 1, q ), 1 )
                                          END IF
                                       END IF
                                    ELSE
                                       IF( D( q ).GE.ONE ) THEN
                                          CALL DAXPY( M, T*APOAQ,
     +                                                A( 1, p ), 1,
     +                                                A( 1, q ), 1 )
                                          CALL DAXPY( M, -CS*SN*AQOAP,
     +                                                A( 1, q ), 1,
     +                                                A( 1, p ), 1 )
                                          D( p ) = D( p ) / CS
                                          D( q ) = D( q )*CS
                                          IF( RSVEC ) THEN
                                             CALL DAXPY( MVL, T*APOAQ,
     +                                                   V( 1, p ), 1,
     +                                                   V( 1, q ), 1 )
                                             CALL DAXPY( MVL,
     +                                                   -CS*SN*AQOAP,
     +                                                   V( 1, q ), 1,
     +                                                   V( 1, p ), 1 )
                                          END IF
                                       ELSE
                                          IF( D( p ).GE.D( q ) ) THEN
                                             CALL DAXPY( M, -T*AQOAP,
     +                                                   A( 1, q ), 1,
     +                                                   A( 1, p ), 1 )
                                             CALL DAXPY( M, CS*SN*APOAQ,
     +                                                   A( 1, p ), 1,
     +                                                   A( 1, q ), 1 )
                                             D( p ) = D( p )*CS
                                             D( q ) = D( q ) / CS
                                             IF( RSVEC ) THEN
                                                CALL DAXPY( MVL,
     +                                               -T*AQOAP,
     +                                               V( 1, q ), 1,
     +                                               V( 1, p ), 1 )
                                                CALL DAXPY( MVL,
     +                                               CS*SN*APOAQ,
     +                                               V( 1, p ), 1,
     +                                               V( 1, q ), 1 )
                                             END IF
                                          ELSE
                                             CALL DAXPY( M, T*APOAQ,
     +                                                   A( 1, p ), 1,
     +                                                   A( 1, q ), 1 )
                                             CALL DAXPY( M,
     +                                                   -CS*SN*AQOAP,
     +                                                   A( 1, q ), 1,
     +                                                   A( 1, p ), 1 )
                                             D( p ) = D( p ) / CS
                                             D( q ) = D( q )*CS
                                             IF( RSVEC ) THEN
                                                CALL DAXPY( MVL,
     +                                               T*APOAQ, V( 1, p ),
     +                                               1, V( 1, q ), 1 )
                                                CALL DAXPY( MVL,
     +                                               -CS*SN*AQOAP,
     +                                               V( 1, q ), 1,
     +                                               V( 1, p ), 1 )
                                             END IF
                                          END IF
                                       END IF
                                    END IF
                                 END IF
*
                              ELSE
*              .. have to use modified Gram-Schmidt like transformation
                                 CALL DCOPY( M, A( 1, p ), 1, WORK, 1 )
                                 CALL DLASCL( 'G', 0, 0, AAPP, ONE, M,
     +                                        1, WORK, LDA, IERR )
                                 CALL DLASCL( 'G', 0, 0, AAQQ, ONE, M,
     +                                        1, A( 1, q ), LDA, IERR )
                                 TEMP1 = -AAPQ*D( p ) / D( q )
                                 CALL DAXPY( M, TEMP1, WORK, 1,
     +                                       A( 1, q ), 1 )
                                 CALL DLASCL( 'G', 0, 0, ONE, AAQQ, M,
     +                                        1, A( 1, q ), LDA, IERR )
                                 SVA( q ) = AAQQ*DSQRT( DMAX1( ZERO,
     +                                      ONE-AAPQ*AAPQ ) )
                                 MXSINJ = DMAX1( MXSINJ, SFMIN )
                              END IF
*           END IF ROTOK THEN ... ELSE
*
*           In the case of cancellation in updating SVA(q), SVA(p)
*           recompute SVA(q), SVA(p).
                              IF( ( SVA( q ) / AAQQ )**2.LE.ROOTEPS )
     +                            THEN
                                 IF( ( AAQQ.LT.ROOTBIG ) .AND.
     +                               ( AAQQ.GT.ROOTSFMIN ) ) THEN
                                    SVA( q ) = DNRM2( M, A( 1, q ), 1 )*
     +                                         D( q )
                                 ELSE
                                    T = ZERO
                                    AAQQ = ONE
                                    CALL DLASSQ( M, A( 1, q ), 1, T,
     +                                           AAQQ )
                                    SVA( q ) = T*DSQRT( AAQQ )*D( q )
                                 END IF
                              END IF
                              IF( ( AAPP / AAPP0 ).LE.ROOTEPS ) THEN
                                 IF( ( AAPP.LT.ROOTBIG ) .AND.
     +                               ( AAPP.GT.ROOTSFMIN ) ) THEN
                                    AAPP = DNRM2( M, A( 1, p ), 1 )*
     +                                     D( p )
                                 ELSE
                                    T = ZERO
                                    AAPP = ONE
                                    CALL DLASSQ( M, A( 1, p ), 1, T,
     +                                           AAPP )
                                    AAPP = T*DSQRT( AAPP )*D( p )
                                 END IF
                                 SVA( p ) = AAPP
                              END IF
*
                           ELSE
*        A(:,p) and A(:,q) already numerically orthogonal
                              IF( ir1.EQ.0 )NOTROT = NOTROT + 1
                              PSKIPPED = PSKIPPED + 1
                           END IF
                        ELSE
*        A(:,q) is zero column
                           IF( ir1.EQ.0 )NOTROT = NOTROT + 1
                           PSKIPPED = PSKIPPED + 1
                        END IF
*
                        IF( ( i.LE.SWBAND ) .AND.
     +                      ( PSKIPPED.GT.ROWSKIP ) ) THEN
                           IF( ir1.EQ.0 )AAPP = -AAPP
                           NOTROT = 0
                           GO TO 2103
                        END IF
*
 2002                CONTINUE
*     END q-LOOP
*
 2103                CONTINUE
*     bailed out of q-loop

                     SVA( p ) = AAPP

                  ELSE
                     SVA( p ) = AAPP
                     IF( ( ir1.EQ.0 ) .AND. ( AAPP.EQ.ZERO ) )
     +                   NOTROT = NOTROT + MIN0( igl+KBL-1, N ) - p
                  END IF
*
 2001          CONTINUE
*     end of the p-loop
*     end of doing the block ( ibr, ibr )
 1002       CONTINUE
*     end of ir1-loop
*
*........................................................
* ... go to the off diagonal blocks
*
            igl = ( ibr-1 )*KBL + 1
*
            DO 2010 jbc = ibr + 1, NBL
*
               jgl = ( jbc-1 )*KBL + 1
*
*        doing the block at ( ibr, jbc )
*
               IJBLSK = 0
               DO 2100 p = igl, MIN0( igl+KBL-1, N )
*
                  AAPP = SVA( p )
*
                  IF( AAPP.GT.ZERO ) THEN
*
                     PSKIPPED = 0
*
                     DO 2200 q = jgl, MIN0( jgl+KBL-1, N )
*
                        AAQQ = SVA( q )
*
                        IF( AAQQ.GT.ZERO ) THEN
                           AAPP0 = AAPP
*
*     -#- M x 2 Jacobi SVD -#-
*
*        -#- Safe Gram matrix computation -#-
*
                           IF( AAQQ.GE.ONE ) THEN
                              IF( AAPP.GE.AAQQ ) THEN
                                 ROTOK = ( SMALL*AAPP ).LE.AAQQ
                              ELSE
                                 ROTOK = ( SMALL*AAQQ ).LE.AAPP
                              END IF
                              IF( AAPP.LT.( BIG / AAQQ ) ) THEN
                                 AAPQ = ( DDOT( M, A( 1, p ), 1, A( 1,
     +                                  q ), 1 )*D( p )*D( q ) / AAQQ )
     +                                  / AAPP
                              ELSE
                                 CALL DCOPY( M, A( 1, p ), 1, WORK, 1 )
                                 CALL DLASCL( 'G', 0, 0, AAPP, D( p ),
     +                                        M, 1, WORK, LDA, IERR )
                                 AAPQ = DDOT( M, WORK, 1, A( 1, q ),
     +                                  1 )*D( q ) / AAQQ
                              END IF
                           ELSE
                              IF( AAPP.GE.AAQQ ) THEN
                                 ROTOK = AAPP.LE.( AAQQ / SMALL )
                              ELSE
                                 ROTOK = AAQQ.LE.( AAPP / SMALL )
                              END IF
                              IF( AAPP.GT.( SMALL / AAQQ ) ) THEN
                                 AAPQ = ( DDOT( M, A( 1, p ), 1, A( 1,
     +                                  q ), 1 )*D( p )*D( q ) / AAQQ )
     +                                  / AAPP
                              ELSE
                                 CALL DCOPY( M, A( 1, q ), 1, WORK, 1 )
                                 CALL DLASCL( 'G', 0, 0, AAQQ, D( q ),
     +                                        M, 1, WORK, LDA, IERR )
                                 AAPQ = DDOT( M, WORK, 1, A( 1, p ),
     +                                  1 )*D( p ) / AAPP
                              END IF
                           END IF
*
                           MXAAPQ = DMAX1( MXAAPQ, DABS( AAPQ ) )
*
*        TO rotate or NOT to rotate, THAT is the question ...
*
                           IF( DABS( AAPQ ).GT.TOL ) THEN
                              NOTROT = 0
*           ROTATED  = ROTATED + 1
                              PSKIPPED = 0
                              ISWROT = ISWROT + 1
*
                              IF( ROTOK ) THEN
*
                                 AQOAP = AAQQ / AAPP
                                 APOAQ = AAPP / AAQQ
                                 THETA = -HALF*DABS( AQOAP-APOAQ )/AAPQ
                                 IF( AAQQ.GT.AAPP0 )THETA = -THETA
*
                                 IF( DABS( THETA ).GT.BIGTHETA ) THEN
                                    T = HALF / THETA
                                    FASTR( 3 ) = T*D( p ) / D( q )
                                    FASTR( 4 ) = -T*D( q ) / D( p )
                                    CALL DROTM( M, A( 1, p ), 1,
     +                                          A( 1, q ), 1, FASTR )
                                    IF( RSVEC )CALL DROTM( MVL,
     +                                              V( 1, p ), 1,
     +                                              V( 1, q ), 1,
     +                                              FASTR )
                                    SVA( q ) = AAQQ*DSQRT( DMAX1( ZERO,
     +                                         ONE+T*APOAQ*AAPQ ) )
                                    AAPP = AAPP*DSQRT( DMAX1( ZERO,
     +                                     ONE-T*AQOAP*AAPQ ) )
                                    MXSINJ = DMAX1( MXSINJ, DABS( T ) )
                                 ELSE
*
*                 .. choose correct signum for THETA and rotate
*
                                    THSIGN = -DSIGN( ONE, AAPQ )
                                    IF( AAQQ.GT.AAPP0 )THSIGN = -THSIGN
                                    T = ONE / ( THETA+THSIGN*
     +                                  DSQRT( ONE+THETA*THETA ) )
                                    CS = DSQRT( ONE / ( ONE+T*T ) )
                                    SN = T*CS
                                    MXSINJ = DMAX1( MXSINJ, DABS( SN ) )
                                    SVA( q ) = AAQQ*DSQRT( DMAX1( ZERO,
     +                                         ONE+T*APOAQ*AAPQ ) )
                                    AAPP = AAPP*DSQRT( DMAX1( ZERO, 
     +                                     ONE-T*AQOAP*AAPQ ) )
*
                                    APOAQ = D( p ) / D( q )
                                    AQOAP = D( q ) / D( p )
                                    IF( D( p ).GE.ONE ) THEN
*
                                       IF( D( q ).GE.ONE ) THEN
                                          FASTR( 3 ) = T*APOAQ
                                          FASTR( 4 ) = -T*AQOAP
                                          D( p ) = D( p )*CS
                                          D( q ) = D( q )*CS
                                          CALL DROTM( M, A( 1, p ), 1,
     +                                                A( 1, q ), 1,
     +                                                FASTR )
                                          IF( RSVEC )CALL DROTM( MVL,
     +                                        V( 1, p ), 1, V( 1, q ),
     +                                        1, FASTR )
                                       ELSE
                                          CALL DAXPY( M, -T*AQOAP,
     +                                                A( 1, q ), 1,
     +                                                A( 1, p ), 1 )
                                          CALL DAXPY( M, CS*SN*APOAQ,
     +                                                A( 1, p ), 1,
     +                                                A( 1, q ), 1 )
                                          IF( RSVEC ) THEN
                                             CALL DAXPY( MVL, -T*AQOAP,
     +                                                   V( 1, q ), 1,
     +                                                   V( 1, p ), 1 )
                                             CALL DAXPY( MVL,
     +                                                   CS*SN*APOAQ,
     +                                                   V( 1, p ), 1,
     +                                                   V( 1, q ), 1 )
                                          END IF
                                          D( p ) = D( p )*CS
                                          D( q ) = D( q ) / CS
                                       END IF
                                    ELSE
                                       IF( D( q ).GE.ONE ) THEN
                                          CALL DAXPY( M, T*APOAQ,
     +                                                A( 1, p ), 1,
     +                                                A( 1, q ), 1 )
                                          CALL DAXPY( M, -CS*SN*AQOAP,
     +                                                A( 1, q ), 1,
     +                                                A( 1, p ), 1 )
                                          IF( RSVEC ) THEN
                                             CALL DAXPY( MVL, T*APOAQ,
     +                                                   V( 1, p ), 1,
     +                                                   V( 1, q ), 1 )
                                             CALL DAXPY( MVL,
     +                                                   -CS*SN*AQOAP,
     +                                                   V( 1, q ), 1,
     +                                                   V( 1, p ), 1 )
                                          END IF
                                          D( p ) = D( p ) / CS
                                          D( q ) = D( q )*CS
                                       ELSE
                                          IF( D( p ).GE.D( q ) ) THEN
                                             CALL DAXPY( M, -T*AQOAP,
     +                                                   A( 1, q ), 1,
     +                                                   A( 1, p ), 1 )
                                             CALL DAXPY( M, CS*SN*APOAQ,
     +                                                   A( 1, p ), 1,
     +                                                   A( 1, q ), 1 )
                                             D( p ) = D( p )*CS
                                             D( q ) = D( q ) / CS
                                             IF( RSVEC ) THEN
                                                CALL DAXPY( MVL,
     +                                               -T*AQOAP,
     +                                               V( 1, q ), 1,
     +                                               V( 1, p ), 1 )
                                                CALL DAXPY( MVL,
     +                                               CS*SN*APOAQ,
     +                                               V( 1, p ), 1,
     +                                               V( 1, q ), 1 )
                                             END IF
                                          ELSE
                                             CALL DAXPY( M, T*APOAQ,
     +                                                   A( 1, p ), 1,
     +                                                   A( 1, q ), 1 )
                                             CALL DAXPY( M,
     +                                                   -CS*SN*AQOAP,
     +                                                   A( 1, q ), 1,
     +                                                   A( 1, p ), 1 )
                                             D( p ) = D( p ) / CS
                                             D( q ) = D( q )*CS
                                             IF( RSVEC ) THEN
                                                CALL DAXPY( MVL,
     +                                               T*APOAQ, V( 1, p ),
     +                                               1, V( 1, q ), 1 )
                                                CALL DAXPY( MVL,
     +                                               -CS*SN*AQOAP,
     +                                               V( 1, q ), 1,
     +                                               V( 1, p ), 1 )
                                             END IF
                                          END IF
                                       END IF
                                    END IF
                                 END IF
*
                              ELSE
                                 IF( AAPP.GT.AAQQ ) THEN
                                    CALL DCOPY( M, A( 1, p ), 1, WORK,
     +                                          1 )
                                    CALL DLASCL( 'G', 0, 0, AAPP, ONE,
     +                                           M, 1, WORK, LDA, IERR )
                                    CALL DLASCL( 'G', 0, 0, AAQQ, ONE,
     +                                           M, 1, A( 1, q ), LDA,
     +                                           IERR )
                                    TEMP1 = -AAPQ*D( p ) / D( q )
                                    CALL DAXPY( M, TEMP1, WORK, 1,
     +                                          A( 1, q ), 1 )
                                    CALL DLASCL( 'G', 0, 0, ONE, AAQQ,
     +                                           M, 1, A( 1, q ), LDA,
     +                                           IERR )
                                    SVA( q ) = AAQQ*DSQRT( DMAX1( ZERO,
     +                                         ONE-AAPQ*AAPQ ) )
                                    MXSINJ = DMAX1( MXSINJ, SFMIN )
                                 ELSE
                                    CALL DCOPY( M, A( 1, q ), 1, WORK,
     +                                          1 )
                                    CALL DLASCL( 'G', 0, 0, AAQQ, ONE,
     +                                           M, 1, WORK, LDA, IERR )
                                    CALL DLASCL( 'G', 0, 0, AAPP, ONE,
     +                                           M, 1, A( 1, p ), LDA,
     +                                           IERR )
                                    TEMP1 = -AAPQ*D( q ) / D( p )
                                    CALL DAXPY( M, TEMP1, WORK, 1,
     +                                          A( 1, p ), 1 )
                                    CALL DLASCL( 'G', 0, 0, ONE, AAPP,
     +                                           M, 1, A( 1, p ), LDA,
     +                                           IERR )
                                    SVA( p ) = AAPP*DSQRT( DMAX1( ZERO,
     +                                         ONE-AAPQ*AAPQ ) )
                                    MXSINJ = DMAX1( MXSINJ, SFMIN )
                                 END IF
                              END IF
*           END IF ROTOK THEN ... ELSE
*
*           In the case of cancellation in updating SVA(q)
*           .. recompute SVA(q)
                              IF( ( SVA( q ) / AAQQ )**2.LE.ROOTEPS )
     +                            THEN
                                 IF( ( AAQQ.LT.ROOTBIG ) .AND.
     +                               ( AAQQ.GT.ROOTSFMIN ) ) THEN
                                    SVA( q ) = DNRM2( M, A( 1, q ), 1 )*
     +                                         D( q )
                                 ELSE
                                    T = ZERO
                                    AAQQ = ONE
                                    CALL DLASSQ( M, A( 1, q ), 1, T,
     +                                           AAQQ )
                                    SVA( q ) = T*DSQRT( AAQQ )*D( q )
                                 END IF
                              END IF
                              IF( ( AAPP / AAPP0 )**2.LE.ROOTEPS ) THEN
                                 IF( ( AAPP.LT.ROOTBIG ) .AND.
     +                               ( AAPP.GT.ROOTSFMIN ) ) THEN
                                    AAPP = DNRM2( M, A( 1, p ), 1 )*
     +                                     D( p )
                                 ELSE
                                    T = ZERO
                                    AAPP = ONE
                                    CALL DLASSQ( M, A( 1, p ), 1, T,
     +                                           AAPP )
                                    AAPP = T*DSQRT( AAPP )*D( p )
                                 END IF
                                 SVA( p ) = AAPP
                              END IF
*              end of OK rotation
                           ELSE
                              NOTROT = NOTROT + 1
                              PSKIPPED = PSKIPPED + 1
                              IJBLSK = IJBLSK + 1
                           END IF
                        ELSE
                           NOTROT = NOTROT + 1
                           PSKIPPED = PSKIPPED + 1
                           IJBLSK = IJBLSK + 1
                        END IF
*
                        IF( ( i.LE.SWBAND ) .AND. ( IJBLSK.GE.BLSKIP ) )
     +                      THEN
                           SVA( p ) = AAPP
                           NOTROT = 0
                           GO TO 2011
                        END IF
                        IF( ( i.LE.SWBAND ) .AND.
     +                      ( PSKIPPED.GT.ROWSKIP ) ) THEN
                           AAPP = -AAPP
                           NOTROT = 0
                           GO TO 2203
                        END IF
*
 2200                CONTINUE
*        end of the q-loop
 2203                CONTINUE
*
                     SVA( p ) = AAPP
*
                  ELSE
                     IF( AAPP.EQ.ZERO )NOTROT = NOTROT +
     +                   MIN0( jgl+KBL-1, N ) - jgl + 1
                     IF( AAPP.LT.ZERO )NOTROT = 0
                  END IF

 2100          CONTINUE
*     end of the p-loop
 2010       CONTINUE
*     end of the jbc-loop
 2011       CONTINUE
*2011 bailed out of the jbc-loop
            DO 2012 p = igl, MIN0( igl+KBL-1, N )
               SVA( p ) = DABS( SVA( p ) )
 2012       CONTINUE
*
 2000    CONTINUE
*2000 :: end of the ibr-loop
*
*     .. update SVA(N)
         IF( ( SVA( N ).LT.ROOTBIG ) .AND. ( SVA( N ).GT.ROOTSFMIN ) )
     +       THEN
            SVA( N ) = DNRM2( M, A( 1, N ), 1 )*D( N )
         ELSE
            T = ZERO
            AAPP = ONE
            CALL DLASSQ( M, A( 1, N ), 1, T, AAPP )
            SVA( N ) = T*DSQRT( AAPP )*D( N )
         END IF
*
*     Additional steering devices
*
         IF( ( i.LT.SWBAND ) .AND. ( ( MXAAPQ.LE.ROOTTOL ) .OR.
     +       ( ISWROT.LE.N ) ) )SWBAND = i
*
         IF( ( i.GT.SWBAND+1 ) .AND. ( MXAAPQ.LT.DBLE( N )*TOL ) .AND.
     +       ( DBLE( N )*MXAAPQ*MXSINJ.LT.TOL ) ) THEN
            GO TO 1994
         END IF
*
         IF( NOTROT.GE.EMPTSW )GO TO 1994

 1993 CONTINUE
*     end i=1:NSWEEP loop
* #:) Reaching this point means that the procedure has comleted the given
*     number of iterations.
      INFO = NSWEEP - 1
      GO TO 1995
 1994 CONTINUE
* #:) Reaching this point means that during the i-th sweep all pivots were
*     below the given tolerance, causing early exit.
*
      INFO = 0
* #:) INFO = 0 confirms successful iterations.
 1995 CONTINUE
*
*     Sort the vector D.
      DO 5991 p = 1, N - 1
         q = IDAMAX( N-p+1, SVA( p ), 1 ) + p - 1
         IF( p.NE.q ) THEN
            TEMP1 = SVA( p )
            SVA( p ) = SVA( q )
            SVA( q ) = TEMP1
            TEMP1 = D( p )
            D( p ) = D( q )
            D( q ) = TEMP1
            CALL DSWAP( M, A( 1, p ), 1, A( 1, q ), 1 )
            IF( RSVEC )CALL DSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 )
         END IF
 5991 CONTINUE
*
      RETURN
*     ..
*     .. END OF DGSVJ0
*     ..
      END
