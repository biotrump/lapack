      SUBROUTINE ZLARRV( N, D, L, ISPLIT, M, W, IBLOCK, GERSCH, TOL, Z,
     $                   LDZ, ISUPPZ, WORK, IWORK, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDZ, M, N
      DOUBLE PRECISION   TOL
*     ..
*     .. Array Arguments ..
      INTEGER            IBLOCK( * ), ISPLIT( * ), ISUPPZ( * ),
     $                   IWORK( * )
      DOUBLE PRECISION   D( * ), GERSCH( * ), L( * ), W( * ), WORK( * )
      COMPLEX*16         Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARRV computes the eigenvectors of the tridiagonal matrix
*  T = L D L^T given L, D and the eigenvalues of L D L^T.
*  The input eigenvalues should have high relative accuracy with
*  respect to the entries of L and D. The desired accuracy of the
*  output can be specified by the input parameter TOL.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the n diagonal elements of the diagonal matrix D.
*          On exit, D may be overwritten.
*
*  L       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the unit
*          bidiagonal matrix L in elements 1 to N-1 of L.
*          On exit, L is overwritten.
*
*  ISPLIT  (input) INTEGER array, dimension (N)
*          The splitting points, at which T breaks up into submatrices.
*          The first submatrix consists of rows/columns 1 to
*          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
*          through ISPLIT( 2 ), etc.
*
*  M       (input) INTEGER
*          The total number of eigenvalues found.  0 <= M <= N.
*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*
*  W       (input) DOUBLE PRECISION array, dimension (N)
*          The first M elements of W contain the eigenvalues for
*          which eigenvectors are to be computed.  The eigenvalues
*          should be grouped by split-off block and ordered from
*          smallest to largest within the block ( The output array
*          W from DLARRE is expected here ).
*          Errors in W must be bounded by TOL (see above).
*
*  IBLOCK  (input) INTEGER array, dimension (N)
*          The submatrix indices associated with the corresponding
*          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to
*          the first submatrix from the top, =2 if W(i) belongs to
*          the second submatrix, etc.
*
*  GERSCH   (input) DOUBLE PRECISION array, dimension (2*N)
*           The n Gerschgorin intervals. These are used to restrict
*           the initial search for R, when R is input as 0.
*
*  TOL     (input) DOUBLE PRECISION
*          The absolute error tolerance for the
*          eigenvalues/eigenvectors.
*          Errors in the input eigenvalues must be bounded by TOL.
*          The eigenvectors output have residual norms
*          bounded by TOL, and the dot products between different
*          eigenvectors are bounded by TOL. TOL must be at least
*          N*EPS*|T|, where EPS is the machine precision and |T| is
*          the 1-norm of the tridiagonal matrix.
*
*  Z       (output) COMPLEX*16 array, dimension (LDZ, max(1,M) )
*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix T
*          corresponding to the selected eigenvalues, with the i-th
*          column of Z holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*          Note: the user must ensure that at least max(1,M) columns are
*          supplied in the array Z; if RANGE = 'V', the exact value of M
*          is not known in advance and an upper bound must be used.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) )
*          The support of the eigenvectors in Z, i.e., the indices
*          indicating the nonzero elements in Z. The i-th eigenvector
*          is nonzero only in elements ISUPPZ( 2*i-1 ) through
*          ISUPPZ( 2*i ).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (13*N)
*
*  IWORK   (workspace) INTEGER array, dimension (6*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = 1, internal error in DLARRB (this cannot be
*                             triggered since DLARRB does not set INFO)
*                if INFO = 2, internal error in ZSTEIN
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Inderjit Dhillon, IBM Almaden, USA
*     Osni Marques, LBNL/NERSC, USA
*     Ken Stanley, Computer Science Division, University of
*       California at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MGSSIZ
      PARAMETER          ( MGSSIZ = 20 )
      DOUBLE PRECISION   ZERO, ONE, FOUR
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            MGSCLS
      INTEGER            I, IBEGIN, IEND, IINDC1, IINDC2, IINDR, IINDWK,
     $                   IINFO, IM, IN, INDERR, INDGAP, INDIN1, INDIN2,
     $                   INDLD, INDLLD, INDWRK, ITER, ITMP1, ITMP2, J,
     $                   JBLK, K, KTOT, LSBDPT, MAXITR, NCLUS, NDEPTH,
     $                   NDONE, NEWCLS, NEWFRS, NEWFTT, NEWLST, NEWSIZ,
     $                   NSPLIT, OLDCLS, OLDFST, OLDIEN, OLDLST, OLDNCL,
     $                   P, Q
      DOUBLE PRECISION   EPS, GAP, LAMBDA, MGSTOL, MINGMA, MINRGP,
     $                   NRMINV, RELGAP, RELTOL, RESID, RQCORR, SIGMA,
     $                   TMP1, ZTZ
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DZNRM2
      COMPLEX*16         ZDOTU
      EXTERNAL           DLAMCH, DZNRM2, ZDOTU
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLARRB, DLARRF, ZAXPY, ZDSCAL, ZLAR1V,
     $                   ZLASET, ZSTEIN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, MAX, MIN, SQRT
*     ..
*     .. Local Arrays ..
      INTEGER            TEMP( 1 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INDERR = N + 1
      INDLD = 2*N
      INDLLD = 3*N
      INDGAP = 4*N
      INDIN1 = 5*N + 1
      INDIN2 = 6*N + 1
      INDWRK = 7*N + 1
*
      IINDR = N
      IINDC1 = 2*N
      IINDC2 = 3*N
      IINDWK = 4*N + 1
*
      EPS = DLAMCH( 'Precision' )
*
      DO 10 I = 1, 2*N
         IWORK( I ) = 0
   10 CONTINUE
      DO 20 I = 1, M
         WORK( INDERR+I-1 ) = EPS*ABS( W( I ) )
   20 CONTINUE
      CALL ZLASET( 'Full', N, N, CZERO, CZERO, Z, LDZ )
      MGSTOL = 5.0D0*EPS
*
      NSPLIT = IBLOCK( M )
      IBEGIN = 1
      DO 190 JBLK = 1, NSPLIT
         IEND = ISPLIT( JBLK )
*
*        Find the eigenvectors of the submatrix indexed IBEGIN
*        through IEND.
*
         IF( IBEGIN.EQ.IEND ) THEN
            Z( IBEGIN, IBEGIN ) = ONE
            ISUPPZ( 2*IBEGIN-1 ) = IBEGIN
            ISUPPZ( 2*IBEGIN ) = IBEGIN
            IBEGIN = IEND + 1
            GO TO 190
         END IF
         OLDIEN = IBEGIN - 1
         IN = IEND - OLDIEN
         RELTOL = MIN( 1.0D-2, ONE / DBLE( IN ) )
         IM = IN
         CALL DCOPY( IM, W( IBEGIN ), 1, WORK, 1 )
         DO 30 I = 1, IN - 1
            WORK( INDGAP+I ) = WORK( I+1 ) - WORK( I )
   30    CONTINUE
         WORK( INDGAP+IN ) = MAX( ABS( WORK( IN ) ), EPS )
         NDONE = 0
*
         NDEPTH = 0
         LSBDPT = 1
         NCLUS = 1
         IWORK( IINDC1+1 ) = 1
         IWORK( IINDC1+2 ) = IN
*
*        While( NDONE.LT.IM ) do
*
   40    CONTINUE
         IF( NDONE.LT.IM ) THEN
            OLDNCL = NCLUS
            NCLUS = 0
            LSBDPT = 1 - LSBDPT
            DO 170 I = 1, OLDNCL
               IF( LSBDPT.EQ.0 ) THEN
                  OLDCLS = IINDC1
                  NEWCLS = IINDC2
               ELSE
                  OLDCLS = IINDC2
                  NEWCLS = IINDC1
               END IF
*
*              If NDEPTH > 1, retrieve the relatively robust
*              representation (RRR) and perform limited bisection
*              (if necessary) to get approximate eigenvalues.
*
               J = OLDCLS + 2*I
               OLDFST = IWORK( J-1 )
               OLDLST = IWORK( J )
               IF( NDEPTH.GT.0 ) THEN
                  J = OLDIEN + OLDFST
                  DO 50 K = 1, IN
                     D( IBEGIN+K-1 ) = DBLE( Z( IBEGIN+K-1,
     $                                 OLDIEN+OLDFST ) )
                     L( IBEGIN+K-1 ) = DBLE( Z( IBEGIN+K-1,
     $                                 OLDIEN+OLDFST+1 ) )
   50             CONTINUE
                  SIGMA = L( IEND )
               END IF
               K = IBEGIN
               DO 60 J = 1, IN - 1
                  WORK( INDLD+J ) = D( K )*L( K )
                  WORK( INDLLD+J ) = WORK( INDLD+J )*L( K )
                  K = K + 1
   60          CONTINUE
               IF( NDEPTH.GT.0 ) THEN
                  IINFO = 0
                  CALL DLARRB( IN, D( IBEGIN ), L( IBEGIN ),
     $                         WORK( INDLD+1 ), WORK( INDLLD+1 ),
     $                         OLDFST, OLDLST, SIGMA, RELTOL, WORK,
     $                         WORK( INDGAP+1 ), WORK( INDERR ),
     $                         WORK( INDWRK ), IWORK( IINDWK ), IINFO )
                  IF( IINFO.NE.0 ) THEN
                     INFO = 1
                     RETURN
                  END IF
               END IF
*
*              Classify eigenvalues of the current representation (RRR)
*              as (i) isolated, (ii) loosely clustered or (iii) tightly
*              clustered
*
               NEWFRS = OLDFST
               DO 160 J = OLDFST, OLDLST
                  IF( J.EQ.OLDLST .OR. WORK( INDGAP+J ).GE.RELTOL*
     $                ABS( WORK( J ) ) ) THEN
                     NEWLST = J
                  ELSE
*
*                    continue (to the next loop)
*
                     RELGAP = WORK( INDGAP+J ) / ABS( WORK( J ) )
                     IF( J.EQ.NEWFRS ) THEN
                        MINRGP = RELGAP
                     ELSE
                        MINRGP = MIN( MINRGP, RELGAP )
                     END IF
                     GO TO 160
                  END IF
                  NEWSIZ = NEWLST - NEWFRS + 1
                  MAXITR = 10
                  NEWFTT = OLDIEN + NEWFRS
                  IF( NEWSIZ.GT.1 ) THEN
                     MGSCLS = NEWSIZ.LE.MGSSIZ .AND. MINRGP.GE.MGSTOL
                     IF( .NOT.MGSCLS ) THEN
                        DO 70 K = 1, IN
                           WORK( INDIN1+K-1 ) = DBLE( Z( IBEGIN+K-1,
     $                                          NEWFTT ) )
                           WORK( INDIN2+K-1 ) = DBLE( Z( IBEGIN+K-1,
     $                                          NEWFTT+1 ) )
   70                   CONTINUE
                        CALL DLARRF( IN, D( IBEGIN ), L( IBEGIN ),
     $                               WORK( INDLD+1 ), WORK( INDLLD+1 ),
     $                               NEWFRS, NEWLST, WORK,
     $                               WORK( INDIN1 ), WORK( INDIN2 ),
     $                               WORK( INDWRK ), IWORK( IINDWK ),
     $                               INFO )
                        IF( INFO.EQ.0 ) THEN
                           NCLUS = NCLUS + 1
                           K = NEWCLS + 2*NCLUS
                           IWORK( K-1 ) = NEWFRS
                           IWORK( K ) = NEWLST
                        ELSE
                           INFO = 0
                           IF( MINRGP.GE.MGSTOL ) THEN
                              MGSCLS = .TRUE.
                           ELSE
*
*                             Call ZSTEIN to process this tight cluster.
*                             This happens only if MINRGP <= MGSTOL
*                             and DLARRF returns INFO = 1. The latter
*                             means that a new RRR to "break" the
*                             cluster could not be found.
*
                              WORK( INDWRK ) = D( IBEGIN )
                              DO 80 K = 1, IN - 1
                                 WORK( INDWRK+K ) = D( IBEGIN+K ) +
     $                                              WORK( INDLLD+K )
   80                         CONTINUE
                              DO 90 K = 1, NEWSIZ
                                 IWORK( IINDWK+K-1 ) = 1
   90                         CONTINUE
                              DO 100 K = NEWFRS, NEWLST
                                 ISUPPZ( 2*( IBEGIN+K )-3 ) = 1
                                 ISUPPZ( 2*( IBEGIN+K )-2 ) = IN
  100                         CONTINUE
                              TEMP( 1 ) = IN
                              CALL ZSTEIN( IN, WORK( INDWRK ),
     $                                     WORK( INDLD+1 ), NEWSIZ,
     $                                     WORK( NEWFRS ),
     $                                     IWORK( IINDWK ), TEMP( 1 ),
     $                                     Z( IBEGIN, NEWFTT ), LDZ,
     $                                     WORK( INDWRK+IN ),
     $                                     IWORK( IINDWK+IN ),
     $                                     IWORK( IINDWK+2*IN ), IINFO )
                              IF( IINFO.NE.0 ) THEN
                                 INFO = 2
                                 RETURN
                              END IF
                              NDONE = NDONE + NEWSIZ
                           END IF
                        END IF
                     END IF
                  ELSE
                     MGSCLS = .FALSE.
                  END IF
                  IF( NEWSIZ.EQ.1 .OR. MGSCLS ) THEN
                     KTOT = NEWFTT
                     DO 120 K = NEWFRS, NEWLST
                        ITER = 0
  110                   CONTINUE
                        LAMBDA = WORK( K )
                        CALL ZLAR1V( IN, 1, IN, LAMBDA, D( IBEGIN ),
     $                               L( IBEGIN ), WORK( INDLD+1 ),
     $                               WORK( INDLLD+1 ),
     $                               GERSCH( 2*OLDIEN+1 ),
     $                               Z( IBEGIN, KTOT ), ZTZ, MINGMA,
     $                               IWORK( IINDR+KTOT ),
     $                               ISUPPZ( 2*KTOT-1 ),
     $                               WORK( INDWRK ) )
                        TMP1 = ONE / ZTZ
                        NRMINV = SQRT( TMP1 )
                        RESID = ABS( MINGMA )*NRMINV
                        RQCORR = MINGMA*TMP1
                        IF( K.EQ.IN ) THEN
                           GAP = WORK( INDGAP+K-1 )
                        ELSE IF( K.EQ.1 ) THEN
                           GAP = WORK( INDGAP+K )
                        ELSE
                           GAP = MIN( WORK( INDGAP+K-1 ),
     $                           WORK( INDGAP+K ) )
                        END IF
                        ITER = ITER + 1
                        IF( RESID.GT.TOL*GAP .AND. ABS( RQCORR ).GT.
     $                      FOUR*EPS*ABS( LAMBDA ) ) THEN
                           WORK( K ) = LAMBDA + RQCORR
                           IF( ITER.LT.MAXITR ) THEN
                              GO TO 110
                           END IF
                        END IF
                        IWORK( KTOT ) = 1
                        IF( NEWSIZ.EQ.1 )
     $                     NDONE = NDONE + 1
                        CALL ZDSCAL( IN, NRMINV, Z( IBEGIN, KTOT ), 1 )
                        KTOT = KTOT + 1
  120                CONTINUE
                     IF( NEWSIZ.GT.1 ) THEN
                        ITMP1 = ISUPPZ( 2*NEWFTT-1 )
                        ITMP2 = ISUPPZ( 2*NEWFTT )
                        KTOT = OLDIEN + NEWLST
                        DO 140 P = NEWFTT + 1, KTOT
                           DO 130 Q = NEWFTT, P - 1
                              TMP1 = -ZDOTU( IN, Z( IBEGIN, P ), 1,
     $                               Z( IBEGIN, Q ), 1 )
                              CALL ZAXPY( IN, DCMPLX( TMP1, ZERO ),
     $                                    Z( IBEGIN, Q ), 1,
     $                                    Z( IBEGIN, P ), 1 )
  130                      CONTINUE
                           TMP1 = ONE / DZNRM2( IN, Z( IBEGIN, P ), 1 )
                           CALL ZDSCAL( IN, TMP1, Z( IBEGIN, P ), 1 )
                           ITMP1 = MIN( ITMP1, ISUPPZ( 2*P-1 ) )
                           ITMP2 = MAX( ITMP2, ISUPPZ( 2*P ) )
  140                   CONTINUE
                        DO 150 P = NEWFTT, KTOT
                           ISUPPZ( 2*P-1 ) = ITMP1
                           ISUPPZ( 2*P ) = ITMP2
  150                   CONTINUE
                        NDONE = NDONE + NEWSIZ
                     END IF
                  END IF
                  NEWFRS = J + 1
  160          CONTINUE
  170       CONTINUE
            NDEPTH = NDEPTH + 1
            GO TO 40
         END IF
         J = 2*IBEGIN
         DO 180 I = IBEGIN, IEND
            ISUPPZ( J-1 ) = ISUPPZ( J-1 ) + OLDIEN
            ISUPPZ( J ) = ISUPPZ( J ) + OLDIEN
            J = J + 2
  180    CONTINUE
         IBEGIN = IEND + 1
  190 CONTINUE
*
      RETURN
*
*     End of ZLARRV
*
      END
