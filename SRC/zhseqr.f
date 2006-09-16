      SUBROUTINE ZHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ,
     $                   WORK, LWORK, INFO )
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
      CHARACTER          COMPZ, JOB
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*     Purpose
*     =======
*
*     ZHSEQR computes the eigenvalues of a Hessenberg matrix H
*     and, optionally, the matrices T and Z from the Schur decomposition
*     H = Z T Z**H, where T is an upper triangular matrix (the
*     Schur form), and Z is the unitary matrix of Schur vectors.
*
*     Optionally Z may be postmultiplied into an input unitary
*     matrix Q so that this routine can give the Schur factorization
*     of a matrix A which has been reduced to the Hessenberg form H
*     by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.
*
*     Arguments
*     =========
*
*     JOB   (input) CHARACTER*1
*           = 'E':  compute eigenvalues only;
*           = 'S':  compute eigenvalues and the Schur form T.
*
*     COMPZ (input) CHARACTER*1
*           = 'N':  no Schur vectors are computed;
*           = 'I':  Z is initialized to the unit matrix and the matrix Z
*                   of Schur vectors of H is returned;
*           = 'V':  Z must contain an unitary matrix Q on entry, and
*                   the product Q*Z is returned.
*
*     N     (input) INTEGER
*           The order of the matrix H.  N .GE. 0.
*
*     ILO   (input) INTEGER
*     IHI   (input) INTEGER
*           It is assumed that H is already upper triangular in rows
*           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
*           set by a previous call to ZGEBAL, and then passed to ZGEHRD
*           when the matrix output by ZGEBAL is reduced to Hessenberg
*           form. Otherwise ILO and IHI should be set to 1 and N
*           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
*           If N = 0, then ILO = 1 and IHI = 0.
*
*     H     (input/output) COMPLEX*16 array, dimension (LDH,N)
*           On entry, the upper Hessenberg matrix H.
*           On exit, if INFO = 0 and JOB = 'S', H contains the upper
*           triangular matrix T from the Schur decomposition (the
*           Schur form). If INFO = 0 and JOB = 'E', the contents of
*           H are unspecified on exit.  (The output value of H when
*           INFO.GT.0 is given under the description of INFO below.)
*
*           Unlike earlier versions of ZHSEQR, this subroutine may
*           explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1
*           or j = IHI+1, IHI+2, ... N.
*
*     LDH   (input) INTEGER
*           The leading dimension of the array H. LDH .GE. max(1,N).
*
*     W        (output) COMPLEX*16 array, dimension (N)
*           The computed eigenvalues. If JOB = 'S', the eigenvalues are
*           stored in the same order as on the diagonal of the Schur
*           form returned in H, with W(i) = H(i,i).
*
*     Z     (input/output) COMPLEX*16 array, dimension (LDZ,N)
*           If COMPZ = 'N', Z is not referenced.
*           If COMPZ = 'I', on entry Z need not be set and on exit,
*           if INFO = 0, Z contains the unitary matrix Z of the Schur
*           vectors of H.  If COMPZ = 'V', on entry Z must contain an
*           N-by-N matrix Q, which is assumed to be equal to the unit
*           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,
*           if INFO = 0, Z contains Q*Z.
*           Normally Q is the unitary matrix generated by ZUNGHR
*           after the call to ZGEHRD which formed the Hessenberg matrix
*           H. (The output value of Z when INFO.GT.0 is given under
*           the description of INFO below.)
*
*     LDZ   (input) INTEGER
*           The leading dimension of the array Z.  if COMPZ = 'I' or
*           COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1.
*
*     WORK  (workspace/output) COMPLEX*16 array, dimension (LWORK)
*           On exit, if INFO = 0, WORK(1) returns an estimate of
*           the optimal value for LWORK.
*
*     LWORK (input) INTEGER
*           The dimension of the array WORK.  LWORK .GE. max(1,N)
*           is sufficient, but LWORK typically as large as 6*N may
*           be required for optimal performance.  A workspace query
*           to determine the optimal workspace size is recommended.
*
*           If LWORK = -1, then ZHSEQR does a workspace query.
*           In this case, ZHSEQR checks the input parameters and
*           estimates the optimal workspace size for the given
*           values of N, ILO and IHI.  The estimate is returned
*           in WORK(1).  No error message related to LWORK is
*           issued by XERBLA.  Neither H nor Z are accessed.
*
*
*     INFO  (output) INTEGER
*             =  0:  successful exit
*           .LT. 0:  if INFO = -i, the i-th argument had an illegal
*                    value
*           .GT. 0:  if INFO = i, ZHSEQR failed to compute all of
*                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
*                and WI contain those eigenvalues which have been
*                successfully computed.  (Failures are rare.)
*
*                If INFO .GT. 0 and JOB = 'E', then on exit, the
*                remaining unconverged eigenvalues are the eigen-
*                values of the upper Hessenberg matrix rows and
*                columns ILO through INFO of the final, output
*                value of H.
*
*                If INFO .GT. 0 and JOB   = 'S', then on exit
*
*           (*)  (initial value of H)*U  = U*(final value of H)
*
*                where U is a unitary matrix.  The final
*                value of  H is upper Hessenberg and triangular in
*                rows and columns INFO+1 through IHI.
*
*                If INFO .GT. 0 and COMPZ = 'V', then on exit
*
*                  (final value of Z)  =  (initial value of Z)*U
*
*                where U is the unitary matrix in (*) (regard-
*                less of the value of JOB.)
*
*                If INFO .GT. 0 and COMPZ = 'I', then on exit
*                      (final value of Z)  = U
*                where U is the unitary matrix in (*) (regard-
*                less of the value of JOB.)
*
*                If INFO .GT. 0 and COMPZ = 'N', then Z is not
*                accessed.
*
*     ================================================================
*             Default values supplied by
*             ILAENV(ISPEC,'ZHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).
*             It is suggested that these defaults be adjusted in order
*             to attain best performance in each particular
*             computational environment.
*
*            ISPEC=1:  The ZLAHQR vs ZLAQR0 crossover point.
*                      Default: 75. (Must be at least 11.)
*
*            ISPEC=2:  Recommended deflation window size.
*                      This depends on ILO, IHI and NS.  NS is the
*                      number of simultaneous shifts returned
*                      by ILAENV(ISPEC=4).  (See ISPEC=4 below.)
*                      The default for (IHI-ILO+1).LE.500 is NS.
*                      The default for (IHI-ILO+1).GT.500 is 3*NS/2.
*
*            ISPEC=3:  Nibble crossover point. (See ILAENV for
*                      details.)  Default: 14% of deflation window
*                      size.
*
*            ISPEC=4:  Number of simultaneous shifts, NS, in
*                      a multi-shift QR iteration.
*
*                      If IHI-ILO+1 is ...
*
*                      greater than      ...but less    ... the
*                      or equal to ...      than        default is
*
*                           1               30          NS -   2(+)
*                          30               60          NS -   4(+)
*                          60              150          NS =  10(+)
*                         150              590          NS =  **
*                         590             3000          NS =  64
*                        3000             6000          NS = 128
*                        6000             infinity      NS = 256
*
*                  (+)  By default some or all matrices of this order 
*                       are passed to the implicit double shift routine
*                       ZLAHQR and NS is ignored.  See ISPEC=1 above 
*                       and comments in IPARM for details.
*
*                       The asterisks (**) indicate an ad-hoc
*                       function of N increasing from 10 to 64.
*
*            ISPEC=5:  Select structured matrix multiply.
*                      (See ILAENV for details.) Default: 3.
*
*     ================================================================
*     References:
*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
*       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
*       929--947, 2002.
*
*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
*       of Matrix Analysis, volume 23, pages 948--973, 2002.
*
*     ================================================================
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    ZLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
*
*     ==== NL allocates some local workspace to help small matrices
*     .    through a rare ZLAHQR failure.  NL .GT. NTINY = 11 is
*     .    required and NL .LE. NMIN = ILAENV(ISPEC=1,...) is recom-
*     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
*     .    allows up to six simultaneous shifts and a 16-by-16
*     .    deflation window.  ====
*
      INTEGER            NTINY
      PARAMETER          ( NTINY = 11 )
      INTEGER            NL
      PARAMETER          ( NL = 49 )
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0d0, 0.0d0 ),
     $                   ONE = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO
      PARAMETER          ( RZERO = 0.0d0 )
*     ..
*     .. Local Arrays ..
      COMPLEX*16         HL( NL, NL ), WORKL( NL )
*     ..
*     .. Local Scalars ..
      INTEGER            KBOT, NMIN
      LOGICAL            INITZ, LQUERY, WANTT, WANTZ
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      LOGICAL            LSAME
      EXTERNAL           ILAENV, LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZLACPY, ZLAHQR, ZLAQR0, ZLASET
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     ==== Decode and check the input parameters. ====
*
      WANTT = LSAME( JOB, 'S' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
      WORK( 1 ) = DCMPLX( DBLE( MAX( 1, N ) ), RZERO )
      LQUERY = LWORK.EQ.-1
*
      INFO = 0
      IF( .NOT.LSAME( JOB, 'E' ) .AND. .NOT.WANTT ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.WANTZ ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -5
      ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.NE.0 ) THEN
*
*        ==== Quick return in case of invalid argument. ====
*
         CALL XERBLA( 'ZHSEQR', -INFO )
         RETURN
*
      ELSE IF( N.EQ.0 ) THEN
*
*        ==== Quick return in case N = 0; nothing to do. ====
*
         RETURN
*
      ELSE IF( LQUERY ) THEN
*
*        ==== Quick return in case of a workspace query ====
*
         CALL ZLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z,
     $                LDZ, WORK, LWORK, INFO )
*        ==== Ensure reported workspace size is backward-compatible with
*        .    previous LAPACK versions. ====
         WORK( 1 ) = DCMPLX( MAX( DBLE( WORK( 1 ) ), DBLE( MAX( 1,
     $               N ) ) ), RZERO )
         RETURN
*
      ELSE
*
*        ==== copy eigenvalues isolated by ZGEBAL ====
*
         IF( ILO.GT.1 )
     $      CALL ZCOPY( ILO-1, H, LDH+1, W, 1 )
         IF( IHI.LT.N )
     $      CALL ZCOPY( N-IHI, H( IHI+1, IHI+1 ), LDH+1, W( IHI+1 ), 1 )
*
*        ==== Initialize Z, if requested ====
*
         IF( INITZ )
     $      CALL ZLASET( 'A', N, N, ZERO, ONE, Z, LDZ )
*
*        ==== Quick return if possible ====
*
         IF( ILO.EQ.IHI ) THEN
            W( ILO ) = H( ILO, ILO )
            RETURN
         END IF
*
*        ==== ZLAHQR/ZLAQR0 crossover point ====
*
         NMIN = ILAENV( 1, 'ZHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, ILO,
     $          IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
*
*        ==== ZLAQR0 for big matrices; ZLAHQR for small ones ====
*
         IF( N.GT.NMIN ) THEN
            CALL ZLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI,
     $                   Z, LDZ, WORK, LWORK, INFO )
         ELSE
*
*           ==== Small matrix ====
*
            CALL ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI,
     $                   Z, LDZ, INFO )
*
            IF( INFO.GT.0 ) THEN
*
*              ==== A rare ZLAHQR failure!  ZLAQR0 sometimes succeeds
*              .    when ZLAHQR fails. ====
*
               KBOT = INFO
*
               IF( N.GE.NL ) THEN
*
*                 ==== Larger matrices have enough subdiagonal scratch
*                 .    space to call ZLAQR0 directly. ====
*
                  CALL ZLAQR0( WANTT, WANTZ, N, ILO, KBOT, H, LDH, W,
     $                         ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
*
               ELSE
*
*                 ==== Tiny matrices don't have enough subdiagonal
*                 .    scratch space to benefit from ZLAQR0.  Hence,
*                 .    tiny matrices must be copied into a larger
*                 .    array before calling ZLAQR0. ====
*
                  CALL ZLACPY( 'A', N, N, H, LDH, HL, NL )
                  HL( N+1, N ) = ZERO
                  CALL ZLASET( 'A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ),
     $                         NL )
                  CALL ZLAQR0( WANTT, WANTZ, NL, ILO, KBOT, HL, NL, W,
     $                         ILO, IHI, Z, LDZ, WORKL, NL, INFO )
                  IF( WANTT .OR. INFO.NE.0 )
     $               CALL ZLACPY( 'A', N, N, HL, NL, H, LDH )
               END IF
            END IF
         END IF
*
*        ==== Clear out the trash, if necessary. ====
*
         IF( ( WANTT .OR. INFO.NE.0 ) .AND. N.GT.2 )
     $      CALL ZLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )
*
*        ==== Ensure reported workspace size is backward-compatible with
*        .    previous LAPACK versions. ====
*
         WORK( 1 ) = DCMPLX( MAX( DBLE( MAX( 1, N ) ),
     $               DBLE( WORK( 1 ) ) ), RZERO )
      END IF
*
*     ==== End of ZHSEQR ====
*
      END
