*> \brief \b ZPTSV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZPTSV + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zptsv.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zptsv.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zptsv.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZPTSV( N, NRHS, D, E, B, LDB, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   D( * )
*       COMPLEX*16         B( LDB, * ), E( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZPTSV computes the solution to a complex system of linear equations
*> A*X = B, where A is an N-by-N Hermitian positive definite tridiagonal
*> matrix, and X and B are N-by-NRHS matrices.
*>
*> A is factored as A = L*D*L**H, and the factored form of A is then
*> used to solve the system of equations.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in,out] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (N)
*>          On entry, the n diagonal elements of the tridiagonal matrix
*>          A.  On exit, the n diagonal elements of the diagonal matrix
*>          D from the factorization A = L*D*L**H.
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is COMPLEX*16 array, dimension (N-1)
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal
*>          matrix A.  On exit, the (n-1) subdiagonal elements of the
*>          unit bidiagonal factor L from the L*D*L**H factorization of
*>          A.  E can also be regarded as the superdiagonal of the unit
*>          bidiagonal factor U from the U**H*D*U factorization of A.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)
*>          On entry, the N-by-NRHS right hand side matrix B.
*>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, the leading minor of order i is not
*>                positive definite, and the solution has not been
*>                computed.  The factorization has not been completed
*>                unless i = N.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup complex16OTHERcomputational
*
*  =====================================================================
      SUBROUTINE ZPTSV( N, NRHS, D, E, B, LDB, INFO )
*
*  -- LAPACK computational routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * )
      COMPLEX*16         B( LDB, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZPTTRF, ZPTTRS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPTSV ', -INFO )
         RETURN
      END IF
*
*     Compute the L*D*L**H (or U**H*D*U) factorization of A.
*
      CALL ZPTTRF( N, D, E, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         CALL ZPTTRS( 'Lower', N, NRHS, D, E, B, LDB, INFO )
      END IF
      RETURN
*
*     End of ZPTSV
*
      END
