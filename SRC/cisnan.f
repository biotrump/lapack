      FUNCTION CISNAN( CIN )
      LOGICAL CISNAN
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      COMPLEX          CIN
*     ..
*
*  Purpose
*  =======
*
*  SISNAN returns .TRUE. if its argument is NaN, and .FALSE.
*  otherwise.  To be replaced by the Fortran 2003 intrinsic in the
*  future.
*
*  Arguments
*  =========
*
*  CIN      (input) COMPLEX
*          Input to test for NaN.
*
*  =====================================================================
*
*  .. External Functions ..
      LOGICAL SLAISNAN
      EXTERNAL SLAISNAN
*  ..
*  .. Intrinsic Functions ..
      INTRINSIC REAL, AIMAG
*  ..
*  .. Executable Statements ..
      CISNAN = SLAISNAN( REAL(CIN), REAL(CIN) )
     $   .OR. SLAISNAN( AIMAG(CIN), AIMAG(CIN) )
      END FUNCTION
