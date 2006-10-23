      FUNCTION ZISNAN( ZIN )
      LOGICAL ZISNAN
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      COMPLEX*16       ZIN
*     ..
*
*  Purpose
*  =======
*
*  DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
*  otherwise.  To be replaced by the Fortran 2003 intrinsic in the
*  future.
*
*  Arguments
*  =========
*
*  ZIN      (input) COMPLEX*16
*          Input to test for NaN.
*
*  =====================================================================
*
*  .. External Functions ..
      LOGICAL DLAISNAN
      EXTERNAL DLAISNAN
*  ..
*  .. Intrinsic Functions ..
      INTRINSIC DBLE, DIMAG
*  ..
*  .. Executable Statements ..
      ZISNAN = DLAISNAN( DBLE(ZIN), DBLE(ZIN) )
     $   .OR. DLAISNAN( DIMAG(ZIN), DIMAG(ZIN) )
      END FUNCTION
