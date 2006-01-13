      REAL FUNCTION SCABS1( Z )
*     .. Scalar Arguments ..
      COMPLEX                           Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC                         ABS, REAL, AIMAG
*     ..
      SCABS1 = ABS(REAL(Z)) + ABS(AIMAG(Z))
      RETURN
      END
