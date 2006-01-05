      DOUBLE PRECISION FUNCTION DCABS1( Z )
*     .. Scalar Arguments ..
      DOUBLE COMPLEX                    Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC                         ABS, DBLE, DIMAG
*     ..
      DCABS1 = ABS( DBLE( Z ) ) + ABS( DIMAG( Z ) )
      RETURN
      END
