      DOUBLE PRECISION FUNCTION DCABS1(Z)
C     .. Scalar Arguments ..
      DOUBLE COMPLEX Z
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,DIMAG
C     ..
      DCABS1 = ABS(DBLE(Z)) + ABS(DIMAG(Z))
      RETURN
      END
