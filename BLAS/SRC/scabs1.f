      REAL FUNCTION SCABS1(Z)
C     .. Scalar Arguments ..
      COMPLEX Z
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,AIMAG,REAL
C     ..
      SCABS1 = ABS(REAL(Z)) + ABS(AIMAG(Z))
      RETURN
      END
