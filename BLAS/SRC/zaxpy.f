      SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
c
c     constant times a vector plus a vector.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c

C     .. Scalar Arguments ..
      DOUBLE COMPLEX ZA
      INTEGER INCX,INCY,N
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX ZX(*),ZY(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,IX,IY
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DCABS1
      EXTERNAL DCABS1
C     ..
      IF (N.LE.0) RETURN
      IF (DCABS1(ZA).EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          ZY(IY) = ZY(IY) + ZA*ZX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
c
c        code for both increments equal to 1
c
   20 DO 30 I = 1,N
          ZY(I) = ZY(I) + ZA*ZX(I)
   30 CONTINUE
      RETURN
      END
