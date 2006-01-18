      SUBROUTINE ZSWAP(N,ZX,INCX,ZY,INCY)
c
c     interchanges two vectors.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c
C     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX ZX(*),ZY(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX ZTEMP
      INTEGER I,IX,IY
C     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          ZTEMP = ZX(IX)
          ZX(IX) = ZY(IY)
          ZY(IY) = ZTEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
c
c       code for both increments equal to 1
   20 DO 30 I = 1,N
          ZTEMP = ZX(I)
          ZX(I) = ZY(I)
          ZY(I) = ZTEMP
   30 CONTINUE
      RETURN
      END
