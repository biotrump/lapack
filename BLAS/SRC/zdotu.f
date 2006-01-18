      DOUBLE COMPLEX FUNCTION ZDOTU(N,ZX,INCX,ZY,INCY)
c
c     forms the dot product of two vectors.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
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
      ZTEMP = (0.0d0,0.0d0)
      ZDOTU = (0.0d0,0.0d0)
      IF (N.LE.0) RETURN
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
          ZTEMP = ZTEMP + ZX(IX)*ZY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      ZDOTU = ZTEMP
      RETURN
c
c        code for both increments equal to 1
c
   20 DO 30 I = 1,N
          ZTEMP = ZTEMP + ZX(I)*ZY(I)
   30 CONTINUE
      ZDOTU = ZTEMP
      RETURN
      END
