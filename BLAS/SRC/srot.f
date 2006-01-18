      SUBROUTINE SROT(N,SX,INCX,SY,INCY,C,S)
c
c     applies a plane rotation.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c
C     .. Scalar Arguments ..
      REAL C,S
      INTEGER INCX,INCY,N
C     ..
C     .. Array Arguments ..
      REAL SX(*),SY(*)
C     ..
C     .. Local Scalars ..
      REAL STEMP
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
          STEMP = C*SX(IX) + S*SY(IY)
          SY(IY) = C*SY(IY) - S*SX(IX)
          SX(IX) = STEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
c
c       code for both increments equal to 1
c
   20 DO 30 I = 1,N
          STEMP = C*SX(I) + S*SY(I)
          SY(I) = C*SY(I) - S*SX(I)
          SX(I) = STEMP
   30 CONTINUE
      RETURN
      END
