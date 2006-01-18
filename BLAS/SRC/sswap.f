      SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c
C     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
C     ..
C     .. Array Arguments ..
      REAL SX(*),SY(*)
C     ..
C     .. Local Scalars ..
      REAL STEMP
      INTEGER I,IX,IY,M,MP1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
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
          STEMP = SX(IX)
          SX(IX) = SY(IY)
          SY(IY) = STEMP
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 M = MOD(N,3)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          STEMP = SX(I)
          SX(I) = SY(I)
          SY(I) = STEMP
   30 CONTINUE
      IF (N.LT.3) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
          STEMP = SX(I)
          SX(I) = SY(I)
          SY(I) = STEMP
          STEMP = SX(I+1)
          SX(I+1) = SY(I+1)
          SY(I+1) = STEMP
          STEMP = SX(I+2)
          SX(I+2) = SY(I+2)
          SY(I+2) = STEMP
   50 CONTINUE
      RETURN
      END
