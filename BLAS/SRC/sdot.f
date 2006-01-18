      REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
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
      STEMP = 0.0e0
      SDOT = 0.0e0
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
          STEMP = STEMP + SX(IX)*SY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      SDOT = STEMP
      RETURN
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          STEMP = STEMP + SX(I)*SY(I)
   30 CONTINUE
      IF (N.LT.5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          STEMP = STEMP + SX(I)*SY(I) + SX(I+1)*SY(I+1) +
     +            SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3) + SX(I+4)*SY(I+4)
   50 CONTINUE
   60 SDOT = STEMP
      RETURN
      END
