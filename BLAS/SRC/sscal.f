      SUBROUTINE SSCAL(N,SA,SX,INCX)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to 1.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c
C     .. Scalar Arguments ..
      REAL SA
      INTEGER INCX,N
C     ..
C     .. Array Arguments ..
      REAL SX(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
c
c        code for increment not equal to 1
c
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          SX(I) = SA*SX(I)
   10 CONTINUE
      RETURN
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          SX(I) = SA*SX(I)
   30 CONTINUE
      IF (N.LT.5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          SX(I) = SA*SX(I)
          SX(I+1) = SA*SX(I+1)
          SX(I+2) = SA*SX(I+2)
          SX(I+3) = SA*SX(I+3)
          SX(I+4) = SA*SX(I+4)
   50 CONTINUE
      RETURN
      END
