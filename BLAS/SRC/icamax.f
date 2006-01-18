      INTEGER FUNCTION ICAMAX(N,CX,INCX)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c
C     .. Scalar Arguments ..
      INTEGER INCX,N
C     ..
C     .. Array Arguments ..
      COMPLEX CX(*)
C     ..
C     .. Local Scalars ..
      REAL SMAX
      INTEGER I,IX
C     ..
C     .. External Functions ..
      REAL SCABS1
      EXTERNAL SCABS1
C     ..
      ICAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      ICAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
c
c        code for increment not equal to 1
c
      IX = 1
      SMAX = SCABS1(CX(1))
      IX = IX + INCX
      DO 10 I = 2,N
          IF (SCABS1(CX(IX)).LE.SMAX) GO TO 5
          ICAMAX = I
          SMAX = SCABS1(CX(IX))
    5     IX = IX + INCX
   10 CONTINUE
      RETURN
c
c        code for increment equal to 1
c
   20 SMAX = SCABS1(CX(1))
      DO 30 I = 2,N
          IF (SCABS1(CX(I)).LE.SMAX) GO TO 30
          ICAMAX = I
          SMAX = SCABS1(CX(I))
   30 CONTINUE
      RETURN
      END
