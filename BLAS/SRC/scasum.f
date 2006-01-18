      REAL FUNCTION SCASUM(N,CX,INCX)
c
c     takes the sum of the absolute values of a complex vector and
c     returns a single precision result.
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
      REAL STEMP
      INTEGER I,NINCX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,AIMAG,REAL
C     ..
      SCASUM = 0.0e0
      STEMP = 0.0e0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
c
c        code for increment not equal to 1
c
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
   10 CONTINUE
      SCASUM = STEMP
      RETURN
c
c        code for increment equal to 1
c
   20 DO 30 I = 1,N
          STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
   30 CONTINUE
      SCASUM = STEMP
      RETURN
      END
