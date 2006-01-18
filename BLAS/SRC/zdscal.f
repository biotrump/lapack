      SUBROUTINE ZDSCAL(N,DA,ZX,INCX)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c
C     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX ZX(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,IX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DCMPLX
C     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
c
c        code for increment not equal to 1
c
      IX = 1
      DO 10 I = 1,N
          ZX(IX) = DCMPLX(DA,0.0d0)*ZX(IX)
          IX = IX + INCX
   10 CONTINUE
      RETURN
c
c        code for increment equal to 1
c
   20 DO 30 I = 1,N
          ZX(I) = DCMPLX(DA,0.0d0)*ZX(I)
   30 CONTINUE
      RETURN
      END
