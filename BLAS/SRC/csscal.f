      SUBROUTINE CSSCAL(N,SA,CX,INCX)
c
c     scales a complex vector by a real constant.
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
      COMPLEX CX(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,NINCX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC AIMAG,CMPLX,REAL
C     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
c
c        code for increment not equal to 1
c
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          CX(I) = CMPLX(SA*REAL(CX(I)),SA*AIMAG(CX(I)))
   10 CONTINUE
      RETURN
c
c        code for increment equal to 1
c
   20 DO 30 I = 1,N
          CX(I) = CMPLX(SA*REAL(CX(I)),SA*AIMAG(CX(I)))
   30 CONTINUE
      RETURN
      END
