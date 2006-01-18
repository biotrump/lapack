      COMPLEX FUNCTION CDOTC(N,CX,INCX,CY,INCY)
c
c     forms the dot product of two vectors, conjugating the first
c     vector.
c     jack dongarra, linpack,  3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
c
C     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
C     ..
C     .. Array Arguments ..
      COMPLEX CX(*),CY(*)
C     ..
C     .. Local Scalars ..
      COMPLEX CTEMP
      INTEGER I,IX,IY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC CONJG
C     ..
      CTEMP = (0.0,0.0)
      CDOTC = (0.0,0.0)
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
          CTEMP = CTEMP + CONJG(CX(IX))*CY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      CDOTC = CTEMP
      RETURN
c
c        code for both increments equal to 1
c
   20 DO 30 I = 1,N
          CTEMP = CTEMP + CONJG(CX(I))*CY(I)
   30 CONTINUE
      CDOTC = CTEMP
      RETURN
      END
