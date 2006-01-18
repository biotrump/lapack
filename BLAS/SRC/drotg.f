      SUBROUTINE DROTG(DA,DB,C,S)
c
c     construct givens plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
c
C     .. Scalar Arguments ..
      DOUBLE PRECISION C,DA,DB,S
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION R,ROE,SCALE,Z
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DSIGN,DSQRT
C     ..
      ROE = DB
      IF (DABS(DA).GT.DABS(DB)) ROE = DA
      SCALE = DABS(DA) + DABS(DB)
      IF (SCALE.NE.0.0d0) GO TO 10
      C = 1.0d0
      S = 0.0d0
      R = 0.0d0
      Z = 0.0d0
      GO TO 20
   10 R = SCALE*DSQRT((DA/SCALE)**2+ (DB/SCALE)**2)
      R = DSIGN(1.0d0,ROE)*R
      C = DA/R
      S = DB/R
      Z = 1.0d0
      IF (DABS(DA).GT.DABS(DB)) Z = S
      IF (DABS(DB).GE.DABS(DA) .AND. C.NE.0.0d0) Z = 1.0d0/C
   20 DA = R
      DB = Z
      RETURN
      END
