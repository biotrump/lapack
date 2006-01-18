      SUBROUTINE ZROTG(CA,CB,C,S)

C     .. Scalar Arguments ..
      DOUBLE COMPLEX CA,CB,S
      DOUBLE PRECISION C
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX ALPHA
      DOUBLE PRECISION NORM,SCALE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC CDABS,DCMPLX,DCONJG,DSQRT
C     ..
      IF (CDABS(CA).NE.0.0d0) GO TO 10
      C = 0.0d0
      S = (1.0d0,0.0d0)
      CA = CB
      GO TO 20
   10 CONTINUE
      SCALE = CDABS(CA) + CDABS(CB)
      NORM = SCALE*DSQRT((CDABS(CA/DCMPLX(SCALE,0.0d0)))**2+
     +       (CDABS(CB/DCMPLX(SCALE,0.0d0)))**2)
      ALPHA = CA/CDABS(CA)
      C = CDABS(CA)/NORM
      S = ALPHA*DCONJG(CB)/NORM
      CA = ALPHA*NORM
   20 CONTINUE
      RETURN
      END
