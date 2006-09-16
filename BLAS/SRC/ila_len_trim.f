      INTEGER FUNCTION ILA_LEN_TRIM(SUBNAM)
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     28 August, 2006
!
!     .. Scalar Arguments ..
      CHARACTER*(*) SUBNAM
!     ..
!
!  Purpose
!  =======
!
!  ILA_LEN_TRIM is called from testing and timing routines to remove
!  trailing spaces from its argument.  It is included in the library
!  for possible use within a user's XERBLA error-handing routine.
!
!  Arguments
!  =========
!
!  SUBNAM  (input) CHARACTER*(*)
!          Provides the string.
!
!  RETURN VALUE:  INTEGER
!          = N > 0 : The location of the last non-blank.
!          = 0     : The entire string is blank.
!
!     .. Local Scalars ..
      INTEGER I
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC LEN
!     ..
!     .. Executable Statements ..

      DO I=LEN(SUBNAM), 1, -1
         IF (SUBNAM(I:I) .NE. ' ') THEN
            ILA_LEN_TRIM = I
            RETURN
         END IF
      END DO
      ILA_LEN_TRIM = 0
      END
