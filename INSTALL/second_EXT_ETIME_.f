*> \brief \b SECOND  Using ETIME_
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition
*  ==========
*
*      REAL FUNCTION SECOND( )
*  
*  Purpose
*  =======
*
*>\details \b Purpose:
*>\verbatim
*>
*>  SECOND returns the user time for a process in seconds.
*>  This version gets the time from the system function ETIME_.
*>
*>\endverbatim
*
*  Authors
*  =======
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      REAL FUNCTION SECOND( )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     February 2007
* =====================================================================
*
*     .. Local Scalars ..
      REAL               T1
*     ..
*     .. Local Arrays ..
      REAL               TARRAY( 2 )
*     ..
*     .. External Functions ..
      REAL               ETIME_
      EXTERNAL           ETIME_
*     ..
*     .. Executable Statements ..
*
      T1 = ETIME_( TARRAY )
      SECOND = TARRAY( 1 )
      RETURN
*
*     End of SECOND
*
      END

