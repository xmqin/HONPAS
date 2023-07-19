************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A, B    (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
***********************************************************
!  TOHW: Fortran-ize ROUNDING_MODE
      INTEGER FUNCTION ROUNDING_MODE( )
#ifdef CRAY
      ROUNDING_MODE = 0
#else
      ROUNDING_MODE = 1
#endif
      END
!**********************************************************
!  This version of the routine slightly adapted from that  
!  found in LAPACK3E
!
!  TOHW: manually rename SLAMCH to DLAMCH & remove unnecessary 
!        include/use (which does the renaming in LAPACK3E
!      #include lapacknames.inc
      FUNCTION DLAMCH( CMACH )
!      USE LA_CONSTANTS, ONLY: WP, ONE
      INTEGER, PARAMETER :: WP=KIND(0.0d0)
      REAL(WP) :: DLAMCH
!
!  -- LAPACK auxiliary routine --
!     E. Anderson, Cray Research Inc.
!     April 30, 1996
!     07-23-02:  LAPACK 3E version (eca)
!
!     .. Scalar Arguments ..
      CHARACTER          CMACH
!     ..
!
!  Purpose
!  =======
!
!  SLAMCH determines single precision machine parameters via inquiry
!  functions in Fortran 90, plus a C function to determine the rounding
!  mode.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by SLAMCH:
!          = 'E' or 'e',   SLAMCH := eps
!          = 'S' or 's ,   SLAMCH := sfmin
!          = 'B' or 'b',   SLAMCH := base
!          = 'P' or 'p',   SLAMCH := eps*base
!          = 'N' or 'n',   SLAMCH := t
!          = 'R' or 'r',   SLAMCH := rnd
!          = 'M' or 'm',   SLAMCH := emin
!          = 'U' or 'u',   SLAMCH := rmin
!          = 'L' or 'l',   SLAMCH := emax
!          = 'O' or 'o',   SLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            IC
      REAL(WP)           RMACH

      real(wp)           ONE
!     ..
!     .. External Functions ..
      INTEGER            ROUNDING_MODE
      EXTERNAL           ROUNDING_MODE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DIGITS, EPSILON, HUGE, ICHAR, MAXEXPONENT,
     &                   MINEXPONENT, RADIX, TINY
!     ..
!     .. Executable Statements ..
!
      IC = ICHAR( CMACH )
!
!     Assume ASCII character set and map to upper case.
!
      IF( IC.GE.97 .AND. IC.LE.122 ) IC = IC - 32
!
!     ICHAR( 'E' ) = 69     ICHAR( 'R' ) = 82
!     ICHAR( 'S' ) = 83     ICHAR( 'M' ) = 77
!     ICHAR( 'B' ) = 66     ICHAR( 'U' ) = 85
!     ICHAR( 'P' ) = 80     ICHAR( 'L' ) = 76
!     ICHAR( 'N' ) = 78     ICHAR( 'O' ) = 79
!
      IF( IC.EQ.69 ) THEN
         RMACH = EPSILON( ONE )
!
      ELSE IF( IC.EQ.83 ) THEN
         RMACH = TINY( ONE )
!
      ELSE IF( IC.EQ.66 ) THEN
         RMACH = RADIX( ONE )
!
      ELSE IF( IC.EQ.80 ) THEN
         RMACH = EPSILON( ONE )*RADIX( ONE )
!
      ELSE IF( IC.EQ.78 ) THEN
         RMACH = DIGITS( ONE )
!
      ELSE IF( IC.EQ.82 ) THEN
         RMACH = ROUNDING_MODE()
!
      ELSE IF( IC.EQ.77 ) THEN
         RMACH = MINEXPONENT( ONE )
!
      ELSE IF( IC.EQ.85 ) THEN
         RMACH = TINY( ONE )
!
      ELSE IF( IC.EQ.76 ) THEN
         RMACH = MAXEXPONENT( ONE )
!
      ELSE IF( IC.EQ.79 ) THEN
         RMACH = HUGE( ONE )
      END IF
!
      DLAMCH = RMACH
      RETURN
!
!     End of SLAMCH
!
      END
