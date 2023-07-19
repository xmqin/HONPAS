
      SUBROUTINE COLINEAR( COORPO, COLIN )

C **********************************************************************
C Checks if three points lye in the same straight line (if they are 
C colinear).
C Coded by J. Junquera November'98
C **********************************************************************

      use precision

      IMPLICIT NONE

      REAL(DP), INTENT(IN) ::
     .  COORPO(3,3)
 
      LOGICAL, INTENT(OUT) ::
     .  COLIN

C **********************************************************************
C REAL*8 COORPO(3,3)   : Coordinates of three points. COORPO(POINT,IX)
C LOGICAL COLIN        : True => Three points are colinear
C                        False=> Three points NOT colinear
C **********************************************************************

C Internal variables ---------------------------------------------------

      INTEGER 
     .  IX
 
      REAL(DP)
     .  VEC1(3), VEC2(3), NORMAL(3), EPS


      DATA EPS /1.0e-12_dp/

      DO IX = 1,3
        VEC1(IX) = COORPO(2,IX) - COORPO(1,IX)
        VEC2(IX) = COORPO(3,IX) - COORPO(1,IX)
      ENDDO

      normal =  cross_product( VEC1, VEC2)
     
      IF( (ABS(NORMAL(1)) .LT. EPS) .AND. 
     .    (ABS(NORMAL(2)) .LT. EPS) .AND.
     .    (ABS(NORMAL(3)) .LT. EPS) ) COLIN = .TRUE.

      CONTAINS

      function cross_product(a,b) result(c)
      real(dp), dimension(3), intent(in) :: a, b
      real(dp), dimension(3)             :: c

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)

      end function cross_product

      END


      
      
