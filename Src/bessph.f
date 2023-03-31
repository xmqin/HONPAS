! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!!@LICENSE

      MODULE m_bessph

      CONTAINS

      FUNCTION BESSPH (L,X)

!  RETURNS THE SPHERICAL BESSEL FUNCTION JL(X).
!  REF: ABRAMOWITZ AND STEGUN, FORMULAS 10.1.2 AND 10.1.19
!  WRITTEN BY J.SOLER. NOV/89.

      use precision, only: dp
      use sys, only: die

      integer, intent(in)  :: L
      real(dp), intent(in) :: X
      real(dp)             :: BESSPH

      integer, parameter :: nterms = 100
      real(dp),parameter ::
     .   zero = 0.0_dp, one = 1.0_dp, tiny=1.0e-15_dp

      integer :: i, n
      real(dp) :: switch, term, x2, sum, y, fnm1, fn, fnp1
      character(len=132):: msg

      SWITCH=MAX(1,2*L-1)
      IF (ABS(X).LT.SWITCH) THEN
!        USE POWER SERIES
         TERM=ONE
         DO 10 I=1,L
            TERM=TERM*X/(2*I+1)
   10    CONTINUE
         X2=X*X
         SUM=ZERO
         DO 20 I=1,NTERMS
            SUM=SUM+TERM
            TERM=(-TERM)*X2/(2*I*(2*I+2*L+1))
            IF (ABS(TERM).LT.TINY) GO TO 30
   20    CONTINUE
            WRITE(msg,*) 'BESSPH: SERIES HAS NOT CONVERGED. L,X=',L,X
            call die(trim(msg))
   30    BESSPH=SUM
      ELSE
!        USE EXPLICIT EXPRESSIONS OR RECURRENCE RELATION
         IF (L.EQ.0) THEN
            BESSPH=SIN(X)/X
         ELSEIF (L.EQ.1) THEN
            BESSPH=(SIN(X)/X-COS(X))/X
         ELSE
            Y=ONE/X
            FNM1=SIN(X)*Y
            FN=(FNM1-COS(X))*Y
            DO 40 N=1,L-1
               FNP1=(2*N+1)*Y*FN-FNM1
               FNM1=FN
               FN=FNP1
   40       CONTINUE
            BESSPH=FN
         ENDIF
      ENDIF
      END FUNCTION BESSPH

      END MODULE m_bessph

