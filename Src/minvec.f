! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module m_minvec

      use precision
      use sorting
      use sys

      implicit         none

      public :: minvec

      CONTAINS

      subroutine minvec(B0,BMIN,C)
C *******************************************************************
C  FINDS THE LATTICE BASIS OF MINIMUM LENGTH, I.E. SUCH TAHT ANY 
C  OTHER BASIS (NOT EQUIVALENT BY SYMMETRY) HAS ONE VECTOR LONGER.
C  WRITTEN BY J.MORENO AND J.SOLER. AUGUST 1989 AND OCTOBER 1997.
C ********* INPUT ***************************************************
C REAL*8 B0 : Cell vectors B0(xyj,vector)
C ********* OUTPUT **************************************************
C REAL*8 BMIN : Minimum cell vectors B0(xyj,vector)
C *******************************************************************

C
C  Modules

      real(dp), intent(in)  ::        B0(3,3)
      real(dp), intent(out) ::        BMIN(3,3)
      real(dp), intent(out) ::        C(3,3)

      real(dp)        :: ddot,VOLCEL
      external         ddot,RECLAT,VOLCEL

      integer          I,I1,I2,I3,ITER,J,NITER,IAUX(3)
      real(dp)         AUX(3,3),B(3,3),B2(1,3),BNEW(3),BNEW2,
     .                 EPS,VNEW,V0

      PARAMETER (EPS=1.D-8,NITER=100)

      V0=ABS(VOLCEL(B0))
      IF (V0.LT.EPS)
     $     call die('MINVEC: BASIS VECTORS ARE LINEARLY DEPENDENT')
      DO I=1,3
        DO J=1,3
          B(J,I)=B0(J,I)
        ENDDO
        B2(1,I)=ddot(3,B(1,I),1,B(1,I),1)
      ENDDO

      DO 50 ITER=1,NITER
        CALL ORDIX(B2,1,3,IAUX)
        CALL ORDER(B2,1,3,IAUX)
        CALL ORDER(B ,3,3,IAUX)
        DO I1=0,1
          DO I2=-1,1
            DO I3=-1,1
              IF (I1.EQ.0.AND.I2.NE.1) GO TO 40
              IF (I2.EQ.0.AND.I3.EQ.0) GO TO 40
              BNEW(1)=B(1,1)*I1+B(1,2)*I2+B(1,3)*I3
              BNEW(2)=B(2,1)*I1+B(2,2)*I2+B(2,3)*I3
              BNEW(3)=B(3,1)*I1+B(3,2)*I2+B(3,3)*I3
              BNEW2=ddot(3,BNEW,1,BNEW,1)
              DO I=3,1,-1
                IF (BNEW2+EPS.GE.B2(1,I)) GO TO 40
                CALL VOLNEW(B,BNEW,I,VNEW)
                IF (ABS((VNEW-V0)/V0).LT.EPS) THEN
                  B(1,I)=BNEW(1)
                  B(2,I)=BNEW(2)
                  B(3,I)=BNEW(3)
                  B2(1,I)=BNEW2
                  GO TO 50
                END IF
              ENDDO
  40          CONTINUE
            ENDDO
          ENDDO
        ENDDO
        GOTO 55
  50  CONTINUE

      call die('MINVEC: ERROR: Iteration has not converged')

  55  CONTINUE

      IF (VOLCEL(B).LT.0.D0) THEN
        B(1,3)=-B(1,3)
        B(2,3)=-B(2,3)
        B(3,3)=-B(3,3)
      ENDIF
      CALL RECLAT(B0,AUX,0)
      DO I=1,3
        DO J=1,3
          C(J,I)=NINT(ddot(3,AUX(1,J),1,B(1,I),1))
        ENDDO
      ENDDO
      DO I=1,3
        DO J=1,3
          B(J,I)=B0(J,1)*C(1,I)+B0(J,2)*C(2,I)+B0(J,3)*C(3,I)
        ENDDO
      ENDDO
      DO I=1,3
        DO J=1,3
          BMIN(J,I)=B(J,I)
        ENDDO
      ENDDO

      CONTAINS
      
      subroutine volnew(A,ANEW,INEW,VOL)
      integer          INEW
      real(dp)         A(3,3),ANEW(3),AUX(3,3),VOL
      AUX(1:3,1:3)=A(1:3,1:3)
      AUX(1:3,INEW)=ANEW(1:3)
      VOL=ABS(VOLCEL(AUX))
      end subroutine volnew

      end subroutine minvec
      end module m_minvec
