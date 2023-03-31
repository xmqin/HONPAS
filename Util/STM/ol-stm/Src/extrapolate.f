! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      SUBROUTINE EXTRAPOLATE(NPX,NPY,NPZ,ZREF,ZMIN,ZMAX,UCELL,V0,
     .                       CW,E,K,CWE)

C ***********************************************************************
C Routine that extrapolates the value of the wave function from that
C in a reference plane into the vacuum region
C
C P. Ordejon and N. Lorente, November 2004
C ************************ INPUT ****************************************
C integer NPX,NPY          : Number of points in x,y grid
C integer NPZ              : Number of z-planes in which psi is calculated
C real*8  ZREF             : Z value of the reference plane
C real*8  ZMIN             : Minimum Z value of z-plane scan
C real*8  ZMAX             : Maximum Z value of z-plane scan
C real*8  UCELL(3,3)       : Cell vectors of unit cell
C real*8  V0               : Value of the potential at vacuum
C complex*8 CW(0:NPX-1,0:NPY-1): Values of wavefunction at reference plane
C real*8  E                : Energy of the wavefunction
C real*8  K(3)             : K-point of wave function
C ************************ OUTPUT ***************************************
C complex*8 CWE(0:NPX-1,0:NPY-1,0:NPZ-1): 
C                            Wave Function extrapolated at all heights
C ***********************************************************************

      use precision, only: dp
      USE fftw3_mymod

      IMPLICIT NONE

      INTEGER
     .  NPX, NPY, NPZ

      REAL(DP)
     .  UCELL(3,3), ZREF, ZMIN, ZMAX, E, K(3), V0

      COMPLEX(DP)
     .  CW(0:NPX-1,0:NPY-1), CWE(0:NPX-1,0:NPY-1,0:NPZ-1)

C----------------------------------------------------------
C INTERNAL VARIABLES

      INTEGER ::  NX, NY, NZ, NX1, NY1, IA1, IA2

      REAL(DP) :: DECAY, GX, GY, Z, VOLCEL, EAU, V0AU

      REAL(DP), SAVE ::
     .  STEPZ, A1(3), A2(3), A3(3), VOL, VEC1(3), VEC2(3)

      COMPLEX(DP) :: EXPSI(0:NPX-1,0:NPY-1)

      LOGICAL, SAVE :: FIRST = .true.

      EXTERNAL  VOLCEL

C CHANGE UNITS OF ENERGY TO A.U.

      EAU = E / 27.21162908D0
      V0AU = V0 / 27.21162908D0

C CHECK THAT STATE IS BOUND (E BELOW VACUUM LEVEL)
      IF (EAU .GE. V0AU) THEN
         write(6,*) 'wf NOT considered:'
         write(6,*) ' ENERGY EIGENVALUE ABOVE VACUUM LEVEL'
         RETURN
      endif

C INITIALIZE VARIABLES ............

! the only problem is the decay
! which is calculated now by an ABS
! so in principle this can be commented out
!     N.L. Fevrier 2005

!! .... but what is the physical meaning of a decay *towards* the surface?? (AG)
         
!       IF (ZMIN .LT. ZREF) THEN
!         WRITE(6,*) 'error: the reference plane can not be above'
!         WRITE(6,*) '       the simulation area'
!         WRITE(6,*) '       STM.minZ must be larger than STM.RefZ'
!         STOP
!       ENDIF

      IF (FIRST) THEN

        IF (NPZ .NE. 1) THEN
          STEPZ = (ZMAX - ZMIN) / (NPZ - 1)
        ELSE 
          STEPZ = 0.0D0
        ENDIF
C
        VOL = VOLCEL(UCELL)

        A1 = UCELL(:,1)
        A2 = UCELL(:,2)
        A3 = UCELL(:,3)
        CALL RECIPROCAL(A1,A2,A3,VOL,VEC1,VEC2)

        FIRST = .FALSE.

      ENDIF

C DO DIRECT FOURIER TRANSFORM TO GET SPATIAL FREQUENCIES OF WF AT 
C REFERENCE PLANE ........

      ! Reverse dimensions for f2003 interface !!!
      plan =  fftw_plan_dft_2d (NPY,NPX,CW,CW,FFTW_FORWARD, 
     .                        FFTW_ESTIMATE)
      call fftw_execute_dft (plan,cw,cw)
      call fftw_destroy_plan(plan)

C LOOP OVER SIMULATION HEIGHTS ........
      plan = fftw_plan_dft_2d (NPY,NPX,EXPSI,EXPSI,FFTW_BACKWARD, 
     .                        FFTW_ESTIMATE)

      DO NZ = 1, NPZ
        Z = ZMIN + (NZ-1)*STEPZ

        ! NOTE: No "inward propagation...", despite what Nicolas says above
        ! Points with Z < ZREF are computed from un-propagated wfs
        IF (Z < ZREF) CYCLE  
        
C LOOP OVER POINTS IN XY PLANE ...
        DO NY = 1, NPY
          NY1 = NY-1
C     INDEX CHOSEN TO CENTER G's
          ! Ambiguous arithmetic here and below (precedence?)
          IA2 = MOD (NY1,NPY/2)-NY1/(NPY/2)*(NPY/2)
          DO NX = 1, NPX
            NX1 = NX-1
            IA1 = MOD (NX1,NPX/2)-NX1/(NPX/2)*(NPX/2)
C G COMPONENTS
            GX = IA1*VEC1(1) + IA2*VEC2(1)
            GY = IA1*VEC1(2) + IA2*VEC2(2)

C CALCULATE EXPONENT OF DECAY
           DECAY = SQRT( 2.0D0*(V0AU-EAU) + 
     .                    (GX+K(1))**2 + (GY+K(2))**2)

C EXTEND WAVE FUNCTION TO CURRENT HEIGHT
            EXPSI(NX1,NY1) = CW(NX1,NY1) * EXP(-DECAY*ABS(Z-ZREF)) / 
     .                       (NPX*NPY)
          ENDDO
        ENDDO        !  LOOP XY

C DO BACK FOURIER TRANSFORM TO GET REAL SPACE WF AT 
C REFERENCE PLANE .....

        call fftw_execute_dft (plan,expsi,expsi)

        CWE(:,:,NZ-1) = EXPSI(:,:)

      ENDDO     !  LOOP Z
      
      call fftw_destroy_plan(plan)

      RETURN
      END


        SUBROUTINE  reciprocal(a1,a2,a3,vol,vec1,vec2)
        implicit none
        real*8    vol
        real*8    pv(3),a1(3),a2(3),a3(3)
        real*8    vec1(3),vec2(3), PI_D

        parameter (PI_D = 3.141592653589793238462643383279502884197d0)

        call pro(a2,a3,pv)
              vec1=2.0D0*PI_D*pv/vol
        call pro(a3,a1,pv)
              vec2=2.0D0*PI_D*pv/vol

        END SUBROUTINE reciprocal

        SUBROUTINE pro(v1,v2,pv)
        implicit none
        real*8   v1(3),v2(3)
        real*8   pv(3)

        pv(1)=v1(2)*v2(3)-v1(3)*v2(2)
        pv(2)=v1(3)*v2(1)-v1(1)*v2(3)
        pv(3)=v1(1)*v2(2)-v1(2)*v2(1)

        END SUBROUTINE pro




