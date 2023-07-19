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

      type(fftw3_plan_t) :: plan

      INTEGER
     .  NPX, NPY, NPZ

      REAL(DP)
     .  UCELL(3,3), ZREF, ZMIN, ZMAX, E, K(3), V0

      COMPLEX(DP)
     .  CW(0:NPX-1,0:NPY-1), CWE(0:NPX-1,0:NPY-1,0:NPZ-1)

C----------------------------------------------------------
C INTERNAL VARIABLES

      INTEGER
     .  NX, NY, NZ, NX1, NY1, IA1, IA2

      REAL(DP)
     .  DECAY, GX, GY, Z, VOLCEL, EAU, V0AU

      REAL(DP), SAVE ::
     .  STEPZ, A1(3), A2(3), A3(3), VOL, VEC1(3), VEC2(3)

      COMPLEX(DP)
     .  EXPSI(0:NPX-1,0:NPY-1)

      LOGICAL, SAVE :: FIRST

      EXTERNAL 
     .  VOLCEL

      EXTERNAL
     .  dfftw_plan_dft_2d, dfftw_execute, dfftw_destroy_plan

      DATA FIRST /.TRUE./


C CHANGE UNITS OF ENERGY TO A.U.

      EAU = E / 27.21162908D0
      V0AU = V0 / 27.21162908D0

C CHECK THAT STATE IS BOUND (E BELOW VACUUM LEVEL)
      IF (EAU .LT. V0AU) THEN

C INITIALIZE VARIABLES ............

      IF (FIRST) THEN

! the only problem is the decay
! which is calculated now by an ABS
! so in principle this can be commented out
! N.L. Fevrier 2005
!       IF (ZMIN .LT. ZREF) THEN
!         WRITE(6,*) 'error: the reference plane can not be above'
!         WRITE(6,*) '       the simulation area'
!         WRITE(6,*) '       STM.minZ must be larger than STM.RefZ'
!         STOP
!       ENDIF

C
        IF (NPZ .NE. 1) THEN
          STEPZ = (ZMAX - ZMIN) / (NPZ - 1)
        ELSE 
          STEPZ = 0.0D0
        ENDIF
C
        VOL = VOLCEL(UCELL)
C
!!      This was wrong !
!!        A1=(/UCELL(1,1), UCELL(1,2), UCELL(1,3)/)
!!        A2=(/UCELL(2,1), UCELL(2,2), UCELL(2,3)/)
!!        A3=(/UCELL(3,1), UCELL(3,2), UCELL(3,3)/)

        A1 = UCELL(:,1)
        A2 = UCELL(:,2)
        A3 = UCELL(:,3)

C
        CALL RECIPROCAL(A1,A2,A3,VOL,VEC1,VEC2)
C
        FIRST = .FALSE.
C
      ENDIF

C .....

C      REWIND(7)
C      DO NX=0,NPX-1
C        DO NY=0,NPY-1
C          WRITE(7,'(4f20.10)') UCELL(1,1)*NX/NPX, 
C     .             UCELL(2,2)*NY/NPY,
C     .             DREAL(CW(NX,NY)*DCONJG(CW(NX,NY)))
C        ENDDO
C        WRITE(7,*)
C      ENDDO


C DO DIRECT FOURIER TRANSFORM TO GET SPATIAL FREQUENCIES OF WF AT 
C REFERENCE PLANE ........

      call dfftw_plan_dft_2d (plan,NPX,NPY,CW,CW,FFTW_FORWARD, 
     .                        FFTW_ESTIMATE)
      call dfftw_execute (plan)
      call dfftw_destroy_plan(plan)

C .....

C LOOP OVER SIMULATION HEIGHTS ........
      call dfftw_plan_dft_2d (plan,NPX,NPY,EXPSI,EXPSI,FFTW_BACKWARD, 
     .                        FFTW_ESTIMATE)
      DO NZ = 1, NPZ
        Z = ZMIN + (NZ-1)*STEPZ

C LOOP OVER POINTS IN XY PLANE ...
        DO NY = 1, NPY
          NY1 = NY-1
C INDEX CHOSEN TO CENTER G's
          IA2 = MOD (NY1,NPY/2)-NY1/(NPY/2)*(NPY/2)
          DO NX = 1, NPX
            NX1 = NX-1
            IA1 = MOD (NX1,NPX/2)-NX1/(NPX/2)*(NPX/2)
C G COMPONENTS
            GX = IA1*VEC1(1) + IA2*VEC2(1)
            GY = IA1*VEC1(2) + IA2*VEC2(2)

C CALCULATE EXPONENT OF DECAY
           DECAY = DSQRT( 2.0D0*(V0AU-EAU) + 
     .                    (GX+K(1))**2 + (GY+K(2))**2)

C EXTEND WAVE FUNCTION TO CURRENT HEIGHT
            EXPSI(NX1,NY1) = CW(NX1,NY1) * DEXP(-DECAY*ABS(Z-ZREF)) / 
     .                       (NPX*NPY)
          ENDDO
        ENDDO
C ... END LOOP XY

C DO BACK FOURIER TRANSFORM TO GET REAL SPACE WF AT 
C REFERENCE PLANE .....

        call dfftw_execute (plan)

C .....

        CWE(:,:,NZ-1) = EXPSI(:,:)

      ENDDO
C ........ END LOOP Z
      call dfftw_destroy_plan(plan)


C      DO NX=0,NPX-1
C        DO NY=0,NPY-1
C          WRITE(7,'(3f20.10)') UCELL(1,1)*NX/NPX, 
C     .             UCELL(2,2)*NY/NPY,
C     .             DREAL(CW(NX,NY)*DCONJG(CW(NX,NY)))
C        ENDDO
C        WRITE(7,*)
C      ENDDO
C
C      REWIND(8)
C
C      DO NX=0,NPX-1
C        DO NY=0,NPY-1
C          WRITE(8,'(3f20.10)') UCELL(1,1)*NX/NPX, 
C     .             UCELL(2,2)*NY/NPY,
C     .             DREAL(CWE(NX,NY,NPZ-1)*DCONJG(CWE(NX,NY,NPZ-1)))
C        ENDDO
C        WRITE(8,*)
C      ENDDO



	else
       write(6,*) 'wf NOT considered:'
       write(6,*) ' ENERGY EIGENVALUE ABOVE VACUUM LEVEL'

	endif


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




