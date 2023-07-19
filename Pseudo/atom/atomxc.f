      SUBROUTINE ATOMXC( FUNCTL, AUTHOR, IREL,
     .                   NR, MAXR, RMESH, NSPIN, DENS,
     .                   EX, EC, DX, DC, VXC )

C *******************************************************************
C Finds total exchange-correlation energy and potential for a
C spherical electron density distribution.
C This version implements the Local (spin) Density Approximation and
C the Generalized-Gradient-Aproximation with the 'explicit mesh 
C functional' method of White & Bird, PRB 50, 4954 (1994).
C Gradients are 'defined' by numerical derivatives, using 2*NN+1 mesh
C   points, where NN is a parameter defined below
C Coded by L.C.Balbas and J.M.Soler. December 1996. Version 0.5.
C ************************* INPUT ***********************************
C CHARACTER*(*) FUNCTL : Functional to be used:
C              'LDA' or 'LSD' => Local (spin) Density Approximation
C                       'GGA' => Generalized Gradient Corrections
C                                Uppercase is optional
C CHARACTER*(*) AUTHOR : Parametrization desired:
C     'CA' or 'PZ' => LSD Perdew & Zunger, PRB 23, 5075 (1981)
C           'PW92' => LSD Perdew & Wang, PRB, 45, 13244 (1992). This is
C                     the local density limit of the next:
C            'PBE' => GGA Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)
C                     Uppercase is optional
C INTEGER IREL         : Relativistic exchange? (0=>no, 1=>yes)
C INTEGER NR           : Number of radial mesh points
C INTEGER MAXR         : Physical first dimension of RMESH, DENS and VXC
C REAL*8  RMESH(MAXR)  : Radial mesh points
C INTEGER NSPIN        : NSPIN=1 => unpolarized; NSPIN=2 => polarized
C REAL*8  DENS(MAXR,NSPIN) : Total (NSPIN=1) or spin (NSPIN=2) electron
C                            density at mesh points
C ************************* OUTPUT **********************************
C REAL*8  EX              : Total exchange energy
C REAL*8  EC              : Total correlation energy
C REAL*8  DX              : IntegralOf( rho * (eps_x - v_x) )
C REAL*8  DC              : IntegralOf( rho * (eps_c - v_c) )
C REAL*8  VXC(MAXR,NSPIN) : (Spin) exch-corr potential
C ************************ UNITS ************************************
C Distances in atomic units (Bohr).
C Densities in atomic units (electrons/Bohr**3)
C Energy unit depending of parameter EUNIT below
C ********* ROUTINES CALLED *****************************************
C GGAXC, LDAXC
C *******************************************************************

C Next line is nonstandard but may be suppressed
      IMPLICIT NONE

C Argument types and dimensions
      CHARACTER*(*)     FUNCTL, AUTHOR
      INTEGER           IREL, MAXR, NR, NSPIN
      DOUBLE PRECISION  DENS(MAXR,NSPIN), RMESH(MAXR), VXC(MAXR,NSPIN)
      DOUBLE PRECISION  DC, DX, EC, EX

C Fix the order of the numerical derivatives: the number of radial 
C points used is 2*NN+1
C MAXAUX must be larger than the number of radial mesh points
      INTEGER NN, MAXAUX
      PARAMETER ( NN     =    5 )
      PARAMETER ( MAXAUX = 2000 )

C Fix energy unit:  EUNIT=1.0 => Hartrees,
C                   EUNIT=0.5 => Rydbergs,
C                   EUNIT=0.03674903 => eV
      DOUBLE PRECISION EUNIT
      PARAMETER ( EUNIT = 0.5D0 )

C DVMIN is added to differential of volume to avoid division by zero
      DOUBLE PRECISION DVMIN
      PARAMETER ( DVMIN = 1.D-12 )

C Local variables and arrays
      LOGICAL
     .  GGA
      INTEGER
     .  IN, IN1, IN2, IR, IS, JN
      DOUBLE PRECISION
     .  AUX(MAXAUX), D(2), DECDD(2), DECDGD(3,2), DEXDD(2), DEXDGD(3,2),
     .  DGDM(-NN:NN), DGIDFJ(-NN:NN), DRDM, DVOL, 
     .  EPSC, EPSX, F1, F2, GD(3,2), PI
      EXTERNAL
     .  GGAXC, LDAXC

C Set GGA switch
      IF ( FUNCTL.EQ.'LDA' .OR. FUNCTL.EQ.'lda' .OR.
     .     FUNCTL.EQ.'LSD' .OR. FUNCTL.EQ.'lsd' ) THEN
        GGA = .FALSE.
      ELSEIF ( FUNCTL.EQ.'GGA' .OR. FUNCTL.EQ.'gga') THEN
        GGA = .TRUE.
      ELSE
        WRITE(6,*) 'CELLXC: Unknown functional ', FUNCTL
        STOP
      ENDIF

C Check size of auxiliary array
      IF (GGA .AND. MAXAUX.LT.NR)
     .  STOP 'ATOMXC: Parameter MAXAUX too small'

C Initialize output
      EX = 0
      EC = 0
      DX = 0
      DC = 0
      DO 20 IS = 1,NSPIN
        DO 10 IR = 1,NR
          VXC(IR,IS) = 0
   10   CONTINUE
   20 CONTINUE

C Get number pi
      PI = 4 * ATAN(1.D0)

C Loop on mesh points
      DO 140 IR = 1,NR

C       Find interval of neighbour points to calculate derivatives
        IN1 = MAX(  1, IR-NN ) - IR
        IN2 = MIN( NR, IR+NN ) - IR

C       Find weights of numerical derivation from Lagrange
C       interpolation formula
        DO 50 IN = IN1,IN2
          IF (IN .EQ. 0) THEN
            DGDM(IN) = 0
            DO 30 JN = IN1,IN2
              IF (JN.NE.0) DGDM(IN) = DGDM(IN) + 1.D0 / (0 - JN)
   30       CONTINUE
          ELSE
            F1 = 1
            F2 = 1
            DO 40 JN = IN1,IN2
              IF (JN.NE.IN .AND. JN.NE.0) F1 = F1 * (0  - JN)
              IF (JN.NE.IN)               F2 = F2 * (IN - JN)
   40       CONTINUE
            DGDM(IN) = F1 / F2
          ENDIF
   50   CONTINUE

C       Find dr/dmesh
        DRDM = 0
        DO 60 IN = IN1,IN2
          DRDM = DRDM + RMESH(IR+IN) * DGDM(IN)
   60   CONTINUE

C       Find differential of volume. Use trapezoidal integration rule
        DVOL = 4 * PI * RMESH(IR)**2 * DRDM
C       DVMIN is a small number added to avoid a division by zero
        DVOL = DVOL + DVMIN
        IF (IR.EQ.1 .OR. IR.EQ.NR) DVOL = DVOL / 2
        IF (GGA) AUX(IR) = DVOL

C       Find the weights for the derivative d(gradF(i))/d(F(j)), of
C       the gradient at point i with respect to the value at point j
        IF (GGA) THEN
          DO 80 IN = IN1,IN2
            DGIDFJ(IN) = DGDM(IN) / DRDM
   80     CONTINUE
        ENDIF

C       Find density and gradient of density at this point
        DO 90 IS = 1,NSPIN
          D(IS) = DENS(IR,IS)
   90   CONTINUE
        IF (GGA) THEN
          DO 110 IS = 1,NSPIN
            GD(1,IS) = 0
            GD(2,IS) = 0
            GD(3,IS) = 0
            DO 100 IN = IN1,IN2
              GD(3,IS) = GD(3,IS) + DGIDFJ(IN) * DENS(IR+IN,IS)
  100       CONTINUE
  110     CONTINUE
        ENDIF

C       Find exchange and correlation energy densities and their 
C       derivatives with respect to density and density gradient
        IF (GGA) THEN
          CALL GGAXC( AUTHOR, IREL, NSPIN, D, GD,
     .                EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
        ELSE
          CALL LDAXC( AUTHOR, IREL, NSPIN, D, EPSX, EPSC, DEXDD, DECDD )
        ENDIF

C       Add contributions to exchange-correlation energy and its
C       derivatives with respect to density at all points
        DO 130 IS = 1,NSPIN
          EX = EX + DVOL * D(IS) * EPSX
          EC = EC + DVOL * D(IS) * EPSC
          DX = DX + DVOL * D(IS) * (EPSX - DEXDD(IS))
          DC = DC + DVOL * D(IS) * (EPSC - DECDD(IS))
          IF (GGA) THEN
            VXC(IR,IS) = VXC(IR,IS) + DVOL * ( DEXDD(IS) + DECDD(IS) )
            DO 120 IN = IN1,IN2
              DX= DX - DVOL * DENS(IR+IN,IS) * DEXDGD(3,IS) * DGIDFJ(IN)
              DC= DC - DVOL * DENS(IR+IN,IS) * DECDGD(3,IS) * DGIDFJ(IN)
              VXC(IR+IN,IS) = VXC(IR+IN,IS) + DVOL *
     .               (DEXDGD(3,IS) + DECDGD(3,IS)) * DGIDFJ(IN)
  120       CONTINUE
          ELSE
            VXC(IR,IS) = DEXDD(IS) + DECDD(IS)
          ENDIF
  130   CONTINUE

  140 CONTINUE

C Divide by volume element to obtain the potential (per electron)
      IF (GGA) THEN
        DO 160 IS = 1,NSPIN
          DO 150 IR = 1,NR
            DVOL = AUX(IR)
            VXC(IR,IS) = VXC(IR,IS) / DVOL
  150     CONTINUE
  160   CONTINUE
      ENDIF

C Divide by energy unit
      EX = EX / EUNIT
      EC = EC / EUNIT
      DX = DX / EUNIT
      DC = DC / EUNIT
      DO 180 IS = 1,NSPIN
        DO 170 IR = 1,NR
          VXC(IR,IS) = VXC(IR,IS) / EUNIT
  170   CONTINUE
  180 CONTINUE

      END
