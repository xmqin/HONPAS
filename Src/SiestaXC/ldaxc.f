!!@LICENSE
!
!******************************************************************************
! MODULE m_ldaxc
! Provides routines for LDA XC functional evaluation
!******************************************************************************
!
!   PUBLIC procedures available from this module:
! ldaxc   ! General subroutine for all coded LDA XC functionals
! exchng  ! Local exchange
! pw92c   ! Perdew & Wang, PRB, 45, 13244 (1992) (Correlation only)
! pw92xc  ! Perdew & Wang, PRB, 45, 13244 (1992)
! pzxc    ! Perdew & Zunger, PRB 23, 5075 (1981)
!
!   PUBLIC parameters, types, and variables available from this module:
! none
!
!******************************************************************************
!
!   USED module procedures:
! use sys,        only: die             ! Termination routine
!
!   USED module parameters:
! use precision,  only: dp              ! Real double precision type
!
!   EXTERNAL procedures used:
! none
!
!******************************************************************************

      MODULE m_ldaxc

      implicit none

      PUBLIC::
     .  ldaxc,   ! General subroutine for all coded LDA XC functionals
     .  exchng,  ! Local exchange
     .  pw92c,   ! Perdew & Wang, PRB, 45, 13244 (1992) (Correlation only)
     .  pw92xc,  ! Perdew & Wang, PRB, 45, 13244 (1992)
     .  pzxc     ! Perdew & Zunger, PRB 23, 5075 (1981)

      PRIVATE ! Nothing is declared public beyond this point

      CONTAINS

      SUBROUTINE LDAXC( AUTHOR, IREL, nspin, D, EPSX, EPSC, VX, VC,
     .                  DVXDN, DVCDN )

C ******************************************************************
C Finds the exchange and correlation energies and potentials, in the
C Local (spin) Density Approximation.
C Written by L.C.Balbas and J.M.Soler, Dec'96.
C Non-collinear spin added by J.M.Soler, May'98
C *********** INPUT ************************************************
C CHARACTER*(*) AUTHOR : Parametrization desired:
C     'CA' or 'PZ' => LSD Perdew & Zunger, PRB 23, 5075 (1981)
C           'PW92' => LSD Perdew & Wang, PRB, 45, 13244 (1992)
C                     Uppercase is optional
C INTEGER IREL     : Relativistic exchange? (0=>no, 1=>yes)
C INTEGER nspin    : nspin=1 => unpolarized; nspin=2 => polarized;
C                    nspin=4 => non-collinear polarization
C REAL*8  D(nspin) : Local (spin) density. For non-collinear
C                    polarization, the density matrix is given by:
C                    D(1)=D11, D(2)=D22, D(3)=Real(D12), D(4)=Im(D12)
C *********** OUTPUT ***********************************************
C REAL*8 EPSX, EPSC : Exchange and correlation energy densities
C REAL*8 VX(nspin), VC(nspin) : Exchange and correlation potentials,
C                               defined as dExc/dD(ispin)
C REAL*8 DVXDN(nspin,nspin)  :  Derivative of exchange potential with
C                               respect the charge density, defined 
C                               as DVx(spin1)/Dn(spin2)
C REAL*8 DVCDN(nspin,nspin)  :  Derivative of correlation potential
C                               respect the charge density, defined 
C                               as DVc(spin1)/Dn(spin2)
C *********** UNITS ************************************************
C Lengths in Bohr, energies in Hartrees
C ******************************************************************

      use precision, only : dp
      use sys,       only : die

      implicit          none

      CHARACTER*(*),intent(in) :: AUTHOR ! LDA flavour ('PZ'|'PW92')
      INTEGER, intent(in) :: IREL        ! Relativistic exchange? 0=>no, 1=>yes
      INTEGER, intent(in) :: nspin       ! Number of spin components
      real(dp),intent(in) :: D(nspin)    ! Electron density (matrix)
      real(dp),intent(out):: EPSX        ! Exchange energy per electron
      real(dp),intent(out):: EPSC        ! Correlation energy per electron
      real(dp),intent(out):: VX(nspin)   ! Exchange potential
      real(dp),intent(out):: VC(nspin)   ! Correlation potential
      real(dp),intent(out):: DVXDN(nspin,nspin) ! dVX(spin1)/dDens(spin2)
      real(dp),intent(out):: DVCDN(nspin,nspin) ! dVC(spin1)/dDens(spin2)

      INTEGER           IS, NS, ISPIN1, ISPIN2
      real(dp)          DD(2), DPOL, DTOT, TINY, VCD(2), VPOL, VXD(2)

      PARAMETER ( TINY = 1.D-12 )

      IF (nspin .EQ. 4) THEN
C Find eigenvalues of density matrix (up and down densities
C along the spin direction)
C Note: D(1)=D11, D(2)=D22, D(3)=Real(D12), D(4)=Im(D12)
        NS = 2
        DTOT = D(1) + D(2)
        DPOL = SQRT( (D(1)-D(2))**2 + 4.D0*(D(3)**2+D(4)**2) )
        DD(1) = 0.5D0 * ( DTOT + DPOL )
        DD(2) = 0.5D0 * ( DTOT - DPOL )
      ELSE
        NS = nspin
        DO 10 IS = 1,nspin
cag       Avoid negative densities
          DD(IS) = max(D(IS),0.0d0)
   10   CONTINUE
      ENDIF


      DO ISPIN2 = 1, nspin
        DO ISPIN1 = 1, nspin
          DVXDN(ISPIN1,ISPIN2) = 0.D0
          DVCDN(ISPIN1,ISPIN2) = 0.D0
        ENDDO
      ENDDO

      IF ( AUTHOR.EQ.'CA' .OR. AUTHOR.EQ.'ca' .OR.
     .     AUTHOR.EQ.'PZ' .OR. AUTHOR.EQ.'pz') THEN
        CALL PZXC( IREL, NS, DD, EPSX, EPSC, VXD, VCD, DVXDN, DVCDN )
      ELSEIF ( AUTHOR.EQ.'PW92' .OR. AUTHOR.EQ.'pw92' ) THEN
        CALL PW92XC( IREL, NS, DD, EPSX, EPSC, VXD, VCD )
      ELSE
        call die('LDAXC: Unknown author ' // trim(AUTHOR))
      ENDIF

      IF (nspin .EQ. 4) THEN
C Find dE/dD(ispin) = dE/dDup * dDup/dD(ispin) +
C                     dE/dDdown * dDown/dD(ispin)
        VPOL  = (VXD(1)-VXD(2)) * (D(1)-D(2)) / (DPOL+TINY)
        VX(1) = 0.5D0 * ( VXD(1) + VXD(2) + VPOL )
        VX(2) = 0.5D0 * ( VXD(1) + VXD(2) - VPOL )
        VX(3) = (VXD(1)-VXD(2)) * D(3) / (DPOL+TINY)
        VX(4) = (VXD(1)-VXD(2)) * D(4) / (DPOL+TINY)
        VPOL  = (VCD(1)-VCD(2)) * (D(1)-D(2)) / (DPOL+TINY)
        VC(1) = 0.5D0 * ( VCD(1) + VCD(2) + VPOL )
        VC(2) = 0.5D0 * ( VCD(1) + VCD(2) - VPOL )
        VC(3) = (VCD(1)-VCD(2)) * D(3) / (DPOL+TINY)
        VC(4) = (VCD(1)-VCD(2)) * D(4) / (DPOL+TINY)
      ELSE
        DO 20 IS = 1,nspin
          VX(IS) = VXD(IS)
          VC(IS) = VCD(IS)
   20   CONTINUE
      ENDIF
      END SUBROUTINE LDAXC



      subroutine exchng( IREL, NSP, DS, EX, VX )

C *****************************************************************
C  Finds local exchange energy density and potential
C  Adapted by J.M.Soler from routine velect of Froyen's 
C    pseudopotential generation program. Madrid, Jan'97. Version 0.5.
C  Relativistic exchange modified by JMS, May.2014
C **** Input ******************************************************
C INTEGER IREL    : relativistic-exchange switch (0=no, 1=yes)
C INTEGER NSP     : spin-polarizations (1=>unpolarized, 2=>polarized)
C REAL*8  DS(NSP) : total (nsp=1) or spin (nsp=2) electron density
C **** Output *****************************************************
C REAL*8  EX      : exchange energy density
C REAL*8  VX(NSP) : (spin-dependent) exchange potential
C **** Units ******************************************************
C Densities in electrons/Bohr**3
C Energies in Hartrees
C *****************************************************************

      use precision, only: dp
      implicit none

      integer, intent(in) :: nsp, irel
      real(dp), intent(in)             :: DS(NSP)
      real(dp), intent(out)            :: VX(NSP)
      real(dp), intent(out)            :: EX

      real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp
      real(dp), parameter :: pfive = 0.5_dp, opf = 1.5_dp
!      real(dp), parameter :: c = 137.035999_dp       ! speed of light in a.u.
!      real(dp), parameter :: C014 = 0.014_dp         ! (9*pi/4)^(1/3)/c
      real(dp), parameter :: C014 = 0.0140047747_dp  ! updated JMS, May.2014

      real(dp) :: a0, alp, sb, rs
      real(dp) :: pi, trd, ftrd, tftm
      real(dp) :: d1, d2, d, z, fz, fzp, vxp, exp_var
      real(dp) :: beta, vxf, exf, alb
      real(dp) :: dBETAdD, dETAdD, dGAMMAdD, dPHIdD, dKFdD
      real(dp) :: ETA, GAMMA, KF, PHI

      PI=4*ATAN(ONE)
      TRD = ONE/3
      FTRD = 4*TRD
      TFTM = 2._dp**FTRD-2

      IF (NSP .EQ. 2) THEN
        D1 = MAX(DS(1),ZERO)
        D2 = MAX(DS(2),ZERO)
        D = D1 + D2
        IF (D .LE. ZERO) THEN
          EX = ZERO
          VX(1) = ZERO
          VX(2) = ZERO
          RETURN
        ENDIF
        Z = (D1 - D2) / D
        FZ = ((1+Z)**FTRD+(1-Z)**FTRD-2)/TFTM
        FZP = FTRD*((1+Z)**TRD-(1-Z)**TRD)/TFTM 
      ELSE
        D = DS(1)
        IF (D .LE. ZERO) THEN
          EX = ZERO
          VX(1) = ZERO
          RETURN
        ENDIF
        Z = ZERO
        FZ = ZERO
        FZP = ZERO
      ENDIF

      A0 = (4/(9*PI))**TRD 
      ALP = 2 * TRD                      ! X-alpha parameter
      RS = (3 / (4*PI*D) )**TRD
      VXP = -(3*ALP/(2*PI*A0*RS))        ! VX=-KF/PI
      EXP_VAR = 3*VXP/4                  ! epsX=-3*KF/4/PI
      IF (IREL .EQ. 1) THEN
        ! Ref: MacDonald and Vosco, J.Phys C 12, 2977 (1979)
        BETA = C014/RS                   ! ratio of Fermi to light speed
        SB = SQRT(1+BETA*BETA)
        ALB = LOG(BETA+SB)
        VXP = VXP * (-PFIVE + OPF * ALB / (BETA*SB))
        EXP_VAR = EXP_VAR * (ONE-OPF*((BETA*SB-ALB)/BETA**2)**2) 
      ENDIF

      IF (NSP .EQ. 2) THEN
        VXF = 2**TRD*VXP
        EXF = 2**TRD*EXP_VAR
        VX(1) = VXP + FZ*(VXF-VXP) + (1-Z)*FZP*(EXF-EXP_VAR)
        VX(2) = VXP + FZ*(VXF-VXP) - (1+Z)*FZP*(EXF-EXP_VAR)
        EX    = EXP_VAR + FZ*(EXF-EXP_VAR)
      ELSE
        VX(1) = VXP
        EX    = EXP_VAR
      ENDIF
      END subroutine exchng



      SUBROUTINE PW92C( nspin, Dens, EC, VC )

C ********************************************************************
C Implements the Perdew-Wang'92 local correlation (beyond RPA).
C Ref: J.P.Perdew & Y.Wang, PRB, 45, 13244 (1992)
C Written by L.C.Balbas and J.M.Soler. Dec'96.  Version 0.5.
C ********* INPUT ****************************************************
C INTEGER nspin       : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin) : Local (spin) density
C ********* OUTPUT ***************************************************
C REAL*8  EC        : Correlation energy density
C REAL*8  VC(nspin) : Correlation (spin) potential
C ********* UNITS ****************************************************
C Densities in electrons per Bohr**3
C Energies in Hartrees
C ********* ROUTINES CALLED ******************************************
C None
C ********************************************************************

      use precision, only : dp

C Next line is nonstandard but may be supressed
      implicit          none

C Argument types and dimensions
      INTEGER           nspin
      real(dp)          Dens(nspin), EC, VC(nspin)

C Internal variable declarations
      INTEGER           IG
      real(dp)          A(0:2), ALPHA1(0:2), B, BETA(0:2,4), C,
     .                  DBDRS, DECDD(2), DECDRS, DECDZ, DENMIN, DFDZ,
     .                  DGDRS(0:2), DCDRS, DRSDD, DTOT, DZDD(2),
     .                  F, FPP0, FOUTHD, G(0:2), HALF, ONE,
     .                  P(0:2), PI, RS, THD, THRHLF, ZETA

C Add tiny numbers to avoid numerical errors
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( ONE    = 1.D0 + 1.D-12 )

C Fix some numerical constants
      PARAMETER ( FOUTHD=4.D0/3.D0, HALF=0.5D0,
     .            THD=1.D0/3.D0, THRHLF=1.5D0 )

C Parameters from Table I of Perdew & Wang, PRB, 45, 13244 (92)
      DATA P      / 1.00d0,     1.00d0,     1.00d0     /
      DATA A      / 0.031091d0, 0.015545d0, 0.016887d0 /
      DATA ALPHA1 / 0.21370d0,  0.20548d0,  0.11125d0  /
      DATA BETA   / 7.5957d0,  14.1189d0,  10.357d0,
     .              3.5876d0,   6.1977d0,   3.6231d0,
     .              1.6382d0,   3.3662d0,   0.88026d0,
     .              0.49294d0,  0.62517d0,  0.49671d0 /

C Find rs and zeta
      PI = 4 * ATAN(1.D0)
      IF (nspin .EQ. 1) THEN
        DTOT = MAX( DENMIN, Dens(1) )
        ZETA = 0
        RS = ( 3 / (4*PI*DTOT) )**THD
C       Find derivatives dRs/dDens and dZeta/dDens
        DRSDD = (- RS) / DTOT / 3
        DZDD(1) = 0
      ELSE
        DTOT = MAX( DENMIN, Dens(1)+Dens(2) )
        ZETA = ( Dens(1) - Dens(2) ) / DTOT
        RS = ( 3 / (4*PI*DTOT) )**THD
        DRSDD = (- RS) / DTOT / 3
        DZDD(1) =   (ONE - ZETA) / DTOT
        DZDD(2) = - (ONE + ZETA) / DTOT
      ENDIF

C Find eps_c(rs,0)=G(0), eps_c(rs,1)=G(1) and -alpha_c(rs)=G(2)
C using eq.(10) of cited reference (Perdew & Wang, PRB, 45, 13244 (92))
      DO 20 IG = 0,2
        B = BETA(IG,1) * RS**HALF   +
     .      BETA(IG,2) * RS         +
     .      BETA(IG,3) * RS**THRHLF +
     .      BETA(IG,4) * RS**(P(IG)+1)
        DBDRS = BETA(IG,1) * HALF      / RS**HALF +
     .          BETA(IG,2)                         +
     .          BETA(IG,3) * THRHLF    * RS**HALF +
     .          BETA(IG,4) * (P(IG)+1) * RS**P(IG)
        C = 1 + 1 / (2 * A(IG) * B)
        DCDRS = - ( (C-1) * DBDRS / B )
        G(IG) = (- 2) * A(IG) * ( 1 + ALPHA1(IG)*RS ) * LOG(C)
        DGDRS(IG) = (- 2) *A(IG) * ( ALPHA1(IG) * LOG(C) +
     .                            (1+ALPHA1(IG)*RS) * DCDRS / C )
   20 CONTINUE

C Find f''(0) and f(zeta) from eq.(9)
      C = 1 / (2**FOUTHD - 2)
      FPP0 = 8 * C / 9
      F = ( (ONE+ZETA)**FOUTHD + (ONE-ZETA)**FOUTHD - 2 ) * C
      DFDZ = FOUTHD * ( (ONE+ZETA)**THD - (ONE-ZETA)**THD ) * C

C Find eps_c(rs,zeta) from eq.(8)
      EC = G(0) - G(2) * F / FPP0 * (ONE-ZETA**4) +
     .    (G(1)-G(0)) * F * ZETA**4
      DECDRS = DGDRS(0) - DGDRS(2) * F / FPP0 * (ONE-ZETA**4) +
     .        (DGDRS(1)-DGDRS(0)) * F * ZETA**4
      DECDZ = (- G(2)) / FPP0 * ( DFDZ*(ONE-ZETA**4) - F*4*ZETA**3 ) +
     .        (G(1)-G(0)) * ( DFDZ*ZETA**4 + F*4*ZETA**3 )
      
C Find correlation potential
      IF (nspin .EQ. 1) THEN
        DECDD(1) = DECDRS * DRSDD
        VC(1) = EC + DTOT * DECDD(1)
      ELSE
        DECDD(1) = DECDRS * DRSDD + DECDZ * DZDD(1)
        DECDD(2) = DECDRS * DRSDD + DECDZ * DZDD(2)
        VC(1) = EC + DTOT * DECDD(1)
        VC(2) = EC + DTOT * DECDD(2)
      ENDIF

      END SUBROUTINE PW92C



      SUBROUTINE PW92XC( IREL, nspin, Dens, EPSX, EPSC, VX, VC )

C ********************************************************************
C Implements the Perdew-Wang'92 LDA/LSD exchange correlation
C Ref: J.P.Perdew & Y.Wang, PRB, 45, 13244 (1992)
C Written by L.C.Balbas and J.M.Soler. Dec'96. Version 0.5.
C ********* INPUT ****************************************************
C INTEGER IREL        : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin       : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin) : Local (spin) density
C ********* OUTPUT ***************************************************
C REAL*8  EPSX       : Exchange energy density
C REAL*8  EPSC       : Correlation energy density
C REAL*8  VX(nspin)  : Exchange (spin) potential
C REAL*8  VC(nspin)  : Correlation (spin) potential
C ********* UNITS ****************************************************
C Densities in electrons per Bohr**3
C Energies in Hartrees
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      use precision, only : dp

      implicit          none
      INTEGER           IREL, nspin
      real(dp)          Dens(nspin), EPSX, EPSC, VC(nspin), VX(nspin)

      CALL EXCHNG( IREL, nspin, Dens, EPSX, VX )
      CALL PW92C( nspin, Dens, EPSC, VC )
      END SUBROUTINE PW92XC



      SUBROUTINE PZXC( IREL, NSP, DS, EX, EC, VX, VC, DVXDN, DVCDN )

C *****************************************************************
C  Perdew-Zunger parameterization of Ceperley-Alder exchange and 
C  correlation. Ref: Perdew & Zunger, Phys. Rev. B 23 5075 (1981).
C  Adapted by J.M.Soler from routine velect of Froyen's 
C    pseudopotential generation program. Madrid, Jan'97.
C **** Input *****************************************************
C INTEGER IREL    : relativistic-exchange switch (0=no, 1=yes)
C INTEGER NSP     : spin-polarizations (1=>unpolarized, 2=>polarized)
C REAL*8  DS(NSP) : total (nsp=1) or spin (nsp=2) electron density
C **** Output *****************************************************
C REAL*8  EX            : exchange energy density
C REAL*8  EC            : correlation energy density
C REAL*8  VX(NSP)       : (spin-dependent) exchange potential
C REAL*8  VC(NSP)       : (spin-dependent) correlation potential
C REAL*8  DVXDN(NSP,NSP): Derivative of the exchange potential
C                         respect the charge density, 
C                         Dvx(spin1)/Dn(spin2)
C REAL*8  DVCDN(NSP,NSP): Derivative of the correlation potential
C                         respect the charge density, 
C                         Dvc(spin1)/Dn(spin2)
C **** Units *******************************************************
C Densities in electrons/Bohr**3
C Energies in Hartrees
C *****************************************************************

      use precision, only: dp

      implicit none
      
       integer  :: nsp, irel, isp1, isp2, isp
       real(dp) :: DS(NSP), VX(NSP), VC(NSP), 
     .           DVXDN(NSP,NSP), DVCDN(NSP,NSP)
       real(dp), parameter ::
     $      ZERO=0.D0,ONE=1.D0,PFIVE=.5D0,OPF=1.5D0,PNN=.99D0,
     $      PTHREE=0.3D0,PSEVF=0.75D0,C0504=0.0504D0,
     $      C0254=0.0254D0,C014=0.014D0,C0406=0.0406D0,
     $      C15P9=15.9D0,C0666=0.0666D0,C11P4=11.4D0,
     $      C045=0.045D0,C7P8=7.8D0,C88=0.88D0,C20P59=20.592D0,
     $      C3P52=3.52D0,C0311=0.0311D0,C0014=0.0014D0,
     $      C0538=0.0538D0,C0096=0.0096D0,C096=0.096D0,
     $      C0622=0.0622D0,C004=0.004D0,C0232=0.0232D0,
     $      C1686=0.1686D0,C1P398=1.3981D0,C2611=0.2611D0,
     $      C2846=0.2846D0,C1P053=1.0529D0,C3334=0.3334D0

C    Ceperly-Alder 'ca' constants. Internal energies in Rydbergs.
       real(dp), parameter ::
     $      CON1=1.D0/6, CON2=0.008D0/3, CON3=0.3502D0/3,
     $      CON4=0.0504D0/3, CON5=0.0028D0/3, CON6=0.1925D0/3,
     $      CON7=0.0206D0/3, CON8=9.7867D0/6, CON9=1.0444D0/3,
     $      CON10=7.3703D0/6, CON11=1.3336D0/3

C      X-alpha parameter:
       real(dp), PARAMETER :: ALP = 2.D0 / 3.D0 

C      Other variables converted into parameters by J.M.Soler
       real(dp), parameter ::
     $       TINY = 1.D-6 ,
     $       PI   = 3.14159265358979312_dp,
     $       TWO  = 2.0D0,
     $       HALF = 0.5D0,
     $       TRD  = 1.D0 / 3.D0,
     $       FTRD = 4.D0 / 3.D0,
     $       TFTM = 0.51984209978974638D0,
     $       A0   = 0.52106176119784808D0,
     $       CRS  = 0.620350490899400087D0,
     $       CXP  = (- 3.D0) * ALP / (PI*A0),
     $       CXF  = 1.25992104989487319D0 

       real(dp)  :: d1, d2, d, z, fz, fzp
       real(dp)  :: ex, ec, dfzpdn, rs, vxp, exp_var
       real(dp)  :: beta, sb, alb, vxf, exf, dvxpdn
       real(dp)  :: dvxfdn, sqrs, te, be, ecp, vcp
       real(dp)  :: dtedn, be2, dbedn, dvcpdn, decpdn
       real(dp)  :: ecf, vcf, dvcfdn, decfdn, rslog


C      Find density and polarization
       IF (NSP .EQ. 2) THEN
         D1 = MAX(DS(1),ZERO)
         D2 = MAX(DS(2),ZERO)
         D = D1 + D2
         IF (D .LE. ZERO) THEN
           EX = ZERO
           EC = ZERO
           VX(1) = ZERO
           VX(2) = ZERO
           VC(1) = ZERO
           VC(2) = ZERO
           RETURN
         ENDIF
c
c        Robustness enhancement by Jose Soler (August 2002)
c
         Z = (D1 - D2) / D
         IF (Z .LE. -ONE) THEN
           FZ = (TWO**FTRD-TWO)/TFTM
           FZP = -FTRD*TWO**TRD/TFTM
           DFZPDN = FTRD*TRD*TWO**(-ALP)/TFTM
         ELSEIF (Z .GE. ONE) THEN
           FZ = (TWO**FTRD-TWO)/TFTM
           FZP = FTRD*TWO**TRD/TFTM
           DFZPDN = FTRD*TRD*TWO**(-ALP)/TFTM
         ELSE
           FZ = ((ONE+Z)**FTRD+(ONE-Z)**FTRD-TWO)/TFTM
           FZP = FTRD*((ONE+Z)**TRD-(ONE-Z)**TRD)/TFTM 
           DFZPDN = FTRD*TRD*((ONE+Z)**(-ALP) + (ONE-Z)**(-ALP))/TFTM
         ENDIF
       ELSE
         D = DS(1)
         IF (D .LE. ZERO) THEN
           EX = ZERO
           EC = ZERO
           VX(1) = ZERO
           VC(1) = ZERO
           RETURN
         ENDIF
         Z = ZERO
         FZ = ZERO
         FZP = ZERO
       ENDIF
       RS = CRS / D**TRD

C      Exchange
       VXP = CXP / RS
       EXP_VAR = 0.75D0 * VXP
       IF (IREL .EQ. 1) THEN
         BETA = C014/RS
         IF (BETA .LT. TINY) THEN
           SB = ONE + HALF*BETA**2
           ALB = BETA
         ELSE
           SB = SQRT(1+BETA*BETA)
           ALB = LOG(BETA+SB)
         ENDIF
         VXP = VXP * (-PFIVE + OPF * ALB / (BETA*SB))
         EXP_VAR = EXP_VAR *(ONE-OPF*((BETA*SB-ALB)/BETA**2)**2) 
       ENDIF
       VXF = CXF * VXP
       EXF = CXF * EXP_VAR
       DVXPDN = TRD * VXP / D
       DVXFDN = TRD * VXF / D

C      Correlation 
       IF (RS .GT. ONE) THEN  
         SQRS=SQRT(RS)
         TE = ONE+CON10*SQRS+CON11*RS
         BE = ONE+C1P053*SQRS+C3334*RS
         ECP = -(C2846/BE)
         VCP = ECP*TE/BE
         DTEDN = ((CON10 * SQRS *HALF) + CON11 * RS)*(-TRD/D)
         BE2 = BE * BE
         DBEDN = ((C1P053 * SQRS *HALF) + C3334 * RS)*(-TRD/D)
         DVCPDN = -(C2846/BE2)*(DTEDN - 2.0D0 * TE * DBEDN/BE)
         DECPDN = (C2846/BE2)*DBEDN
         TE = ONE+CON8*SQRS+CON9*RS
         BE = ONE+C1P398*SQRS+C2611*RS
         ECF = -(C1686/BE)
         VCF = ECF*TE/BE
         DTEDN = ((CON8 * SQRS * HALF) + CON9 * RS)*(-TRD/D)
         BE2 = BE * BE
         DBEDN = ((C1P398 * SQRS * HALF) + C2611 * RS)*(-TRD/D)
         DVCFDN = -(C1686/BE2)*(DTEDN - 2.0D0 * TE * DBEDN/BE)
         DECFDN = (C1686/BE2)*DBEDN
       ELSE
         RSLOG=LOG(RS)
         ECP=(C0622+C004*RS)*RSLOG-C096-C0232*RS
         VCP=(C0622+CON2*RS)*RSLOG-CON3-CON4*RS
         DVCPDN = (CON2*RS*RSLOG + (CON2-CON4)*RS + C0622)*(-TRD/D)
         DECPDN = (C004*RS*RSLOG + (C004-C0232)*RS + C0622)*(-TRD/D)
         ECF=(C0311+C0014*RS)*RSLOG-C0538-C0096*RS
         VCF=(C0311+CON5*RS)*RSLOG-CON6-CON7*RS
         DVCFDN = (CON5*RS*RSLOG + (CON5-CON7)*RS + C0311)*(-TRD/D)
         DECFDN = (C0014*RS*RSLOG + (C0014-C0096)*RS + C0311)*(-TRD/D)
       ENDIF

       ISP1 = 1
       ISP2 = 2

C      Find up and down potentials
       IF (NSP .EQ. 2) THEN
         EX    = EXP_VAR + FZ*(EXF-EXP_VAR)
         EC    = ECP + FZ*(ECF-ECP)
         VX(1) = VXP + FZ*(VXF-VXP) + (ONE-Z)*FZP*(EXF-EXP_VAR)
         VX(2) = VXP + FZ*(VXF-VXP) - (ONE+Z)*FZP*(EXF-EXP_VAR)
         VC(1) = VCP + FZ*(VCF-VCP) + (ONE-Z)*FZP*(ECF-ECP)
         VC(2) = VCP + FZ*(VCF-VCP) - (ONE+Z)*FZP*(ECF-ECP)

C        Derivatives of exchange potential respect the density

         DVXDN(ISP1,ISP1) =
     .             DVXPDN
     .              +  FZP*(VXF-VXP-EXF+EXP_VAR)*( 2.D0*D2/(D*D) )
     .              +  FZ*(DVXFDN-DVXPDN)+(1-Z)*FZP*(VXF-VXP)/(4.D0*D)
     .              +  (1-Z)*DFZPDN*(EXF-EXP_VAR)*( 2.D0*D2/(D*D) )
         DVXDN(ISP1,ISP2) =
     .                 DVXPDN
     .              +  FZP*(VXF-VXP-EXF+EXP_VAR)*(-2.D0*D1/(D*D) )
     .              +  FZ*(DVXFDN-DVXPDN)+(1-Z)*FZP*(VXF-VXP)/(4.D0*D)
     .              +  (1-Z)*DFZPDN*(EXF-EXP_VAR)*( -2.D0*D1/(D*D) )
         DVXDN(ISP2,ISP1) =
     .                 DVXPDN
     .              +  FZP*(VXF-VXP-EXF+EXP_VAR)*( 2.D0*D2/(D*D) )
     .              +  FZ*(DVXFDN-DVXPDN)-(1+Z)*FZP*(VXF-VXP)/(4.D0*D)
     .              -  (1+Z)*DFZPDN*(EXF-EXP_VAR)*( 2.D0*D2/(D*D) )
         DVXDN(ISP2,ISP2) =
     .                 DVXPDN
     .              +  FZP*(VXF-VXP-EXF+EXP_VAR)*(-2.D0*D1/(D*D) )
     .              +  FZ*(DVXFDN-DVXPDN)-(1+Z)*FZP*(VXF-VXP)/(4.D0*D)
     .              -  (1+Z)*DFZPDN*(EXF-EXP_VAR)*( -2.D0*D1/(D*D) )

C        Derivatives of correlation potential respect the density

         DVCDN(ISP1,ISP1) =
     .                DVCPDN
     .              + FZP*(VCF-VCP-ECF+ECP)*( 2.D0*D2/(D*D) )
     .              + FZ*(DVCFDN-DVCPDN)+ (1-Z)*FZP*(DECFDN-DECPDN)
     .              + (1-Z)*DFZPDN*(ECF-ECP)*( 2.D0*D2/(D*D) )
         DVCDN(ISP1,ISP2) =
     .                DVCPDN
     .              + FZP*(VCF-VCP-ECF+ECP)*(-2.D0*D1/(D*D) )
     .              + FZ*(DVCFDN-DVCPDN)+ (1-Z)*FZP*(DECFDN-DECPDN)
     .              + (1-Z)*DFZPDN*(ECF-ECP)*( -2.D0*D1/(D*D) )
         DVCDN(ISP2,ISP1) =
     .                DVCPDN
     .              + FZP*(VCF-VCP-ECF+ECP)*( 2.D0*D2/(D*D) )
     .              + FZ*(DVCFDN-DVCPDN)- (1+Z)*FZP*(DECFDN-DECPDN)
     .              - (1+Z)*DFZPDN*(ECF-ECP)*( 2.D0*D2/(D*D) )
         DVCDN(ISP2,ISP2) =
     .                DVCPDN
     .              + FZP*(VCF-VCP-ECF+ECP)*(-2.D0*D1/(D*D) )
     .              + FZ*(DVCFDN-DVCPDN)- (1+Z)*FZP*(DECFDN-DECPDN)
     .              - (1+Z)*DFZPDN*(ECF-ECP)*( -2.D0*D1/(D*D) )

       ELSE
         EX    = EXP_VAR
         EC    = ECP
         VX(1) = VXP
         VC(1) = VCP
         DVXDN(1,1) = DVXPDN
         DVCDN(1,1) = DVCPDN
       ENDIF

C      Change from Rydbergs to Hartrees
       EX = HALF * EX
       EC = HALF * EC
       DO 10 ISP = 1,NSP
         VX(ISP) = HALF * VX(ISP)
         VC(ISP) = HALF * VC(ISP)
         DO 5 ISP2 = 1,NSP
           DVXDN(ISP,ISP2) = HALF * DVXDN(ISP,ISP2)
           DVCDN(ISP,ISP2) = HALF * DVCDN(ISP,ISP2)
    5    CONTINUE
   10  CONTINUE
      END SUBROUTINE PZXC

      END MODULE m_ldaxc

