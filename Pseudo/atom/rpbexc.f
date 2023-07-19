
      SUBROUTINE RPBEXC( IREL, nspin, Dens, GDens,
     .                   EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )

C *********************************************************************
C Implements Hammer's RPBE Generalized-Gradient-Approximation (GGA).
C A revision of PBE (Perdew-Burke-Ernzerhof)
C Ref: Hammer, Hansen & Norskov, PRB 59, 7413 (1999) and
C J.P.Perdew, K.Burke & M.Ernzerhof, PRL 77, 3865 (1996)
C
C Written by M.V. Fernandez-Serra. March 2004. On the PBE routine of
C L.C.Balbas and J.M.Soler. December 1996. Version 0.5.
C ******** INPUT ******************************************************
C INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER nspin          : Number of spin polarizations (1 or 2)
C REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
C                           spin electron density (if nspin=2)
C REAL*8  GDens(3,nspin) : Total or spin density gradient
C ******** OUTPUT *****************************************************
C REAL*8  EX             : Exchange energy density
C REAL*8  EC             : Correlation energy density
C REAL*8  DEXDD(nspin)   : Partial derivative
C                           d(DensTot*Ex)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          exchange potential
C REAL*8  DECDD(nspin)   : Partial derivative
C                           d(DensTot*Ec)/dDens(ispin),
C                           where DensTot = Sum_ispin( Dens(ispin) )
C                          For a constant density, this is the
C                          correlation potential
C REAL*8  DEXDGD(3,nspin): Partial derivative
C                           d(DensTot*Ex)/d(GradDens(i,ispin))
C REAL*8  DECDGD(3,nspin): Partial derivative
C                           d(DensTot*Ec)/d(GradDens(i,ispin))
C ********* UNITS ****************************************************
C Lengths in Bohr
C Densities in electrons per Bohr**3
C Energies in Hartrees
C Gradient vectors in cartesian coordinates
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      implicit          none
      INTEGER           IREL, nspin
      double precision  Dens(nspin), DECDD(nspin), DECDGD(3,nspin),
     .                  DEXDD(nspin), DEXDGD(3,nspin), GDens(3,nspin)

C Internal variables
      INTEGER
     .  IS, IX

      double precision
     .  A, BETA, D(2), DADD, DECUDD, DENMIN,
     .  DF1DD, DF2DD, DF3DD, DF4DD, DF1DGD, DF3DGD, DF4DGD,
     .  DFCDD(2), DFCDGD(3,2), DFDD, DFDGD, DFXDD(2), DFXDGD(3,2),
     .  DHDD, DHDGD, DKFDD, DKSDD, DPDD, DPDZ, DRSDD,
     .  DS(2), DSDD, DSDGD, DT, DTDD, DTDGD, DZDD(2),
     .  EC, ECUNIF, EX, EXUNIF,
     .  F, F1, F2, F3, F4, FC, FX, FOUTHD,
     .  GAMMA, GD(3,2), GDM(2), GDMIN, GDMS, GDMT, GDS, GDT(3),
     .  H, HALF, KAPPA, KF, KFS, KS, MU, PHI, PI, RS, S,
     .  T, THD, THRHLF, TWO, TWOTHD, VCUNIF(2), VXUNIF(2), ZETA

C Lower bounds of density and its gradient to avoid divisions by zero
      PARAMETER ( DENMIN = 1.D-12 )
      PARAMETER ( GDMIN  = 1.D-12 )

C Fix some numerical parameters
      PARAMETER ( FOUTHD=4.D0/3.D0, HALF=0.5D0,
     .            THD=1.D0/3.D0, THRHLF=1.5D0,
     .            TWO=2.D0, TWOTHD=2.D0/3.D0 )

C Fix some more numerical constants
      PI = 4 * ATAN(1.D0)
      BETA = 0.066725D0
      GAMMA = (1 - LOG(TWO)) / PI**2
      MU = BETA * PI**2 / 3
      KAPPA = 0.804D0

C Translate density and its gradient to new variables
      IF (nspin .EQ. 1) THEN
        D(1) = HALF * Dens(1)
        D(2) = D(1)
        DT = MAX( DENMIN, Dens(1) )
        DO 10 IX = 1,3
          GD(IX,1) = HALF * GDens(IX,1)
          GD(IX,2) = GD(IX,1)
          GDT(IX) = GDens(IX,1)
   10   CONTINUE
      ELSE
        D(1) = Dens(1)
        D(2) = Dens(2)
        DT = MAX( DENMIN, Dens(1)+Dens(2) )
        DO 20 IX = 1,3
          GD(IX,1) = GDens(IX,1)
          GD(IX,2) = GDens(IX,2)
          GDT(IX) = GDens(IX,1) + GDens(IX,2)
   20   CONTINUE
      ENDIF
      GDM(1) = SQRT( GD(1,1)**2 + GD(2,1)**2 + GD(3,1)**2 )
      GDM(2) = SQRT( GD(1,2)**2 + GD(2,2)**2 + GD(3,2)**2 )
      GDMT   = SQRT( GDT(1)**2  + GDT(2)**2  + GDT(3)**2  )
      GDMT = MAX( GDMIN, GDMT )

C Find local correlation energy and potential
      CALL PW92C( 2, D, ECUNIF, VCUNIF )

C Find total correlation energy
      RS = ( 3 / (4*PI*DT) )**THD
      KF = (3 * PI**2 * DT)**THD
      KS = SQRT( 4 * KF / PI )
      ZETA = ( D(1) - D(2) ) / DT
      ZETA = MAX( -1.D0+DENMIN, ZETA )
      ZETA = MIN(  1.D0-DENMIN, ZETA )
      PHI = HALF * ( (1+ZETA)**TWOTHD + (1-ZETA)**TWOTHD )
      T = GDMT / (2 * PHI * KS * DT)
      F1 = ECUNIF / GAMMA / PHI**3
      F2 = EXP(-F1)
      A = BETA / GAMMA / (F2-1)
      F3 = T**2 + A * T**4
      F4 = BETA/GAMMA * F3 / (1 + A*F3)
      H = GAMMA * PHI**3 * LOG( 1 + F4 )
      FC = ECUNIF + H

C Find correlation energy derivatives
      DRSDD = - (THD * RS / DT)
      DKFDD =   THD * KF / DT
      DKSDD = HALF * KS * DKFDD / KF
      DZDD(1) =   1 / DT - ZETA / DT
      DZDD(2) = - (1 / DT) - ZETA / DT
      DPDZ = HALF * TWOTHD * ( 1/(1+ZETA)**THD - 1/(1-ZETA)**THD )
      DO 40 IS = 1,2
        DECUDD = ( VCUNIF(IS) - ECUNIF ) / DT
        DPDD = DPDZ * DZDD(IS)
        DTDD = (- T) * ( DPDD/PHI + DKSDD/KS + 1/DT )
        DF1DD = F1 * ( DECUDD/ECUNIF - 3*DPDD/PHI )
        DF2DD = (- F2) * DF1DD
        DADD = (- A) * DF2DD / (F2-1)
        DF3DD = (2*T + 4*A*T**3) * DTDD + DADD * T**4
        DF4DD = F4 * ( DF3DD/F3 - (DADD*F3+A*DF3DD)/(1+A*F3) )
        DHDD = 3 * H * DPDD / PHI
        DHDD = DHDD + GAMMA * PHI**3 * DF4DD / (1+F4)
        DFCDD(IS) = VCUNIF(IS) + H + DT * DHDD

        DO 30 IX = 1,3
          DTDGD = (T / GDMT) * GDT(IX) / GDMT
          DF3DGD = DTDGD * ( 2 * T + 4 * A * T**3 )
          DF4DGD = F4 * DF3DGD * ( 1/F3 - A/(1+A*F3) )
          DHDGD = GAMMA * PHI**3 * DF4DGD / (1+F4)
          DFCDGD(IX,IS) = DT * DHDGD
   30   CONTINUE
   40 CONTINUE

C Find exchange energy and potential
      FX = 0
      DO 60 IS = 1,2
        DS(IS)   = MAX( DENMIN, 2 * D(IS) )
        GDMS = MAX( GDMIN, 2 * GDM(IS) )
        KFS = (3 * PI**2 * DS(IS))**THD
        S = GDMS / (2 * KFS * DS(IS))
cea Hammer's RPBE (Hammer, Hansen & Norskov PRB 59 7413 (99)
cea     F1 = DEXP( - MU * S**2 / KAPPA)
cea     F = 1 + KAPPA * (1 - F1)
cea Following is standard PBE
cea     F1 = 1 + MU * S**2 / KAPPA
cea     F = 1 + KAPPA - KAPPA / F1
cea (If revPBE Zhang & Yang, PRL 80,890(1998),change PBE's KAPPA to 1.245)
        F1 = DEXP( - MU * S**2 / KAPPA)
        F = 1 + KAPPA * (1 - F1)

c       Note nspin=1 in call to exchng...

        CALL EXCHNG( IREL, 1, DS(IS), EXUNIF, VXUNIF(IS) )
        FX = FX + DS(IS) * EXUNIF * F

cMVFS   The derivatives of F  also need to be changed for Hammer's RPBE.
cMVFS   DF1DD = 2 * F1 * DSDD  * ( - MU * S / KAPPA)
cMVFS   DF1DGD= 2 * F1 * DSDGD * ( - MU * S / KAPPA)
cMVFS   DFDD  = -1 * KAPPA * DF1DD
cMVFS   DFDGD = -1 * KAPPA * DFDGD

        DKFDD = THD * KFS / DS(IS)
        DSDD = S * ( -(DKFDD/KFS) - 1/DS(IS) )
c       DF1DD = 2 * (F1-1) * DSDD / S
c       DFDD = KAPPA * DF1DD / F1**2
        DF1DD = 2* F1 * DSDD * ( - MU * S / KAPPA)
        DFDD = -1 * KAPPA * DF1DD
        DFXDD(IS) = VXUNIF(IS) * F + DS(IS) * EXUNIF * DFDD

        DO 50 IX = 1,3
          GDS = 2 * GD(IX,IS)
          DSDGD = (S / GDMS) * GDS / GDMS
c         DF1DGD = 2 * MU * S * DSDGD / KAPPA
c         DFDGD = KAPPA * DF1DGD / F1**2
          DF1DGD =2*F1 * DSDGD * ( - MU * S / KAPPA)
          DFDGD = -1 * KAPPA * DF1DGD
          DFXDGD(IX,IS) = DS(IS) * EXUNIF * DFDGD
   50   CONTINUE
   60 CONTINUE
      FX = HALF * FX / DT

C Set output arguments
      EX = FX
      EC = FC
      DO 90 IS = 1,nspin
        DEXDD(IS) = DFXDD(IS)
        DECDD(IS) = DFCDD(IS)
        DO 80 IX = 1,3
          DEXDGD(IX,IS) = DFXDGD(IX,IS)
          DECDGD(IX,IS) = DFCDGD(IX,IS)
   80   CONTINUE
   90 CONTINUE

      END
