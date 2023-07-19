      SUBROUTINE PW92XC( IREL, NSPIN, DENS, EPSX, EPSC, VX, VC )

C ********************************************************************
C Implements the Perdew-Wang'92 LDA/LSD exchange correlation
C Ref: J.P.Perdew & Y.Wang, PRB, 45, 13244 (1992)
C Written by L.C.Balbas and J.M.Soler. Dec'96. Version 0.5.
C ********* INPUT ****************************************************
C INTEGER IREL        : Relativistic-exchange switch (0=No, 1=Yes)
C INTEGER NSPIN       : Number of spin polarizations (1 or 2)
C REAL*8  DENS(NSPIN) : Local (spin) density
C ********* OUTPUT ***************************************************
C REAL*8  EPSX       : Exchange energy density
C REAL*8  EPSC       : Correlation energy density
C REAL*8  VX(NSPIN)  : Exchange (spin) potential
C REAL*8  VC(NSPIN)  : Correlation (spin) potential
C ********* UNITS ****************************************************
C Densities in electrons per Bohr**3
C Energies in Hartrees
C ********* ROUTINES CALLED ******************************************
C EXCHNG, PW92C
C ********************************************************************

      IMPLICIT          NONE
      INTEGER           IREL, NSPIN
      DOUBLE PRECISION  DENS(NSPIN), EPSX, EPSC, VC(NSPIN), VX(NSPIN)

      CALL EXCHNG( IREL, NSPIN, DENS, EPSX, VX )
      CALL PW92C( NSPIN, DENS, EPSC, VC )
      END
