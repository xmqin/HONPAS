! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      SUBROUTINE READSTS(VOLUME,EMIN, NE, DELTAE, ETA, XSTS, BFUNC,
     .                   IUNITCD,  ARMUNI)


C **********************************************************************
C Reads from the input file the data to calculate STS spectra
C
C Coded by P. Ordejon, July 2004
C **********************************************************************

      USE FDF

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) ::
     .  VOLUME

      INTEGER, INTENT(OUT) ::
     .  NE, IUNITCD, BFUNC
     
      DOUBLE PRECISION, INTENT(OUT) ::
     .  EMIN, DELTAE, ETA, XSTS(3), ARMUNI


C ****** INPUT  ********************************************************
C REAL*8 VOLUME          : Volume of the unit cell
C ****** OUTPUT ********************************************************
C REAL*8  EMIN           : Minimum energy in STS spectrum
C INTEGER NE             : Number of energy points
C REAL*8  DELTAE         : Step in energy plot
C REAL*8  ETA            : Lorenzian broadening of levels
C REAL*8  XSTS(3)        : Position where LDOS is computed
C INTEGER BFUNC          : Broadening function for STS spectra
C                          1 = Lorentzian, 0 = Gaussian
C INTEGER IUNITCD        : Units for the electron density
C                          IUNITCD = 1 => Ele/(bohr)**3
C                          IUNITCD = 2 => Ele/(Ang)**3
C                          IUNITCD = 3 => Ele/(unitcell)
C REAL*8  ARMUNI         : Conversion factors for the charge density
C **********************************************************************

C Internal variables ---------------------------------------------------
C INTEGER ISCALE         : Units for the atomic positions
C                          (ISCALE = 1 => Bohrs, ISCALE = 2 => Ang)
      CHARACTER 
     .  CPF*22, CPF_DEFECT*22,
     .  UCD*22, UCD_DEFECT*22,
     .  BROAD*22, BROAD_DEFECT*22

      INTEGER
     .  IUNIT, IX, JX, ISCALE

      DOUBLE PRECISION EMAX

      LOGICAL
     .  LEQI

      EXTERNAL
     .  LEQI
C ----

      CPF_DEFECT = 'Bohr'

      CPF = FDF_STRING('Denchar.CoorUnits',CPF_DEFECT)

      IF (LEQI(CPF,'bohr')) then
        ISCALE = 1
      ELSEIF (LEQI(CPF,'ang')) then
        ISCALE = 2
      ELSE
        WRITE(6,*)'readsts: ERROR Denchar.CoorUnits must be Ang or Bohr'
        STOP
      ENDIF

      UCD_DEFECT = 'Ele/bohr**3'
      UCD = FDF_STRING('Denchar.DensityUnits',UCD_DEFECT)
      IF (LEQI(UCD,'ele/bohr**3')) then
        IUNITCD = 1
      ELSEIF (LEQI(UCD,'ele/ang**3')) then
        IUNITCD = 2
      ELSEIF (LEQI(UCD,'ele/unitcell')) then
        IUNITCD = 3
      ELSE
       WRITE(6,'(A)')' readsts: ERROR   Wrong Option in Units of      '
       WRITE(6,'(A)')' readsts:  Charge Density                       '
       WRITE(6,'(A)')' readsts:  You must choose one of the following:' 
       WRITE(6,'(A)')' readsts:                                       '
       WRITE(6,'(A)')' readsts:      - Ele/bohr**3                    '
       WRITE(6,'(A)')' readsts:      - Ele/ang**3                     '
       WRITE(6,'(A)')' readsts:      - Ele/unitcell                   '
       STOP
      ENDIF


      BROAD_DEFECT = 'Lorentzian'

      BROAD = FDF_STRING('Denchar.STSBroadeningFunc',BROAD_DEFECT)

      IF (LEQI(BROAD,'lorentzian')) then
        BFUNC  = 1
      ELSEIF (LEQI(BROAD,'gaussian')) then
        BFUNC  = 2
      ELSE
        WRITE(6,*)'readsts: ERROR Denchar.STSBroadeningFunc must be'
        WRITE(6,*)'readsts:       Lorentzian or Gaussian'
        STOP
      ENDIF

      NE = FDF_INTEGER('Denchar.STSEnergyPoints',100)

      IF ( FDF_BLOCK('Denchar.STSPosition',IUNIT) ) THEN
        READ(IUNIT,*)(XSTS(IX),IX=1,3)
      ELSE
        WRITE(6,*)'To calculate STS spectra, you must'
        WRITE(6,*)'specify a position (Denchar.STSPosition)'
        WRITE(6,*)'in the input file'
        STOP
      ENDIF

      EMIN = FDF_PHYSICAL('Denchar.STSEmin',-1.0D0,'eV')
      EMAX = FDF_PHYSICAL('Denchar.STSEmax',1.0D0,'eV')
      ETA  = FDF_PHYSICAL('Denchar.STSEta',0.1D0,'eV')

      DELTAE = (EMAX-EMIN)/DFLOAT(NE-1)

C Scale points coordinates
C   Iscale = 1 => Do nothing
C   Iscale = 2 => Multiply by 1./0.529177 (Ang --> Bohr)

      IF (ISCALE .EQ. 2) THEN
        DO IX = 1,3
          XSTS(IX) = 1.D0 / 0.529177D0 * XSTS(IX)
        ENDDO
      ENDIF 

C Units of Charge Density
C   Iunitcd = 1 => Do nothing
C   Iunitcd = 2 => Multiply by (1.d0 / 0.529177d0) **3 (bohr**3 --> Ang**3)
C   Iunitcd = 3 => Multiply by volume unit cell (in bohrs**3) 

      IF (IUNITCD .EQ. 1) THEN
        ARMUNI = 1.D0
      ELSEIF( IUNITCD .EQ. 2 ) THEN
        ARMUNI = (1.D0 / 0.529177D0)**3 
      ELSEIF( IUNITCD .EQ. 3 ) THEN
        ARMUNI = VOLUME
      ENDIF

      END

