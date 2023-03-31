! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      SUBROUTINE REINIT( MAXO, MAXA, MAXUO, MAXNH, MAXNA, NSPIN, IDIMEN,
     .                   CHARGE, WAVES )

C **********************************************************************
C Read some variables from SIESTA in order to define
C the dimensions of some arrays in DENCHAR 
C
C Coded by J. Junquera 07/01
C Modifications for 3D and wavefunctions by P. Ordejon, June 2003
C **********************************************************************

      USE FDF
      use files, only: slabel

      IMPLICIT NONE

      INTEGER, INTENT(OUT) ::
     .  MAXO, MAXA, MAXUO, NSPIN, MAXNH, MAXNA, IDIMEN
      LOGICAL, INTENT(OUT) ::  CHARGE, WAVES

C **********************************************************************
C INTEGER MAXO           : Maximum number of atomic orbitals in supercell
C INTEGER MAXA           : Maximum number of atoms in supercell
C INTEGER MAXUO          : Number of atomic orbitals in unit cell.
C INTEGER MAXNH          : Maximum number
C                          of basis orbitals interacting, either directly
C                          or through a KB projector, with any orbital
C INTEGER MAXNA          : Maximum number of neighbour of any atom
C INTEGER NSPIN          : Number of different spin polarizations
C                          Nspin = 1 => Non polarized. Nspin = 2 => Polarized
C INTEGER IDIMEN         : Type of calculation: 2D or 3D
C LOGICAL CHARGE         : Should charge density be computed?
C LOGICAL WAVES          : Should wave functions be computed?
C **********************************************************************

C Internal variables --------------------------------------------------

      CHARACTER*30 FNAME1

      CHARACTER 
     .  TYPRUN*2, TYPRUN_DEFECT*2


      INTEGER
     .  UNIT1 

      EXTERNAL
     .  IO_ASSIGN, IO_CLOSE

C Assign the name of the output file -----------------------------------
      slabel = FDF_STRING('SystemLabel','siesta')
      FNAME1 = TRIM(slabel)// '.DIM'

      TYPRUN_DEFECT = '2D'
      TYPRUN = FDF_STRING('Denchar.TypeOfRun',TYPRUN_DEFECT)
      IF (LEQI(TYPRUN,'2D')) THEN
        IDIMEN  = 2
      ELSE IF (LEQI(TYPRUN,'3D')) THEN
        IDIMEN  = 3
      ELSE
        WRITE(6,'(A)')' readpla:  Wrong type of run; must be 2D or 3D  '
        STOP
      ENDIF

      CHARGE = FDF_BOOLEAN('Denchar.PlotCharge',.FALSE.)
      WAVES  = FDF_BOOLEAN('Denchar.PlotWaveFunctions',.FALSE.)

      IF (.NOT. CHARGE .AND. .NOT. WAVES) THEN
        WRITE(6,*)'Denchar.PlotCharge and Denchar.PlotWaveFunctions'
        WRITE(6,*)' are all .FALSE.'
        WRITE(6,*)'At least one of them should be .TRUE.'
        STOP
      ENDIF

      CALL IO_ASSIGN(UNIT1)
        OPEN ( UNIT = UNIT1, FILE = FNAME1, FORM = 'UNFORMATTED',
     .         STATUS = 'UNKNOWN' )

          READ(UNIT1)MAXA
          READ(UNIT1)MAXO
          READ(UNIT1)MAXUO 
          READ(UNIT1)NSPIN
          READ(UNIT1)MAXNH
          READ(UNIT1)MAXNA

      CALL IO_CLOSE(UNIT1)

      END
