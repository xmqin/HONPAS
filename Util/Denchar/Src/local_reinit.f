
      SUBROUTINE REINIT( MAXO, MAXA, MAXUO, MAXNH, MAXNA, NSPIN, IDIMEN,
     .                   CHARGE, WAVES, STS )

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
      LOGICAL, INTENT(OUT) ::
     .  CHARGE, WAVES, STS

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
C LOGICAL STS            : Should STS simulation be done?
C **********************************************************************

C Internal variables --------------------------------------------------

      CHARACTER*33 PASTE

      CHARACTER*30 FNAME1

      CHARACTER 
     .  TYPRUN*2, TYPRUN_DEFECT*2


      INTEGER
     .  UNIT1 

      LOGICAL
     .  LEQI

      EXTERNAL
     .  IO_ASSIGN, IO_CLOSE, PASTE, LEQI

C Assign the name of the output file -----------------------------------
      slabel = FDF_STRING('SystemLabel','siesta')
      FNAME1 = trim(slabel) // '.DIM'

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
      STS    = FDF_BOOLEAN('Denchar.PlotSTS',.FALSE.)


      IF (.NOT. CHARGE .AND. .NOT. WAVES .AND. .NOT. STS) THEN
        WRITE(6,*)'Denchar.PlotCharge, Denchar.PlotWaveFunctions'
        WRITE(6,*)'and Denchar.PlotSTS are all .FALSE.'
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
