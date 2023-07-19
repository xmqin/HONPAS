
      SUBROUTINE REINITSTM( MAXO, MAXA, MAXUO, MAXNH, MAXNA, NSPIN)

C **********************************************************************
C Read some variables from SIESTA in order to define
C the dimensions of some arrays in STM
C
C Coded by P. Ordejon, November 2004
C Based on reinit by J. Junquera 07/01
C **********************************************************************

      USE FDF

      IMPLICIT NONE

      INTEGER, INTENT(OUT) ::
     .  MAXO, MAXA, MAXUO, NSPIN, MAXNH, MAXNA

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
C **********************************************************************

C Internal variables --------------------------------------------------

      CHARACTER*33 PASTE

      CHARACTER*30
     .  SNAME, FNAME1


      INTEGER
     .  UNIT1 

      LOGICAL
     .  LEQI

      EXTERNAL
     .  IO_ASSIGN, IO_CLOSE, PASTE, LEQI

C Assign the name of the output file -----------------------------------
      SNAME = FDF_STRING('SystemLabel','siesta')
      FNAME1 = PASTE(sname,'.DIM')

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

      RETURN

      END
