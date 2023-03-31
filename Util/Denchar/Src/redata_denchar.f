! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      SUBROUTINE REDATA_DENCHAR( MAXO, MAXA, MAXUO, MAXNH, NSPIN, 
     .                   ISA, IPHORB, INDXUO, LASTO,
     .                   CELL, NSC, XA, RMAXO, DATM )

C **********************************************************************
C Read the data files to plot charge density at the points of a plane 
C or at a 3D grid in real space.
C The information is written by the subroutine plcharge in SIESTA,
C when WriteDenchar option is set up to .true. in the fdf input data 
C file.
C
C Coded by J. Junquera 11/98
C Modified by DSP, July 1999
C Modified by J. Junquera, 7/01
C Modified by P. Ordejon to include 3D capabilities, June 2003
C **********************************************************************

      use precision
      USE FDF

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::
     .  MAXO, MAXA, MAXUO, NSPIN, MAXNH

      INTEGER, INTENT(OUT) ::
     .  LASTO(0:MAXA), ISA(MAXA), IPHORB(MAXO), INDXUO(MAXO), NSC(3)

      real(dp), INTENT(OUT) ::
     .  CELL(3,3), XA(3,MAXA), RMAXO, DATM(MAXO)


C **** INPUT ***********************************************************
C INTEGER MAXO           : Maximum number of atomic orbitals in supercell
C INTEGER MAXA           : Maximum number of atoms in supercell
C INTEGER MAXUO          : Maximum number of atomic orbitals in unit cell.
C INTEGER MAXNH          : Maximum number
C                          of basis orbitals interacting, either directly
C                          or through a KB projector, with any orbital
C INTEGER NSPIN          : Number of different spin polarizations
C                          Nspin = 1 => Non polarized. Nspin = 2 => Polarized
C **** OUTPUT **********************************************************
C INTEGER ISA(MAXA)      : Species index of each atom in the supercell
C INTEGER IPHORB(MAXO)   : Orbital index (within atom) of each orbital
C INTEGER INDXUO(MAXO)   : Equivalent orbital in unit cell
C INTEGER LASTO(0:MAXA)  : Last orbital of each atom in array iphorb
C REAL*8  CELL(3,3)      : Supercell vectors CELL(IXYZ,IVECT)
C                          (in bohrs)
C INTEGER NSC(3)         : Num. of unit cells in each supercell direction
C REAL*8  XA(3,MAXA)     : Atomic positions in cartesian coordinates
C                          (in bohrs)
C REAL*8  RMAXO          : Maximum range of basis orbitals
C REAL*8  DATM(MAXO)     : Occupations of basis orbitals in free atom
C **********************************************************************

C Internal variables ---------------------------------------------------

      CHARACTER*30
     .  SNAME, FNAME

      INTEGER
     .  UNIT1, IL, IA, J

      EXTERNAL
     .  IO_ASSIGN, IO_CLOSE


C Assign the name of the output file -----------------------------------
      SNAME = FDF_STRING('SystemLabel','siesta')
      FNAME = TRIM(sname)//'.PLD'

      CALL IO_ASSIGN(UNIT1)

      OPEN ( UNIT = UNIT1, FILE = FNAME, FORM = 'UNFORMATTED',
     .       STATUS = 'UNKNOWN' )
C Dump the tables into a file ------------------------------------------

        READ(UNIT1)RMAXO

        DO IL = 1, MAXO
          READ(UNIT1)IPHORB(IL), INDXUO(IL), DATM(IL)
        ENDDO

        DO IA = 1, MAXA
          READ(UNIT1)ISA(IA)
        ENDDO

        DO IA = 0, MAXA
          READ(UNIT1)LASTO(IA)
        ENDDO

        DO IA = 1,3
          READ(UNIT1)(CELL(J,IA),J=1,3)
        ENDDO

        READ(UNIT1)(NSC(IA),IA=1,3)

        DO IA = 1, MAXA
          READ(UNIT1)(XA(J,IA),J=1,3)
        ENDDO

      CALL IO_CLOSE(UNIT1)

      END SUBROUTINE REDATA_DENCHAR
