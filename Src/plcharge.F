! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine plcharge( MAXO, MAXA, MAXUO, MAXNH, MAXNA, NSPIN, 
     .                     ISA, IPHORB, INDXUO, LASTO, 
     .                     CELL, NSC, XA, RMAXO, DATM )
C **********************************************************************
C Prepare the data files to plot charge density at the points of a plane 
C in real space.
C The information is to be read by the external DENCHAR
C program, to plot charge density contours in 2D
C
C Coded by J. Junquera 11/98
C Modified by DSP, July 1999
C Modified by J. Junquera 7/01
C Modified by J. Junquera 2/02
C **********************************************************************

      use precision
      use old_atmfuncs
      use files,         only : slabel, label_length
      use parallel,      only : Node, Nodes
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      integer, intent(in) ::
     .  MAXO, MAXA, MAXUO, MAXNH, MAXNA, NSPIN

      integer, intent(in) ::
     .  LASTO(0:MAXA), ISA(MAXA), IPHORB(MAXO), INDXUO(MAXO), NSC(3)

      real(dp), intent(in) ::
     .  CELL(3,3), XA(3,MAXA), RMAXO, DATM(MAXO)

C ****** INPUT *********************************************************
C INTEGER MAXO           : Maximum number of atomic orbitals in supercell
C INTEGER MAXA           : Maximum number of atoms in supercell
C INTEGER MAXUO          : Maximum number of atomic orbitals in unit cell.
C INTEGER MAXNH          : Maximum number
C                          of basis orbitals interacting, either directly
C                          or through a KB projector, with any orbital
C INTEGER MAXNA          : Maximum numbers of neighbour for any atom 
C INTEGER NSPIN          : Number of different spin polarizations
C                          Nspin = 1 => Non polarized. Nspin = 2 => Polarized
C INTEGER LASTO(0:MAXA)  : Last orbital of each atom in array iphorb
C INTEGER ISA(MAXA)      : Species index of each atom in the supercell
C INTEGER IPHORB(MAXO)   : Orbital index (within atom) of each orbital
C INTEGER INDXUO(MAXO)   : Equivalent orbital in unit cell
C INTEGER NSC(3)         : Num. of unit cells in each supercell direction
C REAL*8  CELL(3,3)      : Supercell vectors CELL(IXYZ,IVECT)
C                          (in bohrs)
C REAL*8  XA(3,MAXA)     : Atomic positions in cartesian coordinates
C                          (in bohrs)
C REAL*8  RMAXO          : Maximum range of basis orbitals
C REAL*8  DATM(MAXO)     : Occupations of basis orbitals in free atom
C **********************************************************************

C Internal variables ---------------------------------------------------

      character(len=label_length+4) :: fname1
      character(len=label_length+4) :: fname2

      integer
     .  UNIT1, UNIT2, IL, IA, J

#ifdef MPI
      integer MPIerror, MAXNHG
#endif

      external
     .  IO_ASSIGN, IO_CLOSE

#ifdef MPI
      call MPI_AllReduce(MAXNH,MAXNHG,1,MPI_integer,MPI_sum,
     .                   MPI_Comm_World,MPIerror)
#endif

C Write only if we are in the node 0 -----------------------------------
      if (Node .eq. 0) then
C Assign the name of the output file -----------------------------------
        FNAME1 = trim(slabel)//'.DIM'
        FNAME2 = trim(slabel)//'.PLD'

        call IO_ASSIGN(UNIT1)
        open( UNIT = UNIT1, FILE = FNAME1, FORM = 'UNFORMATTED',
     .        STATUS = 'UNKNOWN' )

        write(UNIT1) MAXA
        write(UNIT1) MAXO
        write(UNIT1) MAXUO 
        write(UNIT1) NSPIN
#ifdef MPI
        write(UNIT1) MAXNHG
#else
        write(UNIT1) MAXNH
#endif
        write(UNIT1) MAXNA

        call IO_CLOSE(UNIT1)

        call IO_ASSIGN(UNIT2)

        open( UNIT = UNIT2, FILE = FNAME2, FORM = 'UNFORMATTED',
     .        STATUS = 'UNKNOWN' )
C Dump the tables into a file ------------------------------------------

        write(UNIT2) RMAXO
C       write(6,*) RMAXO

        do IL = 1, MAXO
          write(UNIT2)IPHORB(IL), INDXUO(IL), DATM(IL)
        enddo

        do IA = 1, MAXA
          write(UNIT2)ISA(IA)
        enddo

        do IA = 0, MAXA
          write(UNIT2)LASTO(IA)
        enddo
        
        do IA = 1,3
          write(UNIT2)(CELL(J,IA),J=1,3)
        enddo

        write(UNIT2)(NSC(IA),IA=1,3)

        do IA = 1, MAXA
          write(UNIT2)(XA(J,IA),J=1,3)
        enddo

        call IO_CLOSE(UNIT2)
 
      endif

      return

      end
