! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine coceri(iza, xa, cell, na, sname )
c *******************************************************************
c Writes coordinates in format to be read by CERIUS
c
c It implies atomic symbols, atomic coordinates (fractional format)
c and lattice parameters (modules (in Ang) and angles)
c
c Written by E. Artacho. December 1997.
c ********* INPUT ***************************************************
c integer   iza(na)   : Atomic numbers of different atoms
c double    xa(3,na)  : Atom coordinates (in Bohr)
c double    cell(3,3) : Lattice vectors (in Bohr)
c integer   na        : Number of atoms
c character slabel*20 : Label for file naming
c character sname*150 : Label for the title  
c ******************************************************************

      use precision,      only: dp
      use periodic_table, only: symbol
      use files,          only: slabel, label_length
      use units, only: Pi, Ang

      implicit          none

      character(len=150)            :: sname
      integer                       :: na
      integer                       :: iza(na)
      real(dp)                      :: cell(3,3)
      real(dp)                      :: xa(3,na)
      external          io_assign, io_close

c Internal variables and arrays
 
      character(len=label_length+4) :: fname
      character(len=2)              :: sym
      integer                       :: unit, ix, iv, i, ia
      real(dp)                      :: celang(3)
      real(dp)                      :: cellm(3)
      real(dp)                      :: recell(3,3)
      real(dp)                      :: xac(3)

!     automatic array

      real(dp)                      :: xap(3,na)

c Find lattice parameters out of lattice vectors: first modules:

      do iv = 1, 3
        cellm(iv) = 0.0d0
        do ix = 1, 3
          cellm(iv) = cellm(iv) + cell(ix,iv)*cell(ix,iv)
        enddo
        cellm(iv) = sqrt(cellm(iv))
      enddo

c and angles

      celang(1) = 0.d0
      do ix = 1, 3
        celang(1) = celang(1) + cell(ix,2)*cell(ix,3)
      enddo
      celang(1) = acos(celang(1)/(cellm(2)*cellm(3)))*180.d0/pi
      celang(2) = 0.d0
      do ix = 1, 3
        celang(2) = celang(2) + cell(ix,1)*cell(ix,3)
      enddo
      celang(2) = acos(celang(2)/(cellm(1)*cellm(3)))*180.d0/pi
      celang(3) = 0.d0
      do ix = 1, 3
        celang(3) = celang(3) + cell(ix,1)*cell(ix,2)
      enddo
      celang(3) = acos(celang(3)/(cellm(1)*cellm(2)))*180.d0/pi

c Obtain fractional coordinates (reclat inverts matrix)

      call reclat(cell, recell, 0)
      do ia = 1,na
        do ix = 1,3
          xac(ix) = xa(ix,ia)
        enddo
        do ix = 1,3
          xap(ix,ia) = recell(1,ix) * xac(1) +
     .                 recell(2,ix) * xac(2) +
     .                 recell(3,ix) * xac(3)
        enddo
      enddo


c Find file name

      fname = trim(slabel) // '.xtl'

      write(6,'(/,2a)')'coceri: Writing CERIUS coordinates into file ',
     .                  fname

      call io_assign(unit)
      open( unit, file=fname, form = 'formatted', status='unknown')
      rewind(unit)

c Write file

      write(unit,'(a,a70)') 'TITLE ', sname
      write(unit,'(a)')  'DIMENSION 3'
      write(unit,'(a,6f11.5)') 
     .          'CELL', (cellm(iv)/Ang,iv=1,3), (celang(i),i=1,3)
      write(unit,'(a)') 'SYMMETRY  NUMBER 1  LABEL P1'
      write(unit,'(3a)') 
     .       'SYM MAT  1.000000  0.000000  0.000000  0.000000',
     .              '  1.000000  0.000000  0.000000  0.000000',
     .              '  1.000000 0.0000 0.0000 0.0000'
      write(unit,'(/,a)') 'ATOMS'
      write(unit,'(2a)') 'NAME       X          Y          Z     ',
     .                             'CHARGE   TEMP    OCCUP   SCAT'
      do ia = 1, na
         sym =  symbol(iza(ia))
         write(unit,'(2x,a2,3f11.5,3f8.4,3x,a2)')
     .     sym, (xap(i,ia),i=1,3), 0.d0, 0.d0, 1.d0, sym
      enddo
      write(unit,'(a)') 'EOF'

      call io_close(unit)
      
      end

