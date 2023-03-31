! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine typecell(cell,ctype,lv) 
c *******************************************************************
c Finds out if cell is cubic, fcc or bcc
c
c Written by P. Ordejon, March 1999.
C (from outcell of E. Artacho)
c ********* INPUT ***************************************************
c double precision cell(3,3): Lattice (supercell) vectors
c ********* OUTPUT **************************************************
c character        ctype:   : Type of cell (sc, fcc, bcc or none)
c double precision lv       : Lattice vector
c *******************************************************************

      use precision, only : dp

      implicit none

      real(dp)         :: cell(3,3), lv

      character(len=4) :: ctype

c Internal variables and arrays

      integer          :: iv, ix
      real(dp)         :: cellm(3), celang(3), tol, pi

      parameter (tol=1.d-4)

      data pi / 3.1415926d0 /

c Cell-vector modules

      do iv = 1, 3
         cellm(iv) = 0.d0
         do ix = 1, 3
            cellm(iv) = cellm(iv) + cell(ix,iv)*cell(ix,iv)
         enddo
         cellm(iv) = sqrt(cellm(iv))
      enddo

c Cell-vector angles

      celang(1) = 0.d0
      do ix = 1, 3
         celang(1) = celang(1) + cell(ix,1)*cell(ix,2)
      enddo
      celang(1) = acos(celang(1)/(cellm(1)*cellm(2)))*180.d0/pi
      celang(2) = 0.d0
      do ix = 1, 3
         celang(2) = celang(2) + cell(ix,1)*cell(ix,3)
      enddo
      celang(2) = acos(celang(2)/(cellm(1)*cellm(3)))*180.d0/pi
      celang(3) = 0.d0
      do ix = 1, 3
         celang(3) = celang(3) + cell(ix,2)*cell(ix,3)
      enddo
      celang(3) = acos(celang(3)/(cellm(2)*cellm(3)))*180.d0/pi

      ctype = 'none'

C Check if lattice vectors have same length, and return if not

      if (abs(cellm(1)-cellm(2)) .gt. tol) return
      if (abs(cellm(1)-cellm(3)) .gt. tol) return

C Check if angles are 90 deg, in which case the cell is cubic
      if ((abs(celang(1) - 90.0) .lt. tol) .and.
     .    (abs(celang(2) - 90.0) .lt. tol) .and.
     .    (abs(celang(3) - 90.0) .lt. tol)) then
        ctype = 'sc'
        lv = cellm(1)
        return
      endif

C Check if angles are 60 deg, in which case the cell is fcc
      if ((abs(celang(1) - 60.0) .lt. tol) .and.
     .    (abs(celang(2) - 60.0) .lt. tol) .and.
     .    (abs(celang(3) - 60.0) .lt. tol)) then
        ctype = 'fcc'
        lv = sqrt(2.0d0)*cellm(1)
        return
      endif
   
C Check if angles are 109.47122 deg, in which case the cell is bcc
      if ((abs(celang(1) - acos(-1./3.)*180.0/pi) .lt. tol) .and.
     .    (abs(celang(2) - acos(-1./3.)*180.0/pi) .lt. tol) .and.
     .    (abs(celang(3) - acos(-1./3.)*180.0/pi) .lt. tol)) then
        ctype = 'bcc'
        lv = (2.0d0/sqrt(3.0d0))*cellm(1)
        return
      endif
      
      return
      end
