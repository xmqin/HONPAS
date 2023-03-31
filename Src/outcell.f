! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine outcell(cell) 

c *******************************************************************
c Writes lattice vectors and lattice parameters
c
c Written by E. Artacho, December 1997.
c ********* INPUT ***************************************************
c double precision cell(3,3): Lattice (supercell) vectors
c *******************************************************************

      use precision, only : dp
      use siesta_cml
      use units, only : Ang, pi

      implicit none

      real(dp), intent(in)   :: cell(3,3)

c Internal variables and arrays

      integer    :: iv, ix
      real(dp)   :: cellm(3), celang(3), volume, volcel

      external volcel

c Writing cell vectors

      write(6,'(/,a,3(/,a,3f12.6))')
     . 'outcell: Unit cell vectors (Ang):',
     . ('    ', (cell(ix,iv)/Ang,ix=1,3), iv =1,3)

c Cell-vector modules

      do iv = 1,3
        cellm(iv) = dot_product(cell(:,iv),cell(:,iv))
        cellm(iv) = sqrt(cellm(iv))
      enddo

c Cell-vector angles

      celang(1) = dot_product(cell(:,1),cell(:,2))
      celang(1) = acos(celang(1)/(cellm(1)*cellm(2)))*180._dp/pi
      celang(2) = dot_product(cell(:,1),cell(:,3))
      celang(2) = acos(celang(2)/(cellm(1)*cellm(3)))*180._dp/pi
      celang(3) = dot_product(cell(:,2),cell(:,3))
      celang(3) = acos(celang(3)/(cellm(2)*cellm(3)))*180._dp/pi

      write(6,'(/,a,3f12.6)')
     . 'outcell: Cell vector modules (Ang)   :',
     .          (cellm(iv)/Ang,iv=1,3)
      write(6,'(a,3f12.4)')
     . 'outcell: Cell angles (23,13,12) (deg):',
     .          (celang(iv),iv=3,1,-1)

c Cell volume

      volume = volcel( cell )
      write(6,'(a,f12.4)')
     .  'outcell: Cell volume (Ang**3)        :', volume/Ang**3

      if (cml_p) then
        call cmlAddCrystal(xf=mainXML, 
     .       title='Lattice Parameters', fmt='r6',
     .       alpha=celang(1), beta=celang(2), gamma=celang(3),
     .       a=cellm(1)/Ang, b=cellm(2)/Ang, c=cellm(3)/Ang )
      endif

      end subroutine outcell
