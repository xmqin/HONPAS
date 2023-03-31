! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine pixmol(iza, xa, na )
c *******************************************************************
c Writes and accumulates coordinates to be animated by Xmol
c Written by E. Artacho. February 1999. Modified to open/close 2003
c ********* INPUT ***************************************************
c integer   iza(na)   : Atomic numbers of different atoms
c double    xa(3,na)  : Atom coordinates
c integer   na        : Number of atoms
c character slabel*20 : Label for file naming
c *******************************************************************

      use precision,      only: dp
      use periodic_table, only: symbol
      use files,          only: slabel, label_length
      use units, only: Ang

      implicit          none

      integer                             :: na
      integer                             :: iza(na)
      real(dp)                            :: xa(3,na)
      external          io_assign, io_close

c Internal variables and arrays
 
      character(len=label_length+4) :: fname
      integer :: i, ia
      integer :: unit
c -------------------------------------------------------------------

      character(len=2) :: sym

      fname = trim(slabel) //'.ANI'

      call io_assign(unit)
      open( unit, file=fname, form = 'formatted', position='append',
     .      status='unknown')

      write(unit,'(i5)') na
      write(unit,*)
      do ia = 1, na
         sym =  symbol(iza(ia))
         write(unit,'(a2,2x,3f12.6)') sym , (xa(i,ia)/Ang,i=1,3)
      enddo

      call io_close(unit)
      
      return
      end

