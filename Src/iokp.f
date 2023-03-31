! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine iokp( nk, points, weight )
c *******************************************************************
c Saves k-points (only writing) Bohr^-1
c Emilio Artacho, Feb. 1999
c ********** INPUT **************************************************
c integer nk           : Number k points
c real*8  points(3,nk) : k-point coordinates
c real*8  weight(3,nk) : k-point weight
c *******************************************************************

      use files,     only : slabel, label_length
      use precision, only : dp

      implicit none

      integer, intent(in) :: nk
      real(dp), intent(in) :: points(3,nk), weight(nk)
      external :: io_assign, io_close

c Internal 
      character(len=label_length+3) :: fname
      integer :: ik, iu

      fname = trim(slabel) // '.KP'

      call io_assign( iu )
      open( iu, file=fname, form='formatted', status='unknown' )      

      write(iu,'(i10)') nk
      do ik = 1, nk
        write(iu,'(i10,3(tr1,e13.6),tr3,e12.6)')
     &      ik, points(:,ik), weight(ik)
      end do
      
      call io_close( iu )

      return
      end
