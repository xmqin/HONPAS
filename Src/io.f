! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
c
c Copyright Alberto Garcia, 1996, 1997, 1998
c
c This module implements an interface to the FORTRAN logical unit
c system. Based on code by Richard Maine.
c
c
c Alberto Garcia, December 30, 1996
c Rewritten as a single subroutine 
c with multiple entry points, March 7, 1998
c Rewritten to use the module m_io. All calls are deferred to m_io,
c    Nick Papior, July, 2018
c---------------------------------------------------------------
c
      subroutine io
      use m_io, only: seterr_ => io_seterr
      use m_io, only: setout_ => io_setout
      use m_io, only: geterr_ => io_geterr
      use m_io, only: getout_ => io_getout
      use m_io, only: assign_ => io_assign
      use m_io, only: reserve_ => io_reserve
      use m_io, only: close_ => io_close
      use m_io, only: status_ => io_status

      implicit none

      integer :: unit
      
      entry io_seterr(unit)
      call seterr_(unit)
      return
      entry io_setout(unit)
      call setout_(unit)
      return

      entry io_geterr(unit)
      call geterr_(unit)
      return
      entry io_getout(unit)
      call getout_(unit)
      return

      entry io_assign(unit)
      call assign_(unit)
      return

      entry io_reserve(unit)
      call reserve_(unit)
      return

      entry io_close(unit)
      call close_(unit)
      return

      entry io_status
      call status_()
      return

      end subroutine
