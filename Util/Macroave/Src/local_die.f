! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!--------------------------------------------------
! Stand-alone 'die' routine for use by libraries and
! low-level modules.
!
! Each program using the module or library needs to
! provide a routine with the proper interface, but
! accomodating the needs and conventions of the program.
!
! Routines using this functionality should include
! the following
!
!     interface
!      subroutine die(str)
!      character(len=*), intent(in)  :: str
!      end subroutine die
!     end interface
!
!------------------------------------------------------

      subroutine die(str)

      character(len=*), intent(in)  :: str

      write(6,'(a)') trim(str)
      write(0,'(a)') trim(str)

      stop

      end subroutine die
