! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

      module sys
!
!     Termination and messaging routines, MPI aware
!
      implicit none

      public :: die      ! Prints an error message and calls MPI_Abort
      private

      CONTAINS

!--------------------------------------------------
      subroutine die(str)

      character(len=*), intent(in), optional   :: str

      if (present(str)) then
         write(6,'(a)') trim(str)
         write(0,'(a)') trim(str)
      endif

      stop
      end subroutine die

      end module sys

! Stand-alone copy      
      subroutine die(str)

      character(len=*), intent(in)  :: str

         write(6,'(a)') trim(str)
         write(0,'(a)') trim(str)

      stop
      end subroutine die
      
