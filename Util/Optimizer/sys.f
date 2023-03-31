! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      module sys
      public :: die
      CONTAINS
      subroutine die(str)

      character(len=*), intent(in), optional   :: str

         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
      stop
      end subroutine die
      end module sys

