! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module sys
!
!     Termination and messaging routines, MPI aware
!
      public :: die, bye, message

      CONTAINS

      subroutine message(str)

      character(len=*), intent(in), optional   :: str


         if (present(str)) then
            write(6,'(a)') trim(str)
         endif

      end subroutine message
!
!--------------------------------------------------
      subroutine die(str)

      character(len=*), intent(in), optional   :: str

         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
      STOP
      end subroutine die

!---------------------------------------------------------
      subroutine bye(str)

      character(len=*), intent(in), optional   :: str

         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
         write(6,'(a)') 'Requested End of Run. Bye!!'
!!                                       endif

      stop

      end subroutine bye

      end module sys

