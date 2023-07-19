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
      module files
!
!     Contains the short system label, used to generate file names
!     slabel is currently set in reinit.
!
      integer, parameter, public                  :: label_length = 60
      character(len=label_length), save, public   :: slabel

      private

      end module files
