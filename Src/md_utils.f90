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
module md_utils

  use precision, only:  dp

  integer, save, public  :: md_io_unit = 99
  logical, save, public  :: restart_with_mdinfo  = .false.
  logical, save,  :: file_opened  = .false.

  CONTAINS

    subroutine add_to_md_file(natoms,xa,va,cell,vcell,nose,nosedot)

      integer, intent(in)                                  :: natoms
      real(dp), dimension(3,natoms), intent(in)            :: xa, va
      real(dp), dimension(3,3), intent(in), optional       :: cell, vcell
      real(dp), intent(in), optional                       :: nose, nosedot

!     Here we can save x, xa, va for MD  (experimental)
!
        write(99,*) natoms
        do ia = 1,natoms
          write(99,'(i4,3f14.8,3x,3f14.8)')
     .      iza(ia),(xa(i,ia),i=1,3),(va(i,ia),i=1,3)
        enddo

        have_nose = present(nose)
        have_nosedot = present(nosedot)

        if (have_nose) then
           write(99,*) nose
        endif

        end subroutine add_to_md_file
