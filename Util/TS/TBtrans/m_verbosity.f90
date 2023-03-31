! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_verbosity

  ! Default verbosity
  integer, save :: verbosity = 5

contains

  subroutine init_verbosity(fdf_name,default,verb)
    use fdf, only : fdf_get
    character(len=*), intent(in) :: fdf_name
    integer, intent(in) :: default
    integer, intent(out), optional :: verb
    
    integer :: v

    v = fdf_get(fdf_name,default)
    if ( present(verb) ) then
       verb = v
    else
       verbosity = v
    end if

  end subroutine init_verbosity

end module m_verbosity

    
    
