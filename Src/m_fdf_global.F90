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
module m_fdf_global

use precision, only: sp, dp
use fdf
use parallel, only: ionode
use m_mpi_utils, only: broadcast

implicit none

private

public :: fdf_global_get
interface fdf_global_get
   module procedure get_dp, get_int, get_bool
   module procedure get_sp, get_phys, get_str
end interface

private :: get_dp, get_int, get_bool
private :: get_sp, get_phys, get_str

CONTAINS

subroutine get_dp(x,name,default)
real(dp), intent(out)      :: x
real(dp), intent(in)       :: default
character(len=*), intent(in) :: name

if (ionode) x = fdf_double(name,default)
call Broadcast(x)

end subroutine get_dp

subroutine get_phys(x,name,default,unit)
real(dp), intent(out)      :: x
real(dp), intent(in)       :: default
character(len=*), intent(in) :: name
character(len=*), intent(in) :: unit

if (ionode) x = fdf_physical(name,default,unit)
call Broadcast(x)

end subroutine get_phys

subroutine get_sp(x,name,default)
real(sp), intent(out)      :: x
real(sp), intent(in)       :: default
character(len=*), intent(in) :: name

if (ionode) x = fdf_single(name,default)
call Broadcast(x)

end subroutine get_sp

subroutine get_str(x,name,default)
character(len=*), intent(out)      :: x
character(len=*), intent(in)       :: default
character(len=*), intent(in) :: name

if (ionode) x = fdf_string(name,default)
call Broadcast(x)

end subroutine get_str

subroutine get_int(i,name,default)
integer, intent(out)       :: i
integer, intent(in)        :: default
character(len=*), intent(in) :: name

if (ionode) i = fdf_integer(name,default)
call Broadcast(i)

end subroutine get_int

subroutine get_bool(b,name,default)
logical, intent(out)       :: b
logical, intent(in)        :: default
character(len=*), intent(in) :: name

if (ionode) b = fdf_boolean(name,default)
call Broadcast(b)

end subroutine get_bool

end module m_fdf_global
