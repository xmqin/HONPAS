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
module m_wallclock
!
! A simple wall-time-stamper
! It should be called only from the master node
!

use m_walltime, only: wall_time

public :: wallclock

integer, parameter, private :: dp = selected_real_kind(14,200)

integer, private,  save    :: wt
logical, private,  save    :: first = .true.

private

CONTAINS

subroutine wallclock(str)
character(len=*), intent(in) :: str
!
real(dp) :: t

if (first) then
   call io_assign(wt)
   open(wt,file='CLOCK',form='formatted',status='unknown')
   first = .false.
endif

call wall_time(t)
write(wt,"(a,f18.3)") str, t

end subroutine wallclock

end module m_wallclock

