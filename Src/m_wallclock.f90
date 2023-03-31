! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_wallclock
!
! A simple wall-time-stamper
! It should be called only from the master node
!

use m_walltime, only: wall_time

public :: wallclock
public :: wallclock_shutdown

integer, parameter, private :: dp = selected_real_kind(14,200)

integer, private,  save    :: wt = -1
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

subroutine wallclock_shutdown()

  if ( wt >= 0 ) then
     call io_close(wt)
     ! Signal it has to be opened again.
     first = .true.
  end if

end subroutine wallclock_shutdown

end module m_wallclock

