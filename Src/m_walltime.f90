! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!!@LICENSE
!
module m_walltime
!
! Implements a wall_time routine to return the wall time, instead of
! the cpu time. It tries to deal with possible wraparounds of the system
! clock's counter.
!
implicit none

public :: wall_time

integer, parameter, private :: dp = selected_real_kind(14,200)
integer, parameter, private :: i8 = selected_int_kind(12)

integer(i8), private, save     :: last_count
integer(i8), private, save     :: max_count
real(dp), private, save    :: last_time
real(dp), private, save    :: count_rate
logical, private, save     :: first = .true.


private

CONTAINS

subroutine wall_time(t)
real(dp), intent(out)  :: t

real(dp)  :: elapsed_time

integer(i8)       :: count_rate_int
integer(i8)       :: count

      if (first) then
         CALL system_clock (count_rate=count_rate_int)
         CALL system_clock (count_max=max_count)
         count_rate = real(count_rate_int,kind=dp)
         first = .false.
         CALL system_clock (last_count)
         t = 0.0_dp
         last_time = t
         RETURN
      endif

CALL system_clock (count)

! Watch out for wrap-around...
! This works if the system counter wrapped around just once...
! ... be liberal in your use of this routine.

if (count < last_count) then
   elapsed_time = (max_count - last_count + count)/count_rate
else
   elapsed_time = (count-last_count)/count_rate
endif

t = last_time + elapsed_time
last_time = t
last_count = count

end subroutine wall_time

end module m_walltime

