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
module m_walltime
!
! Implements a wall_time routine to return the wall time, instead of
! the cpu time. It tries to deal with possible wraparounds of the system
! clock's counter.
!
public :: wall_time

integer, parameter, private :: dp = selected_real_kind(14,200)

integer, private, save     :: last_count
integer, private, save     :: max_count
real(dp), private, save    :: last_time
real(dp), private, save    :: count_rate
logical, private, save     :: first = .true.


private

CONTAINS

subroutine wall_time(t)
real(dp), intent(out)  :: t

real(dp)  :: elapsed_time

integer       :: count_rate_int
integer       :: count

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

