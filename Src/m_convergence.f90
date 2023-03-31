! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_convergence

! A very simple derived type to monitor convergence of 
! a magnitude. 
! Just the last value added is kept, but a convergence check is
! performed first. Note that only absolute convergence is
! monitored. 
!
integer, parameter, private  :: dp = selected_real_kind(14,100)

type, public :: converger_t
 private
 integer     :: counter = 0
 logical     :: converged = .false.
 real(dp)    :: tolerance = 0.0_dp
 real(dp)    :: value
end type converger_t

public :: set_tolerance, is_converged, add_value, reset

CONTAINS

!--------------------------------------
subroutine reset(conv)
type(converger_t), intent(inout) :: conv

conv%converged = .false.
conv%counter   = 0

end subroutine reset

!--------------------------------------
subroutine set_tolerance(conv,tolerance)
type(converger_t), intent(inout) :: conv
real(dp), intent(in)             :: tolerance

conv%tolerance = tolerance
conv%counter   = 0

end subroutine set_tolerance

!--------------------------------------
function tolerance(conv)
type(converger_t), intent(in) :: conv
real(dp)                      :: tolerance

tolerance = conv%tolerance

end function tolerance

!--------------------------------------
function is_converged(conv) result(converged)
type(converger_t), intent(in) :: conv
logical                       :: converged

converged = conv%converged
end function is_converged

!--------------------------------------
subroutine add_value(conv,value)
type(converger_t), intent(inout) :: conv
real(dp), intent(in)             :: value

if (conv%counter > 0) then
    if (abs(conv%value-value) < conv%tolerance) then
       conv%converged = .true.
    endif
endif
conv%counter = conv%counter + 1
conv%value = value

end subroutine add_value
!--------------------------------------

end module m_convergence
