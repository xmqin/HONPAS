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
module m_history

!
! Implements a simple cyclic structure to hold a set of 
! values, recycling the oldest as new ones are added.
!

implicit none

integer, private, parameter :: dp=selected_real_kind(14,100)

type, public :: history_t
private
    integer :: size
    logical :: full
    integer :: n
    real(dp), dimension(:), pointer :: val
end type history_t

public :: initialize_history, add_to_history, reset_history
public :: n_history_elements, history_values, destroy_history

CONTAINS

!---------------------------------------------------
subroutine initialize_history(h,max_size)
type(history_t), intent(inout)  :: h
integer, intent(in)             :: max_size

h%size = max_size
allocate(h%val(max_size))
h%full = .false.
h%n = 0

end subroutine initialize_history

!---------------------------------------------------
subroutine reset_history(h)
type(history_t), intent(inout)  :: h

h%full = .false.
h%n = 0

end subroutine reset_history
!---------------------------------------------------

subroutine add_to_history(h,x)
type(history_t), intent(inout)   :: h
real(dp), intent(in)             :: x

integer :: i

if ( .not. h%full) then
   i = h%n + 1
   h%val(i) = x
   h%n = i
   if (i == h%size) h%full = .true.
else
   h%val = cshift(h%val,1)
   h%val(h%size) = x
endif

end subroutine add_to_history
!---------------------------------------------------

function n_history_elements(h) result (n)
type(history_t), intent(inout)       :: h
integer                              :: n

n = h%n

end function n_history_elements
!---------------------------------------------------

function history_values(h,nlast) result (a)
type(history_t), intent(in)          :: h
integer, intent(in)                  :: nlast
real(dp), dimension(nlast)           :: a

if (nlast > h%n) then ! STOP "requested more values than stored"
   print *, "requested, stored: ", nlast, h%n
   stop
endif
a(1:nlast) = h%val(h%n-nlast+1:h%n)

end function history_values
!
!---------------------------------------------------
subroutine destroy_history(h)
type(history_t), intent(inout)  :: h

if (associated(h%val)) deallocate(h%val)

end subroutine destroy_history

end module m_history
