module m_wxml_elstack

use m_wxml_error

implicit none

private

!
! Simple stack to keep track of which elements have appeared so far
!
! Initial stack size:
integer, parameter, private            :: STACK_SIZE_INIT = 10
! Multiplier when stack is exceeded:
real, parameter, private            :: STACK_SIZE_MULT = 1.5

type, private :: elstack_item
  character(len=200) :: data
end type

type, public :: elstack_t
private
      integer                                   :: n_items
      type(elstack_item), pointer, dimension(:) :: stack
end type elstack_t

public  :: push_elstack, pop_elstack, init_elstack, reset_elstack, print_elstack
public  :: get_top_elstack, is_empty, get_elstack_signature
public  :: len

interface len
   module procedure number_of_items
end interface
private :: number_of_items

interface is_empty
      module procedure is_empty_elstack
end interface
private :: is_empty_elstack

CONTAINS

!-----------------------------------------------------------------
subroutine init_elstack(elstack)
  type(elstack_t), intent(inout)  :: elstack

  allocate(elstack%stack(STACK_SIZE_INIT))
  elstack%n_items = 0

end subroutine init_elstack

!-----------------------------------------------------------------
subroutine reset_elstack(elstack)
  type(elstack_t), intent(inout)  :: elstack

  deallocate(elstack%stack)
  call init_elstack(elstack)

end subroutine reset_elstack

!-----------------------------------------------------------------
subroutine resize_elstack(elstack)
  type(elstack_t), intent(inout)  :: elstack
  type(elstack_item), pointer, dimension(:) :: temp
  integer :: s

  s = size(elstack%stack)

  temp=>elstack%stack
  allocate(elstack%stack(nint(s*STACK_SIZE_MULT)))
  elstack%stack(:s) = temp
  deallocate(temp)

end subroutine resize_elstack

!-----------------------------------------------------------------
function is_empty_elstack(elstack) result(answer)
type(elstack_t), intent(in)  :: elstack
logical                    :: answer

answer = (elstack%n_items == 0)
end function is_empty_elstack

!-----------------------------------------------------------------
function number_of_items(elstack) result(n)
type(elstack_t), intent(in)  :: elstack
integer                      :: n

n = elstack%n_items
end function number_of_items

!-----------------------------------------------------------------
subroutine push_elstack(item,elstack)
character(len=*), intent(in)      :: item
type(elstack_t), intent(inout)  :: elstack

integer   :: n

n = elstack%n_items
if (n == size(elstack%stack)) then
  call resize_elstack(elstack)
endif
n = n + 1
elstack%stack(n)%data = item
elstack%n_items = n

end subroutine push_elstack

!-----------------------------------------------------------------
subroutine pop_elstack(elstack,item)
type(elstack_t), intent(inout)     :: elstack
character(len=*), intent(out)        :: item

integer   :: n

n = elstack%n_items
if (n == 0) then
      call wxml_error("Element stack empty")
endif
item = elstack%stack(n)%data
elstack%n_items = n - 1

end subroutine pop_elstack

!-----------------------------------------------------------------
subroutine get_top_elstack(elstack,item)
!
! Get the top element of the stack, *without popping it*.
!
type(elstack_t), intent(in)        :: elstack
character(len=*), intent(out)        :: item

integer   :: n

n = elstack%n_items
if (n == 0) then
  call wxml_error("Element stack empty")
endif
item = elstack%stack(n)%data

end subroutine get_top_elstack

!-----------------------------------------------------------------
subroutine print_elstack(elstack,unit)
type(elstack_t), intent(in)   :: elstack
integer, intent(in)           :: unit
integer   :: i

do i = elstack%n_items, 1, -1
  write(unit=unit,fmt=*) trim(elstack%stack(i)%data)
enddo

end subroutine print_elstack

!-------------------------------------------------------------
subroutine get_elstack_signature(elstack,string)
type(elstack_t), intent(in)   :: elstack
character(len=*), intent(out) :: string
integer   :: i, length, j

string = ' '
j = 0
do i = 1, elstack%n_items
   length = len_trim(elstack%stack(i)%data)
   string(j+1:j+1) = "/"
   j = j+1
   string(j+1:j+length) = trim(elstack%stack(i)%data)
   j = j + length
enddo

end subroutine get_elstack_signature

end module m_wxml_elstack
