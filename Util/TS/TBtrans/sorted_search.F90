! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code segment has been fully created by:
! Nick R. Papior, 2018, nickpapior@gmail.com
! Please conctact the author, prior to re-using this code.

! A module to implement reduced search functions.
!
! In TranSiesta and TBtrans we are often searching for unique
! values in a sorted array.
! I.e. we will quite often search the same array for numerous data values
! either in consecutive order or not.
module sorted_search_m

  use intrinsic_missing, only: SFIND
  implicit none

  private

  public :: ssearch_t

  type ssearch_t
    !< The list of items to search in, this list *MUST* be 
    integer, pointer :: list(:)

    !< Index where we are to start from
    integer :: prev = 1
    !< Size of array `list`
    integer :: n = 0

  end type ssearch_t

  public :: ssearch_init
  public :: ssearch_find

contains

  subroutine ssearch_init(this, array)
    type(ssearch_t), intent(inout) :: this
    integer, target :: array(:)

    this%list => array
    this%prev = 1
    this%n = size(array)

  end subroutine ssearch_init

  function ssearch_find(this, val) result(idx)
    type(ssearch_t), intent(inout) :: this
    integer, intent(in) :: val
    integer :: idx

    if ( val == this%list(this%prev) ) then
      idx = this%prev
      this%prev = min(this%prev+1, this%n)
    else if ( val < this%list(this%prev) ) then
      idx = SFIND(this%list(1:this%prev), val)
      this%prev = max(1, idx)
    else if ( this%list(this%prev) < val ) then
      idx = SFIND(this%list(this%prev:), val)
      if ( idx > 0 ) then
        idx = this%prev - 1 + idx
        this%prev = min(this%n, idx+1)
      end if
    else
      idx = 0
    end if

  end function ssearch_find

end module sorted_search_m
    
