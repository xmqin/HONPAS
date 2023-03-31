module m_object_debug

  public :: set_object_debug_level, get_object_debug_level
  
  integer, public :: debug_level  = 0
  logical, public :: ok_to_print_debug_info = .false.

  CONTAINS

    subroutine set_object_debug_level(level)
      integer, intent(in) :: level
      
      debug_level = level
      if (debug_level > 0) then
         ok_to_print_debug_info = .true.
      endif

    end subroutine set_object_debug_level

!--------------------
    subroutine get_object_debug_level(level)
      integer, intent(out) :: level
      
      level = debug_level

    end subroutine get_object_debug_level

  end module m_object_debug
!
! Stand-alone function
!
  function print_debug_object_info() result(ok)
    use m_object_debug, only: ok_to_print_debug_info 
    logical :: ok

    ok = ok_to_print_debug_info

  end function print_debug_object_info
