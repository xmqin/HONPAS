! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_timer_tree

  !
  ! Timing framework that knows about the tree structure among
  ! timed sections: an enclosing section has sub-sections as
  ! 'children', and this relationship is preserved in the report.
  ! Future extensions include the generation of a (possibly sorted)
  ! flat report, and the accumulation of timings for specific classes
  ! of sections.
  !
  ! This module keeps a "global_section" which is the parent of all
  ! first-level sections. Users should *not* handle this section.
  !
  ! Exported routines, all accepting a string 'secname' (section name)
  ! as argument:
  !
  !  timer_on    :   Flags the section as 'active' and (resumes) timing it.
  !
  !  timer_off   :   Removes the 'active' flag and computes the elapsed time.
  !                  The argument 'secname' is optional. If not present, or if
  !                  its value is 'all', all active sections are closed.
  !
  !  timer_report:   Produces a report. The argument 'secname' is optional. 
  !                  If not present, or if its value is 'all', a full report
  !                  rooted at the 'global_section' is produced.
  !                  A section-specific report covers the section itself and
  !                  its (first-level) children.
  !
  !  (It can be argued that the special meaning of 'all' is not desirable.
  !   The same functionality can be achieved by the presence or absence of
  !   the optional 'secname' argument)
  !  
  !  The timing is based by default on calls to wall_time()
  !  This behavior can be changed by setting the public module variable 
  !  'use_walltime' to .false.
  ! 
  !  64-bit integers are used in the system_clock counters to maximize
  !  precision and avoid wrap-around issues.
  !
  !  This version of the library can only deal with a single timing framework,
  !  as it uses global variables. This restriction could be removed in the
  !  future, if needed.
  !
  !  The user must provide an external routine 'die' with the interface
  !  specified below to properly abort the calling program in the event
  !  of an error (the routine could simply ignore timing errors...)
  !
  ! Alberto Garcia, January 2013-, re-using some pieces of code by Jose Soler
  ! 
  implicit none

  integer, parameter  :: dp = selected_real_kind(10,100)
  ! initial size of children array
  integer, parameter  :: INITIAL_N_CHILDREN = 3
  ! size increase in each extension of the children array
  integer, parameter  :: N_INC_CHILDREN = 3
  ! maximum length of the section name
  integer, parameter  :: NAME_LENGTH = 40

  ! Derived type to hold time data
  type times_t
     private
     character(len=NAME_LENGTH) :: name=' '  ! Name of timed section
     integer :: nCalls=0           ! Number of calls made to section
     real(dp):: totTime=0          ! Total time for this section
     real(dp):: startTime=0         ! Time upon start of current cycle
  end type times_t

  type section_t
     type(times_t) :: data
     logical       :: active=.false.     ! Is program time being counted?
     integer       :: nchildren = 0
     type(section_t), dimension(:), pointer :: child
     type(section_t), pointer               :: caller => null()
  end type section_t

  type(section_t), pointer :: global_section => null()
  type(section_t), pointer :: last_active => null()

  real(dp)       :: globaltime
  logical,public :: use_walltime = .true.

  type(section_t), pointer :: p
  type(times_t), pointer   :: pd

  real(dp)                 :: t_current
  real(dp)                 :: deltaTime
  character(len=256) :: msg


  public :: timer_on, timer_off, timer_report
  private

  interface
     subroutine die(str)
       character(len=*), intent(in), optional :: str
     end subroutine die
  end interface

CONTAINS

  !------------------------------------------------
  subroutine timer_on(secname)
    character(len=*)    :: secname

    integer :: loc

    ! Use an extra enclosing level for everything,
    ! so that multiple user "trees" can be supported

    if (.not. associated(global_section)) then
       allocate(global_section)
       p => global_section
       p%active = .true.
       pd => p%data
       pd%name = "global_section"
       pd%nCalls = pd%nCalls + 1
       last_active => p
       call current_time( t_current )
       pd%startTime = t_current
    endif

    ! Find proper place
    p => last_active
    if (p%nchildren==0) then
       allocate(p%child(INITIAL_N_CHILDREN))
    endif
    call child_index(secname,p%child,loc)
    if (loc == 0) then
       ! New child
       if (p%nchildren == size(p%child)) then
          call expand_array(p%child)
       endif
       p%nchildren = p%nchildren + 1
       !print *, "New child: " // trim(secname) // " of " // trim(p%data%name)
       p=>p%child(p%nchildren)
       p%caller => last_active
    else
       ! Existing child
       p => p%child(loc)
       p%caller => last_active
    endif

    p%active = .true.
    pd => p%data
    pd%name = secname
    pd%nCalls = pd%nCalls + 1
    last_active => p

    call current_time( t_current )
    pd%startTime = t_current

  end subroutine timer_on

  !------------------------------------------------
  recursive subroutine timer_off(secname)
    character(len=*), optional    :: secname

    logical :: close_all

    close_all = .true.
    if (present(secname)) then
       if (trim(secname) /= "all") then
          close_all = .false.
       endif
    endif

    if (close_all) then
       call timer_all_off()
    else
       p => last_active
       if (trim(p%data%name) /= trim(secname)) then
          write(msg,"(a,a,a,a)") "Wrong sequence in 'timer_off'. Last: ", &
               trim(p%data%name), ". This: ", trim(secname)
          call die(msg)
       endif

       if (.not. p%active) then
          write(msg,"(a,a)") trim(p%data%name), " not active!"
          call die(msg)
       endif
       p%active = .false.
       pd => p%data

       call current_time( t_current )  
       deltaTime = t_current - pd%startTime
       pd%totTime  = pd%totTime + deltaTime

       if (associated(p%caller)) then
          last_active => p%caller
       else
          call die('timer_tree: attempt to close global_section?')
       endif
    endif
  end subroutine timer_off

  !------------------------------------------------
  subroutine timer_all_off()
    !
    !  Closes all outstanding active sections
    !

    if (.not. associated(last_active)) RETURN

    p => last_active
    do
       if (associated(p,global_section)) exit
       call timer_off(p%data%name)
       p => p%caller
    enddo
  end subroutine timer_all_off

  !------------------------------------------------
  subroutine timer_report(secname,file_unit)
    character(len=*), optional    :: secname
    integer, intent(in), optional :: file_unit
    
    integer :: i, loc
    logical :: full

    full = .true.
    if (present(secname)) then
       if (trim(secname) /= "all") then
          full = .false.
       endif
    endif

    if (full) then
       ! Do a full report
       call timer_report_global()
    else
       p => last_active
       call child_index(secname,p%child,loc)
       if (loc == 0) then
          write(msg,"(a,a)") "Timing report requested for stale section: ", &
               trim(secname)
          call die(msg)
       endif
       p => p%child(loc)
 
       globaltime = p%data%totTime
       write(6,"(a20,T30,a6,a12,a8)") "Section","Calls","Walltime","% sect."
       call walk_tree(p,0,maxlevel=1)
    endif

  end subroutine timer_report
  !------------------------------------------------
  subroutine timer_report_global()

    integer :: i, js_lun
    type(times_t), pointer :: qd

    p => global_section
    ! Assign to the global section the sum of the times
    ! of its children
    globaltime = 0
    do i = 1, p%nchildren
       qd => p%child(i)%data
       globaltime = globaltime + qd%totTime
    enddo
    p%data%totTime = globaltime + 1.0e-6_dp

    ! Open JSON
    call get_unit_number(js_lun)
    open(unit=js_lun, file="time.json", form="formatted", &
        action="write",position="rewind")
    
    write(6,"(/,a20,T30,a6,a12,a8)") "Section","Calls","Walltime","%"
    ! Due to border logic we need this top-level wrapping
    write(js_lun,"(a)") "{"
    call walk_tree(p,0, js_lun=js_lun)
    write(js_lun,"(a)") "}"
    
    close(js_lun)
    
  end subroutine timer_report_global

  !------------------------------------------------
  recursive subroutine walk_tree(p,level,maxlevel,js_lun)
    type(section_t), intent(in),target  :: p
    integer, intent(in)          :: level
    integer, intent(in), optional:: maxlevel
    integer, intent(in), optional:: js_lun

    integer :: i
    character(len=64) fmtstr, fmt_json, fmt_json_head
    logical :: json_output

    ! Determine whether we should output JSON
    ! If the unit is there, we output to JSON as well
    json_output = present(js_lun)
    
    if (present(maxlevel)) then
       if (level > maxlevel) RETURN
    end if

    pd => p%data
    write(fmtstr,"(a,i0,a1,a)") "(", level+1, "x", ",a20,T30,i6,f12.3,f8.2)"
    
    write(6,fmtstr) pd%name, pd%nCalls, pd%totTime, 100*pd%totTime/globaltime
    if (json_output) then
       write(fmt_json,"(a,i0,a1,a)") "(",2*level+1,"x",",a,i0,a,f12.3,a,f8.2)"
       write(fmt_json_head,"(a,i0,a)") "(", 2*level+1, "x,a)"
       write(js_lun,fmt=fmt_json,advance="no") &
           '"' // trim(pd%name) // '": { "_calls": ', pd%nCalls,   &
                                      ', "_time": ', pd%totTime,   &
                                      ', "_%": ', 100*pd%totTime/globaltime
    end if
    if (p%nchildren /= 0) then
       if (json_output) write(js_lun,"(a)") ","
       do i=1,p%nchildren
          call walk_tree(p%child(i),level+1,maxlevel, js_lun)
          if (json_output) then
             if (i < p%nchildren) then
                write(js_lun,fmt="(a)") ','
             else
                write(js_lun,*)  ! just new line
             endif
          endif
       enddo
       ! End the parent section with an indented '}'
       if (json_output) write(js_lun,fmt=fmt_json_head,advance="no") "}"
    else
       if (json_output) write(js_lun,fmt="(a)",advance="no") "}"
    endif

  end subroutine walk_tree

  !------------------------------------------------
  subroutine child_index( name, childData, ichild ) 

    implicit none
    character(len=*),intent(in) :: name  ! Child name
    type(section_t), dimension(:) ,intent(in) :: childData
    integer,intent(out)                       :: ichild     ! 0 if not found

    integer :: kChild, nsize

    ichild = 0
    nsize = size(childData)

    if (nsize>0) then
       search: do kChild = 1,nsize
          if (trim(childData(kChild)%data%name) == trim(name)) then
             iChild = kChild
             exit search
          end if
       end do search
    end if

  end subroutine child_index

  subroutine expand_array(a)
    type(section_t), pointer :: a(:)

    integer :: n
    type(section_t), pointer :: tmp(:) => null()

    n = size(a)
    allocate(tmp(n))
    tmp(:) = a(:)
    deallocate(a)
    allocate(a(n+N_INC_CHILDREN))
    a(1:n) = tmp(1:n)
    deallocate(tmp)

  end subroutine expand_array

  !------------------------------------------------
  subroutine current_time(t)
    !
    ! CPU or walltime, depending on the setting of 'use_walltime'
    !
    use m_walltime, only: wall_time
    real(dp), intent(out) :: t

    real  :: treal   ! for call to cpu_time

    if (use_walltime) then
       call wall_time(t)
    else
       call cpu_time(treal)
       t = treal
    endif
  end subroutine current_time

subroutine get_unit_number(lun)
integer, intent(out) :: lun

logical :: used
integer :: iostat

do lun= 10, 99
   inquire(unit=lun, opened=used, iostat=iostat)
   if (iostat .ne. 0) used = .true.
   if (.not. used) return
enddo
call die("Cannot get unit for timer output")
end subroutine get_unit_number

end module m_timer_tree
