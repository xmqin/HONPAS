! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This is now a wrapper for the functionality in modules m_timer_tree
! and m_timer
!
! subroutine timer( prog, iOpt )
!
!  FINDS AND PRINTS THE CPU TIME SPENT IN DIFFERENT ROUTINES AND/OR
!   PROGRAM SECTIONS. IT MUST BE CALLED WITH IOPT=1 AT THE BEGINNING
!   OF EACH ROUTINE AND WITH IOPT=2 AT THE END OF IT.
!  ARGUMENTS:
!    PROG: INPUT ARBITRARY NAME FOR THE ROUTINE AND/OR PROGRAM SECTION
!    IOPT: INPUT OPTION PARAMETER:
!      IOPT = 0  => SET ALL TIMES TO ZERO (OPTIONAL)
!      IOPT = 1  => BEGIN COUNTING TIME FOR A ROUTINE
!      IOPT = 2  => STOP  COUNTING TIME FOR A ROUTINE
!      IOPT = 3  => PRINT TIME FOR A ROUTINE OR FOR ALL (IF PROG='ALL')
!  ROUTINE TIMES INCLUDE THAT SPENT IN THE ROUTINES THEY CALL
!  WRITTEN BY J.SOLER. JUL.2009
!=============================================================================
module timer_options
 logical, public :: use_tree_timer = .true.
 logical, public :: use_parallel_timer = .false.
end module timer_options

recursive subroutine timer( prog, iOpt )

! Module procedures used
  use sys,     only: die           ! Termination routine
  use parallel,only: node
  use timer_options, only: use_tree_timer, use_parallel_timer

  ! New 'tree-based' module by A. Garcia
  use m_timer_tree, only: timer_on   ! Start counting time
  use m_timer_tree, only: timer_off    ! Stop counting time
  use m_timer_tree, only: timer_report  ! Write all times

  ! Standard module by J. Soler
  use m_timer, only: timer_init    ! Initialize all times
  use m_timer, only: timer_start   ! Start counting time
  use m_timer, only: timer_stop    ! Stop counting time
  use m_timer, only: jms_timer_report=>timer_report  ! Write all times

#ifdef TRACING
  use extrae_module
  use extrae_eventllist
#endif

! Arguments
  implicit none
  character(len=*),intent(in):: prog   ! Name of program to time
  integer,         intent(in):: iOpt   ! Action option

#ifdef TRACING
  integer*8                           :: extrae_eventnumber

  if (iOpt==0) then
    extrae_maxEventNumber = 0

  else if (iOpt==1) then
    ! check if prog is in list, find eventnumber or add to list with new number
    extrae_eventnumber = getNumber(eventlist, prog)
    if (extrae_eventnumber == NOT_FOUND) then
      extrae_eventnumber = addToList(eventlist, prog)
    end if
    call extrae_eventandcounters(1000, extrae_eventnumber)

  else if (iOpt==2) then
    call extrae_eventandcounters(1000, 0)

  else if (iOpt==3) then
    ! write file with eventnumber-map
! if (Node == 0) write (*,*) 'extraeLIST write list '
!     if (Node == 0) then
!       call writeList(eventlist)
!     end if
!     call deleteList(eventlist)
  else
    call die('timer: ERROR: invalid iOpt value')
  end if

#endif

if (use_tree_timer) then

  ! New method, based on walltime in the master node only,
  ! and with a tree structure for the report

 if (Node == 0) then
  if (iOpt==0) then
    ! The timer is initialized in the first call to timer_on...
    !     call timer_init()
  else if (iOpt==1) then
    call timer_on( prog )
  else if (iOpt==2) then
    call timer_off( prog )
  else if (iOpt==3) then
    call timer_report(prog) 
  else
    call die('timer: ERROR: invalid iOpt value')
  end if
 end if

endif

if (use_parallel_timer) then
! Select action
  if (iOpt==0) then
     call timer_init()
  else if (iOpt==1) then
    call timer_start( prog )
  else if (iOpt==2) then
    call timer_stop( prog )
  else if (iOpt==3) then
    call jms_timer_report( prog, printNow=.true. )
  else
    call die('timer: ERROR: invalid iOpt value')
  end if

end if

end subroutine timer

