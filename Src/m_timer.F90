! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!!@LICENSE
!
!===============================================================================
! MODULE m_timer
! Provides utility routines for CPU timing
! Written by J.M.Soler. July 2009
!===============================================================================
! Used MPI routines and parameters
!  use mpi_siesta, only: MPI_AllGather
!  use mpi_siesta, only: MPI_Bcast
!  use mpi_siesta, only: MPI_Integer
!  use mpi_siesta, only: MPI_Double_Precision
!  use mpi_siesta, only: MPI_Character
!  use mpi_siesta, only: MPI_COMM_WORLD
!
! Other used module procedures
!  use sys,        only: die             ! Termination routine
!  use m_walltime, only: wall_time       ! Wall time routine
!  use m_io,       only: io_assign, io_close  ! Get and reserve an available IO unit
!  use m_io,       only: io_close        ! Free up IO unit after use
!  use moreParallelSubs, only: copyFile  ! Copies a file across nodes
!
! Used module parameters and variables
!  use precision,  only: dp              ! Double precision real kind
!  use parallel,   only: nNodes => nodes ! Number of processors
!  use parallel,   only: myNode => node  ! My processor's index
!===============================================================================
! Public procedures provided by this module
!  timer_init,  &! Initialize timing
!  timer_start, &! Start counting time
!  timer_stop,  &! Stop counting time
!  timer_get,   &! Get counted times
!  timer_report  ! Write times report on a file
!
! Public types provided by this module
!  times_t       ! Derived data type to hold counted times
!
! Public parameters, variables and arrays provided by this module:
!  none
!===============================================================================
! GENERAL USAGE:
!   program myProg
!     use m_timer, only: timer_init, timer_start, timer_stop, timer_report
!     call timer_init()
!     call timer_start('myProg')
!     call mySub1()
!     do i=1,100
!       call mySub2(i)
!     end do
!     call timer_stop('myProg')
!     call timer_report(file='myProg.times',printNow=.true.)
!   end program myProg
!   subroutine mySub1()
!     use m_timer, only: timer_get, timer_start, timer_stop
!     real*8:: myTime
!     call timer_start('mySub1')
!     do something
!     call timer_stop('mySub1')
!     call timer_get('mySub1',totTime=myTime)
!     print*, 'mySub1: my CPU time was ', myTime
!   end subroutine mySub1
!   subroutine mySub2(i)
!     use m_timer, only: timer_start, timer_stop
!     call timer_start('mySub2')
!     do something
!     call timer_stop('mySub2')
!   end subroutine mySub2
! - To obtain communication times in parallel execution, 
!   compile with -DMPI_TIMING
!===============================================================================
! SUBROUTINE timer_init()   
!   Initializes timing
! USAGE:
!   See GENERAL USAGE section
!   Note: Call timer_init only once in your entire program, unless 
!         you want to re-start all times
! ALGORITHMS:
!   Sets to zero the initial time and the number of timed programs
!===============================================================================
! SUBROUTINE timer_start( prog )   
!   Start counting time for a program or code section
! INPUT:
!   character(len=*) prog  ! Name of program of code section to be timed
!                          ! prog ideally should be at least of length 4.
! USAGE:
!   See GENERAL USAGE section
! BEHAVIOUR:
!   Stops with an error message if timer_start is called twice with the same
!   argument, without an intermediate call to timer_stop with that argument
! ALGORITHMS:
!   Calls intrinsic routine cpu_time and stores the present time, associated
!   to the prog name, for future use by timer_stop
!===============================================================================
! SUBROUTINE timer_stop( prog )   
!   Stops counting time for a program or code section
! INPUT:
!   character(len=*) prog  ! Name of program of code section to be timed
!                          ! prog ideally should be at least of length 4.
! USAGE:
!   See GENERAL USAGE section
! BEHAVIOUR:
!   Stops with an error message if timer_stop is called without a previous
!   call to timer_start with the same argument
! ALGORITHMS:
!   Calls intrinsic routine cpu_time and subtracts, from the present time,  
!   that previously stored by timer_start, and associated to the prog name.
!   Then it adds the difference to the total time associated to that name.
!   All those times are stored in an internal array of type(times_t), and
!   size equal to the number of timed programs.
!===============================================================================
! SUBROUTINE timer_get( prog, active, nCalls, &
!                       startTime, stopTime, totTime, commTime )  
!   Returns CPU times of a program or code section
! INPUT:
!   character(len=*):: prog  ! Name of program or code section
! OPTIONAL OUTPUT:
!   logical  active ! Is timing active for this prog?
!   integer  nCalls ! Number of calls to this prog?
!   real(dp) startTime ! Last start time for this prog
!   real(dp) stopTime  ! Last stop time for this prog
!   real(dp) totTime   ! Tot. accumulated time for prog
!   real(dp) commTime  ! Tot. communication time for prog
! USAGE:
!   See GENERAL USAGE section
!   Except with the 'active' argument, it must be called only after timer_stop
!   and before next call to timer_start, i.e. for non-active programs
! BEHAVIOUR:
!   Stops with an error message if called with arguments other than 'active'
!   for an active program or piece of code.
! ALGORITHMS:
!   Just copies the data stored for the prog name into the times argument
!===============================================================================
! SUBROUTINE timer_report( prog, unit, file, printNow, threshold )
!   Writes a report file of CPU times stored for one prog, or for all progs if
!   the prog argument is not present. Program times include those spent in the
!   subroutines that they call.
! OPTIONAL INPUT:
!   character(len=*):: prog      ! Name of program or code section
!   integer,        :: unit      ! IO file unit (used only in parallel exec.)
!   character(len=*):: file      ! IO file name (used only in parallel exec.)
!   logical         :: printNow  ! Print report now?
!   real(dp)        :: threshold ! Min. fractional time to be reported
! USAGE:
!   See GENERAL USAGE section
!   It can be called several times, with different arguments. In this case, 
!   the last values prevail and get stored for future calls.
! BEHAVIOUR:
! - If prog is not present, or prog=='all', it prints a full report, in the 
!   specified unit or file, of all the CPU times that have been profiled by 
!   timer_start--timer_stop. Otherwise, it prints a single line, in the
!   standard output, of the specified program or code section.
! - In serial execution, to keep backwards compatibility, the report is written
!   in the standard output, and arguments unit and file are not used.
! - In parallel execution, successive reports are overwritten, i.e. only the
!   last report remains written, unless different files are specified.
! - If unit is present and unit>0, argument file is not used neither in that
!   nor in future calls. If unit==0, the present or stored file is used.
!   If unit is not present and file is present, that file is used in that and
!   future calls.
! - If neither unit nor file are present, the (parallel) report is written on
!   file 'timer_report'
! - If printNow is not present, the report is NOT written.
! - If prog name is not found, it stops with an error message.
! - In the full report, program times are written in the order of the first
!   call to timer_start
! - In parallel execution, the reported times are those spent in the node with
!   the largest total CPU time, excluding communications.
! ALGORITHMS:
! - In parallel execution, the total calculation time (excluding commun.) of 
!   all nodes is first found using MPI_All_Gather, and the node with the 
!   largest value is designed the busyNode that will write the report.
! - Since the order in which prog times are stored may be different in 
!   different nodes, the busyNode broadcasts the name(s) of the prog(s) 
!   whose time(s) it wants to write. Each node then finds its time for that
!   prog and sends it to the busyNode, so that it can determine and print
!   the load balancing for that prog (specifically the min/max ratio of the
!   calculation time spent in that prog by the different nodes).
! - After the report is written by the busyNode, it is sent using copyFile to 
!   the root node, that writes it in its file system.
!===============================================================================

MODULE m_timer

! Used module procedures
  use sys,        only: die             ! Termination routine
  use m_io,       only: io_assign, io_close ! Get and reserve an available IO unit
  use m_walltime, only: wall_time       ! Wall time routine
  use moreParallelSubs, only: copyFile  ! Copies a file across nodes
  use parallel,   only: parallel_init   ! Initialize parallel variables
#ifdef MPI
  use mpi_siesta
#endif

! Used module parameters and variables
  use precision,  only: dp              ! Double precision real kind
  use parallel,   only: nNodes => nodes ! Number of processors
  use parallel,   only: myNode => node  ! My processor's index

  implicit none

! Public procedures provided by this module
PUBLIC:: &
  timer_init,  &! Initialize timing
  timer_start, &! Start counting time
  timer_stop,  &! Stop counting time
  timer_get,   &! Get counted times
  timer_report  ! Write times report on a file

! Public types, parameters, variables and arrays provided by this module:
! none

PRIVATE ! Nothing is declared public beyond this point

! Parameters
  character(len=*),parameter:: myName  = 'timer '
  character(len=*),parameter:: errHead = myName//'ERROR: '
  integer, parameter:: maxLength  = 32       ! Max prog name length
  integer, parameter:: maxProgs   = 500      ! Max number of timed progs
  real(dp),parameter:: minTotTime = 1.e-6_dp ! Min total CPU time in sec

! Derived type to hold time data
  type times_t
    private
    character(len=maxLength):: name=' '  ! Name of program or code section
    logical :: active=.false.     ! Is program time being counted?
    integer :: nCalls=0           ! Number of calls made to the program
    real(dp):: totTime=0          ! Total time for this program
    real(dp):: comTime=0          ! Communications time
    real(dp):: lastTime=0         ! Total time in last call
    real(dp):: lastComTime=0      ! Communications time in last call
  end type times_t

! Module variables and arrays
  real(dp),save :: minRepTime = 0.0_dp  ! Min reported CPU time fraction
  integer, save :: nProgs=0             ! Number of timed programs
  type(times_t), target, save:: progData(maxProgs) ! Holds data of timed progs
  character(len=150),save:: reportFile = 'timer_report'  ! File name for report
  integer, save :: reportUnit=0         ! IO unit for report file
  real(dp),save :: time0=0.0_dp         ! CPU time at timer initialization
  real(dp),save :: wallTime0=0.0_dp     ! Wall-time at timer initialization
  logical, save :: writingTimes=.false. ! Are we writing a times report?

CONTAINS

!===============================================================================

subroutine print_report( prog )   ! Write a report of counted times 
#ifdef _OPENMP
  use omp_lib
#endif
! Arguments
  implicit none
  character(len=*),intent(in):: prog   ! Name of program or code section

! Internal variables and arrays progsWriterNode
  character(len=*),parameter:: myReportFile = 'TIMES'
  real(dp),allocatable:: nodeCalTime(:), nodeComTime(:), nodeTotTime(:)
  character(len=maxLength):: progsWriterNode(maxProgs)=' ' ! Prog. names in
                                                           ! writer node
  character(len=maxLength):: progName
#ifndef _OPENMP
  real    :: treal                 ! Single precision to call cpu_time
#endif
  real(dp):: myCalTime, progCalTime, progComTime, progTotTime
  real(dp):: timeNow, totalCalTime, totalComTime, totalTime
  real(dp):: wallTime
  integer :: busyNode, totalComCalls, iProg, iu, jProg
  integer :: node, nProgsWriterNode, progCalls, writerNode
  logical :: found, opened
#ifdef MPI
  integer:: MPIerror, MPIstatus(MPI_STATUS_SIZE), MPItag
#endif

! If already writing a report, return to avoid an infinite recursion
  if (writingTimes) return
  writingTimes = .true.

! Find present CPU time and convert it to double precision
#ifdef _OPENMP
  timeNow = omp_get_wtime( )
#else
  call cpu_time( treal )
  timeNow = treal
#endif
  totalTime = timeNow - time0
  call wall_time( wallTime )
  wallTime = wallTime - wallTime0

! Do nothing if total time is too small
  if (totalTime<minTotTime) go to 999  ! Go to exit point

! Make sure that parallel variables are set
  call parallel_init()

! Select what program's time is to be written
  if (prog=='ALL' .or. prog=='all') then  ! Write all times

    ! Write different things for serial and parallel execution
    if (nNodes==1) then  ! Serial execution

      ! Write header
      write(6,'(/,a,f12.3)') 'timer: Elapsed wall time (sec) =', wallTime
      write(6,'(a)') 'timer: CPU execution times (sec):'
      write(6,'(/,a15,a9,2a12,a9)') &
        'Routine        ', 'Calls', 'Time/call', 'Tot.time', '%'

      ! Write program times
      do iProg = 1,nProgs
        progName    = progData(iProg)%name
        progCalls   = progData(iProg)%nCalls
        progTotTime = progData(iProg)%totTime
        write(6,'(a15,i9,2f12.3,f9.2)') &
          progName, progCalls, progTotTime/progCalls, &
                    progTotTime, progTotTime/totalTime*100
      end do ! iProg
      write(6,*) ' '

    else ! Parallel execution

      ! Find total communication calls and time
      totalComCalls = 0
      totalComTime = 0
      do iProg = 1,nProgs
        if (progData(iProg)%name(1:4)=='MPI_') then
          totalComCalls = totalComCalls + progData(iProg)%nCalls
          totalComTime  = totalComTime  + progData(iProg)%totTime
        end if
      end do ! iProg
      totalCalTime = totalTime - totalComTime

      ! Find time spent by each node
      allocate( nodeCalTime(0:nNodes-1), &
                nodeComTime(0:nNodes-1), &
                nodeTotTime(0:nNodes-1) )
#ifdef MPI
      call MPI_AllGather( totalTime,     1, MPI_Double_Precision, &
                          nodeTotTime, 1, MPI_Double_Precision, &
                          MPI_Comm_World, MPIerror )
      call MPI_AllGather( totalComTime,     1, MPI_Double_Precision, &
                          nodeComTime, 1, MPI_Double_Precision, &
                          MPI_Comm_World, MPIerror )
#else
      nodeTotTime(myNode) = totalTime
      nodeComTime(myNode) = totalComTime
#endif
      nodeCalTime(:) = nodeTotTime(:) - nodeComTime(:)

      ! Find which processor required most calculation time
      busyNode = 0
      do node = 1,nNodes-1
        if (nodeCalTime(node) > nodeCalTime(busyNode)) busyNode = node
      end do

      ! Set writing node
      writerNode = busyNode
! DEBUG
!      writerNode = 0
!      writerNode = 1
! END DEBUG

      ! Broadcast program names from writing node
      if (myNode==writerNode) then
        nProgsWriterNode = nProgs
        progsWriterNode(1:nProgs) = progData(1:nProgs)%name
      end if
#ifdef MPI
      call MPI_Bcast( nProgsWriterNode, 1, MPI_Integer, &
                      writerNode, MPI_Comm_World, MPIerror )
      call MPI_Bcast( progsWriterNode, nProgsWriterNode*maxLength, &
                      MPI_Character, writerNode, MPI_Comm_World, MPIerror )
#endif
! DEBUG
!      if (myNode==1) then
!        print'(a,i4,2a18)', &
!          ('timer: iProg,progs=',iProg,progs(iProg), &
!            progsWriterNode(iProg),iProg=1,nProgsWriterNode)
!      end if
! END DEBUG

      ! Write header
      if (myNode==writerNode) then
        ! Open local report file
        call io_assign( iu )
        open( unit=iu, file=myReportFile, form='formatted', status='unknown' )
         
        write(iu,'(/,a,i6)') &
         'timer: Number of nodes = ', nNodes
        write(iu,'(a,i6)') &
         'timer: Busiest calculating node was node = ', busyNode
        write(iu,'(a,i6)') &
         'timer: Times refer to node = ', writerNode
        write(iu,'(a,f12.3)') &
         'timer: Total elapsed wall-clock time (sec) = ', wallTime
        write(iu,'(a)') &
         'timer: CPU times (sec):'
        write(iu,'(a,3f12.3,f8.3)') &
          'Calc: Sum, Avge, myNode, Avg/Max =', &
          sum(nodeCalTime), sum(nodeCalTime)/nNodes, totalCalTime, &
          sum(nodeCalTime)/nNodes / maxval(nodeCalTime)
#ifdef MPI_TIMING
        write(iu,'(a,3f12.3,f8.3)') &
          'Comm: Sum, Avge, myNode, Avg/Max =', &
          sum(nodeComTime), sum(nodeComTime)/nNodes, totalComTime, &
          sum(nodeComTime)/nNodes / maxval(nodeComTime)
#else
	totalComTime = huge(1.0_dp) ! Avoid division by zero in prog table output
	write(iu,'(a)') 'No communications time available. Compile with -DMPI_TIMING'
#endif
        write(iu,'(a,3f12.3,f8.3)') &
          'Tot:  Sum, Avge, myNode, Avg/Max =', &
          sum(nodeTotTime), sum(nodeTotTime)/nNodes, totalTime, &
          sum(nodeTotTime)/nNodes / maxval(nodeTotTime)
        write(iu,'(a,e10.3,a)') &
          'Program times printed if > ', minRepTime,' of total'
        write(iu,'(/,(a15,a9,2(a12,a9),a8))') &
          'Program        ', 'Calls', 'Prg.com', ' Prg.com', &
          'Prg.tot', ' Prg.tot', ' Nod.avg',                  &
          '               ', '     ', '       ', '/Tot.com', &
          '       ', '/Tot.tot', '/Nod.max'
      endif ! (myNode==writerNode)

      ! Write program times
      do iProg = 1,nProgsWriterNode

        ! Find my jProg match to writing-node's iProg
        jProg = prog_index( progsWriterNode(iProg), found )

        ! Send program times to writing node
        if (found) then
          myCalTime = progData(jProg)%totTime - progData(jProg)%comTime
        else
          myCalTime = 0
          nProgs = nProgs-1  ! Delete new prog made by prog_index
        end if
#ifdef MPI
        ! For some unknown reason AllGather seems to be faster than Gather
        call MPI_AllGather( myCalTime,      1, MPI_Double_Precision, &
                            nodeCalTime(:), 1, MPI_Double_Precision, &
                            MPI_Comm_World, MPIerror)
#else
        nodeCalTime(writerNode) = myCalTime
#endif

        ! Write times for this program
        progName    = progData(iProg)%name
        progCalls   = progData(iProg)%nCalls
        progTotTime = progData(iProg)%totTime
        progComTime = progData(iProg)%comTime
        progCalTime = progTotTime - progComTime
        if (myNode==writerNode .and. progTotTime/totalTime>minRepTime) then
          if (maxval(nodeCalTime(:))==0._dp) then ! Only comm. time
            write(iu,'(a15,i9,2(f12.3,f9.4),f8.3)') &
              progName, progCalls, progComTime, progComTime/totalComTime, &
                        progTotTime, progTotTime/totalTime
          else                                    ! Some calc. time
            write(iu,'(a15,i9,2(f12.3,f9.4),f8.3)') &
              progName, progCalls, progComTime, progComTime/totalComTime, &
                        progTotTime, progTotTime/totalTime, &
              sum(nodeCalTime(:))/nNodes / maxval(nodeCalTime(:))
          endif ! (maxval(nodeCalTime(:))==0._dp)
        endif ! (myNode==writerNode)

      end do ! iProg

      ! Write total communications time
      if (myNode==writerNode) then
#ifdef MPI_TIMING
        write(iu,'(a15,i9,2(f12.3,f9.4))') &
         'MPI total      ', totalComCalls, &
          totalComTime, 1., totalComTime, totalComTime/totalTime
#else
	write(iu,'(a)') 'No communications time available. Compile with -DMPI_TIMING'
#endif
      endif ! (myNode==writerNode)


! DEBUG
!      ! Write time spent within final call to timer
!      call wall_time( wallTime1 )
!      call cpu_time( treal )
!      if (myNode==writerNode) write(iu,'(a,2f12.3,/)') &
!        'timer: Wall-clock and CPU times in final timer call =', &
!        wallTime1-wallTime, treal-time
! END DEBUG

      ! Write definitions
      if (myNode==writerNode) then
        write(iu,'(/,a)') 'Definitions:'
        write(iu,'(a)') &
          'Program: Name of a time-profiled routine (or code section)',     &
          'Calls:   Number of profiled calls to a program',                 &
          'Prg.com: CPU time in communications in one program in one node', &
          'Tot.com: CPU time in communications in all progs. in one node',  &
          'Prg.tot: total CPU time in one program in one node',             &
          'Tot.tot: total CPU time in all programs in one node',            &
          'Nod.avg: average calculation time in one program across nodes',  &
          'Nod.max: maximum calculation time in one program across nodes',  &
          'Calculation time: CPU time excluding communications', ' '

        ! Copy report file to the file system of root node
        call io_close( iu )
        
      endif ! (myNode==writerNode)

      if (reportUnit>0) then ! Append report to file already open
#ifdef MPI
        ! Find report file name in node 0 and send it to writer node
        if (myNode==0) then
          inquire( unit=reportUnit, opened=opened, name=reportFile )
          if (.not.opened) call die(errHead//' IO unit not open')
          call MPI_Send( reportFile, len(reportFile), MPI_Character, &
                         writerNode, 0, MPI_Comm_World, MPIerror )
        else if (myNode==writerNode) then
          call MPI_Recv( reportFile, len(reportFile), MPI_Character, &
                         0, MPItag, MPI_Comm_World, MPIstatus, MPIerror )
        end if ! (myNode==0)
#endif
        call copyFile( srcNode=writerNode, srcFile=myReportFile, &
                       dstNode=0, dstFile=reportFile, &
                       writeOption='append' )
      else ! (reportUnit==0) => use report file name
        call copyFile( srcNode=writerNode, srcFile=myReportFile, &
                       dstNode=0, dstFile=reportFile, &
                       writeOption='overwrite' )
      end if ! (reportUnit>0)

      deallocate( nodeCalTime, nodeComTime, nodeTotTime )

    endif ! (nNodes==1)

  else ! (prog/='all')  =>  write only one prog's time (only in root node)

    if (myNode==0) then
      ! Look for program name among those stored
      iProg = prog_index( prog, found )
      if (.not.found) call die(errHead//'not found prog = '//trim(prog))
      ! Write times
      write(6,'(a,a10,i6,f12.3,f7.2)') 'timer: Routine,Calls,Time,% = ', &
        prog, progData(iProg)%nCalls, progData(iProg)%totTime, &
        progData(iProg)%totTime/totalTime*100
    endif ! (myNode==0)

  end if ! (prog=='all')

! Exit point
999 continue
  writingTimes = .false.

end subroutine print_report

!===============================================================================

integer function prog_index( prog, found )   ! Get index of program in my data

  implicit none
  character(len=*),intent(in) :: prog  ! Program name
  logical,optional,intent(out):: found ! Was program already in my data?

! Internal variables
  integer,save :: iProg=0
  integer :: is, jProg, kProg
  logical :: progFound

! Check name length
  if (len(trim(prog)) > maxLength) then
    call die(errHead//'maxLength too small for prog = '//trim(prog))
  end if

! Look for program name among those stored. Look back and forth, from last 
! found, to find it quickly in sequential searches
  progFound = .false.
  if (nProgs>0) then
! DEBUG
!        if (myNode==0) print*, 'prog_index:'
! END DEBUG
    search: do jProg = 0,nProgs/2
      do is = -1,+1,2
        kProg = iProg + is*jProg
        kProg = modulo( kProg-1, nProgs ) + 1
! DEBUG
!        if (myNode==0) &
!          print'(a,i4,2(2x,a))', 'prog_index: iProg, prog(iProg), prog =', &
!          kProg, trim(progData(kProg)%name), trim(prog)
! END DEBUG
!        if (progData(kProg)%name == prog) then
        if (trim(progData(kProg)%name) == trim(prog)) then
          progFound = .true.
          iProg = kProg
          exit search
        end if
      end do ! is
    end do search
  end if ! (nProgs>0)

! Store new prog name
  if (.not.progFound) then
    nProgs = nProgs+1
    if (nProgs>maxProgs) call die(errHead//'parameter maxProgs too small')
    iProg = nProgs
    progData(iProg)%name = prog
  end if ! (.not.progFound)

! Set output
  prog_index = iProg
  if (present(found)) found = progFound

end function prog_index

!===============================================================================

subroutine timer_get( prog, active, nCalls, &    ! Returns times of a program
                      totTime, commTime, lastTime, lastCommTime )

  implicit none
  character(len=*), intent(in) :: prog   ! Name of program or code section
  logical, optional,intent(out):: active ! Is timing active for this prog?
  integer, optional,intent(out):: nCalls ! Number of calls to this prog?
  real(dp),optional,intent(out):: totTime   ! Tot. accumulated time for prog
  real(dp),optional,intent(out):: commTime  ! Tot. communication time for prog
  real(dp),optional,intent(out):: lastTime  ! Time in last call of this prog
  real(dp),optional,intent(out):: lastCommTime  ! Communic. time in last call

  integer:: iProg
  logical:: progActive

! Find program index
  iProg = prog_index( prog )

! Copy times data
  progActive = progData(iProg)%active
  if (present(active))         active       = progActive
  if (progActive) then
    call die(errHead//'timer_get: program '//trim(prog)//' is active')
  else
    if (present(nCalls))       nCalls       = progData(iProg)%nCalls
    if (present(totTime))      totTime      = progData(iProg)%totTime
    if (present(commTime))     commTime     = progData(iProg)%comTime
    if (present(lastTime))     lastTime     = progData(iProg)%lastTime
    if (present(lastCommTime)) lastCommTime = progData(iProg)%lastComTime
  end if

end subroutine timer_get

!===============================================================================

subroutine timer_init()   ! Initialize timing
#ifdef _OPENMP
  use omp_lib
#else
! Internal variables
  real    :: treal
#endif

  call wall_time( wallTime0 )
#ifdef _OPENMP
  time0 = omp_get_wtime( )
#else
  call cpu_time( treal )       ! Notice single precision
  time0 = treal
#endif
  nProgs = 0

! (Re)initialize data array
  progData(:)%name=' '           ! Name of program or code section
  progData(:)%active=.false.     ! Is program time being counted?
  progData(:)%nCalls=0           ! Number of calls made to the program
  progData(:)%totTime=0          ! Total time for this program
  progData(:)%comTime=0          ! Communications time
  progData(:)%lastTime=0         ! Total time in last call
  progData(:)%lastComTime=0      ! Communications time in last call

end subroutine timer_init

!===============================================================================

SUBROUTINE timer_report( prog, unit, file, printNow, threshold )

implicit none

character(len=*), optional, intent(in) :: prog      ! Program or code section
integer,          optional, intent(in) :: unit      ! IO file unit
character(len=*), optional, intent(in) :: file      ! IO file name
logical,          optional, intent(in) :: printNow  ! Print report now?
real(dp),         optional, intent(in) :: threshold ! Min. fract. time to report

if (present(unit)) then
  reportUnit = unit
else if (present(file)) then 
  reportUnit = 0
  reportFile = file
end if

if (present(threshold)) minRepTime = threshold

if (present(printNow)) then
  if (printNow) then
    if (present(prog)) then
      call print_report( prog )
    else
      call print_report( 'all' )
    end if
  end if
end if

END SUBROUTINE timer_report

! ==================================================================

subroutine timer_start( prog )   ! Start counting time for a program
#ifdef _OPENMP
  use omp_lib
#endif

  implicit none
  character(len=*),intent(in):: prog  ! Name of program of code section

! Internal variables
  integer :: iProg
#ifndef _OPENMP
  real    :: treal
#endif
  real(dp):: timeNow

! Do not change data if writing a report
  if (writingTimes) return

! Find present CPU time and convert it to double precision
#ifdef _OPENMP
  timeNow = omp_get_wtime( )
#else
  call cpu_time( treal )       ! Notice single precision
  timeNow = treal
#endif

! Find program index
  iProg = prog_index( prog )

! Check that program was not already active
  if (progData(iProg)%active) &
    call die(errHead//'timer_start: already active prog = '//trim(prog))

! Start timer for iProg. Use lastTime to store start time provisionally.
  progData(iProg)%active = .true.
  progData(iProg)%nCalls = progData(iProg)%nCalls + 1
  progData(iprog)%lastTime = timeNow
  progData(iprog)%lastComTime = 0

end subroutine timer_start

!===============================================================================

recursive subroutine timer_stop( prog )   ! Stop counting time for a program
#ifdef _OPENMP
  use omp_lib
#endif

  implicit none
  character(len=*),intent(in):: prog     ! Name of program of code section

! Internal variables
  integer :: iProg, jProg
#ifndef _OPENMP
  real    :: treal
#endif
  real(dp):: deltaTime, timeNow
  logical :: found

  if (trim(prog) == 'all') then
     call timer_all_stop()
     RETURN
  endif

! Do not change data if writing a report
  if (writingTimes) return

! Find present CPU time and convert it to double precision
#ifdef _OPENMP
  timeNow = omp_get_wtime( )
#else
  call cpu_time( treal )       ! Notice single precision
  timeNow = treal
#endif

! Find program index
  iProg = prog_index( prog, found )

! DEBUG
  if (.not.found) then
    call parallel_init()  ! Make sure that myNode is defined
    if (myNode==0) then
      do jProg = 1,nProgs
        print'(a,i6,2x,a)', 'timer_stop: iProg, prog =', &
        jProg, trim(progData(jProg)%name)
      end do
    end if
    call die(errHead//'timer_stop: not found prog = '//trim(prog))
  end if
! END DEBUG

! Check that program was already active
  if (.not.progData(iProg)%active) &
    call die(errHead//'timer_stop: not active prog = '//trim(prog))

! Stop timer for iProg
  progData(iProg)%active = .false.
  deltaTime = timeNow - progData(iprog)%lastTime
  progData(iprog)%lastTime = deltaTime
  progData(iProg)%totTime  = progData(iProg)%totTime + deltaTime

! Add communication time to all active programs
  if (len(prog) < 4) then
     ! Do nothing -- maybe warn the user that prog is not very meaningful
  else if (prog(1:4)=='MPI_' .or. prog(1:4)=='mpi_') then  
    ! We are dealing with a communication routine
    do jProg = 1,nProgs
      if (progData(jProg)%active) then
        progData(jProg)%comTime = progData(jProg)%comTime + deltaTime
        progData(jProg)%lastComTime = progData(jProg)%lastComTime + deltaTime
      end if
    end do
    ! Add comm. time also to iProg, since progData(iProg)%active=.false. by now
    progData(iProg)%comTime = progData(iProg)%comTime + deltaTime
    progData(iProg)%lastComTime = progData(iProg)%lastComTime + deltaTime
  end if

end subroutine timer_stop

subroutine timer_all_stop() 

integer :: i
type(times_t), pointer :: p

do i=1, maxProgs
   p => progData(i)
   if (p%active) then
      call timer_stop(p%name)
   endif
enddo
end subroutine timer_all_stop
   

END MODULE m_timer

