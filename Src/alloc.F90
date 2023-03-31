! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!!@LICENSE
!
! ==================================================================
! Allocation, reallocation, and deallocation utility routines
! Written by J.M.Soler. May 2000.
! ==================================================================
! SUBROUTINE alloc_default( old, new, restore, &
!                           copy, shrink, imin, routine )
!   Sets defaults for allocation
! INPUT (optional):
!   type(allocDefaults) restore : default settings to be restored
!   logical             copy    : Copy old array to new array?
!   logical             shrink  : Reduce array size?
!   integer             imin    : First index (typically 1 in Fortan,
!                                              0 in C)
!   character(len=*)    routine : Name of calling routine
! OUTPUT (optional):
!   type(allocDefaults) old     : default settings before the call
!   type(allocDefaults) new     : default settings after the call
! BEHAVIOR:
!   All these defaults can be superseeded by optional arguments in
!     each call to re_alloc.
!   Initial default values: copy    = .true.
!                           shrink  = .true.
!                           imin    = 1
!                           routine = 'unknown'
!   If restore is present together with any of copy, shrink, imin, or
!   routine, these are applied AFTER resetting the restore defaults.
! USAGE:
!   In order to restore the allocation defaults possibly set by the
! calling routine, the suggested construction is:
!   use alloc_module
!   type(allocDefaults) oldDefaults
!   call alloc_default( old=oldDefaults, routine=..., &
!                       copy=..., shrink=... )
!   call re_alloc(...)
!   call alloc_default( restore=oldDefaults )
! Notice that, if the restore call is skipped, the new defaults will
! stay in effect until a new call to alloc_dafault is made.
! ==================================================================
! SUBROUTINE alloc_report( level, unit, file, printNow, threshold )
!   Sets the output file for the allocation report
! INPUT (optional):
!   integer      :: level     : Level (detail) of report
!   integer      :: unit      : Output file unit
!   character*(*):: file      : Output file name
!   logical      :: printNow  : If present & true => print report now
!   real(dp)     :: threshold : Memory threshold (in bytes) to print
!                               the memory use of any given array 
! BEHAVIOR:
!   The detail/extent of the report increses with the value of level:
! level=0 : no report at all (the default)
! level=1 : only total memory peak and where it occurred
! level=2 : detailed report created but printed only upon request
! level=3 : detailed report printed at every new memory peak
! level=4 : print every individual reallocation or deallocation
!   If unit is present, alloc_report merely takes note of it for
! future use, assuming that it has been already open outside.
! In this case, file is not used.
!   If unit is absent, and file is present, a file with that
! name is open for future use.
!   If both arguments are absent, a file named 'alloc_report'
! is open for future use.
!   If alloc_report is called with printNow=.true. several times in
! a program, with the same unit or file argument, the subsequent 
! reports are written consecutively in the same file, each with a 
! time stamp header.
!   If threshold is not present, threshold=0 is assumed.
!   In parallel execution, the report sections that involve every
! reallocation (levels 1, 3, and 4) are written only by node 0.
! The section that is written upon request (level 2) is written
! only by the node with the highest peak of memory up to that time,
! but it contains a summary of the memory used by all other nodes.
!   In parallel execution, the nodes that share the same file
! system (e.g. different chip cores or NFS-connected nodes) write
! on the same file. Otherwise they write on files with the same name 
! in their local disks.
! ==================================================================
! SUBROUTINE re_alloc( array, [i1min,] i1max,
!                      [[i2min,] i2max, [[i3min,] i3max]],
!                      name, routine, copy, shrink )
! INPUT:
!   integer       :: i1min        : Lower bound of first dimension
!                                   If not present, it is fixed by 
!                                   the last call to alloc_default.
!                                   If present and the rank is 2(3),
!                                   then i2min(&i3min) must also be
!                                   present
!   integer       :: i1max        : Upper bound of first dimension
!   integer       :: i2min,i2max  : Bounds of second dimension, if
!                                   applicable
!   integer       :: i3min,i3max  : Bounds of third dimension, if appl.
!
! INPUT (optional):
!   character*(*) :: name         : Actual array name or a label for it
!   character*(*) :: routine      : Name of the calling routine
!                                     or routine section
!   logical       :: copy         : Save (copy) contents of old array 
!                                   to new array?
!   logical       :: shrink       : Reallocate if the new array bounds 
!                                   are contained within the old ones? 
!                                   If not present, copy and/or shrink
!                                      are fixed by the last call to
!                                      alloc_default. 
! INPUT/OUTPUT:
!   TYPE, pointer :: array : Array to be allocated or reallocated.
!                            Implemented types and ranks are:
!                              integer,    rank 1, 2, 3
!                              integer*8,  rank 1
!                              real*4,     rank 1, 2, 3, 4
!                              real*8,     rank 1, 2, 3, 4
!                              complex*16, rank 1, 2
!                              logical,    rank 1, 2, 3
!                              character(len=*), rank 1
! BEHAVIOR:
!   Pointers MUST NOT enter in an undefined state. Before using them
! for the first time, they must be nullified explicitly. Alternatively,
! in f95, they can be initialized as null() upon declaration.
!   If argument array is not associated on input, it is just allocated.
!   If array is associated and has the same bounds (or smaller bonds
! and shrink is false) nothing is done. Thus, it is perfectly safe and
! efficient to call re_alloc repeatedly without deallocating the array.
! However, subroutine dealloc is provided to eliminate large arrays
! when they are not needed.
!   In order to save (copy) the contents of the old array, the new array
! needs to be allocated before deallocating the old one. Thus, if the
! contents are not needed, or if reducing memory is a must, calling
! re_alloc with copy=.false. makes it to deallocate before allocating.
!   The elements that are not copied (because copy=.false. or because
! they are outside the bounds of the input array) return with value
! zero (integer and real), .false. (logical), or blank (character).
!   If imin>imax for any dimension, the array pointer returns
! associated to a zero-size array.
!   Besides allocating or reallocating the array, re_alloc keeps a count
! of memory usage (and prints a report in a file previously fixed by a
! call to alloc_report). Thus, it is not recommended to call re_alloc
! to reallocate an array not allocated by it the first time.
!   If name is not present, the memory count associated to the
! allocation/deallocation is still made, but the allocation report
! will account the array under the generic name 'unknown'.
!   The routine argument is NOT used to classify the counting of
! memory usage. The classification uses only the name argument.
! This is because a pointer may be allocated in one routine and
! deallocated in a different one (i.e. when it is used to return an
! array whose size is not known in advance). However, the routine
! argument is reported when alloc_report=4, and it is also used to
! report in which routine the memory peack occurs. If you want the
! routine name to be used for classification, you should include it
! as part of the name argument, like in name='matInvert '//'aux'.
! ==================================================================---
! SUBROUTINE de_alloc( array, name, routine )
! INPUT (optional):
!   character*(*) :: name    : Actual array name or a label for it
!   character*(*) :: routine : Name of the calling routine
!                                or routine section
! INPUT/OUTPUT:
!   TYPE, pointer :: array : Array be deallocated (same types and
!                            kinds as in re_alloc).
! BEHAVIOR:
!   Besides deallocating the array, re_alloc decreases the count of
! memory usage previously counted by re_alloc. Thus, dealloc should 
! not be called to deallocate an array not allocated by re_alloc.
! Equally, arrays allocated or reallocated by re_alloc should be 
! deallocated by dealloc.
! ==================================================================---
! SUBROUTINE alloc_count( delta_size, type, name, routine )
! INPUT:
!   integer         :: delta_size : +/-size(array)
!                                   + => allocation
!                                   - => deallocation
!   character       :: type       : 'I' => integer
!                                   'E' => integer*8
!                                   'R' => real*4
!                                   'D' => real*8
!                                   'L' => logical
!                                   'S' => character (string)
! INPUT (optional):
!   character(len=*) :: name      : Actual array name or a label for it
!   character(len=*) :: routine   : Name of the calling routine
!                                   or routine section
! USAGE:
!   integer,         allocatable:: intArray(:)
!   double precision,allocatable:: doubleArray(:)
!   complex,         allocatable:: complexArray(:)
!   allocate( intArray(n), doubleArray(n), complexArray(n) )
!   call alloc_count(   +n, 'I', 'intArray',  programName )
!   call alloc_count(   +n, 'D', 'doubleArray', programName )
!   call alloc_count( +2*n, 'R', 'complexArray', programName )
!   deallocate( intArray, doubleArray, complexArray )
!   call alloc_count(   -n, 'I', 'intArray',  programName )
!   call alloc_count(   -n, 'D', 'doubleArray', programName )
!   call alloc_count( -2*n, 'R', 'complexArray', programName )
!
! ==================================================================---

MODULE alloc

  use precision, only: sp        ! Single precision real type
  use precision, only: dp        ! Double precision real type
  use precision, only: i8b       ! 8-byte integer
  use parallel,  only: Node      ! My processor node index
  use parallel,  only: Nodes     ! Number of parallel processors
  use parallel,  only: ionode    ! Am I the I/O processor?
  use parallel,  only: parallel_init  ! Initialize parallel variables
  use sys,       only: die       ! Termination routine
  use m_io,      only: io_assign ! Get and reserve an available IO unit
#ifdef MPI
  use mpi_siesta
#endif

  implicit none

PUBLIC ::             &
  alloc_default,      &! Sets allocation defaults
  alloc_report,       &! Sets log report defaults
  re_alloc,           &! Allocation/reallocation
  de_alloc,           &! Deallocation
  alloc_count,        &! Memory counting for external allocs
  allocDefaults        ! Derived type to hold allocation defaults

PRIVATE      ! Nothing is declared public beyond this point

  interface de_alloc
    module procedure &
      dealloc_i1, dealloc_i2, dealloc_i3,             &
      dealloc_E1,                                     &
      dealloc_r1, dealloc_r2, dealloc_r3, dealloc_r4, &
      dealloc_d1, dealloc_d2, dealloc_d3, dealloc_d4, &
      dealloc_d5,                                     &
      dealloc_c1, dealloc_c2,                         &
      dealloc_z1, dealloc_z2, dealloc_z3, dealloc_z4, &
      dealloc_l1, dealloc_l2, dealloc_l3,             &
      dealloc_s1
  end interface

  interface re_alloc
    module procedure &
      realloc_i1, realloc_i2, realloc_i3,             &
      realloc_E1,                                     &
      realloc_r1, realloc_r2, realloc_r3, realloc_r4, &
      realloc_d1, realloc_d2, realloc_d3, realloc_d4, &
      realloc_d5,                                     &
      realloc_c1, realloc_c2,                         &
      realloc_z1, realloc_z2, realloc_z3, realloc_z4, & 
      realloc_l1, realloc_l2, realloc_l3,             &
      realloc_s1
!    module procedure & ! AG: Dangerous!!!
!      realloc_i1s, realloc_i2s, realloc_i3s,              &
!      realloc_r1s, realloc_r2s, realloc_r3s, realloc_r4s, &
!      realloc_d1s, realloc_d2s, realloc_d3s, realloc_d4s, &
!      realloc_l1s, realloc_l2s, realloc_l3s
  end interface

  ! Initial default values
  character(len=*), parameter :: &
    DEFAULT_NAME = 'unknown_name'         ! Array name default
  character(len=*), parameter :: &
    DEFAULT_ROUTINE = 'unknown_routine'   ! Routine name default
  integer, save ::               &
    REPORT_LEVEL = 0,            &! Level (detail) of allocation report
    REPORT_UNIT  = 0              ! Output file unit for report
  character(len=50), save ::     &
    REPORT_FILE = 'alloc_report'  ! Output file name for report
  real(dp), save ::              &
    REPORT_THRESHOLD = 0          ! Memory threshold (in bytes) to print
                                  ! the memory use of any given array 

  ! Derived type to hold allocation default options
  type allocDefaults
    private
    logical ::          copy = .true.              ! Copy default
    logical ::          shrink = .true.            ! Shrink default
    integer ::          imin = 1                   ! Imin default
    character(len=32):: routine = DEFAULT_ROUTINE  ! Routine name default
  end type allocDefaults

  ! Object to hold present allocation default options
  type(allocDefaults), save :: DEFAULT

  ! Internal auxiliary type for a binary tree
  type TREE
    character(len=80)   :: name  ! Name of an allocated array
    real(DP)            :: mem   ! Present memory use of the array
    real(DP)            :: max   ! Maximum memory use of the array
    real(DP)            :: peak  ! Memory use of the array during
                                 !   peak of total memory
    type(TREE), pointer :: left  ! Pointer to data of allocated arrays 
                                 !   preceeding in alphabetical order
    type(TREE), pointer :: right ! Pointer to data of allocated arrays 
                                 !   trailing in alphabetical order
  end type TREE

  ! Global variables used to store allocation data
  real(DP),   parameter     :: MBYTE = 1.e6_dp
  type(TREE), pointer, save :: REPORT_TREE
  real(DP),            save :: TOT_MEM  = 0._dp
  real(DP),            save :: PEAK_MEM = 0._dp
  character(len=80),   save :: PEAK_ARRAY = ' '
  character(len=32),   save :: PEAK_ROUTINE = ' '
  integer,             save :: MAX_LEN  = 0
  
  ! Other common variables
  integer :: IERR
  logical :: ASSOCIATED_ARRAY, NEEDS_ALLOC, NEEDS_COPY, NEEDS_DEALLOC

CONTAINS

! ==================================================================

SUBROUTINE alloc_default( old, new, restore,          &
                          routine, copy, shrink, imin )
implicit none
type(allocDefaults), optional, intent(out) :: old, new
type(allocDefaults), optional, intent(in)  :: restore
character(len=*),    optional, intent(in)  :: routine
logical,             optional, intent(in)  :: copy, shrink
integer,             optional, intent(in)  :: imin

if (present(old))     old = DEFAULT
if (present(restore)) DEFAULT = restore
if (present(copy))    DEFAULT%copy   = copy
if (present(shrink))  DEFAULT%shrink = shrink
if (present(imin))    DEFAULT%imin   = imin
if (present(routine)) DEFAULT%routine = routine
if (present(new))     new = DEFAULT

END SUBROUTINE alloc_default

! ==================================================================

SUBROUTINE alloc_report( level, unit, file, printNow, threshold, &
  shutdown)

implicit none

integer,          optional, intent(in) :: level, unit
character(len=*), optional, intent(in) :: file
logical,          optional, intent(in) :: printNow
real(dp),         optional, intent(in) :: threshold
logical,          optional, intent(in) :: shutdown

logical :: is_open

#ifdef MPI
integer :: MPIerror
#endif

if (present(level)) then
   REPORT_LEVEL = level
end if

if (REPORT_LEVEL <= 0) RETURN

if (node == 0) then
  if (present(unit)) then  ! Assume that unit has been open outside
    if (unit > 0) then
      REPORT_UNIT = unit
      if (present(file)) then
        REPORT_FILE = file
      else
        REPORT_FILE = 'unknown'
      end if
    end if
  else if (present(file)) then    ! If file is the same, do nothing
    if (file /= REPORT_FILE) then ! Check if file was open outside
      REPORT_FILE = file
      inquire( file=REPORT_FILE, opened=is_open, number=REPORT_UNIT )
      if (.not.is_open) then         ! Open new file
        call io_assign(REPORT_UNIT)
        open( REPORT_UNIT, file=REPORT_FILE, status='unknown')
        write(REPORT_UNIT,*) ' '  ! Overwrite previous reports
      end if
    end if
  else if (REPORT_UNIT==0) then   ! No unit has been open yet
    REPORT_FILE = 'alloc_report'
    call io_assign(REPORT_UNIT)
    open( REPORT_UNIT, file=REPORT_FILE, status='unknown')
    write(REPORT_UNIT,*) ' '      ! Overwrite previous reports
  end if
end if

#ifdef MPI
! Distribute information to other nodes and open REPORT_UNIT
! NP:
!   I am not too happy about this
!   This forces a certain unit to be Bcasted to the others.
!   If that unit is already open on the other nodes then you will have
!   this module and some other module write to the same file.
!   Perhaps we should just open the file from this node with
!   the node id appended?
call MPI_Bcast(REPORT_FILE,50,MPI_character,0,MPI_Comm_World,MPIerror)
! JMS: open file only in node 0
!if (node > 0) then
!  open( REPORT_UNIT, file=REPORT_FILE )
!end if
#endif

if (present(threshold)) REPORT_THRESHOLD = threshold

if (present(printNow)) then
  if (printNow) call print_report( )
end if

if (present(shutdown)) then
   if ( shutdown .and. REPORT_UNIT /= 0 ) then
      inquire( REPORT_UNIT, opened=is_open )
      if ( is_open ) call io_close(REPORT_UNIT)
   end if
end if

END SUBROUTINE alloc_report

! ==================================================================
! Integer array reallocs
! ==================================================================

SUBROUTINE realloc_i1( array, i1min, i1max, &
                       name, routine, copy, shrink )
! Arguments
implicit none
integer, dimension(:),      pointer    :: array
integer,                    intent(in) :: i1min
integer,                    intent(in) :: i1max
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
logical,          optional, intent(in) :: copy
logical,          optional, intent(in) :: shrink

! Internal variables and arrays
character, parameter           :: type='I'
integer, parameter             :: rank=1
integer, dimension(:), pointer :: old_array
integer, dimension(2,rank)     :: b, c, new_bounds, old_bounds

! Get old array bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array          ! Keep pointer to old array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if

! Copy new requested array bounds
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)

! Find if it is a new allocation or a true reallocation,
! and if the contents need to be copied (saved)
! Argument b returns common bounds
! Options routine also reads common variable ASSOCIATED_ARRAY,
! and it sets NEEDS_ALLOC, NEEDS_DEALLOC, and NEEDS_COPY
call options( b, c, old_bounds, new_bounds, copy, shrink )
! Deallocate old space
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine )
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if

! Allocate new space
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0
end if

! Copy contents and deallocate old space
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -size(old_array), type, name, routine )
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_i1

! ==================================================================
SUBROUTINE realloc_i2( array, i1min,i1max, i2min,i2max,       &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='I'
integer, parameter                     :: rank=2
integer, dimension(:,:),    pointer    :: array, old_array
integer,                    intent(in) :: i1min, i1max, i2min, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min /)
new_bounds(2,:) = (/ i1max, i2max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2)) =  &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_i2
! ==================================================================

SUBROUTINE realloc_i3( array, i1min,i1max, i2min,i2max, i3min,i3max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='I'
integer, parameter                     :: rank=3
integer, dimension(:,:,:),  pointer    :: array, old_array
integer,                    intent(in) :: i1min,i1max, i2min,i2max, &
                                          i3min,i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min /)
new_bounds(2,:) = (/ i1max, i2max, i3max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3)) =  &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_i3
! ==================================================================
SUBROUTINE realloc_E1( array, i1min, i1max, &
                       name, routine, copy, shrink )
! Arguments
implicit none
integer(i8b), dimension(:),    pointer    :: array
integer,                    intent(in) :: i1min
integer,                    intent(in) :: i1max
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
logical,          optional, intent(in) :: copy
logical,          optional, intent(in) :: shrink

! Internal variables and arrays
character, parameter                   :: type='I'
integer, parameter                     :: rank=1
integer(i8b), dimension(:), pointer       :: old_array
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds

! Get old array bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array          ! Keep pointer to old array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if

! Copy new requested array bounds
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)

! Find if it is a new allocation or a true reallocation,
! and if the contents need to be copied (saved)
! Argument b returns common bounds
! Options routine also reads common variable ASSOCIATED_ARRAY,
! and it sets NEEDS_ALLOC, NEEDS_DEALLOC, and NEEDS_COPY
call options( b, c, old_bounds, new_bounds, copy, shrink )

! Deallocate old space
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if

! Allocate new space
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0
end if

! Copy contents and deallocate old space
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if

END SUBROUTINE realloc_E1

! ==================================================================
! Single precision real array reallocs
! ==================================================================
SUBROUTINE realloc_r1( array, i1min, i1max,        &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='R'
integer, parameter                     :: rank=1
real(SP), dimension(:),     pointer    :: array, old_array
integer,                    intent(in) :: i1min, i1max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._sp
end if
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_r1
! ==================================================================
SUBROUTINE realloc_r2( array, i1min,i1max, i2min,i2max, &
                       name, routine, copy, shrink )
implicit none
character, parameter             :: type='R'
integer, parameter               :: rank=2
real(SP), dimension(:,:),   pointer    :: array, old_array
integer,                    intent(in) :: i1min, i1max, i2min, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min /)
new_bounds(2,:) = (/ i1max, i2max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._sp
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2)) =  &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_r2
! ==================================================================
SUBROUTINE realloc_r3( array, i1min,i1max, i2min,i2max, i3min,i3max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='R'
integer, parameter                     :: rank=3
real(SP), dimension(:,:,:), pointer    :: array, old_array
integer,                    intent(in) :: i1min,i1max, i2min,i2max, &
                                          i3min,i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array          ! Keep pointer to old array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min /)
new_bounds(2,:) = (/ i1max, i2max, i3max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._sp
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3)) =  &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_r3
! ==================================================================
SUBROUTINE realloc_r4( array, i1min,i1max, i2min,i2max, &
                              i3min,i3max, i4min,i4max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='R'
integer, parameter                     :: rank=4
real(SP), dimension(:,:,:,:), pointer  :: array, old_array
integer,                    intent(in) :: i1min,i1max, i2min,i2max, &
                                          i3min,i3max, i4min,i4max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min, i4min /)
new_bounds(2,:) = (/ i1max, i2max, i3max, i4max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2), &
                 b(1,3):b(2,3),b(1,4):b(2,4)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._sp
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4))= &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_r4
! ==================================================================
! Double precision real array reallocs
! ==================================================================
SUBROUTINE realloc_d1( array, i1min, i1max,        &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='D'
integer, parameter                     :: rank=1
real(DP), dimension(:),     pointer    :: array, old_array
integer,                    intent(in) :: i1min, i1max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_d1
! ==================================================================
SUBROUTINE realloc_d2( array, i1min,i1max, i2min,i2max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='D'
integer, parameter                     :: rank=2
real(DP), dimension(:,:),   pointer    :: array, old_array
integer,                    intent(in) :: i1min, i1max, i2min, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
integer                                :: i1, i2
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min /)
new_bounds(2,:) = (/ i1max, i2max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
!      array(c(1,1):c(2,1),c(1,2):c(2,2)) =  &
!  old_array(c(1,1):c(2,1),c(1,2):c(2,2))
  do i2 = c(1,2),c(2,2)
  do i1 = c(1,1),c(2,1)
    array(i1,i2) = old_array(i1,i2)
  end do
  end do
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_d2
! ==================================================================
SUBROUTINE realloc_d3( array, i1min,i1max, i2min,i2max, i3min,i3max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='D'
integer, parameter                     :: rank=3
real(DP), dimension(:,:,:), pointer    :: array, old_array
integer,                    intent(in) :: i1min,i1max, i2min,i2max, &
                                          i3min,i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
integer                                :: i1, i2, i3
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min /)
new_bounds(2,:) = (/ i1max, i2max, i3max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
!      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3)) =  &
!  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3))
  do i3 = c(1,3),c(2,3)
  do i2 = c(1,2),c(2,2)
  do i1 = c(1,1),c(2,1)
    array(i1,i2,i3) = old_array(i1,i2,i3)
  end do
  end do
  end do
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_d3
! ==================================================================
SUBROUTINE realloc_d4( array, i1min,i1max, i2min,i2max, &
                              i3min,i3max, i4min,i4max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='D'
integer, parameter                     :: rank=4
real(DP), dimension(:,:,:,:), pointer  :: array, old_array
integer,                    intent(in) :: i1min,i1max, i2min,i2max, &
                                          i3min,i3max, i4min,i4max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min, i4min /)
new_bounds(2,:) = (/ i1max, i2max, i3max, i4max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2), &
                 b(1,3):b(2,3),b(1,4):b(2,4)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4))= &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_d4
! ==================================================================
SUBROUTINE realloc_d5( array, i1min,i1max, i2min,i2max, &
                              i3min,i3max, i4min,i4max, &
                              i5min,i5max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                    :: type='D'
integer, parameter                      :: rank=5
real(DP), dimension(:,:,:,:,:), pointer :: array, old_array
integer,                    intent(in)  :: i1min,i1max, i2min,i2max, &
                                           i3min,i3max, i4min,i4max, &
                                           i5min,i5max
character(len=*), optional, intent(in)  :: name, routine
logical,          optional, intent(in)  :: copy, shrink
integer, dimension(2,rank)              :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min, i4min, i5min /)
new_bounds(2,:) = (/ i1max, i2max, i3max, i4max, i5max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine )
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2), &
                 b(1,3):b(2,3),b(1,4):b(2,4), &
                 b(1,5):b(2,5)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
  array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4), &
        c(1,5):c(2,5))= &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4), &
            c(1,5):c(2,5))
  call alloc_count( -size(old_array), type, name, routine )
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_d5

! ==================================================================
! Single precision complex array reallocs
! ==================================================================
SUBROUTINE realloc_c1( array, i1min, i1max,        &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='S'
integer, parameter                     :: rank=1
complex(SP), dimension(:),  pointer    :: array, old_array
integer,                    intent(in) :: i1min, i1max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -2*size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( 2*size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -2*size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_c1
! ==================================================================
SUBROUTINE realloc_c2( array, i1min,i1max, i2min,i2max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='S'
integer, parameter                     :: rank=2
complex(SP), dimension(:,:),  pointer  :: array, old_array
integer,                    intent(in) :: i1min, i1max, i2min, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
integer                                :: i1, i2
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min /)
new_bounds(2,:) = (/ i1max, i2max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -2*size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( 2*size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
!      array(c(1,1):c(2,1),c(1,2):c(2,2)) =  &
!  old_array(c(1,1):c(2,1),c(1,2):c(2,2))
  do i2 = c(1,2),c(2,2)
  do i1 = c(1,1),c(2,1)
    array(i1,i2) = old_array(i1,i2)
  end do
  end do
  call alloc_count( -2*size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_c2

! ==================================================================
! Double precision complex array reallocs
! ==================================================================
SUBROUTINE realloc_z1( array, i1min, i1max,        &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='D'
integer, parameter                     :: rank=1
complex(DP), dimension(:),  pointer    :: array, old_array
integer,                    intent(in) :: i1min, i1max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -2*size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( 2*size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -2*size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_z1
! ==================================================================
SUBROUTINE realloc_z2( array, i1min,i1max, i2min,i2max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='D'
integer, parameter                     :: rank=2
complex(DP), dimension(:,:),  pointer  :: array, old_array
integer,                    intent(in) :: i1min, i1max, i2min, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
integer                                :: i1, i2
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min /)
new_bounds(2,:) = (/ i1max, i2max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -2*size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( 2*size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
!      array(c(1,1):c(2,1),c(1,2):c(2,2)) =  &
!  old_array(c(1,1):c(2,1),c(1,2):c(2,2))
  do i2 = c(1,2),c(2,2)
  do i1 = c(1,1),c(2,1)
    array(i1,i2) = old_array(i1,i2)
  end do
  end do
  call alloc_count( -2*size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_z2
SUBROUTINE realloc_z3( array, i1min,i1max, i2min,i2max, i3min,i3max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='D'
integer, parameter                     :: rank=3
complex(DP), dimension(:,:,:), pointer :: array, old_array
integer,                    intent(in) :: i1min, i1max, i2min, i2max
integer,                    intent(in) :: i3min, i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
integer                                :: i1, i2, i3
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min /)
new_bounds(2,:) = (/ i1max, i2max, i3max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -2*size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( 2*size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
!      array(c(1,1):c(2,1),c(1,2):c(2,2)) =  &
!  old_array(c(1,1):c(2,1),c(1,2):c(2,2))
  do i3 = c(1,3),c(2,3)
  do i2 = c(1,2),c(2,2)
  do i1 = c(1,1),c(2,1)
    array(i1,i2,i3) = old_array(i1,i2,i3)
  end do
  end do
  end do
  call alloc_count( -2*size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_z3
SUBROUTINE realloc_z4( array, i1min,i1max, i2min,i2max, &
                              i3min,i3max, i4min,i4max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                       :: type='D'
integer, parameter                         :: rank=4
complex(DP), dimension(:,:,:,:),  pointer  :: array, old_array
integer,                    intent(in)     :: i1min, i1max, i2min, i2max, &
                                              i3min, i3max, i4min, i4max
character(len=*), optional, intent(in)     :: name, routine
logical,          optional, intent(in)     :: copy, shrink
integer, dimension(2,rank)                 :: b, c, new_bounds, old_bounds
integer                                    :: i1, i2, i3, i4
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min, i4min /)
new_bounds(2,:) = (/ i1max, i2max, i3max, i4max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -2*size(old_array), type, name, routine )
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3),b(1,4):b(2,4)), &
            stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( 2*size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
!      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4)) =  &
!  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4))
  do i4 = c(1,4),c(2,4)
  do i3 = c(1,3),c(2,3)
  do i2 = c(1,2),c(2,2)
  do i1 = c(1,1),c(2,1)
    array(i1,i2,i3,i4) = old_array(i1,i2,i3,i4)
  end do
  end do
  end do
  end do
  call alloc_count( -2*size(old_array), type, name, routine )
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_z4
! ==================================================================
! Logical array reallocs
! ==================================================================
SUBROUTINE realloc_l1( array, i1min,i1max,  &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='L'
integer, parameter                     :: rank=1
logical, dimension(:),      pointer    :: array, old_array
integer,                    intent(in) :: i1min,i1max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = .false.
end if
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_l1
! ==================================================================
SUBROUTINE realloc_l2( array, i1min,i1max, i2min,i2max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='L'
integer, parameter                     :: rank=2
logical, dimension(:,:),    pointer    :: array, old_array
integer,                    intent(in) :: i1min,i1max, i2min,i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min /)
new_bounds(2,:) = (/ i1max, i2max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = .false.
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2)) = &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_l2
! ==================================================================
SUBROUTINE realloc_l3( array, i1min,i1max, i2min,i2max, i3min,i3max, &
                       name, routine, copy, shrink )
implicit none
character, parameter                   :: type='L'
integer, parameter                     :: rank=3
logical, dimension(:,:,:),  pointer    :: array, old_array
integer,                    intent(in) :: i1min,i1max, i2min,i2max, &
                                          i3min,i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
integer, dimension(2,rank)             :: b, c, new_bounds, old_bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array 
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if
new_bounds(1,:) = (/ i1min, i2min, i3min /)
new_bounds(2,:) = (/ i1max, i2max, i3max /)
call options( b, c, old_bounds, new_bounds, copy, shrink )
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = .false.
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3)) = &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if
END SUBROUTINE realloc_l3
! ==================================================================
! Realloc routines with assumed lower bound = 1
!AG: Extremely dangerous -- do not use.
! ==================================================================
SUBROUTINE realloc_i1s( array, i1max, &
                        name, routine, copy, shrink )
! Arguments
implicit none
integer, dimension(:),      pointer    :: array
integer,                    intent(in) :: i1max
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
logical,          optional, intent(in) :: copy
logical,          optional, intent(in) :: shrink

call realloc_i1( array, DEFAULT%imin, i1max, &
                 name, routine, copy, shrink )

END SUBROUTINE realloc_i1s
! ==================================================================
SUBROUTINE realloc_i2s( array, i1max, i2max,  &
                        name, routine, copy, shrink )
implicit none
integer, dimension(:,:),    pointer    :: array
integer,                    intent(in) :: i1max, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_i2( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_i2s
! ==================================================================
SUBROUTINE realloc_i3s( array, i1max, i2max, i3max,  &
                        name, routine, copy, shrink )
implicit none
integer, dimension(:,:,:),  pointer    :: array
integer,                    intent(in) :: i1max, i2max, i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_i3( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 DEFAULT%imin, i3max,                             &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_i3s
! ==================================================================
SUBROUTINE realloc_r1s( array, i1max, &
                        name, routine, copy, shrink )
implicit none
real(SP), dimension(:),     pointer    :: array
integer,                    intent(in) :: i1max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_r1( array, DEFAULT%imin, i1max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_r1s
! ==================================================================
SUBROUTINE realloc_r2s( array, i1max, i2max, &
                        name, routine, copy, shrink )
implicit none
real(SP), dimension(:,:),   pointer    :: array
integer,                    intent(in) :: i1max, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_r2( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_r2s
! ==================================================================
SUBROUTINE realloc_r3s( array, i1max, i2max, i3max, &
                        name, routine, copy, shrink )
implicit none
real(SP), dimension(:,:,:), pointer    :: array
integer,                    intent(in) :: i1max, i2max, i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_r3( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 DEFAULT%imin, i3max,                             &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_r3s
! ==================================================================
SUBROUTINE realloc_r4s( array, i1max, i2max, i3max, i4max, &
                        name, routine, copy, shrink )
implicit none
real(SP), dimension(:,:,:,:), pointer  :: array
integer,                    intent(in) :: i1max, i2max, i3max, i4max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_r4( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                        DEFAULT%imin, i3max, DEFAULT%imin, i4max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_r4s
! ==================================================================
SUBROUTINE realloc_d1s( array, i1max, &
                        name, routine, copy, shrink )
implicit none
real(DP), dimension(:),     pointer    :: array
integer,                    intent(in) :: i1max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_d1( array, DEFAULT%imin, i1max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_d1s
! ==================================================================
SUBROUTINE realloc_d2s( array, i1max, i2max, &
                        name, routine, copy, shrink )
implicit none
real(DP), dimension(:,:),   pointer    :: array
integer,                    intent(in) :: i1max, i2max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_d2( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_d2s
! ==================================================================
SUBROUTINE realloc_d3s( array, i1max, i2max, i3max, &
                        name, routine, copy, shrink )
implicit none
real(DP), dimension(:,:,:), pointer    :: array
integer,                    intent(in) :: i1max, i2max, i3max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_d3( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 DEFAULT%imin, i3max,                             &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_d3s
! ==================================================================
SUBROUTINE realloc_d4s( array, i1max, i2max, i3max, i4max, &
                        name, routine, copy, shrink )
implicit none
real(DP), dimension(:,:,:,:), pointer  :: array
integer,                    intent(in) :: i1max, i2max, i3max, i4max
character(len=*), optional, intent(in) :: name, routine
logical,          optional, intent(in) :: copy, shrink
call realloc_d4( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                        DEFAULT%imin, i3max, DEFAULT%imin, i4max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_d4s
! ==================================================================
SUBROUTINE realloc_l1s( array, i1max, &
                        name, routine, copy, shrink )
implicit none
logical, dimension(:),      pointer    :: array
integer,                    intent(in) :: i1max
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
logical,          optional, intent(in) :: copy
logical,          optional, intent(in) :: shrink
call realloc_l1( array, DEFAULT%imin, i1max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_l1s
! ==================================================================
SUBROUTINE realloc_l2s( array, i1max, i2max, &
                        name, routine, copy, shrink )
implicit none
logical, dimension(:,:),    pointer    :: array
integer,                    intent(in) :: i1max, i2max
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
logical,          optional, intent(in) :: copy
logical,          optional, intent(in) :: shrink
call realloc_l2( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 name, routine, copy, shrink )
END SUBROUTINE realloc_l2s
! ==================================================================
SUBROUTINE realloc_l3s( array, i1max, i2max, i3max, &
                        name, routine, copy, shrink )
implicit none
logical, dimension(:,:,:),  pointer    :: array
integer,                    intent(in) :: i1max, i2max, i3max
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
logical,          optional, intent(in) :: copy
logical,          optional, intent(in) :: shrink
call realloc_l3( array, DEFAULT%imin, i1max, DEFAULT%imin, i2max, &
                 DEFAULT%imin, i3max, name, routine, copy, shrink )
END SUBROUTINE realloc_l3s
! ==================================================================
! Character vector realloc
! ==================================================================
SUBROUTINE realloc_s1( array, i1min, i1max, &
                       name, routine, copy, shrink )
! Arguments
implicit none
character(len=*), dimension(:),      pointer    :: array
integer,                    intent(in) :: i1min
integer,                    intent(in) :: i1max
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
logical,          optional, intent(in) :: copy
logical,          optional, intent(in) :: shrink

! Internal variables and arrays
character, parameter           :: type='S'
integer, parameter             :: rank=1
character(len=len(array)), dimension(:), pointer :: old_array
integer, dimension(2,rank)     :: b, c, new_bounds, old_bounds

! Get old array bounds
ASSOCIATED_ARRAY = associated(array)
if (ASSOCIATED_ARRAY) then
  old_array => array          ! Keep pointer to old array
  old_bounds(1,:) = lbound(old_array)
  old_bounds(2,:) = ubound(old_array)
end if

! Copy new requested array bounds
new_bounds(1,:) = (/ i1min /)
new_bounds(2,:) = (/ i1max /)

! Find if it is a new allocation or a true reallocation,
! and if the contents need to be copied (saved)
! Argument b returns common bounds
! Options routine also reads common variable ASSOCIATED_ARRAY,
! and it sets NEEDS_ALLOC, NEEDS_DEALLOC, and NEEDS_COPY
call options( b, c, old_bounds, new_bounds, copy, shrink )

! Deallocate old space
if (NEEDS_DEALLOC .and. .not.NEEDS_COPY) then
  call alloc_count( -size(old_array)*len(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if

! Allocate new space
if (NEEDS_ALLOC) then
  allocate( array(b(1,1):b(2,1)), stat=IERR )
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array)*len(array), type, name, routine )
  array = ''
end if

! Copy contents and deallocate old space
if (NEEDS_COPY) then
  array(c(1,1):c(2,1)) = old_array(c(1,1):c(2,1))
  call alloc_count( -size(old_array)*len(old_array), type, name, routine ) 
  deallocate(old_array,stat=IERR)
  call alloc_err( IERR, name, routine, old_bounds )
end if

END SUBROUTINE realloc_s1
! ==================================================================
! Dealloc routines
! ==================================================================
SUBROUTINE dealloc_i1( array, name, routine )

! Arguments
implicit none
integer, dimension(:),      pointer    :: array
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine

if (associated(array)) then
  call alloc_count( -size(array), 'I', name, routine )
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if

END SUBROUTINE dealloc_i1

! ==================================================================
SUBROUTINE dealloc_i2( array, name, routine )
implicit none
integer, dimension(:,:),    pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'I', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )

end if
END SUBROUTINE dealloc_i2
! ==================================================================
SUBROUTINE dealloc_i3( array, name, routine )
implicit none
integer, dimension(:,:,:),  pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'I', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )

end if
END SUBROUTINE dealloc_i3
! ==================================================================
SUBROUTINE dealloc_E1( array, name, routine )

! Arguments
implicit none
integer(i8b), dimension(:),      pointer    :: array
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine

if (associated(array)) then
  call alloc_count( -size(array), 'I', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )

end if

END SUBROUTINE dealloc_E1

! ==================================================================
SUBROUTINE dealloc_r1( array, name, routine )
implicit none
real(SP), dimension(:),     pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'R', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_r1
! ==================================================================
SUBROUTINE dealloc_r2( array, name, routine )
implicit none
real(SP), dimension(:,:),   pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'R', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_r2
! ==================================================================
SUBROUTINE dealloc_r3( array, name, routine )
implicit none
real(SP), dimension(:,:,:), pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'R', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_r3
! ==================================================================
SUBROUTINE dealloc_r4( array, name, routine )
implicit none
real(SP), dimension(:,:,:,:), pointer  :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'R', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_r4
! ==================================================================
SUBROUTINE dealloc_d1( array, name, routine )
implicit none
real(DP), dimension(:),     pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'D', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_d1
! ==================================================================
SUBROUTINE dealloc_d2( array, name, routine )
implicit none
real(DP), dimension(:,:),   pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'D', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_d2
! ==================================================================
SUBROUTINE dealloc_d3( array, name, routine )
implicit none
real(DP), dimension(:,:,:), pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'D', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_d3
! ==================================================================
SUBROUTINE dealloc_d4( array, name, routine )
implicit none
real(DP), dimension(:,:,:,:), pointer  :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'D', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_d4
! ==================================================================
SUBROUTINE dealloc_d5( array, name, routine )
implicit none
real(DP), dimension(:,:,:,:,:), pointer :: array
character(len=*), optional, intent(in)  :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'D', name, routine )
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_d5
! ==================================================================
! COMPLEX versions
!
SUBROUTINE dealloc_c1( array, name, routine )
implicit none
complex(SP), dimension(:),   pointer   :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -2*size(array), 'S', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_c1
! ==================================================================
SUBROUTINE dealloc_c2( array, name, routine )
implicit none
complex(SP), dimension(:,:),  pointer  :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -2*size(array), 'S', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_c2
! ==================================================================
SUBROUTINE dealloc_z1( array, name, routine )
implicit none
complex(DP), dimension(:),   pointer   :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -2*size(array), 'D', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_z1
! ==================================================================
SUBROUTINE dealloc_z2( array, name, routine )
implicit none
complex(DP), dimension(:,:),  pointer  :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -2*size(array), 'D', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_z2
! ==================================================================
SUBROUTINE dealloc_z3( array, name, routine )
implicit none
complex(DP), dimension(:,:,:), pointer :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -2*size(array), 'D', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_z3
! ==================================================================
SUBROUTINE dealloc_z4( array, name, routine )
implicit none
complex(DP), dimension(:,:,:,:),  pointer  :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -2*size(array), 'D', name, routine )
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_z4
! ==================================================================
SUBROUTINE dealloc_l1( array, name, routine )
implicit none
logical, dimension(:),      pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'L', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_l1
! ==================================================================
SUBROUTINE dealloc_l2( array, name, routine )
implicit none
logical, dimension(:,:),    pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'L', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_l2
! ==================================================================
SUBROUTINE dealloc_l3( array, name, routine )
implicit none
logical, dimension(:,:,:),  pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'L', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_l3
! ==================================================================
SUBROUTINE dealloc_s1( array, name, routine )
implicit none
character(len=*), dimension(:), pointer :: array
character(len=*), optional, intent(in)  :: name, routine
if (associated(array)) then
  call alloc_count( -size(array)*len(array), 'S', name, routine ) 
  deallocate(array,stat=IERR)
  call alloc_err( IERR, name, routine )
end if
END SUBROUTINE dealloc_s1

! ==================================================================
! Internal subroutines
! ==================================================================

SUBROUTINE options( final_bounds, common_bounds, &
                    old_bounds, new_bounds, copy, shrink )
! Arguments
integer, dimension(:,:), intent(out) :: final_bounds
integer, dimension(:,:), intent(out) :: common_bounds
integer, dimension(:,:), intent(in)  :: old_bounds
integer, dimension(:,:), intent(in)  :: new_bounds
logical,       optional, intent(in)  :: copy
logical,       optional, intent(in)  :: shrink

! Internal variables and arrays
logical want_shrink


!! AG*****
!  It might be worthwhile to check whether the user
!  atttemps to use bounds which do not make sense,
!  such as zero, or with upper<lower...
!!***

! Find if it is a new allocation or a true reallocation,
! and if the contents need to be copied (saved)
if (ASSOCIATED_ARRAY) then

  ! Check if array bounds have changed
  if ( all(new_bounds==old_bounds) ) then
    ! Old and new arrays are equal. Nothing needs to be done
    NEEDS_ALLOC   = .false. 
    NEEDS_DEALLOC = .false.
    NEEDS_COPY    = .false.
  else 

    ! Want to shrink?
    if (present(shrink)) then
      want_shrink = shrink
    else
      want_shrink = DEFAULT%shrink
    end if

    if (.not. want_shrink  &
        .and. all(new_bounds(1,:)>=old_bounds(1,:)) &
        .and. all(new_bounds(2,:)<=old_bounds(2,:)) ) then
      ! Old array is already fine. Nothing needs to be done
      NEEDS_ALLOC   = .false. 
      NEEDS_DEALLOC = .false.
      NEEDS_COPY    = .false.
    else
      ! Old array needs to be substituted by a new array
      NEEDS_ALLOC   = .true.
      NEEDS_DEALLOC = .true.
      if (present(copy)) then
        NEEDS_COPY = copy
      else
        NEEDS_COPY = DEFAULT%copy
      end if

      ! Ensure that bounds shrink only if desired
      if (want_shrink) then
        final_bounds(1,:) = new_bounds(1,:)
        final_bounds(2,:) = new_bounds(2,:)
      else
        final_bounds(1,:) = min( old_bounds(1,:), new_bounds(1,:) )
        final_bounds(2,:) = max( old_bounds(2,:), new_bounds(2,:) )
      end if

      ! Find common section of old and new arrays
      common_bounds(1,:) = max( old_bounds(1,:), final_bounds(1,:) )
      common_bounds(2,:) = min( old_bounds(2,:), final_bounds(2,:) )
    end if

  end if

else
  ! Old array does not exist. Allocate new one
  NEEDS_ALLOC   = .true. 
  NEEDS_DEALLOC = .false.
  NEEDS_COPY    = .false.
  final_bounds(1,:) = new_bounds(1,:)
  final_bounds(2,:) = new_bounds(2,:)
end if

END SUBROUTINE options

! ==================================================================

SUBROUTINE alloc_count( delta_size, type, name, routine )

implicit none

integer, intent(in)          :: delta_size  ! +/-size(array)
character, intent(in)        :: type        ! 'I' => integer
                                            ! 'E' => integer*8
                                            ! 'R' => real*4
                                            ! 'D' => real*8
                                            ! 'L' => logical
                                            ! 'S' => character (string)
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine

character(len=32)   :: aname, rname
character(len=1)    :: memType, task
real(DP)            :: delta_mem
logical             :: newPeak
logical,  save      :: header_written = .false.
logical,  save      :: tree_nullified = .false.
integer             :: memSize

! Set routine name
if (present(routine)) then
  rname = routine
else
  rname = DEFAULT%routine
end if

! Call siesta counting routine 'memory'
! Switched off in Aug.2009, and made 'memory' call alloc_count instead
! if (delta_size > 0) then
!   task = 'A'
! else
!   task = 'D'
! end if
! select case( type )
! case ('R')         ! Real --> single
!   memType = 'S'
!   memSize = abs(delta_size)
! case ('S')         ! Character (string) --> integer/4
!   memType = 'I'
!   memSize = abs(delta_size) / type_mem('I')
! case default
!   memType = type
!   memSize = abs(delta_size)
! end select
! call memory( task, memType, memSize, trim(rname) )


if (REPORT_LEVEL <= 0) return

! Compound routine+array name
if (present(name) .and. present(routine)) then
  aname = trim(routine)//' '//name
else if (present(name) .and. DEFAULT%routine/=DEFAULT_ROUTINE) then
  aname = trim(DEFAULT%routine)//' '//name
else if (present(name)) then
  aname = name
else if (present(routine)) then
  aname = trim(routine)//' '//DEFAULT_NAME
else if (DEFAULT%routine/=DEFAULT_ROUTINE) then
  aname = trim(DEFAULT%routine)//' '//DEFAULT_NAME
else
  aname = DEFAULT_ROUTINE//' '//DEFAULT_NAME
end if

MAX_LEN = max( MAX_LEN, len(trim(aname)) )

! Find memory increment and total allocated memory
delta_mem = delta_size * type_mem(type)
TOT_MEM = TOT_MEM + delta_mem
if (TOT_MEM > PEAK_MEM+0.5_dp) then
  newPeak = .true.
  PEAK_MEM = TOT_MEM
  PEAK_ARRAY = aname
  PEAK_ROUTINE = rname
!  print'(/,a,f18.6),a,/)',
!    'alloc_count: Memory peak =', PEAK_MEM/MBYTE, ' Mbytes'
else
  newPeak = .false.
end if

! Add/subtract/classify array memory
if (REPORT_LEVEL > 1) then
  if (.not.tree_nullified) then
    nullify(report_tree)
    tree_nullified = .true.
  end if
  call tree_add( report_tree, aname, delta_mem )
  if (newPeak) call tree_peak( report_tree )
end if

! Print report, but only in node 0, as not all 
! processors may follow the same route here
!   The detail/extent of the report increses with the value of level:
! level=0 : no report at all (the default)
! level=1 : only total memory peak and where it occurred
! level=2 : detailed report created but printed only upon request
! level=3 : detailed report printed at every new memory peak
! level=4 : print every individual reallocation or deallocation


if (newPeak .and. (REPORT_LEVEL==1 .or. REPORT_LEVEL==3) .and. &
    node == 0) then
  call print_report(.false.)
end if

if (REPORT_LEVEL == 4 .and. node == 0) then
  if (.not.header_written) then
    write(REPORT_UNIT,'(/,a7,9x,1x,a4,28x,1x,2a15)') &
     'Routine', 'Name', 'Incr. (MB)', 'Total (MB)'
    header_written = .true.
  end if
  write(REPORT_UNIT,'(a16,1x,a32,1x,2f15.6)') &
     rname, aname, delta_mem/MBYTE, TOT_MEM/MBYTE
end if
END SUBROUTINE alloc_count

! ==================================================================

INTEGER FUNCTION type_mem( var_type )
!
! It is not clear that the sizes assumed are universal for
! non-Cray machines...
!
implicit none
character, intent(in) :: var_type
character(len=40)     :: message

select case( var_type )
#ifdef OLD_CRAY
  case('I')
    type_mem = 8
  case('R')
    type_mem = 8
  case('L')
    type_mem = 8
#else
  case('I')
    type_mem = 4
  case('R')
    type_mem = 4
  case('L')
    type_mem = 4
#endif
case('E')
  type_mem = 8
case('D')
  type_mem = 8
case('S')
  type_mem = 1
case default
  write(message,"(2a)") &
    'alloc_count: ERROR: unknown type = ', var_type
  call die(trim(message))
end select

END FUNCTION type_mem

! ==================================================================

RECURSIVE SUBROUTINE tree_add( t, name, delta_mem )

implicit none
type(TREE),       pointer    :: t
character(len=*), intent(in) :: name
real(DP),         intent(in) :: delta_mem

logical, save :: warn_negative = .true.

if (.not.associated(t)) then
  allocate( t )
  t%name = name
  t%mem  = delta_mem
  t%max  = delta_mem
  t%peak = 0._dp
  nullify( t%left, t%right )
else if (name == t%name) then
  t%mem = t%mem + delta_mem
  ! The abs is to handle the case of apparent de_allocs without re_allocs,
  ! caused by routine/name argument mismatches
  if (abs(t%mem) > abs(t%max)) t%max = t%mem
else if ( llt(name,t%name) ) then
  call tree_add( t%left, name, delta_mem )
else
  call tree_add( t%right, name, delta_mem )
end if

if (warn_negative .and. t%mem<0._dp) then
  call parallel_init()   ! Make sure that node and Nodes are initialized
  if (Node==0) then
   write(6,'(/,a,/,2a,/,a,f18.0,a)')  &
      'WARNING: alloc-realloc-dealloc name mismatch',  &
      '         Name: ', trim(name),                   &
      '         Size: ', t%mem, ' Bytes'
    if (Nodes>1) write(6,'(9x,a,i6)') 'Node:', Node
    write(6,'(9x,a)') 'Subsequent mismatches will not be reported'
    warn_negative = .false.  ! Print this warning only once
  end if
end if

END SUBROUTINE tree_add

! ==================================================================

RECURSIVE SUBROUTINE tree_peak( t )

implicit none
type(TREE), pointer :: t

if (.not.associated(t)) return

t%peak = t%mem
call tree_peak( t%left )
call tree_peak( t%right )

END SUBROUTINE tree_peak

! ==================================================================

RECURSIVE SUBROUTINE tree_print( t )

implicit none
type(TREE), pointer :: t

if (.not.associated(t)) return

call tree_print( t%left )

if (abs(t%max) >= REPORT_THRESHOLD) then
  write(REPORT_UNIT,'(a,1x,3f15.6,f9.2)') &
    t%name(1:MAX_LEN), t%mem/MBYTE, t%max/MBYTE, t%peak/MBYTE, &
    100._dp * t%peak / (PEAK_MEM + tiny(PEAK_MEM) )
end if

call tree_print( t%right )

END SUBROUTINE tree_print

! ==================================================================

SUBROUTINE print_report(all)

implicit none

! Whether MPI reductions should be performed
! If, not then ensure that no MPI calls are performed
logical, intent(in), optional :: all 

character(len=80)   :: string = 'Name'
character           :: date*8, time*10, zone*5
integer             :: iNode, peakNode
real(dp)            :: maxPeak
real(dp),allocatable:: nodeMem(:), nodePeak(:)
logical             :: lall

#ifdef MPI
integer           :: MPIerror
#endif

! Enables parallel call of print-report
lall = .true.
if ( present(all) ) lall = all

! Only if MPI should all be used
if ( lall ) then
! Make sure that variables node and Nodes are initialized
call parallel_init()
end if

! Allocate and initialize two small arrays
allocate( nodeMem(0:Nodes-1), nodePeak(0:Nodes-1) )
! initialize (in case all == .false.)
nodeMem = 0._dp
nodePeak = 0._dp

! Initializations for Nodes=1 (serial case)
nodeMem(node) = TOT_MEM
nodePeak(node) = PEAK_MEM
peakNode = node

! In parallel, find the memory values of all nodes
#ifdef MPI
if (Nodes > 1 .and. lall ) then
  ! Gather the present and peak memories of all nodes
  call MPI_AllGather( TOT_MEM, 1, MPI_double_precision, &
                      nodeMem, 1, MPI_double_precision, &
                      MPI_COMM_WORLD, MPIerror )
  call MPI_AllGather( PEAK_MEM, 1, MPI_double_precision, &
                      nodePeak, 1, MPI_double_precision, &
                      MPI_COMM_WORLD, MPIerror )
  ! Find the node with the highest peak of memory
  maxPeak = 0
  do iNode = 0,Nodes-1
    if (nodePeak(iNode) > maxPeak) then
      peakNode = iNode
      maxPeak = nodePeak(iNode)
    end if
  end do ! iNode
  ! Change the writing node for the peak-node information
  if (node==0 .and. peakNode/=0) close( REPORT_UNIT )
  call MPI_Barrier( MPI_COMM_WORLD, MPIerror )
  if (node==peakNode .and. peakNode/=0) then
     call io_assign(REPORT_UNIT)
     open( unit=REPORT_UNIT, file=REPORT_FILE, &
          status='unknown', position='append' )
  end if
end if ! (Nodes>1)
#endif

! The report is printed by the highest-peak node
if (node == peakNode) then

  ! AG: Commented out to allow multiple batches of information
  ! if (REPORT_LEVEL < 4) rewind(REPORT_UNIT)

  call date_and_time( date, time, zone )

  write(REPORT_UNIT,'(/,a,16a)')                &
    'Allocation summary at ',                   &
    date(1:4),'/',date(5:6),'/',date(7:8),' ',  &
    time(1:2),':',time(3:4),':',time(5:10),' ', &
    zone(1:3),':',zone(4:5)

  if (Nodes > 1) then
    write(REPORT_UNIT,'(/,(a,f18.6,a))')            &
      'Present memory all nodes : ', sum(nodeMem)/MBYTE,  ' MB', &
      'Added peak mem all nodes : ', sum(nodePeak)/MBYTE, ' MB', &
      'Min peak memory in a node: ', minval(nodePeak)/MBYTE, ' MB', &
      'Max peak memory in a node: ', maxval(nodePeak)/MBYTE, ' MB'
!   Impractical for many nodes:
!    write(REPORT_UNIT,'(/,a,/,(i6,f12.6))') &
!      'Memory peaks of individual nodes (Mb):', &
!      (iNode,nodePeak(iNode)/MBYTE,iNode=0,Nodes-1)
    write(REPORT_UNIT,'(/,a,i6)') &
      'Maximum peak of memory occurred in node:', peakNode
  end if

  write(REPORT_UNIT,'(2(/,a,f18.6,a),/,2a,/,2a)')            &
    'Present memory allocation: ', TOT_MEM/MBYTE,  ' MB', &
    'Maximum memory allocation: ', PEAK_MEM/MBYTE, ' MB', &
    'Occurred after allocating: ', trim(PEAK_ARRAY),         &
    'In routine:                ', trim(PEAK_ROUTINE)

  if (REPORT_LEVEL > 1) then
    if (REPORT_THRESHOLD > 0._dp) then
      write(REPORT_UNIT,'(/,a,f12.6,a,/,a,1x,3a15,a9)') &
        'Allocated sizes (in MByte) of arrays larger than ', &
        REPORT_THRESHOLD/MBYTE, ' MB:', &
        string(1:MAX_LEN), 'Present', 'Maximum', 'At peak', '%'
    else
      write(REPORT_UNIT,'(/,a,/,a,1x,3a15,a9)') &
        'Allocated array sizes (in MByte):', &
        string(1:MAX_LEN), 'Present', 'Maximum', 'At peak', '%'
    end if
    call tree_print( report_tree )
  end if

  ! Close file if not common IO node
  if ( node /= 0 ) call io_close(REPORT_UNIT)
end if ! (node == peakNode)

! Change again the writing node for the rest of the report
#ifdef MPI
if (lall) call MPI_Barrier( MPI_COMM_WORLD, MPIerror )
if (node==0 .and. peakNode/=0) &
     open( unit=REPORT_UNIT, file=REPORT_FILE, &
     status='unknown', position='append' )
#endif

deallocate( nodeMem, nodePeak )

END SUBROUTINE print_report

! ==================================================================

SUBROUTINE alloc_err( ierr, name, routine, bounds )

implicit none

integer,                    intent(in) :: ierr
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
integer, dimension(:,:), optional, intent(in) :: bounds

integer i

if (ierr/=0) then
  if (ionode) print*, 'alloc_err: allocate status error', ierr
  if (present(name).and.present(routine)) then
    if (ionode) print*, 'alloc_err: array ', name, &
               ' requested by ', routine
  elseif (present(name)) then
    if (ionode) print*, 'alloc_err: array ', name, &
               ' requested by unknown'
  elseif (present(routine)) then
    if (ionode) print* , 'alloc_err: array unknown', &
               ' requested by ', routine
  endif
  if (ionode.and.present(bounds))                            &
    print '(a,i3,2i10)', ('alloc_err: dim, lbound, ubound:',  &
          i,bounds(1,i),bounds(2,i),                         &
          i=1,size(bounds,dim=2))            

  call die('alloc_err: allocate error')
end if

END SUBROUTINE alloc_err

! ==================================================================

END MODULE alloc
