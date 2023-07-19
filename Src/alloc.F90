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
MODULE alloc

      use precision, only: sp, dp
      use parallel,   only : Node, Nodes, ionode
      use sys, only: die
#ifdef MPI
      use mpi_siesta
#endif

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
!     each call to realloc.
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
!   call realloc(...)
!   call alloc_default( restore=oldDefaults )
! Notice that, if the restore call is skipped, the new defaults will
! stay in effect until a new call to alloc_dafault is made.
! ==================================================================
! SUBROUTINE alloc_report( level, unit, file, printNow )
!   Sets the output file for the allocation report
! INPUT (optional):
!   integer       :: level    : Level (detail) of report
!   integer       :: unit     : Output unit
!   character*(*) :: file     : Output file name
!   logical       :: printNow : If present & true => print report now
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
! ==================================================================
! SUBROUTINE realloc( array, [i1min,] i1max,
!                     [[i2min,] i2max, [[i3min,] i3max]],
!                     name, routine, copy, shrink )
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
!                              integer, rank 1, 2, or 3
!                              real*4,  rank 1, 2, or 3
!                              real*8,  rank 1, 2, or 3
!                              logical, rank 1, 2, or 3
!                              character(len=*), rank 1
! BEHAVIOR:
!   Pointers MUST NOT enter in an undefined state. Before using them
! for the first time, they must be nullified explicitly. Alternatively,
! in f95, the can be initialized as null() upon declaration.
!   If argument array is not associated on input, it is just allocated.
!   If array is associated and has the same bounds (or smaller bonds
! and shrink is false) nothing is done. Thus, it is perfectly safe and
! efficient to call realloc repeatedly without deallocating the array.
! However, subroutine dealloc is provided to eliminate large arrays
! when they are not needed.
!   In order to save (copy) the contents of the old array, the new array
! needs to be allocated before deallocating the old one. Thus, if the
! contents are not needed, or if reducing memory is a must, calling
! realloc with copy=.false. makes it to deallocate before allocating.
!   The elements that are not copied (because copy=.false. or because
! they are outside the bounds of the input array) return with value
! zero (integer and real), .false. (logical), or blank (character).
!   If imin>imax for any dimension, the array pointer returns
! associated to a zero-size array.
!   Besides allocating or reallocating the array, realloc keeps a count
! of memory usage (and prints a report in a file previously fixed by a
! call to alloc_report). Thus, it is not recommended to call realloc
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
! SUBROUTINE dealloc( array, name, routine )
! INPUT (optional):
!   character*(*) :: name    : Actual array name or a label for it
!   character*(*) :: routine : Name of the calling routine
!                                or routine section
! INPUT/OUTPUT:
!   TYPE, pointer :: array : Array be deallocated (same types and
!                            kinds as in realloc).
! BEHAVIOR:
!   Besides deallocating the array, realloc decreases the count of
! memory usage previously counted by realloc. Thus, dealloc should 
! not be called to deallocate an array not allocated by realloc.
! Equally, arrays allocated or reallocated by realloc should be 
! deallocated by dealloc.
! ==================================================================---

implicit none

PUBLIC ::             &
  alloc_default,      &! Sets allocation defaults
  alloc_report,       &! Sets log report defaults
  re_alloc,           &! Allocation/reallocation
  de_alloc,           &! Deallocation
  allocDefaults        ! Derived type to hold allocation defaults

PRIVATE      ! Nothing is declared public beyond this point


interface de_alloc
  module procedure dealloc_d1, dealloc_d2, dealloc_d3, dealloc_d4, &
                   dealloc_i1, dealloc_i2, dealloc_i3,             &
                   dealloc_l1, dealloc_l2, dealloc_l3,             &
                   dealloc_s1,                         &
                   dealloc_r1, dealloc_r2, dealloc_r3, dealloc_r4, &
                   dealloc_z1, dealloc_z2
end interface

interface re_alloc
  module procedure &
    realloc_d1,  realloc_i1,  realloc_l1,  realloc_r1,  realloc_z1, &
    realloc_d2,  realloc_i2,  realloc_l2,  realloc_r2,  realloc_z2, &
    realloc_d3,  realloc_i3,  realloc_l3,  realloc_r3,  realloc_r4, &
    realloc_d4,                                         &
    realloc_s1                                         
!AG: Dangerous!!!    realloc_d2s, realloc_i2s, realloc_l2s, realloc_r2s, &
!   realloc_d3s, realloc_i3s, realloc_l3s, realloc_r3s, &
!    realloc_d4s, &

end interface

! Derived type to hold allocation default options
type allocDefaults
  private
  logical           copy
  logical           shrink
  integer           imin
  character(len=32) routine
end type allocDefaults

! Initial default values
type(allocDefaults), save ::   &
  DEFAULT = allocDefaults(     &
    .true.,                    &! Copy default
    .true.,                    &! Shrink default
     1,                        &! Imin default
    'unknown' )                 ! Routine name default
character(len=*), parameter :: &
  DEFAULT_NAME = 'unknown'      ! Array name default
integer, save ::               &
  REPORT_LEVEL = 0,            &! Level (detail) of allocation report
  REPORT_UNIT  = 0              ! Output file unit for report
character(len=50), save ::     &
  REPORT_FILE = 'alloc_report'  ! Output file name for report
  
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

SUBROUTINE alloc_report( level, unit, file, printNow )

implicit none

integer,          optional, intent(in) :: level, unit
character(len=*), optional, intent(in) :: file
logical,          optional, intent(in) :: printNow

logical open

#ifdef MPI
integer MPIerror
#endif

if (present(level)) then
  REPORT_LEVEL = level
end if

if (node == 0) then
  if (present(unit)) then
    if (unit > 0) then
      REPORT_UNIT = unit
      if (present(file)) then
        REPORT_FILE = file
      else
        REPORT_FILE = 'unknown'
      end if
    end if
  else if (present(file)) then
    if (file /= REPORT_FILE) then
      REPORT_FILE = file
      inquire( file=REPORT_FILE, opened=open, number=REPORT_UNIT )
      if (.not.open) then
        call io_assign(REPORT_UNIT)
        open( REPORT_UNIT, file=REPORT_FILE, status='unknown')
      end if
    end if
  else if (REPORT_UNIT==0) then
    REPORT_FILE = 'alloc_report'
    call io_assign(REPORT_UNIT)
    open( REPORT_UNIT, file=REPORT_FILE, status='unknown')
  end if
end if

#ifdef MPI
! Distribute information to other nodes and open REPORT_UNIT
call MPI_Bcast(REPORT_LEVEL,1,MPI_integer,0,MPI_Comm_World,MPIerror)
call MPI_Bcast(REPORT_UNIT,1,MPI_integer,0,MPI_Comm_World,MPIerror)
call MPI_Bcast(REPORT_FILE,50,MPI_character,0,MPI_Comm_World,MPIerror)

if (node > 0) then
  open( REPORT_UNIT, file=REPORT_FILE )
end if
#endif

if (present(printNow)) then
  if (printNow) call print_report
end if

END SUBROUTINE alloc_report

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
  deallocate(old_array)
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
  deallocate(old_array)
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
  deallocate(old_array)
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
  deallocate(old_array)
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
  deallocate(old_array)
end if
if (NEEDS_ALLOC) then
  allocate(array(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3)) =  &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array)
end if
END SUBROUTINE realloc_i3
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
  deallocate(old_array)
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
  deallocate(old_array)
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
  deallocate(old_array)
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
  deallocate(old_array)
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
  deallocate(old_array)
end if
if (NEEDS_ALLOC) then
  allocate(array(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._sp
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3)) =  &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array)
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
  deallocate(old_array)
end if
if (NEEDS_ALLOC) then
  allocate(array(b(1,1):b(2,1),b(1,2):b(2,2), &
                 b(1,3):b(2,3),b(1,4):b(2,4)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._sp
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4))= &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array)
end if
END SUBROUTINE realloc_r4
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
  deallocate(old_array)
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
  deallocate(old_array)
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
  deallocate(old_array)
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
  deallocate(old_array)
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
  deallocate(old_array)
end if
if (NEEDS_ALLOC) then
  allocate(array(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)),stat=IERR)
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
  deallocate(old_array)
end if
END SUBROUTINE realloc_d3
! ==================================================================
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
  deallocate(old_array)
end if
if (NEEDS_ALLOC) then
  allocate(array(b(1,1):b(2,1),b(1,2):b(2,2), &
                 b(1,3):b(2,3),b(1,4):b(2,4)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = 0._dp
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4))= &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3),c(1,4):c(2,4))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array)
end if
END SUBROUTINE realloc_d4
! ==================================================================
! COMPLEX versions
!
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
  deallocate(old_array)
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
  deallocate(old_array)
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
  deallocate(old_array)
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
  deallocate(old_array)
end if
END SUBROUTINE realloc_z2
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
  deallocate(old_array)
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
  deallocate(old_array)
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
  deallocate(old_array)
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
  deallocate(old_array)
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
  deallocate(old_array)
end if
if (NEEDS_ALLOC) then
  allocate(array(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)),stat=IERR)
  call alloc_err( IERR, name, routine, new_bounds )
  call alloc_count( size(array), type, name, routine )
  array = .false.
end if
if (NEEDS_COPY) then
      array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3)) = &
  old_array(c(1,1):c(2,1),c(1,2):c(2,2),c(1,3):c(2,3))
  call alloc_count( -size(old_array), type, name, routine ) 
  deallocate(old_array)
end if
END SUBROUTINE realloc_l3
! ==================================================================
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
  deallocate(old_array)
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
  deallocate(old_array)
end if

END SUBROUTINE realloc_s1

! ==================================================================
! ==================================================================
SUBROUTINE dealloc_i1( array, name, routine )

! Arguments
implicit none
integer, dimension(:),      pointer    :: array
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine

if (associated(array)) then
  call alloc_count( -size(array), 'I', name, routine ) 
  deallocate(array)
end if

END SUBROUTINE dealloc_i1

! ==================================================================
SUBROUTINE dealloc_i2( array, name, routine )
implicit none
integer, dimension(:,:),    pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'I', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_i2
! ==================================================================
SUBROUTINE dealloc_i3( array, name, routine )
implicit none
integer, dimension(:,:,:),  pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'I', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_i3
! ==================================================================
SUBROUTINE dealloc_r1( array, name, routine )
implicit none
real(SP), dimension(:),     pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'R', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_r1
! ==================================================================
SUBROUTINE dealloc_r2( array, name, routine )
implicit none
real(SP), dimension(:,:),   pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'R', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_r2
! ==================================================================
SUBROUTINE dealloc_r3( array, name, routine )
implicit none
real(SP), dimension(:,:,:), pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'R', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_r3
! ==================================================================
SUBROUTINE dealloc_r4( array, name, routine )
implicit none
real(SP), dimension(:,:,:,:), pointer  :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'R', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_r4
! ==================================================================
SUBROUTINE dealloc_d1( array, name, routine )
implicit none
real(DP), dimension(:),     pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'D', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_d1
! ==================================================================
SUBROUTINE dealloc_d2( array, name, routine )
implicit none
real(DP), dimension(:,:),   pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'D', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_d2
! ==================================================================
SUBROUTINE dealloc_d3( array, name, routine )
implicit none
real(DP), dimension(:,:,:), pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'D', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_d3
! ==================================================================
SUBROUTINE dealloc_d4( array, name, routine )
implicit none
real(DP), dimension(:,:,:,:), pointer  :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'D', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_d4
! ==================================================================
! COMPLEX versions
!
SUBROUTINE dealloc_z1( array, name, routine )
implicit none
complex(DP), dimension(:),   pointer   :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -2*size(array), 'D', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_z1
! ==================================================================
SUBROUTINE dealloc_z2( array, name, routine )
implicit none
complex(DP), dimension(:,:),  pointer  :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -2*size(array), 'D', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_z2
! ==================================================================
SUBROUTINE dealloc_l1( array, name, routine )
implicit none
logical, dimension(:),      pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'L', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_l1
! ==================================================================
SUBROUTINE dealloc_l2( array, name, routine )
implicit none
logical, dimension(:,:),    pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'L', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_l2
! ==================================================================
SUBROUTINE dealloc_l3( array, name, routine )
implicit none
logical, dimension(:,:,:),  pointer    :: array
character(len=*), optional, intent(in) :: name, routine
if (associated(array)) then
  call alloc_count( -size(array), 'L', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_l3
! ==================================================================
SUBROUTINE dealloc_s1( array, name, routine )
implicit none
character(len=*), dimension(:), pointer :: array
character(len=*), optional, intent(in)  :: name, routine
if (associated(array)) then
  call alloc_count( -size(array)*len(array), 'S', name, routine ) 
  deallocate(array)
end if
END SUBROUTINE dealloc_s1
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
if (delta_size > 0) then
  task = 'A'
else
  task = 'D'
end if
select case( type )
case ('R')         ! Real --> single
  memType = 'S'
  memSize = abs(delta_size)
case ('S')         ! Character (string) --> integer/4
  memType = 'I'
  memSize = abs(delta_size) / type_mem('I')
case default
  memType = type
  memSize = abs(delta_size)
end select
call memory( task, memType, memSize, trim(rname) )

if (REPORT_LEVEL <= 0) return

! Compound routine+array name
if (present(name)) then
  aname = name
else
  aname = DEFAULT_NAME
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

! Print report - only do this if number of nodes is 1 as
! not all processors made follow the same route here
if (newPeak .and. (REPORT_LEVEL==1 .or. REPORT_LEVEL==3) .and. &
    node == 0) then
  call print_report
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
  t%max = max( t%max, t%mem )
else if ( llt(name,t%name) ) then
  call tree_add( t%left, name, delta_mem )
else
  call tree_add( t%right, name, delta_mem )
end if

if (warn_negative .and. t%mem<0._dp) then
  write(6,'(/,a,/,2a,/,a,f18.0,a,/,a,/)')  &
    'WARNING: alloc-realloc-dealloc name mismatch',  &
    '         Name: ', trim(name),                   &
    '         Size: ', t%mem, ' Bytes',              &
    'Subsequent mismatches will not be reported'
!!!!!!!  warn_negative = .false.
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

write(REPORT_UNIT,'(a,1x,3f15.6,f9.2)') &
  t%name(1:MAX_LEN), t%mem/MBYTE, t%max/MBYTE, t%peak/MBYTE, &
  100._dp * t%peak / (PEAK_MEM + tiny(PEAK_MEM) )

call tree_print( t%right )

END SUBROUTINE tree_print

! ==================================================================

SUBROUTINE print_report

implicit none

character(len=80) :: string = 'Name'
character         :: date*8, time*10, zone*5

#ifdef MPI
integer           :: MPIerror, nodePEAK
integer           :: Status(MPI_Status_Size)
real(DP)          :: G_TOT_MEM, G_PEAK_MEM(2), L_PEAK_MEM(2)
#endif

#ifdef MPI
if (nodes > 1) then
!
! For parallel version need to find the values across all nodes
!

! Create combination of local peak memory value and node number
  L_PEAK_MEM(1) = PEAK_MEM
  L_PEAK_MEM(2) = node

! Find maximum memory values and their nodes
  call MPI_Reduce(TOT_MEM,G_TOT_MEM,1,MPI_double_precision, &
    MPI_Sum,0,MPI_Comm_World,MPIerror)
  call MPI_Reduce(L_PEAK_MEM,G_PEAK_MEM,1,MPI_2double_precision, &
    MPI_MaxLoc,0,MPI_Comm_World,MPIerror)

! Tell the rest of the nodes where the peak is
  if (node == 0) nodePEAK = nint(G_PEAK_MEM(2))
  call MPI_Bcast(nodePEAK,1,MPI_integer,0,MPI_Comm_World,MPIerror)

! Get details of peak use - if the maximum is on node 0 then just copy
  if (nodePEAK > 0) then

    if (node == 0) then
      call MPI_Send(G_TOT_MEM,1,MPI_integer,nodePEAK, &
        1,MPI_Comm_World,MPIerror)
    endif

    if (node == nodePEAK) then
      call MPI_Recv(G_TOT_MEM,1,MPI_integer,0, &
        1,MPI_Comm_World,Status,MPIerror)
    endif

  end if

else

  nodePEAK  = 0
  G_TOT_MEM = TOT_MEM

end if
#endif

if (node == 0) then

! AG: Allow multiple batches of information
!!!      if (REPORT_LEVEL < 4) rewind(REPORT_UNIT)

  call date_and_time( date, time, zone )

  write(REPORT_UNIT,'(/,a,16a)')                &
    'Allocation summary at ',                   &
    date(1:4),'/',date(5:6),'/',date(7:8),' ',  &
    time(1:2),':',time(3:4),':',time(5:10),' ', &
    zone(1:3),':',zone(4:5)

  write(REPORT_UNIT,'(2(/,a,f18.6,a),/,2a,/,2a)')            &
    'Present memory allocation: ', TOT_MEM/MBYTE,  ' Mbyte', &
    'Maximum memory allocation: ', PEAK_MEM/MBYTE, ' Mbyte', &
    'Occurred after allocating: ', trim(PEAK_ARRAY),         &
    'In routine:                ', trim(PEAK_ROUTINE)
#ifdef MPI
  if (nodes > 1) then
    write(REPORT_UNIT,'(/,a,f18.6,a)')            &
      'Present memory all nodes : ', G_TOT_MEM/MBYTE,  ' Mbyte'
  end if
#endif

  if (REPORT_LEVEL > 1) then
    write(REPORT_UNIT,'(/,a,/,a,1x,3a15,a9)') &
      'Allocated array sizes (in Mbyte):', &
      string(1:MAX_LEN), 'Present', 'Maximum', 'At peak', '%'
    call tree_print( report_tree )
  end if

end if

END SUBROUTINE print_report

! ==================================================================

SUBROUTINE alloc_err( ierr, name, routine, bounds )

implicit none

integer,                    intent(in) :: ierr
character(len=*), optional, intent(in) :: name
character(len=*), optional, intent(in) :: routine
integer, dimension(:,:),    intent(in) :: bounds

integer i
#ifdef MPI
integer mpierror
#endif

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
  if (ionode) print '(a,i3,2i8)', ('alloc_err: dim, lbound, ubound:', &
                      i,bounds(1,i),bounds(2,i),         &
                      i=1,size(bounds,dim=2))

  call die('alloc_err: allocate error')
end if

END SUBROUTINE alloc_err

! ==================================================================

END MODULE alloc
