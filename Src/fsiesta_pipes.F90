! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module fsiesta

! Support routines for siesta-as-a-subroutine in Unix/Linux.
! The routines that handle the other side of the communication are
! in module iopipes of siesta program.
! Usage:
!   call siesta_launch( label, nnodes, mpi_comm, launcher, localhost )
!     character(len=*),intent(in) :: label  : Name of siesta process
!                                             (prefix of its .fdf file)
!     integer,optional,intent(in) :: nnodes : Number of MPI nodes
!     integer,optional,intent(in) :: mpi_comm : not used in this version
!     character(len=*),intent(in),optional:: launcher : full launch command
!     logical,optional,intent(in) :: localhost : will siesta run at localhost?
!                                                (not used in this version)
!
!   call siesta_units( length, energy )
!     character(len=*),intent(in) :: length : Physical unit of length
!     character(len=*),intent(in) :: energy : Physical unit of energy
!
!   call siesta_forces( label, na, xa, cell, energy, fa, stress )
!     character(len=*), intent(in) :: label      : Name of siesta process
!     integer,          intent(in) :: na         : Number of atoms
!     real(dp),         intent(in) :: xa(3,na)   : Cartesian coords
!     real(dp),optional,intent(in) :: cell(3,3)  : Unit cell vectors
!     real(dp),optional,intent(out):: energy     : Total energy
!     real(dp),optional,intent(out):: fa(3,na)   : Atomic forces
!     real(dp),optional,intent(out):: stress(3,3): Stress tensor
!   call siesta_quit( label )
!     character(len=*),intent(in) :: label  : Name of one siesta process,
!                                             or 'all' to stop all procs.
! Behaviour:
! - If nnodes is not present among siesta_launch arguments, or nnodes<2, 
!   a serial siesta process will be launched. Otherwise, a parallel 
!   mpirun process will be launched. In this case, the mpi launching
!   command (e.g., "mpiexec <options> -n ") can be specified in the 
!   optional argument mpi_launcher
! - If siesta_units is not called, length='Ang', energy='eV' are
!   used by default. If it is called more than once, the units in the
!   last call become in effect.
! - The physical units set by siesta_units are used for all the siesta
!   processes launched
! - If siesta_forces is called without a previous call to siesta_launch
!   for that label, it assumes that the siesta process has been launched
!   (and the communication pipes created) externally in the shell.
!   In this case, siesta_forces only opens its end of the pipes and begins
!   communication through them.
! - If argument cell is not present in the call to siesta_forces, or if
!   the cell has zero volume, it is assumed that the system is a molecule,
!   and a supercell is generated automatically by siesta so that the 
!   different images do not overlap. In this case the stress returned
!   has no physical meaning.
! - The stress is defined as dE/d(strain)/Volume, with a positive sign
!   when the system tends to contract (negative pressure)
! - The following events result in a stopping error message:
!   - siesta_launch is called twice with the same label
!   - siesta_forces finds a communication error trough the pipes
!   - siesta_quit is called without a prior call to siesta_launch or
!     siesta_forces for that label
! - If siesta_quit is not called for a launched siesta process, that
!   process will stay listening indefinitedly to the pipe and will need
!   to be killed in the shell.
! - siesta_units may be called either before or after siesta_launch
! J.M.Soler and A.Garcia. Nov.2003

! *** Note ***
! Make sure that you have a working "flush" subroutine in your system,
! otherwise the process might hang.

#ifdef __NAG__
  use f90_unix_proc, only: system
#endif

  implicit none

PUBLIC :: siesta_launch, siesta_units, siesta_forces, siesta_quit

PRIVATE ! Nothing is declared public beyond this point

! Holds data on siesta processes and their communication pipes
  type proc
    private
    character(len=80) :: label     ! Name of process
    integer           :: iuc, iuf  ! I/O units for coords/forces commun.
  end type proc

! Global module variables
  integer, parameter :: max_procs = 100
  integer, parameter :: dp = kind(1.d0)
  type(proc),   save :: p(max_procs)
  integer,      save :: np=0
  character(len=32), save :: xunit = 'Ang'
  character(len=32), save :: eunit = 'eV'
  character(len=32), save :: funit = 'eV/Ang'
  character(len=32), save :: sunit = 'eV/Ang**3'

CONTAINS

!---------------------------------------------------

subroutine siesta_launch( label, nnodes, mpi_comm, launcher, localhost )
  implicit none
  character(len=*),         intent(in) :: label
  integer,         optional,intent(in) :: nnodes
  integer,         optional,intent(in) :: mpi_comm
  character(len=*),optional,intent(in) :: launcher
  logical,         optional,intent(in) :: localhost ! Not used in this version

  character(len=32) :: cpipe, fpipe
  character(len=80) :: task
  integer           :: ip, iu

!  print*, 'siesta_launch: launching process ', trim(label)

! Check that pipe does not exist already
  if (idx(label) /= 0) &
    print*, 'siesta_launch: ERROR: process for label ', trim(label), &
            ' already launched'

! Create pipes
  cpipe = trim(label)//'.coords'
  fpipe = trim(label)//'.forces'
  task = 'mkfifo '//trim(cpipe)//' '//trim(fpipe)
  call system(task)

! Open this side of pipes
  call open_pipes( label )

! Send wait message to coordinates pipe
  ip = idx( label )
  iu = p(ip)%iuc
  write(iu,*) 'wait'
  call pxfflush(iu)

! Start siesta process
  if (present(launcher)) then
    task = launcher
  elseif (present(nnodes) .and. nnodes>1) then
    write(task,*) ' mpirun -np ', nnodes, &
                  ' siesta < ', trim(label)//'.fdf > ', trim(label)//'.out &'
  else
    write(task,*) ' siesta < ', trim(label)//'.fdf > ', trim(label)//'.out &'
  endif
  print*,'siesta_launch: task = ',trim(task)
  call system(task)

end subroutine siesta_launch

!---------------------------------------------------

subroutine siesta_units( length, energy )
  implicit none
  character(len=*), intent(in) :: length, energy
  xunit = length
  eunit = energy
  funit = trim(eunit)//'/'//trim(xunit)
  sunit = trim(eunit)//'/'//trim(xunit)//'**3'
end subroutine siesta_units

!---------------------------------------------------

subroutine siesta_forces( label, na, xa, cell, energy, fa, stress )
  implicit none
  character(len=*),   intent(in) :: label
  integer,            intent(in) :: na
  real(dp),           intent(in) :: xa(3,na)
  real(dp), optional, intent(in) :: cell(3,3)
  real(dp), optional, intent(out):: energy
  real(dp), optional, intent(out):: fa(3,na)
  real(dp), optional, intent(out):: stress(3,3)

  integer           :: i, ia, ip, iu, n
  character(len=80) :: message
  real(dp)          :: e, f(3,na), s(3,3), c(3,3)

! Find system index
  ip = idx( label )
  if (ip==0) then
    call open_pipes( label )
    ip = idx( label )
  end if

! Copy unit cell
  if (present(cell)) then
    c = cell
  else
!**AG: Careful with this ...
    c = 0._dp
  end if

! Print coords for debugging
  print'(/,2a)',         'siesta_forces: label = ', trim(label)
  print'(3a,/,(3f12.6))','siesta_forces: cell (',trim(xunit),') =',c
  print'(3a,/,(3f12.6))','siesta_forces: xa (',trim(xunit),') =', xa

! Write coordinates to pipe
  iu = p(ip)%iuc
  write(iu,*) 'begin_coords'
  write(iu,*) trim(xunit)
  write(iu,*) trim(eunit)
  do i = 1,3
    write(iu,*) c(:,i)
  end do
  write(iu,*) na
  do ia = 1,na
    write(iu,*) xa(:,ia)
  end do
  write(iu,*) 'end_coords'
  call pxfflush(iu)

! Read forces from pipe
  iu = p(ip)%iuf
  read(iu,*) message
  if (message=='error') then
    read(iu,*) message
    call die( 'siesta_forces: siesta ERROR:' // trim(message))
  else if (message/='begin_forces') then
    call die('siesta_forces: ERROR: unexpected header:' // trim(message))
  end if
  read(iu,*) e
  read(iu,*) s
  read(iu,*) n
  if (n /= na) then
    print*, 'siesta_forces: ERROR: na mismatch: na, n =', na, n
    call die()
  end if
  do ia = 1,na
    read(iu,*) f(:,ia)
  end do
  read(iu,*) message
  if (message/='end_forces') then
    call die('siesta_forces: ERROR: unexpected trailer:' // trim(message))
  end if

! Print forces for debugging
  print'(3a,f12.6)',      'siesta_forces: energy (',trim(eunit),') =', e
  print'(3a,/,(3f12.6))', 'siesta_forces: stress (',trim(sunit),') =', s
  print'(3a,/,(3f12.6))', 'siesta_forces: forces (',trim(funit),') =', f

! Copy results to output arguments
  if (present(energy)) energy = e
  if (present(fa))     fa     = f
  if (present(stress)) stress = s

end subroutine siesta_forces

!---------------------------------------------------
subroutine siesta_quit( label )
  implicit none
  character(len=*), intent(in) :: label

  integer :: ip

  if (label == 'all') then
    ! Stop all siesta processes
    do ip = 1,np
      call siesta_quit_process( p(ip)%label )
    end do
  else
    ! Stop one siesta process
     call siesta_quit_process( label )
  endif

end subroutine siesta_quit

subroutine siesta_quit_process(label)
  implicit none
  character(len=*), intent(in) :: label

  integer :: ip, iuc, iuf
  character(len=20) message

    ip = idx(label)      ! Find process index
    if (ip==0) &
      call die('siesta_quit: ERROR: unknown label: ' // trim(label))

    iuc = p(ip)%iuc      ! Find cooordinates pipe unit
    write(iuc,*) 'quit'  ! Send quit signal through pipe
    call pxfflush(iuc)
    iuf = p(ip)%iuf                  ! Find forces pipe unit
    read(iuf,*) message              ! Receive response from pipe
    if (message == 'quitting') then  ! Check answer
      close(iuc,status="delete")     ! Close coordinates pipe
      close(iuf,status="delete")     ! Close forces pipe
      if (ip < np) then              ! Move last process to this slot
        p(ip)%label = p(np)%label
        p(ip)%iuc   = p(np)%iuc
        p(ip)%iuf   = p(np)%iuf
      end if
      np = np - 1                    ! Remove process
    else
      call die('siesta_quit: ERROR: unexpected response: ' // trim(message))
    end if ! (message)

end subroutine siesta_quit_process

!---------------------------------------------------

subroutine open_pipes( label )
  implicit none
  character(len=*), intent(in) :: label

  integer           :: iuc, iuf
  character(len=80) :: cpipe, fpipe

! Check that pipe does not exist already
  if (idx(label) /= 0) &
    print *, 'open_pipes: ERROR: pipes for ', trim(label), &
            ' already opened'

! Get io units for pipes
  call get_io_units( iuc, iuf )

! Open pipes
  cpipe = trim(label)//'.coords'
  fpipe = trim(label)//'.forces'
  open( unit=iuc, file=cpipe, form="formatted", &
        status="old", position="asis")
  open( unit=iuf, file=fpipe, form="formatted", &
        status="old", position="asis")

! Store data of process
  np = np + 1
  if (np > max_procs) then
    stop 'siesta_launch: ERROR: parameter max_procs too small'
  else
    p(np)%label = label
    p(np)%iuc = iuc
    p(np)%iuf = iuf
  end if

end subroutine open_pipes

!---------------------------------------------------

subroutine get_io_units( iuc, iuf)
! Finds two available I/O unit numbers
  implicit none
  integer, intent(out)  :: iuc, iuf

  integer :: i, ip
  logical :: unit_used

  iuc = 0
  iuf = 0
  do i = 10, 99
    inquire(unit=i,opened=unit_used)  ! This does not work with pipes
    do ip = 1,np
      unit_used = unit_used .or. i==p(ip)%iuc .or. i==p(ip)%iuf
    end do
    if (.not. unit_used) then
      if (iuc==0) then
        iuc = i
      else
        iuf = i
        return
      end if
    endif
  enddo
  stop 'fsiesta:get_io_units: ERROR: cannot find free I/O unit'
end subroutine get_io_units

!---------------------------------------------------

integer function idx( label )
! Finds which of the stored labels is equal to the input label
  implicit none
  character(len=*), intent(in) :: label
  integer :: i
  do i = 1,np
    if (label==p(i)%label) then
      idx = i
      return
    end if
  end do
  idx = 0
end function idx

!---------------------------------------------------
!
! Private copy of die
!
subroutine die(msg)
character(len=*), intent(in), optional :: msg

if (present(msg)) then
   write(6,*) msg
endif
STOP
end subroutine die

end module fsiesta
