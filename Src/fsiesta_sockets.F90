! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! module fsiesta
!
! Support routines for siesta-as-a-subroutine.
! The routines that handle the other side of the communication are
! in module iosockets of siesta program.
! Used modules
!   f90sockets : support for socket communications (in file fsockets.f90)
! Usage:
!   call siesta_launch( label, nnodes, mpi_comm, launcher, localhost )
!     character(len=*),intent(in) :: label  : Name of siesta process
!                                             (prefix of its .fdf file)
!     integer,optional,intent(in) :: nnodes : Number of MPI nodes
!     integer,optional,intent(in) :: mpi_comm : not used in this version
!     character(len=*),intent(in),optional:: launcher : full launch command
!     logical,optional,intent(in) :: localhost : will siesta run at localhost?
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
! - localhost=.true. is assumed by default, if not present at siesta_launch,
!   or if siesta_launch is not called. 
! - If localhost=.false., the IP address of the driver program's host must 
!   be given as Master.address in siesta .fdf file.
! - If siesta_units is not called, length='Ang', energy='eV' are
!   used by default. If it is called more than once, the units in the
!   last call become in effect.
! - The physical units set by siesta_units are used for all the siesta
!   processes launched
! - If siesta_forces is called without a previous call to siesta_launch
!   for that label, it assumes that the siesta process will be launched
!   externally in the shell (AFTER the driver program calls siesta_forces).
!   In this case, siesta_forces creates the socket and waits until siesta
!   connects to it.
! - If argument cell is not present in the call to siesta_forces, or if
!   the cell has zero volume, it is assumed that the system is a molecule,
!   and a supercell is generated automatically by siesta so that the 
!   different images do not overlap. In this case the stress returned
!   has no physical meaning.
! - The stress is defined as dE/d(strain)/Volume, with a positive sign
!   when the system tends to contract (negative pressure)
! - The following events result in a stopping error message:
!   - siesta_launch is called twice with the same label
!   - siesta_forces finds a communication error trough the socket
!   - siesta_quit is called without a prior call to siesta_launch or
!     siesta_forces for that label
! - siesta_units may be called either before or after siesta_launch
! J.M.Soler, A.Garcia, and M.Ceriotti. Nov.2003, Mar.2015

MODULE fsiesta

use f90sockets, only: create_socket, writebuffer, readbuffer
#ifdef __NAG__
  use f90_unix_proc, only: system
#endif

  implicit none

PUBLIC :: siesta_launch, siesta_units, siesta_forces, siesta_quit

PRIVATE ! Nothing is declared public beyond this point

! Derived type to hold data on siesta processes and their communication sockets
  type proc
    private
    character(len=256) :: label     ! name of siesta process
    character(len=1024):: host      ! host IP address of the siesta process
    integer            :: inet      ! communication socket type
    integer            :: port      ! port number for the commun. socket
    integer            :: socket    ! socket id
  end type proc

! Global module variables
! WARNING: MSGLEN must be equal in siesta module iosockets
  integer, parameter :: MSGLEN = 80      ! string length of socket messages
  integer, parameter :: unit_len = 32    ! string length of physical units
  integer, parameter :: max_procs = 100  ! max. simultaneous siesta processes
  integer, parameter :: dp = kind(1.d0)  ! double precision real kind
  type(proc),   save :: p(max_procs)     ! data of siesta processes
  integer,      save :: np = 0           ! present number of siesta processes    
  integer,      save :: totp = 0         ! total siesta processes ever started

! Default driver's physical units (reset by siesta_units)
  character(len=unit_len), save :: xunit = 'Ang'       ! length unit
  character(len=unit_len), save :: eunit = 'eV'        ! energy unit
  character(len=unit_len), save :: funit = 'eV/Ang'    ! force unit
  character(len=unit_len), save :: sunit = 'eV/Ang**3' ! stress unit

CONTAINS

!---------------------------------------------------

subroutine siesta_launch( label, nnodes, mpi_comm, launcher, localhost )
  implicit none
  character(len=*),          intent(in) :: label
  integer,         optional, intent(in) :: nnodes
  integer,         optional, intent(in) :: mpi_comm
  character(len=*),optional, intent(in) :: launcher
  logical,         optional, intent(in) :: localhost

  character(len=1024):: task
  integer            :: ip, isocket

  print*, 'siesta_launch: launching process ', trim(label)

! Check that siesta process does not exist already
  if (idx(label) /= 0) &
    print*, 'siesta_launch: ERROR: process for label ', trim(label), &
            ' already launched'

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

! Create and open new socket for communication with siesta process
  call open_new_socket(label,localhost)

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

  integer              :: i, ia, ip, n, socket
  character(len=MSGLEN):: message
  real(dp)             :: e, f(3*na), s(9), c(9)

! Find process index
  ip = idx( label )
  if (ip==0) then
    call open_new_socket( label )
    ip = idx( label )
  end if

! Copy unit cell
  if (present(cell)) then
    c = reshape(cell,(/9/))
  else
!**AG: Careful with this ...
    c = 0._dp
  end if

! Print coords for debugging
  print'(/,2a)',         'siesta_forces: label = ', trim(label)
  print'(3a,/,(3f12.6))','siesta_forces: cell (',trim(xunit),') =',c
  print'(3a,/,(3f12.6))','siesta_forces: xa (',trim(xunit),') =', xa

! Write coordinates to socket
  socket = p(ip)%socket
  call writebuffer( socket, 'begin_coords' )
  call writebuffer( socket, xunit )
  call writebuffer( socket, eunit )
  call writebuffer( socket, c, 9 )
  call writebuffer( socket, na )
  call writebuffer( socket, reshape(xa,(/3*na/)), 3*na )
  call writebuffer( socket, 'end_coords' )

! Read forces from socket
  call readbuffer( socket, message )
  if (message=='error') then
    call readbuffer( socket, message )
    call die( 'siesta_forces: siesta ERROR:' // trim(message))
  else if (message/='begin_forces') then
    call die('siesta_forces: ERROR: unexpected header:' // trim(message))
  end if
  call readbuffer( socket, e )
  call readbuffer( socket, s, 9 )
  call readbuffer( socket, n )
  if (n /= na) then
    print*, 'siesta_forces: ERROR: na mismatch: na, n =', na, n
    call die()
  end if
  call readbuffer( socket, f, 3*na )
  call readbuffer( socket, message )
  if (message/='end_forces') then
    call die('siesta_forces: ERROR: unexpected trailer:' // trim(message))
  end if

! Print forces for debugging
  print'(3a,f12.6)',      'siesta_forces: energy (',trim(eunit),') =', e
  print'(3a,/,(3f12.6))', 'siesta_forces: stress (',trim(sunit),') =', s
  print'(3a,/,(3f12.6))', 'siesta_forces: forces (',trim(funit),') =', f

! Copy results to output arguments
  if (present(energy)) energy = e
  if (present(fa))     fa     = reshape(f,(/3,na/))
  if (present(stress)) stress = reshape(s,(/3,3/))

end subroutine siesta_forces

!---------------------------------------------------

subroutine siesta_quit( label )
  implicit none
  character(len=*), intent(in) :: label

  integer :: ip

  if (label == 'all') then
    ! Stop all siesta processes
    do ip = np,1,-1
      call siesta_quit_process( p(ip)%label )
    end do
  else
    ! Stop one siesta process
    call siesta_quit_process( label )
  endif

end subroutine siesta_quit

!---------------------------------------------------

subroutine siesta_quit_process(label)
  implicit none
  character(len=*), intent(in) :: label

  integer :: ip, socket
  character(len=MSGLEN):: message

    ip = idx(label)      ! Find process index
    if (ip==0) &
      call die('siesta_quit: ERROR: unknown label: ' // trim(label))

    print*,'siesta_quit: stopping siesta process ',trim(label)
    socket = p(ip)%socket
    call writebuffer(socket,'quit')  ! Send quit signal to server
    call readbuffer(socket,message)  ! Receive response
    if (message == 'quitting') then  ! Check answer
      if (ip < np) then              ! Move last process to this slot
        p(ip)%label  = p(np)%label
        p(ip)%host   = p(np)%host
        p(ip)%inet   = p(np)%inet
        p(ip)%port   = p(np)%port
        p(ip)%socket = p(np)%socket
      end if
      np = np - 1                    ! Remove process
    else
      call die('siesta_quit: ERROR: unexpected response: ' // trim(message))
    end if ! (message)

end subroutine siesta_quit_process

!---------------------------------------------------

subroutine open_new_socket( label, localhost )
  implicit none
  character(len=*),intent(in) :: label
  logical,optional,intent(in) :: localhost

  integer,parameter:: iu=87
  character(len=32):: host
  logical          :: local

! Check that process does not exist already
  if (idx(label) /= 0) &
    print*, 'fsiesta ERROR: siesta process ', trim(label),' exists already'

! Get my IP address, unless socket communication is local
  if (present(localhost)) then
    local = localhost
  else
    local = .true.
  endif
  if (local) then
    host = 'localhost'
  else
    call system("ifconfig eth0 | awk '/inet / { print $2 }'" // &
                "| sed 's/addr://' > my_ip_addr")
    open(iu,file='my_ip_addr')
    read(iu,*) host
    close(iu)
!    call system("rm -f my_ip_addr")
  endif

! Store data of process
  np = np+1              ! present number of processes
  totp = totp+1          ! total number of processes that were created
  if (np > max_procs) then
    stop 'fsiesta ERROR: parameter max_procs too small'
  else
    p(np)%label = label
    p(np)%host = host
    p(np)%inet = 1
    p(np)%port = 10000+totp
  end if

! Create socket
  call create_socket( p(np)%socket, p(np)%inet, p(np)%port, p(np)%host )

end subroutine open_new_socket

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

subroutine die(msg)
  character(len=*), intent(in), optional :: msg
  if (present(msg)) write(6,*) trim(msg)
  stop
end subroutine die

END MODULE fsiesta

