! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

!*****************************************************************************
! module iosockets
!*****************************************************************************
! Handles communications between siesta-as-a-subroutine and a driver 
! program, transfering coordinates and forces trough POSIX sockets.
! The routines that handle the other side of the communication are
! in file fsiesta_sockets.f90 and/or in the i-PI code.
! J.M.Soler, A.Garcia, and M.Ceriotti, 2014
!*****************************************************************************
!
!   Modules used:
! use precision,    only: dp
! use parallel,     only: IOnode
! use fdf
! use f90sockets,   only: open_socket, writebuffer, readbuffer
! use sys,          only: die, bye
! use m_mpi_utils,  only: broadcast
! use cellSubs,     only: volcel
! use mpi_siesta
!
!   Public derived types:
! none
!
!   Public variables and arrays:
! none
!
!   Public procedures:
! coordsFromSocket, &! receive coords and cell vectors from master program
! forcesToSocket     ! send energy, forces, and stress to master program
!
!*****************************************************************************
! subroutine coordsFromSocket( na, xa, cell )
!   integer,  intent(in) :: na            ! Number of atoms
!   real(dp), intent(out):: xa(3,na)      ! Atomic coordinates (bohr)
!   real(dp), intent(out):: cell(3,3)     ! Lattice vectors (bohr)
!*****************************************************************************
! subroutine forcesToSocket( na, energy, forces, stress )
!   integer,  intent(in) :: na            ! Number of atoms
!   real(dp), intent(in) :: energy        ! Total energy (Ry)
!   real(dp), intent(in) :: forces(3,na)  ! Atomic forces (Ry/bohr)
!   real(dp), intent(in) :: stress(3,3)   ! Stress tensor (Ry/bohr^3)
!*****************************************************************************

module iosockets

! Used modules
  use precision,    only: dp
  use parallel,     only: IOnode
  use fdf
  use f90sockets,   only: open_socket, writebuffer, readbuffer, close_socket
  use sys,          only: die, bye
  use m_mpi_utils,  only: broadcast
  use cellSubs,     only: volcel
#ifdef MPI
  use mpi_siesta
#endif

  implicit none

! Public procedures
PUBLIC :: &
  coordsFromSocket, &! receive coords and cell vectors from master program
  forcesToSocket     ! send energy, forces, and stress to master program

PRIVATE  ! Nothing is declared public beyond this point

! Private module parameters
! WARNING: IPI_MSGLEN and MSGLEN must be equal in i-PI and fsiesta, respectively
  integer,                 parameter ::   IPI_MSGLEN = 12
  integer,                 parameter ::       MSGLEN = 80
  integer,                 parameter ::     unit_len = 32
  character(len=unit_len), parameter :: siesta_xunit = 'bohr'
  character(len=unit_len), parameter :: siesta_eunit = 'ry'
  character(len=unit_len), parameter ::    ipi_xunit = 'bohr'
  character(len=unit_len), parameter ::    ipi_eunit = 'hartree'

! Private module variables
  integer, save :: socket
  real(dp),save :: cellv
  character(len=32),      save:: master='unknown', stype='unknown'
  character(len=unit_len),save:: master_eunit='unknown', &
                                 master_xunit='unknown'
#ifdef MPI
  integer :: MPIerror
#endif

CONTAINS

!*****************************************************************************

subroutine coordsFromSocket( na, xa, cell )
  ! Reads coordinates from socket
  implicit none
  integer,  intent(in)  :: na         ! Number of atoms
  real(dp), intent(out) :: xa(3,na)   ! Atomic coordinates (bohr)
  real(dp), intent(out) :: cell(3,3)  ! Lattice vectors (bohr)

! Local parameters
  character(len=*),parameter:: myName='coordsFromSocket '

! Local variables and arrays
  logical, save            :: firstTime = .true.
  CHARACTER(len=MSGLEN)    :: header, message
  CHARACTER(len=1024)      :: host
  INTEGER                  :: inet, port, n
  REAL(dp)                 :: aux(9), c(9)
  REAL(dp),ALLOCATABLE     :: x(:)

! Open the socket, shared by the receive and send sides of communication
  if (firstTime) then
    master = fdf_get( "Master.code",  "fsiesta")
    host = fdf_get(   "Master.address",  "localhost")
    port = fdf_get(   "Master.port",  10001)
    stype = fdf_get(  "Master.socketType",  "inet")

    if (leqi(stype,'unix')) then
      inet = 0
    elseif (leqi(stype,'inet')) then
      inet = 1
    else
      call die(myName//'ERROR: unknown socket type '//trim(stype))
    endif

    if (IOnode) then
      write(*,'(/,a,i4,i8,2x,a)') &
          myName//'opening socket: inet,port,host=',inet,port,trim(host)
      call open_socket(socket, inet, port, host)
    endif

    firstTime = .false.
  end if ! (firstTime .and. IOnode)

! Read header from socket
  header = ''
  if (IOnode) then
    do
      if (leqi(master,'i-pi')) then
        call readbuffer(socket, header, IPI_MSGLEN)

        ! Immediately stop if requested
        if (trim(header)=='EXIT') then ! we are done!
          call close_socket(socket)
          call bye(myName//'STOP requested by driver')
        end if
        
        if (trim(header)/='STATUS') exit ! do loop
        message = 'READY'
        call writebuffer(socket, message, IPI_MSGLEN)
      elseif (leqi(master,'fsiesta')) then
        call readbuffer(socket, header)

        ! Immediately stop if requested
        if (trim(header)=='quit') then
          call writebuffer(socket,'quitting')
          call close_socket(socket)
          call bye(myName//'STOP requested by driver')
        elseif (trim(header)=='wait') then
          call writebuffer(socket,'ready')
          cycle ! do loop
        else
          exit ! do loop
        endif
      else
        call die(myName//'ERROR: unknown master')
      endif ! leqi(master)
    enddo
  endif ! IOnode
#ifdef MPI
  call MPI_Bcast(header, MSGLEN, MPI_Character, 0, MPI_Comm_World, MPIerror)
#endif

! Read cell vectors from socket, as a single buffer vector
  if (IOnode) then
    if (leqi(master,"i-pi") .and. trim(header)=="POSDATA") then 
      call readbuffer(socket, c, 9)
      ! Please see below for the cell array
      ! If this is to be used, it should *ALSO* be transposed!
      call readbuffer(socket, aux, 9)    ! not used in siesta
      master_xunit = ipi_xunit
      master_eunit = ipi_eunit
    elseif (leqi(master,"fsiesta") .and. trim(header)=="begin_coords") then
      call readbuffer(socket, master_xunit)
      call readbuffer(socket, master_eunit)
      call readbuffer(socket, c, 9)
    else
      call die(myName//'ERROR: unexpected header: '//trim(header))
    end if
  endif

! Broadcast and copy cell from buffer
#ifdef MPI
  call MPI_Bcast(master_xunit, unit_len, MPI_Character,0,  &
                 MPI_Comm_World, MPIerror)
  call MPI_Bcast(master_eunit, unit_len, MPI_Character,0,  &
                 MPI_Comm_World, MPIerror)
  call MPI_Bcast(c,9, MPI_Double_Precision,0, MPI_Comm_World, MPIerror)
#endif
  ! I-Pi assumes row-major order of *only* the cell parameters
  ! So we need to transpose the cell
  ! This is a very bad practice since the same is happening in
  ! ASE I-pi implementation.
  ! The cell is transposed when transfering, but not the coordinates.
  ! Essentially one could leave out *any* transposes in all implementations
  ! and it would work! However, this is now a legacy implementation
  ! that should not be changed, neither in I-pi, nor ASE.
  cell = TRANSPOSE( RESHAPE( c, (/3,3/) ) )

! Read and check number of atoms
  if (ionode) then
    call readbuffer(socket, n)
    if (n/=na) call die(myName//'ERROR: unexpected number of atoms')
  endif

! Read and broadcast atomic coordinates
  allocate(x(3*na))
  if (IOnode) call readbuffer(socket, x, 3*na)
#ifdef MPI
  call MPI_Bcast(x,3*na, MPI_Double_Precision,0, MPI_Comm_World, MPIerror)
#endif
  xa = RESHAPE( x, (/3,na/) )
  deallocate(x)

! Read trailing message
  if (IOnode .and. leqi(master,'fsiesta')) then
    call readbuffer(socket, message)
    if (message/='end_coords') then
      call die(myName//'ERROR: unexpected trailer:'//trim(message))
    end if
  endif

  if ( IONode ) then
    ! Print coordinates and cell vectors received
    write(*,'(/,4a,/,(3f12.6))') myName,'cell (',trim(master_xunit),') =', cell
    write(*,'(4a,/,(3f12.6))') myName,'coords (',trim(master_xunit),') =', xa
  end if

! Convert physical units
  cell = cell * fdf_convfac( master_xunit, siesta_xunit )
  xa   = xa   * fdf_convfac( master_xunit, siesta_xunit )

! Compute and save the cell volume, to be used by forcesToSocket
  cellv = volcel(cell)
   
end subroutine coordsFromSocket

!*****************************************************************************

subroutine forcesToSocket( na, energy, forces, stress )
! Writes stress and forces to socket
  implicit none
  integer,  intent(in) :: na            ! Number of atoms
  real(dp), intent(in) :: energy        ! Total energy (Ry)
  real(dp), intent(in) :: forces(3,na)  ! Atomic forces (Ry/bohr)
  real(dp), intent(in) :: stress(3,3)   ! Stress tensor (Ry/bohr^3)

! Local parameters
  character(len=*),parameter:: myName='forcesToSocket '

! Local variables and arrays
  character(len=MSGLEN):: header, message
  real(dp)             :: e, s(9), vir(9)
  real(dp),allocatable :: f(:)

! Copy input to local variables
  allocate(f(3*na))
  e = energy
  f(:) = reshape( forces, (/3*na/) )
  s(:) = reshape( stress, (/9/) )
  
! Convert physical units
  e = e * fdf_convfac( siesta_eunit, master_eunit )
  f(:) = f(:) * fdf_convfac( siesta_eunit, master_eunit ) &
        / fdf_convfac( siesta_xunit, master_xunit )
  s(:) = s(:) * fdf_convfac( siesta_eunit, master_eunit ) &
        / fdf_convfac( siesta_xunit, master_xunit )**3

! Find virial tensor for i-pi, in master's units
  vir(:) = -s * cellv * fdf_convfac( siesta_xunit, master_xunit )**3

! Write forces to socket
  if (IOnode) then
    if (leqi(master,'i-pi')) then
      do
        call readbuffer(socket, header, IPI_MSGLEN)
        if (trim(header)=='STATUS') then        ! inform i-pi of my status
          message = 'HAVEDATA'
          call writebuffer(socket,message,IPI_MSGLEN)
          cycle ! do loop
        elseif (trim(header)=='GETFORCE') then  ! proceed to send forces
          exit ! do loop
        else
          call die(myName//'ERROR: unexpected header from i-pi')
        endif
      enddo
      message = 'FORCEREADY'
      call writebuffer(socket,message,IPI_MSGLEN)
      call writebuffer(socket,e)
      call writebuffer(socket,na)
      call writebuffer(socket,f,3*na)
      call writebuffer(socket,vir,9)       
      ! i-pi can also receive an arbitrary string, that will be printed 
      ! out to the "extra" trajectory file. This is useful if you want 
      ! to return additional information, e.g. atomic charges, wannier 
      ! centres, etc. one must return the number of characters, then
      ! the string. here we just send back zero characters.
      call writebuffer(socket,0)
    elseif (leqi(master,'fsiesta')) then
!      do
!        call readbuffer(socket, header)
!        if (trim(header)/='wait') exit
!      enddo
      call writebuffer(socket,'begin_forces')
      call writebuffer(socket,e)
      call writebuffer(socket,s,9)
      call writebuffer(socket,na)
      call writebuffer(socket,f,3*na)
      call writebuffer(socket,'end_forces')
    else
      call die(myName//'ERROR: unknown master')
    endif ! leqi(master)
  end if ! IOnode

  if ( IONode ) then
    ! Print energy, forces, and stress tensor sent to master
    write(*,'(/,a,f12.6)') myName// 'energy ('//trim(master_eunit)//') =', e
    write(*,'(a,/,(3f12.6))') myName//'stress ('//trim(master_eunit)//'/'//trim(master_xunit)//'^3) =', s
    write(*,'(a,/,(3f12.6))') myName//'forces ('//trim(master_eunit)//'/'//trim(master_xunit)//') =', f
  end if
  deallocate(f)

end subroutine forcesToSocket

end module iosockets
