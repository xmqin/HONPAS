! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module iopipes
! Handles communications between siesta-as-a-subroutine and a driver 
! program, transfering coordinates and forces trough Unix pipes.
! The routines that handle the other side of the communication are
! in module fsiesta.
! J.M.Soler and A.Garcia. Nov.2003

! *** Note ***
! Make sure that you have a working "flush" subroutine in your system,
! and that it is compiled-in in file pxf.F90 through the appropriate
! preprocessor symbols. Otherwise the process might hang.

use precision, only: dp
use parallel, only: IOnode
use fdf
use sys, only: die, bye
#ifdef MPI
      use mpi_siesta
#endif

  implicit none

  external :: io_assign, io_close

PUBLIC :: coordsFromPipe, forcesToPipe

PRIVATE  ! Nothing is declared public beyond this point

  character(len=*), parameter :: siesta_xunit = 'Bohr'
  character(len=*), parameter :: siesta_eunit = 'Ry'
  character(len=*), parameter :: siesta_funit = 'Ry/Bohr'
  character(len=*), parameter :: siesta_sunit = 'Ry/Bohr**3'

  integer,           save :: iuc, iuf
  character(len=32), save :: xunit, eunit

#ifdef MPI
  integer :: MPIerror
#endif

CONTAINS

subroutine coordsFromPipe( na, xa, cell )
! Reads coordinates from pipe
  implicit none
  integer,  intent(in)  :: na         ! Number of atoms
  real(dp), intent(out) :: xa(3,na)   ! Atomic coordinates (bohr)
  real(dp), intent(out) :: cell(3,3)  ! Lattice vectors (bohr)

  logical, save     :: firstTime = .true.
  integer           :: n
  character(len=80) :: fname, task

! Open this side of pipe
  if (firstTime .and. IOnode) then

    ! Get pipe name (SystemLabel MUST be equal to the fdf file prefix)
    fname = fdf_string('SystemLabel',' ')
    fname = trim(fname)//'.coords'

    ! Get pipe unit
    call io_assign( iuc )

    ! Open pipe
    open( unit=iuc, file=fname, form='formatted', status='old', &
          position='asis' )

    firstTime = .false.
  end if ! (firstTime .and. IOnode)

! Read coordinates from pipe
  if (IOnode) then
    do 
      read(iuc,*) task
      if (trim(task)/='wait') exit
    end do
  endif

#ifdef MPI
  call MPI_Bcast(task,80, MPI_Character,0,  &
                 MPI_Comm_World, MPIerror)
#endif

   if (trim(task)=='quit') then
      if (IOnode) then
         write(iuf,*) 'quitting'
         call pxfflush(iuf)
      endif
      call bye('coordsFromPipe: STOP requested by driver')

   else if (trim(task)=='begin_coords') then
      if (IONode) then
         read(iuc,*) xunit
         read(iuc,*) eunit
         read(iuc,*) cell
         read(iuc,*) n
         if (n /= na) call die('coordsFromPipe: ERROR: na mismatch')
         read(iuc,*) xa
         read(iuc,*) task
      endif
#ifdef MPI
      call MPI_Bcast(task,80, MPI_Character,0,  &
                     MPI_Comm_World, MPIerror)
#endif
      if (trim(task)=='end_coords') then
        if (IOnode) then
           print '(/,3a,/,(3f12.6))', &
          'coordsFromPipe: cell (',trim(xunit),') =', cell
           print '(  3a,/,(3f12.6))', &
          'coordsFromPipe: xa (',trim(xunit),') =', xa
          ! Convert coordinate units
          cell = cell * fdf_convfac( xunit, siesta_xunit )
          xa   = xa   * fdf_convfac( xunit, siesta_xunit )
        endif
#ifdef MPI
        call MPI_Bcast( cell(1,1), 9   , MPI_Double_Precision, 0, &
             MPI_Comm_World, MPIerror)
        call MPI_Bcast( xa(1,1)  , 3*na, MPI_Double_Precision, 0, &
             MPI_Comm_World, MPIerror)
#endif
      else
        call die('coordsFromPipe: ERROR: coords not complete')
      end if
    else
      call die('coordsFromPipe: ERROR: unknown task')
    end if ! task

end subroutine coordsFromPipe


subroutine forcesToPipe( na, energy, forces, stress )
! Writes stress and forces to pipe
  implicit none
  integer,  intent(in) :: na            ! Number of atoms
  real(dp), intent(in) :: energy        ! Total energy (Ry)
  real(dp), intent(in) :: forces(3,na)  ! Atomic forces (Ry/bohr)
  real(dp), intent(in) :: stress(3,3)   ! Stress tensor (Ry/bohr^3)

  logical, save     :: firstTime = .true.
  integer           :: i, ia
  character(len=80) :: fname, funit, sunit
  real(dp)          :: e, f(3,na), s(3,3)

! Open this side of pipe
  if (firstTime .and. IOnode) then

    ! Get pipe name
    fname = fdf_string('ThisFilePrefix',' ')
    if (fname==' ') fname = fdf_string('SystemLabel',' ')
    fname = trim(fname)//'.forces'

    ! Get pipe unit
    call io_assign( iuf )

    ! Open pipe
    open( unit=iuf, file=fname, form='formatted', status='old', &
          position='asis' )

    firstTime = .false.
  end if ! (firstTime .and. IOnode)

  if (IOnode) then

! Convert physical units
    funit = trim(eunit)//'/'//trim(xunit)
    sunit = trim(eunit)//'/'//trim(xunit)//'**3'
    e = energy * fdf_convfac( siesta_eunit, eunit )
    s = stress * fdf_convfac( siesta_sunit, sunit )
    f = forces * fdf_convfac( siesta_funit, funit )

! Print forces in output file
    print '(/,3a,f12.6)',    'forcesToPipe: energy (',trim(eunit),') =', e
    print '(3a,/,(3f12.6))', 'forcesToPipe: stress (',trim(sunit),') =', s
    print '(3a,/,(3f12.6))', 'forcesToPipe: forces (',trim(funit),') =', f
  endif

! Write forces to pipe
  if (IOnode) then
    write(iuf,*) 'begin_forces'
    write(iuf,*) e
    do i = 1,3
      write(iuf,*) s(:,i)
    end do
    write(iuf,*) na
    do ia = 1,na
      write(iuf,*) f(:,ia)
    end do
    write(iuf,*) 'end_forces'
    call pxfflush(iuf)
  end if ! IOnode

end subroutine forcesToPipe

end module iopipes
