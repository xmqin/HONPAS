! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!-----------------------------------------------------------------------------
!
! module siesta_master
!
! Keeps a repository of coordinates, forces, and other properties, and handles
! their transmission to/from a master program that uses siesta as a subroutine.
! Since siesta was written as a main program, not as a subroutine, the flow of
! data to/from the master program is not straighforward as subroutine arguments
! and it was simpler to adress it with this module. Furthermore, there are two
! communication modes available with the master program, through unix pipes and
! through MPI, what justifies the following flow organization:
! Using unix pipes (master+fsiesta_pipes compiled together):
!   master <--> fsiesta_pipes <--pipe--> iopipes <--> siesta_master <--> siesta
! Using MPI (master+siesta compiled together):
!   master <--> fsiesta_mpi <--> siesta_master <--> siesta
!
!-----------------------------------------------------------------------------
!
! Public procedures provided by this module:
!   coordsFromMaster     : Get atomic coordinates from master program
!   forcesToMaster       : Send atomic forces to master program
!   getForcesForMaster   : Get atomic forces to be sent to master program
!   getPropertyForMaster : Get other properties to be sent to master program
!   setPropertyForMaster : Set other properties to be sent to master program
!   setCoordsFromMaster  : Set atomic positions sent by master program
!   setMasterUnits       : Set physical units used by master program
!
! Public variables:
!   siesta_server        : Is siesta a server?
!   siesta_subroutine    : Is siesta a subroutine?
!   input_file           : Name of fdf data file
!
!-----------------------------------------------------------------------------
!
! Interfaces of public procedures:
!
! subroutine coordsFromMaster( na, xa, cell )
!   integer, intent(in) :: na        ! Number of atoms
!   real(dp),intent(out):: xa(3,na)  ! Atomic coordinates
!   real(dp),intent(out):: cell(3,3) ! Unit cell vectors
! end subroutine coordsFromMaster
!
! subroutine forcesToMaster( na, Etot, fa, stress )
!   integer, intent(in):: na          ! Number of atoms
!   real(dp),intent(in):: Etot        ! Total energy
!   real(dp),intent(in):: fa(3,na)    ! Atomic forces
!   real(dp),intent(in):: stress(3,3) ! Stress tensor
! end subroutine coordsFromMaster
!
! subroutine getForcesForMaster( na, fa, stress )
!   integer, intent(in) :: na          ! Number of atoms
!   real(dp),intent(out):: Etot        ! Total energy
!   real(dp),intent(out):: fa(3,na)    ! Atomic forces
!   real(dp),intent(out):: stress(3,3) ! Stress tensor
! end subroutine getForcesForMaster
!
! subroutine setCoordsFromMaster( na, xa, cell )
!   integer, intent(in):: na        ! Number of atoms
!   real(dp),intent(in):: xa(3,na)  ! Atomic coordinates
!   real(dp),intent(in):: cell(3,3) ! Unit cell vectors
! end subroutine setCoordsFromMaster
!
! subroutine setMasterUnits( xunit, eunit )
!   character(len=*),intent(in):: xunit  ! Physical unit of length
!   character(len=*),intent(in):: eunit  ! Physical unit of energy
! end subroutine setMasterUnits
!
! subroutine getPropertyForMaster( property, valueSize, value, units, error )
!   character(len=*),intent(in) :: property   ! Property name
!   integer,         intent(in) :: valueSize  ! Size of value array
!   real(dp),        intent(out):: value(:)   ! Property value(s)
!   character(len=*),intent(out):: units      ! Physical units
!   character(len=*),intent(out):: error      ! Error name
! end subroutine getPropertyForMaster
!
! subroutine setPropertyForMaster( property, valueSize, value, units )
!   character(len=*),intent(in) :: property   ! Property name
!   integer,         intent(in) :: valueSize  ! Size of value array
!   real(dp),        intent(in) :: value(*)   ! Property value(s)
!   character(len=*),intent(in) :: units      ! Physical units
! end subroutine setPropertyForMaster
!
!-----------------------------------------------------------------------------

MODULE siesta_master

! Used module parameters and procedures
  use precision, only: dp              ! Double precision real kind
  use sys,       only: die             ! Termination routine
  use fdf,       only: fdf_get         ! Reading fdf-options
  use fdf,       only: fdf_convfac     ! Conversion of physical units

  use iosockets, only: coordsFromSocket, forcesToSocket
  use iopipes,   only: coordsFromPipe  ! Read coordinates from pipe
  use iopipes,   only: forcesToPipe    ! Write forces to pipe

  implicit none

! Public procedures
PUBLIC :: &
  coordsFromMaster,     &! Get atomic coordinates from master program
  forcesToMaster,       &! Send atomic forces to master program
  getForcesForMaster,   &! Get atomic forces to be sent to master program
  getPropertyForMaster, &! Get other properties to be sent to master program
  setPropertyForMaster, &! Set other properties to be sent to master program
  setCoordsFromMaster,  &! Set atomic positions sent by master program
  setMasterUnits         ! Set physical units used by master program

! Public variables
  logical,public,save:: siesta_server     = .false. ! Is siesta a server?
  logical,public,save:: siesta_subroutine = .false. ! Is siesta a subroutine?
  character(len=132),public,save:: input_file = ' ' ! fdf data file

PRIVATE ! Nothing is declared public beyond this point

! Derived type to hold a physical property
  type propType
    private
    character(len=32):: name = ' '
    character(len=32):: units = ' '
    integer          :: size = 0
    real(dp),pointer :: value(:)
  end type propType

! Physical units used by siesta
  character(len=*), parameter :: siesta_xunit = 'Bohr'
  character(len=*), parameter :: siesta_eunit = 'Ry'
  character(len=*), parameter :: siesta_funit = 'Ry/Bohr'
  character(len=*), parameter :: siesta_sunit = 'Ry/Bohr**3'

! Global module variables and default values
  integer,     parameter:: maxProps = 100      ! Max number of properties
  integer,          save:: nProps   = 0        ! Number of stored properties
  integer,          save:: nAtoms   = 0        ! Number of atoms
  real(dp),         save:: ucell(3,3) = 0      ! Unit cell vectors
  real(dp),         save:: stressT(3,3) = 0    ! Stress tensor
  real(dp),         save:: energy = 0          ! Total energy
  real(dp),pointer, save:: coords(:,:)         ! Atomic coordinates
  real(dp),pointer, save:: forces(:,:)         ! Atomic forces
  type(propType),   save:: prop(maxProps)      ! Stored physical properties
  character(len=32),save:: master_xunit = 'Ang'       ! Master's unit of length 
  character(len=32),save:: master_eunit = 'eV'        ! Master's unit of energy
  character(len=32),save:: master_funit = 'eV/Ang'    ! Master's unit of force
  character(len=32),save:: master_sunit = 'eV/Ang**3' ! Master's unit of stress

CONTAINS

!-----------------------------------------------------------------------------

subroutine coordsFromMaster( na, xa, cell )

  implicit none
  integer, intent(in) :: na        ! Number of atoms
  real(dp),intent(out):: xa(3,na)  ! Atomic coordinates
  real(dp),intent(out):: cell(3,3) ! Unit cell vectors
  character*32 :: iface            ! Interface mode

! In the pipes version, this is where we first know that we are a server
  siesta_server = .true.

  if (siesta_subroutine) then
    if (na==nAtoms) then
      xa = coords
      cell = ucell
    else
      call die('coordsFromMaster: ERROR: number-of-atoms mismatch')
    endif
  else  
    iface = fdf_get( "Master.interface",  "pipe")   
    if ( iface == "pipe") call coordsFromPipe( na, xa, cell )
    if ( iface == "socket") call coordsFromSocket (na, xa, cell )
  end if ! (siesta_subroutine)
  
end subroutine coordsFromMaster

!-----------------------------------------------------------------------------

subroutine forcesToMaster( na, Etot, fa, stress )

  implicit none
  integer, intent(in):: na          ! Number of atoms
  real(dp),intent(in):: Etot        ! Total energy
  real(dp),intent(in):: fa(3,na)    ! Atomic forces
  real(dp),intent(in):: stress(3,3) ! Stress tensor
  character*32 :: iface            ! Interface mode

  if (siesta_subroutine) then
    if (na==nAtoms) then
      energy = Etot
      forces = fa
      stressT = stress
    else ! (na/=nAtoms)
      call die('coordsFromMaster: ERROR: number-of-atoms mismatch')
    endif ! (na==nAtoms)
  else  
    iface = fdf_get( "Master.interface",  "pipe")     
    if ( iface == "pipe") call forcesToPipe( na, Etot, fa, stress )
    if ( iface == "socket") call forcesToSocket( na, Etot, fa, stress )
  end if ! (siesta_subroutine)

end subroutine forcesToMaster

!-----------------------------------------------------------------------------

subroutine getForcesForMaster( na, Etot, fa, stress )

  implicit none
  integer, intent(in) :: na          ! Number of atoms
  real(dp),intent(out):: Etot        ! Total energy
  real(dp),intent(out):: fa(3,na)    ! Atomic forces
  real(dp),intent(out):: stress(3,3) ! Stress tensor

! Check that number of atoms is right, then copy data from repository
! and convert them to master's physical units
  if (na==nAtoms) then
    Etot   = energy  * fdf_convfac( siesta_eunit, master_eunit )
    fa     = forces  * fdf_convfac( siesta_funit, master_funit )
    stress = stressT * fdf_convfac( siesta_sunit, master_sunit )
  else
    call die('getForcesForMaster: ERROR: number-of-atoms mismatch')
  end if

end subroutine getForcesForMaster

!-----------------------------------------------------------------------------

subroutine setCoordsFromMaster( na, xa, cell )

  implicit none
  integer, intent(in):: na        ! Number of atoms
  real(dp),intent(in):: xa(3,na)  ! Atomic coordinates
  real(dp),intent(in):: cell(3,3) ! Unit cell vectors

! Allocate arrays for repository the first time
  if (nAtoms==0) then
    nAtoms = na
    allocate( coords(3,na), forces(3,na) )
  end if

! Check that number of atoms is right, then copy data to repository
! after converting them to siesta physical units
  if (na==nAtoms) then
    coords = xa   * fdf_convfac( master_xunit, siesta_xunit )
    ucell  = cell * fdf_convfac( master_xunit, siesta_xunit )
  else
    call die('setCoordsFromMaster: ERROR: number-of-atoms mismatch')
  end if

end subroutine setCoordsFromMaster

!-----------------------------------------------------------------------------

subroutine setMasterUnits( xunit, eunit )

  implicit none
  character(len=*),intent(in):: xunit  ! Physical unit of length
  character(len=*),intent(in):: eunit  ! Physical unit of energy

  master_xunit = xunit
  master_eunit = eunit
  master_funit = trim(eunit)//'/'//trim(xunit)         ! Unit of force
  master_sunit = trim(eunit)//'/'//trim(xunit)//'**3'  ! Unit of stress

end subroutine setMasterUnits

!-----------------------------------------------------------------------------

subroutine getPropertyForMaster( property, valueSize, value, units, error )

  implicit none
  character(len=*),intent(in) :: property   ! Property name
  integer,         intent(in) :: valueSize  ! Size of value array
  real(dp),        intent(out):: value(:)   ! Property value(s)
  character(len=*),intent(out):: units      ! Physical units
  character(len=*),intent(out):: error      ! Error name

  logical:: found
  integer:: iProp

  found = .false.
  do iProp = 1,nProps
    if (prop(iProp)%name == property) then
      if (prop(iProp)%size == valueSize) then
        value = prop(iProp)%value
        units = prop(iProp)%units
        error = ' '
      else
        error = 'wrong_size'
      end if
      found = .true.
      exit ! iProp loop
    end if
  end do ! iProp
  if (.not.found) error = 'unknown_property'

end subroutine getPropertyForMaster

!-----------------------------------------------------------------------------

subroutine setPropertyForMaster( property, valueSize, value, units )

  implicit none
  character(len=*),intent(in) :: property   ! Property name
  integer,         intent(in) :: valueSize  ! Size of value array
  real(dp),        intent(in) :: value(*)   ! Property value(s)
  character(len=*),intent(in) :: units      ! Physical units

  logical:: found
  integer:: iProp

! Check if property is already stored
  found = .false.
  do iProp = 1,nProps
    if (prop(iProp)%name == property) then
      found = .true.
      exit ! iProp loop
    end if
  end do ! iProp

! Increase the number of stored properties, if necessary
  if (.not.found) then
    nProps = nProps + 1
    iProp = nProps
  end if

! (Re)allocate the array to store value(s)
  if (prop(iProp)%size /= valueSize) then
    if (prop(iProp)%size /= 0) deallocate( prop(iProp)%value )
    allocate( prop(iProp)%value(valueSize) )
  end if

! Store data
  prop(iProp)%name  = property
  prop(iProp)%size  = valueSize
  prop(iProp)%value = value(1:valueSize)
  prop(iProp)%units = units

end subroutine setPropertyForMaster

END MODULE siesta_master

