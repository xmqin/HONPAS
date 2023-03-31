! *** Module: nao2gto_dumper ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Utilities to monitor sensitive NAO2GTO variables
!!
!! This module dumps snapshots of sensitive NAO2GTO variables to a netCDF
!! file named 'nao2gto_dump.nc'.
!!
!! \author Yann Pouillon
!!
!! \par History
!!      - 06.2018 Created [Yann Pouillon]
! *****************************************************************************
module nao2gto_dumper

  use, intrinsic :: iso_fortran_env, only: error_unit

#ifdef CDF
  use netcdf
#endif

  implicit none

  private

  ! Global constants --------------------------------------------------------

  !> Set default real precision (avoids dependencies on SIESTA)
  integer, parameter :: dp = kind(1.0d0)

  ! Global netCDF data ------------------------------------------------------

  !> netCDF ID of the dumper file 
  integer, private :: ncid = -1

  ! netCDF dimensions -------------------------------------------------------

  !> Shared storage space for dimension names
  character(len=5), parameter, private :: nc_dkeys(7) = &
&   ["xyz  ", "nx   ", "ny   ", "nz   ", "norb ", "nspin", "time "]

  !> Shared storage space for dimension values
  integer, target, private :: nc_dvals(7) = [-1, -1, -1, -1, -1, -1, -1]

  !> Shared storage space for dimension IDs
  integer, target, private :: nc_dimids(7) = [-1, -1, -1, -1, -1, -1, -1]

  ! Aliases for dimension values
  integer, pointer :: xyz => nc_dvals(1)     !< Number of spatial dimensions
  integer, pointer :: nx => nc_dvals(2)      !< Number of points along X
  integer, pointer :: ny => nc_dvals(3)      !< Number of points along Y
  integer, pointer :: nz => nc_dvals(4)      !< Number of points along Z
  integer, pointer :: norb => nc_dvals(5)    !< Number of orbitals
  integer, pointer :: nspin => nc_dvals(6)   !< Number of spin components
  integer, pointer :: ntime => nc_dvals(7)   !< Time steps (unlimited)

  ! Aliases for dimension IDs
  integer, pointer, private :: xyz_id => nc_dimids(1)     !< ID of xyz
  integer, pointer, private :: nx_id => nc_dimids(2)      !< ID of nx
  integer, pointer, private :: ny_id => nc_dimids(3)      !< ID of ny
  integer, pointer, private :: nz_id => nc_dimids(4)      !< ID of nz
  integer, pointer, private :: norb_id => nc_dimids(5)    !< ID of norb
  integer, pointer, private :: nspin_id => nc_dimids(6)   !< ID of nspin
  integer, pointer, private :: ntime_id => nc_dimids(7)   !< ID of ntime

  ! netCDF variables --------------------------------------------------------

  !> netCDF ID of the electronic density variable
  integer, private :: dens_varid = -1

  !> netCDF ID of the electronic density gradient variable
  integer, private :: grad_varid = -1

  !> netCDF ID of the density matrix variable
  integer, private :: dm_varid = -1

  ! Internal counters -------------------------------------------------------

  !> Current index of the time series for the electronic density
  integer, private :: istep_dens = 1

  !> Current index of the time series for the density matrix
  integer, private :: istep_dm = 1

  ! -------------------------------------------------------------------------

  public :: &
&   nao2gto_dumper_init, &
&   nao2gto_dumper_free, &
&   nao2gto_dump_density, &
&   nao2gto_dump_dm

contains

  ! ***************************************************************************
  !> \brief Constructor: initializes the NAO2GTO data dumper
  !!
  !! \note
  !!     If called multiple times, this routine will silently ignore
  !!     the excess calls.
  !!
  !! \author Yann Pouillon
  !!
  !! \param[in] nx0: number of grid points along X
  !! \param[in] ny0: number of grid points along Y
  !! \param[in] nz0: number of grid points along Z
  !! \param[in] norb0: number of orbitals
  !! \param[in] nspin0: number of spin degrees of freedom
  ! ***************************************************************************
  subroutine nao2gto_dumper_init(nx0, ny0, nz0, norb0, nspin0)

    implicit none

    ! Arguments
    integer, intent(in) :: nx0, ny0, nz0, norb0, nspin0

    ! Local variables
    integer :: idim

    ! -------------------------------------------------------------------------

#ifdef CDF
    ! Check that dimensions are unused
    if ( any(nc_dimids /= -1) ) then
      return
    end if

    ! Set dimension values
    nc_dvals(:) = [3, nx0, ny0, nz0, norb0, nspin0, NF90_UNLIMITED]

    ! Open file
    call netcdf_check( nf90_create("nao2gto_dump.nc", NF90_CLOBBER, ncid) )

    ! Define dimensions
    do idim = 1,size(nc_dkeys, 1)
      call netcdf_check( nf90_def_dim(ncid, trim(nc_dkeys(idim)), &
&       nc_dvals(idim), nc_dimids(idim)) )
    end do

    ! Define electronic density variables
    call netcdf_check( nf90_def_var(ncid, "dens_data", NF90_DOUBLE, &
&     [nx_id, ny_id, nz_id, nspin_id, ntime_id], dens_varid) )
    call netcdf_check( nf90_def_var(ncid, "grad_data", NF90_DOUBLE, &
&     [xyz_id, nx_id, ny_id, nz_id, nspin_id, ntime_id], grad_varid) )

    ! Define density matrix variables
    call netcdf_check( nf90_def_var(ncid, "dm_data", NF90_DOUBLE, &
&     [norb_id, norb_id, nspin_id, ntime_id], &
&     dm_varid) )

    ! Get ready to dump data
    call netcdf_check( nf90_enddef(ncid) )
#else
    write(unit=error_unit, fmt='(A)') &
&     "nao2gto_dumper_init: netCDF disabled => no dump"
#endif

  end subroutine nao2gto_dumper_init

  ! ***************************************************************************
  !> \brief Destructor: terminates the NAO2GTO data dumper
  !!
  !! \author Yann Pouillon
  !!
  ! ***************************************************************************
  subroutine nao2gto_dumper_free()

    implicit none

    ! -------------------------------------------------------------------------

#ifdef CDF
    ! Close netCDF file
    if ( ncid /= -1 ) then
      call netcdf_check( nf90_close(ncid) )
    end if
#else
    write(unit=error_unit, fmt='(A)') &
&     "nao2gto_dumper_free: netCDF disabled => no dump"
#endif

    ! Reset netCDF variables
    ncid = -1
    nc_dvals(:) = -1
    nc_dimids(:) = -1

  end subroutine nao2gto_dumper_free

  ! ***************************************************************************
  !> \brief Dumps the NAO2GTO density and its gradient, one snapshot at a time
  !!
  !! \author Yann Pouillon
  !!
  !! \param[in] dens_data: array(nx,ny,nz,nspin) containing the
  !!                       electronic density
  !! \param[in] grad_data: array(3,nx,ny,nz,nspin) containing the gradient
  !!                       of the electronic density
  ! ***************************************************************************
  subroutine nao2gto_dump_density(dens_data, grad_data)

    implicit none

    ! Arguments
    real(dp), intent(in) :: dens_data(nx,ny,nz,nspin)
    real(dp), intent(in) :: grad_data(xyz,nx,ny,nz,nspin)

    ! Local variables
    integer :: dens_count(6), dens_start(6)

    ! -------------------------------------------------------------------------

    dens_start = [1, 1, 1, 1, 1, istep_dens]
    dens_count = [xyz, nx, ny, nz, nspin, 1]

#ifdef CDF
    call netcdf_check( nf90_put_var(ncid, dens_varid, dens_data, &
&     start=dens_start(2:6), count=dens_count(2:6)) )
    call netcdf_check( nf90_put_var(ncid, grad_varid, grad_data, &
&     start=dens_start(1:6), count=dens_count(1:6)) )
#else
    write(unit=error_unit, fmt='(A)') &
&     "nao2gto_dump_density: netCDF disabled => no dump"
#endif

    istep_dens = istep_dens + 1

  end subroutine nao2gto_dump_density

  ! ***************************************************************************
  !> \brief Dumps the NAO2GTO density matrix, one snapshot at a time
  !!
  !! \author Yann Pouillon
  !!
  !! \param[in] dm_data: array(norb,norb,nspin) representing the density matrix
  ! ***************************************************************************
  subroutine nao2gto_dump_dm(dm_data)

    implicit none

    ! Arguments
    real(dp), intent(in) :: dm_data(norb, norb, nspin)

    ! Local variables
    integer :: dm_count(4), dm_start(4)

    ! -------------------------------------------------------------------------

    dm_start = [1, 1, 1, istep_dm]
    dm_count = [norb, norb, nspin, 1]

#ifdef CDF
    call netcdf_check( nf90_put_var(ncid, dm_varid, dm_data, &
&     start=dm_start, count=dm_count) )
#else
    write(unit=error_unit, fmt='(A)') &
&     "nao2gto_dump_dm: netCDF disabled => no dump"
#endif

    istep_dm = istep_dm + 1

  end subroutine nao2gto_dump_dm

  ! ***************************************************************************
  ! *** Utility routines                                                    ***
  ! ***************************************************************************

  ! ***************************************************************************
  !> \brief Checks the return value of netCDF calls
  !!
  !! \author Yann Pouillon
  !!
  !! \param[in] exitcode: return value of a netCDF routine
  ! ***************************************************************************
  subroutine netcdf_check(exitcode)

    implicit none

    ! Arguments
    integer, intent(in) :: exitcode

    ! -------------------------------------------------------------------------

#ifdef CDF
    if ( exitcode /= NF90_NOERR ) then
      write(unit=error_unit, fmt='(A/4X,A)') &
&       "netcdf_check: Error:", trim(nf90_strerror(exitcode))
      stop 1
    end if
#else
    write(unit=error_unit, fmt='(A)') &
&     "netcdf_check: netCDF disabled => no error checking"
#endif

  end subroutine netcdf_check

end module nao2gto_dumper
