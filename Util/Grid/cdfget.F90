! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_grid
#ifdef CDF

implicit none

integer, parameter  :: dp = selected_real_kind(14,100)
integer, parameter  :: sp = selected_real_kind(6,20)

public :: get_cdf_grid, put_cdf_grid

private
type, public :: grid_t
  integer   :: n(3)
  real(dp)  :: cell(3,3)
  real(dp), pointer  :: grid(:,:,:) => null()
end type

CONTAINS

subroutine put_cdf_grid(gp,filename)
use netcdf

implicit none

character(len=*), intent(in) :: filename
type(grid_t), intent(in)    :: gp

integer  :: ncid 
integer  :: xyz_id, step_id, abc_id, spin_id
integer  :: cell_id, n1_id, n2_id, n3_id, gridfunc_id

integer   ::    nspin  ! Number of spins 
integer   ::    n(3), ig, jg, kg, i, j, k, ip
integer   ::   ispin, iostat, ix, iy, iz, iret

!-----------------------------------------------------

nspin = 1
call check( nf90_create(filename,NF90_CLOBBER,ncid))
       iret = nf90_def_dim(ncid,'xyz',3,xyz_id)
       call check(iret)
       iret = nf90_def_dim(ncid,'abc',3,abc_id)
       call check(iret)
       iret = nf90_def_dim(ncid,'spin',nspin,spin_id)
       call check(iret)

       n(:) = gp%n(:)

       iret = nf90_def_dim(ncid,'n1',n(1),n1_id)
       call check(iret)
       iret = nf90_def_dim(ncid,'n2',n(2),n2_id)
       call check(iret)
       iret = nf90_def_dim(ncid,'n3',n(3),n3_id)
       call check(iret)

       iret = nf90_def_var(ncid,'cell',nf90_float,(/xyz_id,abc_id/),cell_id)
       call check(iret)
       iret = nf90_put_att(ncid,cell_id,'Description', &
               "Cell vectors in Bohr: xyz, abc")
       call check(iret)

       iret = nf90_def_var(ncid,'gridfunc',nf90_float,(/n1_id,n2_id,n3_id,spin_id/),gridfunc_id)
       call check(iret)
       iret = nf90_put_att(ncid,gridfunc_id,'Description', &
               "Grid function -- ")
       call check(iret)

       iret = nf90_enddef(ncid)
       call check(iret)
!
       iret = nf90_put_var(ncid, cell_id, gp%cell, start = (/1, 1 /), count = (/3, 3/) )
       call check(iret)

   ispin = 1
      iret = nf90_put_var(ncid, gridfunc_id, gp%grid, start = (/1, 1, 1, ispin /), &
           count = (/n(1), n(2), n(3), 1/) )
      call check(iret)
   
     call check( nf90_close(ncid) )

end subroutine put_cdf_grid
!--------
subroutine get_cdf_grid(filename,gp)
use netcdf

implicit none

character(len=*), intent(in) :: filename
type(grid_t), intent(out)    :: gp




integer  :: ncid 
integer  :: xyz_id, step_id, abc_id, spin_id
integer  :: cell_id, n1_id, n2_id, n3_id, gridfunc_id

integer   ::    nspin  ! Number of spins 
integer   ::    n(3), ig, jg, kg, i, j, k, ip
integer   ::   ispin, iostat, ix, iy, iz, iret

real(dp)  ::    cell(3,3), cell_volume

!-----------------------------------------------------

call check( nf90_open(filename,NF90_NOWRITE,ncid))
       call check( nf90_inq_dimid(ncid,'spin',spin_id) )
       call check( nf90_inquire_dimension(ncid, dimid=spin_id, len=nspin) )

       call check( nf90_inq_dimid(ncid,'n1',n1_id) )
       call check( nf90_inquire_dimension(ncid, dimid=n1_id, len=n(1)) )
       call check( nf90_inq_dimid(ncid,'n2',n2_id) )
       call check( nf90_inquire_dimension(ncid, dimid=n2_id, len=n(2)) )
       call check( nf90_inq_dimid(ncid,'n3',n3_id) )
       call check( nf90_inquire_dimension(ncid, dimid=n3_id, len=n(3)) )

       call check( nf90_inq_varid(ncid, "cell", cell_id) )
       call check( nf90_inq_varid(ncid, "gridfunc", gridfunc_id) )

       iret = nf90_get_var(ncid, cell_id, cell, start = (/1, 1 /), &
                        count = (/3, 3/) )
       call check(iret)

   gp%n(:) = n(:)
   gp%cell = cell

   allocate(gp%grid(n(1),n(2),n(3)))

   ispin = 1
      iret = nf90_get_var(ncid, gridfunc_id, gp%grid, start = (/1, 1, 1, ispin /), &
           count = (/n(1), n(2), n(3), 1/) )
      call check(iret)
   
     call check( nf90_close(ncid) )

end subroutine get_cdf_grid

subroutine check(code)
use netcdf, only: nf90_noerr, nf90_strerror
integer, intent(in) :: code
if (code /= nf90_noerr) then
  print *, "netCDF error: " // NF90_STRERROR(code)
  STOP
endif
end subroutine check

#endif /* CDF */
end module m_grid

!--------------------
program test
#ifndef CDF
  print *, "No netCDF support at compile time"
#else

use m_grid

type(grid_t)   :: gp

call get_cdf_grid("Rho.grid.nc",gp)

print *, gp%cell
print *, gp%n
print *, gp%grid(1,1,1)

call put_cdf_grid(gp,"Copy.Of.Rho.grid.nc")

#endif /* CDF */
end program test

