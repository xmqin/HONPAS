! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
program cdf2grid
#ifndef CDF
  print *, "No netCDF support at compile time"
#else

!
! Converts a Grid (binary) file to netCDF format
!
! Usage: cdf2grid
!
!        The input file is called "GridFunc.nc" (a symbolic link would do)
!        The output file is named "GRIDFUNC" 
!
use netcdf

implicit none

integer, parameter  :: dp = selected_real_kind(14,100)
integer, parameter  :: sp = selected_real_kind(6,20)

integer  :: ncid 
integer  :: xyz_id, step_id, abc_id, nspin_id
integer  :: cell_id, n1_id, n2_id, n3_id, gridfunc_id

integer   ::    nspin  ! Number of spins 
integer   ::    mesh(3)

integer   ::   ispin, iostat, ix, iy, iz, iret

real(dp)  ::    cell(3,3)
real(sp), dimension(:), allocatable :: gridfunc

!-----------------------------------------------------

open(unit=1,file="GRIDFUNC",form="unformatted",status="unknown",action="write", &
            position="rewind",iostat=iostat)
if (iostat /= 0) then
  print *, "File GRIDFUNC cannot be opened"
  STOP
endif


call check( nf90_open('GridFunc.nc',NF90_NOWRITE,ncid))
     
  call check( nf90_inq_dimid(ncid,'spin',nspin_id) )
       call check( nf90_inquire_dimension(ncid, dimid=nspin_id, len=nspin) )
       call check( nf90_inq_dimid(ncid,'n1',n1_id) )
       call check( nf90_inquire_dimension(ncid, dimid=n1_id, len=mesh(1)) )
       call check( nf90_inq_dimid(ncid,'n2',n2_id) )
       call check( nf90_inquire_dimension(ncid, dimid=n2_id, len=mesh(2)) )
       call check( nf90_inq_dimid(ncid,'n3',n3_id) )
       call check( nf90_inquire_dimension(ncid, dimid=n3_id, len=mesh(3)) )

       call check( nf90_inq_varid(ncid, "cell", cell_id) )
       call check( nf90_inq_varid(ncid, "gridfunc", gridfunc_id) )

       iret = nf90_get_var(ncid, cell_id, cell, start = (/1, 1 /), count = (/3, 3/) )
       call check(iret)

write(1) cell(1:3,1:3)
write(1) mesh(1:3), nspin

allocate(gridfunc(1:mesh(1)))
!
       do ispin=1,nspin
          do iz=1,mesh(3)
             do iy=1,mesh(2)
                iret = nf90_get_var(ncid, gridfunc_id, gridfunc(1:mesh(1)), start = (/1, iy, iz, ispin /), &
                                                                            count = (/mesh(1), 1, 1, 1/) )
                call check(iret)
                write(1) (gridfunc(ix),ix=1,mesh(1))
             enddo
          enddo
       enddo

       close(1)
       call check( nf90_close(ncid) )

CONTAINS

subroutine check(code)
use netcdf, only: nf90_noerr, nf90_strerror
integer, intent(in) :: code
if (code /= nf90_noerr) then
  print *, "netCDF error: " // NF90_STRERROR(code)
  STOP
endif
end subroutine check

#endif /* CDF */
end program cdf2grid
