! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
! 
! This file is part of the SIESTA package.
!
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
program cdf_laplacian
#ifndef CDF
  print *, "No netCDF support at compile time"
#else

!
! Computes the (-)Laplacian of a grid function in netCDF format
!
!        The input file is called "GridFunc.nc"
!
!  It computes the laplacian of the original function by
!  multiplying the fourier coefficients by -G^2 and performing
!  the inverse fourier transform.  
!
!  If used on a density file, the first coefficient times the cell 
!  volume should give the number of electrons in the unit cell.
!  (As a convenience, this is written to standard error).
!
!  Only the first spin component is processed in this prototype version.
!
!  Alberto Garcia, November 2015
!
use netcdf
use m_getopts

implicit none

integer, parameter  :: dp = selected_real_kind(14,100)
integer, parameter  :: sp = selected_real_kind(6,20)

integer  :: ncid 
integer  :: xyz_id, step_id, abc_id, spin_id
integer  :: cell_id, n1_id, n2_id, n3_id, gridfunc_id

integer   ::    nspin  ! Number of spins 
integer   ::    n(3), nc(3), ig, jg, kg, i, j, k, ip
integer   ::    n1, n2, n3, mesh(3)
integer   ::    i1, i2, i3
integer   ::   ispin, iostat, ix, iy, iz, iret

real(dp)  ::    cell(3,3), cell_volume, g(3), g2, b(3,3)
real(sp), dimension(:), allocatable :: gridfunc
complex(dp), dimension(:), allocatable :: aa, lap_aa
character(len=256) :: densfile, lapfile

external fft3d

      character(len=200) :: opt_arg
      character(len=10)  :: opt_name 
      integer :: nargs, n_opts

            logical :: supercell = .false.
      logical :: debug = .false.
      integer :: filter_flag = 0
      logical :: densfile_given = .false.
      logical :: lapfile_given = .false.

      !     
!     Process options
!     
      n_opts = 0
      do
      call getopts('hd:i:o:f:',opt_name,opt_arg,n_opts,iostat)
      if (iostat /= 0) exit
      select case(opt_name)
      case ('d')
         debug = .true.
      case ('f')
         read(opt_arg,*) filter_flag
      case ('i')
         read(opt_arg,*) densfile
         densfile_given = .true.
      case ('o')
         read(opt_arg,*) lapfile
         lapfile_given = .true.

      case ('h')
         write(0,*) "Usage: cdf_laplacian [-dh] [-f filter] -i DENSITY_FILE -o LAPFILE"
         write(0,*) " -f FLAG : filter according to flag (+1/-1)"
         write(0,*) " -d : debug"
         STOP
      case ('?',':')
         write(0,*) "Invalid option: ", opt_arg(1:1)
         write(0,*) "Usage: g2c_ng [-dht] -s STRUCT_FILE -g LAPFILE"
         write(0,*) "Use -h option for manual"
         STOP
      end select
      enddo

      if (.not. densfile_given) then
         write(0,*) "You must give the density file ( -i option)"
         write(0,*) "Usage: cdf_laplacian [-dh] [-f filter] -i DENSITY_FILE -o LAPFILE"
         STOP
      endif
      if (.not. lapfile_given) then
         write(0,*) "You must give the grid file ( -g option)"
         write(0,*) "Usage: g2c_ng [-dht] -s DENSITY_FILE -g LAPFILE"
         STOP
      endif

!-----------------------------------------------------

call check( nf90_open(trim(densfile),NF90_NOWRITE,ncid))
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
       cell_volume = volcel(cell)
!  Find reciprocal lattice vectors
   call reclat(CELL, B, 1 )

   allocate(gridfunc(product(n(1:3))))
   allocate(aa(product(n(1:3))))
   
   ! Laplacian (with minus sign)
   allocate(lap_aa(product(n(1:3))))
   lap_aa(:) = 0.0_dp
   
   ispin = 1
   iret = nf90_get_var(ncid, gridfunc_id, gridfunc, start = (/1, 1, 1, ispin /), &
           count = (/n(1), n(2), n(3), 1/) )
      call check(iret)

   n1 = n(1)
   n2 = n(2)
   n3 = n(3)

   nc(1:3) = (n(1:3) - 1) /2
   
   aa = cmplx(gridfunc,kind=dp)
   call fft3d(aa,(/n1, n2, n3/),+1)
   do ip = 1, size(aa)
      i = modulo(ip-1,n1) + 1
      j = modulo( (ip - i)/n1, n2) + 1
      k = (ip - i - (j-1)*n1)/(n1*n2)  + 1

      ig = i - n1 - 1
      if (ig < -nc(1)) ig = ig + n1

      jg = j - n2 - 1
      if (jg < -nc(2)) jg = jg + n2

      kg = k - n3 - 1
      if (kg < -nc(3)) kg = kg + n3

      i1 = ig; i2 = jg; i3 = kg
      
      G(1)= B(1,1) * I1 + B(1,2) * I2 + B(1,3) * I3
      G(2)= B(2,1) * I1 + B(2,2) * I2 + B(2,3) * I3
      G(3)= B(3,1) * I1 + B(3,2) * I2 + B(3,3) * I3
      G2  = G(1)**2 + G(2)**2 + G(3)**2

      lap_aa(ip) = g2*aa(ip)
   enddo

   call check( nf90_close(ncid) )

   ! To convert to number of electrons
   write(0,"(a,f20.6)") "Cell volume * aa(1) (number of electrons): ", &
           real(aa(1),kind=dp)*cell_volume

   ! Inverse Fourier transform to get the laplacian
   call fft3d(lap_aa,(/n1, n2, n3/),-1)
   ! Normalization !!!
        
   gridfunc = real(lap_aa,kind=sp)

   ! Filter values
   ! If flag>0, we will keep the positive values, and so on
   where (filter_flag*gridfunc < 0.0_dp)
      gridfunc = 0.0_dp
   end where
   
   mesh(1) = n1; mesh(2) = n2; mesh(3) = n3

   call check( nf90_create(trim(lapfile),NF90_CLOBBER,ncid))
       iret = nf90_def_dim(ncid,'xyz',3,xyz_id)
       call check(iret)
       iret = nf90_def_dim(ncid,'abc',3,abc_id)
       call check(iret)
       iret = nf90_def_dim(ncid,'spin',nspin,spin_id)
       call check(iret)
       iret = nf90_def_dim(ncid,'n1',mesh(1),n1_id)
       call check(iret)
       iret = nf90_def_dim(ncid,'n2',mesh(2),n2_id)
       call check(iret)
       iret = nf90_def_dim(ncid,'n3',mesh(3),n3_id)
       call check(iret)

       iret = nf90_def_var(ncid,'cell',nf90_float,(/xyz_id,abc_id/),cell_id)
       call check(iret)
       iret = nf90_put_att(ncid,cell_id,'Description', &
               "Cell vectors in Bohr: xyz, abc")
       call check(iret)

       iret = nf90_def_var(ncid,'gridfunc',nf90_float,(/n1_id,n2_id,n3_id,spin_id/),gridfunc_id)
       call check(iret)
       iret = nf90_put_att(ncid,gridfunc_id,'Description', &
            "(-) Laplacian of Charge Density ")
       if (filter_flag /= 0) then
          iret = nf90_put_att(ncid,gridfunc_id,'Filtered with flag: ', filter_flag)
       endif
       call check(iret)

       iret = nf90_enddef(ncid)
       call check(iret)

!
       iret = nf90_put_var(ncid, cell_id, cell, start = (/1, 1 /), count = (/3, 3/) )
       call check(iret)

       iret = nf90_put_var(ncid, gridfunc_id, gridfunc, start = (/1, 1, 1, ispin /), &
       count = (/n(1), n(2), n(3), 1/) )
       call check(iret)

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

FUNCTION VOLCEL( C )
real(dp) :: volcel

!  CALCULATES THE VOLUME OF THE UNIT CELL

      real(dp) ::  C(3,3)
      VOLCEL = ( C(2,1)*C(3,2) - C(3,1)*C(2,2) ) * C(1,3) + &
               ( C(3,1)*C(1,2) - C(1,1)*C(3,2) ) * C(2,3) + &
               ( C(1,1)*C(2,2) - C(2,1)*C(1,2) ) * C(3,3)
      VOLCEL = ABS( VOLCEL )
END function volcel

#endif /* CDF */
end program cdf_laplacian
