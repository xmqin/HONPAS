! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!
!> Utilities to dea with real-space grid functions
module m_gridfunc

  implicit          none

  integer, parameter, private :: dp = selected_real_kind(10,100)
  integer, parameter, private :: sp = selected_real_kind(5,10)
  integer, parameter          :: grid_p = sp  ! NOTE this

  type, public :: gridfunc_t
     !> Lattice vectors (by columns, in bohr)
     real(dp)              ::  cell(3,3)   
     !> Possible shift for the origin of the box, useful in some
     !> situations (sub-boxes, for example) (in bohr)
     real(dp)              ::  origin(3) = [0.0_dp, 0.0_dp, 0.0_dp]
     !> The 'grid' format can in principle hold functions that
     !> are not periodic (for example, values in a sub-box)
     logical               ::  is_periodic(3) = [ .true., .true., .true.]
     !> Number of points along each cell direction
     integer               ::  n(3) = [ 0,0,0 ]
     !> Extra dimension, typically the spin index
     integer               ::  nspin
     !> Array holding the values of the grid function
     real(grid_p), allocatable ::  val(:,:,:,:)
  end type gridfunc_t

      public :: write_gridfunc, read_gridfunc, clean_gridfunc
#ifdef CDF
      public :: read_gridfunc_netcdf, write_gridfunc_netcdf
#endif
      public :: get_planar_average, grid_p, monoclinic_z
      public :: get_values_at_point
      public :: monoclinic
    
      private

 CONTAINS

   subroutine clean_gridfunc(gf)
     type(gridfunc_t), intent(inout) :: gf

     if (allocated(gf%val)) then
        deallocate(gf%val)
     endif
     gf%n(:) = 0
     gf%cell(:,:) = 0.0_dp
     gf%origin(:) = 0.0_dp
     gf%is_periodic(:) = .true.

   end subroutine clean_gridfunc

   subroutine write_gridfunc(gf,fname)
     type(gridfunc_t), intent(in) :: gf
     character(len=*), intent(in) :: fname

     integer :: n(3)
     integer :: isp, ix, iy, iz
     integer :: iu

     n = gf%n

     call get_lun(iu)
     open(unit=iu,file=fname,form="unformatted",status="unknown", &
          position="rewind",action="write")
     write(iu) gf%cell, gf%origin, gf%is_periodic
     write(iu) n, gf%nspin
     do isp=1,gf%nspin
        do iz=1,n(3)
           do iy=1,n(2)
              write(iu) (gf%val(ix,iy,iz,isp),ix=1,n(1))
           enddo
        enddo
     enddo
     close( iu )

   end subroutine write_gridfunc

   !> Evaluates a 'grid function' at an arbitrary point in the box.
   !> The function is computed using linear interpolation.
   subroutine get_values_at_point(gf, xfrac, vals)
     !> The object representing the grid function
     type(gridfunc_t), intent(in) :: gf
     !> The fractional coordinates of the point
     real(dp), intent(in) :: xfrac(3)
     !> Array holding the values (with the 'spin' dimension)
     real(grid_p), allocatable :: vals(:)

     integer  ::  n(3), lo(3), hi(3)
     real(dp) ::  r(3), x(3), y(3)
     integer  ::  i, j, k
     real(dp) ::  nk

     n(:) = gf%n(:)

!           Find the right 3D "grid cube" and the reduced coordinates
!           of the point in it. The double mod assures that negative
!           numbers are well treated (the idea is to bring the coordinates
!           to the [0,n(k)) interval)
! 
     do k = 1, 3
        nk = real(n(k),kind=dp)
        r(k) =  modulo(n(k)*xfrac(k),nk)
        lo(k) = int(r(k))  
        hi(k) = mod ( lo(k)+1, n(k) ) 
        x(k) = r(k) - lo(k)
        y(k) = 1 - x(k)
     enddo

     ! Switch to 1-based array convention
     lo(:) = lo(:) + 1
     hi(:) = hi(:) + 1
     
!      compute charge density by linear interpolation

     vals(:) = gf%val(lo(1),lo(2),lo(3),:) * y(1) * y(2) * y(3) + &
          gf%val(lo(1),lo(2),hi(3),:) * y(1) * y(2) * x(3) + &
          gf%val(lo(1),hi(2),lo(3),:) * y(1) * x(2) * y(3) + &
          gf%val(lo(1),hi(2),hi(3),:) * y(1) * x(2) * x(3) + &
          gf%val(hi(1),lo(2),lo(3),:) * x(1) * y(2) * y(3) + &
          gf%val(hi(1),lo(2),hi(3),:) * x(1) * y(2) * x(3) + &
          gf%val(hi(1),hi(2),lo(3),:) * x(1) * x(2) * y(3) + &
          gf%val(hi(1),hi(2),hi(3),:) * x(1) * x(2) * x(3) 

   end subroutine get_values_at_point

   
   !> Computes a simple-minded 'planar average' of a grid
   !> function. This is a mathematical average. For a more
   !> physical average, see the 'macroave' utility
   subroutine get_planar_average(gf, dim, val)
     !> The object representing the grid function
     type(gridfunc_t), intent(in) :: gf
     !> The dimension that indexes the average values.
     !> These are obtained by summing over the two other dimensions.
     integer, intent(in) :: dim
     !> Array holding the averages (second dimension is the 'spin'))
     real(grid_p), allocatable :: val(:,:)

     integer :: isp, ix, iy, iz, n(3), spin_dim
     real(grid_p) :: sum

     n = gf%n
     spin_dim = size(gf%val,dim=4)

     allocate(val(n(dim),spin_dim))

     do isp=1,gf%nspin
        select case (dim)
        case (3)
           do iz=1,n(3)
              sum = 0.0_grid_p
              do iy=1,n(2)
                 do ix=1,n(1)
                    sum = sum + gf%val(ix,iy,iz,isp)
                 enddo
              enddo
              val(iz,isp) = sum/(n(1)*n(2))
           enddo
        case (2)
           do iy=1,n(2)
              sum = 0.0_grid_p
              do iz=1,n(3)
                 do ix=1,n(1)
                    sum = sum + gf%val(ix,iy,iz,isp)
                 enddo
              enddo
              val(iy,isp) = sum/(n(1)*n(3))
           enddo
        case (1)
           do ix=1,n(1)
              sum = 0.0_grid_p
              do iz=1,n(3)
                 do iy=1,n(2)
                    sum = sum + gf%val(ix,iy,iz,isp)
                 enddo
              enddo
              val(ix,isp) = sum/(n(2)*n(3))
           enddo
        end select
        
     enddo

   end subroutine get_planar_average

   subroutine read_gridfunc(fname,gf)
     type(gridfunc_t), intent(inout) :: gf
     character(len=*), intent(in) :: fname

     integer :: n(3) 
     integer :: isp, ix, iy, iz
     integer :: iu, iostat

     call clean_gridfunc(gf)

     call get_lun(iu)
     open(unit=iu,file=fname,form="unformatted",status="old", &
          position="rewind",action="read")

     read(iu,iostat=iostat) gf%cell, gf%origin, gf%is_periodic
     if (iostat /= 0) then
        backspace(iu)
        read(iu) gf%cell
        gf%origin(:) = 0.0_dp
        gf%is_periodic(:) = .true.
     endif

     read(iu) gf%n, gf%nspin
     n = gf%n
     allocate(gf%val(n(1),n(2),n(3),gf%nspin))

     do isp=1,gf%nspin
        do iz=1,n(3)
           do iy=1,n(2)
              read(iu) (gf%val(ix,iy,iz,isp),ix=1,n(1))
           enddo
        enddo
     enddo
     close( iu )

   end subroutine read_gridfunc

#ifdef CDF

subroutine write_gridfunc_netcdf(gf,filename,description)
use netcdf

implicit none

character(len=*), intent(in) :: filename
type(gridfunc_t), intent(in)    :: gf
character(len=*), intent(in), optional :: description

integer  :: ncid 
integer  :: xyz_id, step_id, abc_id, spin_id
integer  :: cell_id, n1_id, n2_id, n3_id, gridfunc_id
integer  :: origin_id, is_periodic_id

integer   ::    nspin  ! Number of spins 
integer   ::    n(3), int_from_bool(3), i
integer   ::   ispin, iostat, ix, iy, iz, iret

!-----------------------------------------------------

call check( nf90_create(filename,NF90_CLOBBER,ncid))
       iret = nf90_def_dim(ncid,'xyz',3,xyz_id)
       call check(iret)
       iret = nf90_def_dim(ncid,'abc',3,abc_id)
       call check(iret)
       iret = nf90_def_dim(ncid,'spin',gf%nspin,spin_id)
       call check(iret)

       n(:) = gf%n(:)

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

       iret = nf90_def_var(ncid,'origin',nf90_float,(/xyz_id/),origin_id)
       call check(iret)
       iret = nf90_put_att(ncid,origin_id,'Description', &
               "Origin of coordinates: xyz")
       call check(iret)

       iret = nf90_def_var(ncid,'is_periodic',nf90_int,(/abc_id/),is_periodic_id)
       call check(iret)
       iret = nf90_put_att(ncid,cell_id,'Description', &
               "Periodic along cell vectors (1:yes, 0:no): abc")
       call check(iret)

       iret = nf90_def_var(ncid,'gridfunc',nf90_float,(/n1_id,n2_id,n3_id,spin_id/),gridfunc_id)
       call check(iret)
       if (present(description)) then
          iret = nf90_put_att(ncid,gridfunc_id,'Description', &
               trim(description))
       else
          iret = nf90_put_att(ncid,gridfunc_id,'Description', &
               "Grid function -- ")
       endif
       call check(iret)

       iret = nf90_enddef(ncid)
       call check(iret)
!
       iret = nf90_put_var(ncid, cell_id, gf%cell, start = (/1, 1 /), count = (/3, 3/) )
       call check(iret)
       iret = nf90_put_var(ncid, origin_id, gf%origin, start = [1], count = [3] )
       call check(iret)
       int_from_bool(:) = 0
       do i = 1, 3
          if (gf%is_periodic(i)) int_from_bool(i) = 1
       enddo
       iret = nf90_put_var(ncid, is_periodic_id, int_from_bool, start = [1], count = [3] )
       call check(iret)

      iret = nf90_put_var(ncid, gridfunc_id, gf%val, start = (/1, 1, 1, 1 /), &
           count = (/n(1), n(2), n(3), nspin/) )
      call check(iret)
   
     call check( nf90_close(ncid) )

end subroutine write_gridfunc_netcdf
!--------

subroutine read_gridfunc_netcdf(filename,gf)
use netcdf

implicit none

character(len=*), intent(in) :: filename
type(gridfunc_t), intent(inout)    :: gf

integer  :: ncid 
integer  :: xyz_id, step_id, abc_id, spin_id
integer  :: cell_id, n1_id, n2_id, n3_id, gridfunc_id
integer  :: origin_id, is_periodic_id

integer   ::    nspin  ! Number of spins 
integer   ::    n(3), int_from_bool(3), i
integer   ::   ispin, iostat, ix, iy, iz, iret

logical   ::   has_origin, has_is_periodic
real(dp)  ::    cell(3,3)
!-----------------------------------------------------

call clean_gridfunc(gf)

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

       iret =  nf90_inq_varid(ncid, "origin", origin_id)
       has_origin = (iret == nf90_noerr)
       iret =  nf90_inq_varid(ncid, "is_periodic", is_periodic_id)
       has_is_periodic = (iret == nf90_noerr)
          
       call check( nf90_inq_varid(ncid, "gridfunc", gridfunc_id) )

       iret = nf90_get_var(ncid, cell_id, cell, start = (/1, 1 /), &
            count = (/3, 3/) )

       if (has_origin) then
          iret = nf90_get_var(ncid, origin_id, gf%origin, start = [1], count = [3] )
          call check(iret)
       else
          gf%origin(:) = 0.0_dp
       endif
       
       if (has_is_periodic) then
          iret = nf90_get_var(ncid, is_periodic_id, int_from_bool, start = [1], count = [3] )
          call check(iret)
          do i = 1, 3
             gf%is_periodic(i) = (int_from_bool(i) == 1)
          enddo
       else
          gf%is_periodic(:) = .true.
       endif
       
   gf%n(:) = n(:)
   gf%cell = cell

   allocate(gf%val(n(1),n(2),n(3),nspin))

   iret = nf90_get_var(ncid, gridfunc_id, gf%val, start = (/1, 1, 1, 1 /), &
        count = (/n(1), n(2), n(3), nspin/) )
   call check(iret)
   
   call check( nf90_close(ncid) )

end subroutine read_gridfunc_netcdf

subroutine check(code)
use netcdf, only: nf90_noerr, nf90_strerror
integer, intent(in) :: code
if (code /= nf90_noerr) then
  print *, "netCDF error: " // NF90_STRERROR(code)
  STOP
endif
end subroutine check

#endif /* CDF */

 subroutine get_lun(lun)
   integer, intent(out) :: lun

   logical :: busy
   do lun = 1, 99
      inquire(unit=lun,opened=busy)
      if (.not. busy) RETURN
   enddo
   call die("Cannot get free unit")
 end subroutine get_lun

 subroutine die(str)
   character(len=*), intent(in) :: str
   write(0,"(a)") str
   stop
 end subroutine die

 !> Physically, it should not matter whether the 'c' axis points along
 !> z, since rotations are irrelevant. However, some programs are
 !> hard-wired for the z-direction. Hence the name. See the
 !> 'monoclinic' function for the rotationally invariant
 !> version.
 
 function monoclinic_z(cell)
   real(dp), intent(in) :: cell(3,3)
   logical monoclinic_z

   real(dp), parameter :: tol = 1.0e-8_dp

   ! Check that the 'c' vector is along z, and 'a' and 'b' in the XY plane

   monoclinic_z =  (abs(CELL(3,1)) < tol   &
               .and. abs(CELL(3,2)) < tol  &
               .and. abs(CELL(1,3)) < tol  &
               .and. abs(CELL(2,3)) < tol )

   end function monoclinic_z

   !> Checks that the lattice vector of index 'dim' is orthogonal to
   !> the other two.

   function monoclinic(cell,dim)
     real(dp), intent(in) :: cell(3,3)
     integer, intent(in)  :: dim
     logical monoclinic

     real(dp), parameter :: tol = 1.0e-8_dp

   ! Check that the 'dim' vector is orthogonal to the other two
     monoclinic = .true.
     if (abs(dot_product(cell(:,dim),cell(:,next(dim,1)))) > tol) then
        monoclinic = .false.
     endif
     if (abs(dot_product(cell(:,dim),cell(:,next(dim,2)))) > tol) then
        monoclinic = .false.
     endif

   CONTAINS
     !> generates cyclic neighbors of dimension 'i' with shift 'm'
     !> For example:
     !>   x --> y:   next(1,1) = 2
     !>   y --> x:   next(2,2) = 1
     function next(i,m) result(ind)
       integer, intent(in) :: i, m
       integer :: ind
       
       ind = i + m
       if (ind > 3) then
          ind = ind - 3
       endif
     end function next
   end function monoclinic

 end module m_gridfunc

