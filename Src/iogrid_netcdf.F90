! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module iogrid_netcdf

use parallel, only: Node, Nodes, ProcessorY
use alloc, only: re_alloc, de_alloc
use sys, only: die
#ifdef MPI
use mpi_siesta
#endif

implicit none

!     Private type to hold mesh distribution data

      TYPE meshDisType
      PRIVATE
      integer          :: nMesh(3)   ! Number of total mesh div. in each axis
      integer, pointer :: box(:,:,:) ! Mesh box bounds of each node:
                                     ! box(1,iAxis,iNode)=lower bounds
                                     ! box(2,iAxis,iNode)=upper bounds
      END TYPE meshDisType

      type(meshDisType),  target, save :: distr
      logical, save                    :: initialized = .false.

      public :: set_box_limits, write_grid_netcdf, read_grid_netcdf
      private

CONTAINS

subroutine set_box_limits(mesh,nsm)
!
! Set up easier-to-handle per-processor sub-mesh pointers
! Inspired by a routine written by Rogeli Grima (BSC)
!
  integer, intent(in)  :: mesh(3)     ! Mesh divisions (fine points) 
  integer, intent(in)  :: nsm         ! Number of fine points per big point

 integer :: nm(3)
 integer, allocatable :: npt_node(:)

 integer :: ProcessorZ, PP
 integer :: blocY, blocZ, nremY, nremZ
 integer :: dimY, dimZ, dimX, iniY, iniZ
 integer :: iNode, npt_total, npt_mesh
 integer :: PY, PZ, j


! Allocate memory for the current distribution
 if (.not. initialized) then
    nullify( distr%box )
    call re_alloc( distr%box, 1, 2, 1, 3, 1, Nodes, name='distr%box',  &
                             routine='set_box_limits' )
    initialized = .true.
 endif

 allocate(npt_node(1:Nodes))

 nm(1:3) = mesh(1:3) / nsm
 npt_mesh = mesh(1)*mesh(2)*mesh(3)
 distr%nMesh(1:3) = mesh(1:3)

     ProcessorZ = Nodes/ProcessorY

     blocY = (nm(2)/ProcessorY) 
     blocZ = (nm(3)/ProcessorZ) 
     nremY = nm(2) - blocY*ProcessorY
     nremZ = nm(3) - blocZ*ProcessorZ

     dimX = nm(1)*nsm

     npt_total = 0

     PP   = 1
     iniY = 1
     do PY = 1, ProcessorY

        dimY = blocY
        if (PY.LE.nremY) dimY = dimY + 1  ! Add extra points starting from the first nodes
        dimY = dimY * nsm                 ! For fine points

        iniZ = 1
        do PZ = 1, ProcessorZ
           dimZ = blocZ
           if (PZ.LE.nremZ) dimZ = dimZ + 1
           dimZ = dimZ*nsm                 ! For fine points
           
           distr%box(1,1,PP) = 1
           distr%box(2,1,PP) = dimX
           distr%box(1,2,PP) = iniY
           distr%box(2,2,PP) = iniY + dimY - 1
           distr%box(1,3,PP) = iniZ
           distr%box(2,3,PP) = iniZ + dimZ - 1

           npt_node(PP) = dimX * dimY * dimZ
           npt_total = npt_total + npt_node(PP)

           iniZ = iniZ + dimZ
           PP   = PP + 1
           
        enddo
        iniY = iniY + dimY
     enddo

     if (npt_total /= npt_mesh) then
        if (Node == 0) then
           write(6,*) "Nominal npt: ", npt_mesh, " /= assigned npt:", npt_total
        endif
        call die()
     endif

! JMS: commented out. 2009/02/06
!     if (Node == 0) then
!        do iNode= 1, Nodes
!           write(6,"(a,i4,a,i12,3x,3(i4,a1,i4))") "iogrid_netcdf: -- Node ", iNode, " :", npt_node(iNode), &
!                           (distr%box(1,j,iNode), ":", distr%box(2,j,iNode), j=1,3)
!        enddo
!     endif

     deallocate(npt_node)

end subroutine set_box_limits

subroutine write_grid_netcdf(cell,mesh,nspin,npt_l,gridfunc,name)
use precision, only: dp, grid_p
#ifdef CDF
use netcdf
#endif

      real(dp), intent(in)         ::     cell(3,3)    ! Lattice vectors
      integer, intent(in)          ::     mesh(3)      ! Number of mesh divisions of each lattice vector
      integer, intent(in)          ::     nspin
      integer, intent(in)          ::     npt_l         ! Total number of local points in rho, per spin
      real(grid_p), intent(in)     ::     gridfunc(npt_l,nspin)   ! Grid function
      character(len=*), intent(in) ::     name         ! Identifier

#ifdef CDF
integer  :: ncid , iret, ispin
integer  :: xyz_id, step_id, abc_id, spin_id
integer  :: cell_id, n1_id, n2_id, n3_id, gridfunc_id

#ifdef MPI
      integer, dimension(:), allocatable  :: npt_node 
      real(grid_p), dimension(:), pointer :: grid_buf  => null()
      integer  :: MPIerror, stat(MPI_STATUS_SIZE), count, BNode
      integer  :: iNode, max_npt
      integer  :: nel(3), lb(3)
#endif

    if (Node == 0) then
       iret = nf90_create(trim(name)//".grid.nc",NF90_SHARE+NF90_WRITE,ncid)
       call check(iret)
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
               "Grid function -- " // trim(name))
       call check(iret)

       iret = nf90_enddef(ncid)
       call check(iret)

!
       iret = nf90_put_var(ncid, cell_id, cell, start = (/1, 1 /), &
                        count = (/3, 3/) )
       call check(iret)
    endif

#ifdef MPI
    ! Get number of grid points in each of the processors and allocate buffer
    if (Node == 0) then
       allocate (npt_node(0:Nodes-1))
       do iNode = 0, Nodes-1
          nel(1:3) = distr%box(2,1:3,iNode+1) - distr%box(1,1:3,iNode+1) + 1
          npt_node(iNode) = product(nel(1:3))
       enddo
       max_npt = maxval(npt_node(0:Nodes-1))
       nullify(grid_buf)
       call re_alloc(grid_buf,1,max_npt,name="grid_buf",routine="write_grid_netcdf")
    endif

   do ispin = 1, nspin               ! Outer loop to simplify the logic, as we cannot
                                     ! send non-contiguous arrays
      do BNode = 0, Nodes - 1
         if (Node == 0) then
            if (BNode == 0) then
               grid_buf(1:npt_node(0)) = gridfunc(1:npt_l,ispin)     ! Could do with pointer
            else
               call MPI_Recv(grid_buf,npt_node(BNode),MPI_grid_real,BNode,ispin,  &
                                                        MPI_Comm_World,stat,mpierror)
               !!call MPI_GET_COUNT(stat, MPI_CHARACTER, count, mpierror)
               !! Note that we get bytes, and we need elements... need some king of "sizeof"
               !! for this test to work correctly
               !!if (count /= npt_node(BNode)) then
               !!   print *, 'Task ', Node ,': Received', count, 'char(s) from task',        &
               !!         stat(MPI_SOURCE), 'with tag',stat(MPI_TAG)
               !!   print *, "** But needed ", npt_node(BNode)
               !!   call die()
               !!endif
            endif
            !
            ! Fill in the information at the proper place
            lb(1:3) = distr%box(1,1:3,BNode+1) 
            nel(1:3) = distr%box(2,1:3,BNode+1) - distr%box(1,1:3,BNode+1) + 1
            call check( nf90_put_var(ncid,gridfunc_id,grid_buf(1:npt_node(BNode)),           &
                                                    start=(/ lb(1), lb(2), lb(3), ispin/),  &
                                                    count=(/nel(1), nel(2), nel(3), 1 /) ) )

         else      ! Non-master nodes

            if (Node == BNode) then
               call MPI_Send(gridfunc(1:npt_l,ispin),npt_l,MPI_grid_real,0,ispin,MPI_Comm_World,mpierror)
            endif

         endif

      enddo           ! Bnode
   enddo        ! ispin

   if (Node == 0) then
      call de_alloc(grid_buf,name="grid_buf", routine="write_grid_netcdf")
      deallocate(npt_node)
   endif

#else

   do ispin = 1, nspin
      iret = nf90_put_var(ncid, gridfunc_id, gridfunc, start = (/1, 1, 1, ispin /), &
           count = (/mesh(1), mesh(2), mesh(3), 1/) )
      call check(iret)
   enddo

#endif

       if (Node == 0) then
          iret = nf90_close(ncid)
          call check(iret)
       endif

#endif  /* CDF */
end subroutine write_grid_netcdf

!----------------------------------------
subroutine read_grid_netcdf(mesh,nspin,npt_l,gridfunc,name)
use precision, only: dp, grid_p

#ifdef CDF
use netcdf
#endif

implicit none

      integer, intent(in)          ::     mesh(3)      ! Number of mesh divisions of each lattice vector
      integer, intent(in)          ::     nspin
      integer, intent(in)          ::     npt_l         ! Total number of local points in rho, per spin
      real(grid_p), intent(out)    ::     gridfunc(npt_l,nspin)   ! Grid function
      character(len=*), intent(in) ::     name         ! Identifier

#ifdef CDF
integer  :: ncid , iret, ispin
integer  :: spin_id
integer  :: n1_id, n2_id, n3_id, gridfunc_id

integer  :: nspin_file, mesh_file(3)
character(len=40)  :: filename

#ifdef MPI
      integer, dimension(:), allocatable  :: npt_node 
      real(grid_p), dimension(:), pointer :: grid_buf  => null()
      integer  :: MPIerror, stat(MPI_STATUS_SIZE), count, BNode
      integer  :: iNode, max_npt
      integer  :: nel(3), lb(3)
#endif

    if (Node == 0) then
       filename = trim(name)//".IN.grid.nc"
       write(6,"(a)") "Reading " // trim(name) // " from file: " // trim(filename)
       iret = nf90_open(filename,NF90_NOWRITE,ncid)
       call check(iret)

       call check( nf90_inq_dimid(ncid,'spin',spin_id) )
       call check( nf90_inquire_dimension(ncid, dimid=spin_id, len=nspin_file) )
       if (nspin_file /= nspin) call die("Nspin mismatch")

       call check( nf90_inq_dimid(ncid,'n1',n1_id) )
       call check( nf90_inquire_dimension(ncid, dimid=n1_id, len=mesh_file(1)) )
       call check( nf90_inq_dimid(ncid,'n2',n2_id) )
       call check( nf90_inquire_dimension(ncid, dimid=n2_id, len=mesh_file(2)) )
       call check( nf90_inq_dimid(ncid,'n3',n3_id) )
       call check( nf90_inquire_dimension(ncid, dimid=n3_id, len=mesh_file(3)) )

       if (any(mesh(1:3) /= mesh_file(1:3)) ) then
          print *, "Mesh in file: ", mesh_file
          print *, "Mesh in program: ", mesh
          call die("Mesh mismatch")
       endif

       call check( nf90_inq_varid(ncid, "gridfunc", gridfunc_id) )

    endif

#ifdef MPI
    ! Get number of grid points in each of the processors and allocate buffer
    if (Node == 0) then
       allocate (npt_node(0:Nodes-1))
       do iNode = 0, Nodes-1
          nel(1:3) = distr%box(2,1:3,iNode+1) - distr%box(1,1:3,iNode+1) + 1
          npt_node(iNode) = product(nel(1:3))
       enddo
       max_npt = maxval(npt_node(0:Nodes-1))
       nullify(grid_buf)
       call re_alloc(grid_buf,1,max_npt,name="grid_buf",routine="write_grid_netcdf")
    endif

   do ispin = 1, nspin               ! Outer loop to simplify the logic, as we cannot
                                     ! send non-contiguous arrays
      do BNode = 0, Nodes - 1
         if (Node == 0) then
            lb(1:3) = distr%box(1,1:3,BNode+1) 
            nel(1:3) = distr%box(2,1:3,BNode+1) - distr%box(1,1:3,BNode+1) + 1
            call check( nf90_get_var(ncid,gridfunc_id,grid_buf(1:npt_node(BNode)),           &
                                                    start=(/ lb(1), lb(2), lb(3), ispin/),  &
                                                    count=(/nel(1), nel(2), nel(3), 1 /) ) )
            if (BNode == 0) then
               gridfunc(1:npt_l,ispin) = grid_buf(1:npt_node(0))  ! Could do with pointer
            else
               call MPI_Send(grid_buf,npt_node(BNode),MPI_grid_real,BNode,ispin,  &
                                                        MPI_Comm_World,mpierror)
            endif

         else      ! Non-master nodes

            if (Node == BNode) then
               call MPI_Recv(gridfunc(1:npt_l,ispin),npt_l,MPI_grid_real,0,ispin,MPI_Comm_World,stat,mpierror)
            endif

         endif

      enddo           ! Bnode
   enddo        ! ispin

   if (Node == 0) then
      call de_alloc(grid_buf,name="grid_buf", routine="read_grid_netcdf")
      deallocate(npt_node)
   endif

#else

   do ispin = 1, nspin
      iret = nf90_get_var(ncid, gridfunc_id, gridfunc, start = (/1, 1, 1, ispin /), &
           count = (/mesh(1), mesh(2), mesh(3), 1/) )
      call check(iret)
   enddo

#endif

       if (Node == 0) then
          iret = nf90_close(ncid)
          call check(iret)
       endif

#endif  /* CDF */
end subroutine read_grid_netcdf

#ifdef CDF
subroutine check(code)
use netcdf, only: nf90_noerr, nf90_strerror
integer, intent(in) :: code
if (code /= nf90_noerr) call die("netCDF error: " // NF90_STRERROR(code))
end subroutine check
#endif

end module iogrid_netcdf


