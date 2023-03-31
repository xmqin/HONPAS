! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module iowfs_netcdf

#ifdef CDF
use netcdf
use sys, only: die
use precision, only: dp
use parallel, only: Node, Nodes
use alloc, only: re_alloc, de_alloc

#ifdef MPI
use mpi_siesta
use parallelsubs, only : GLobalToLocalOrb, WhichNodeOrb
#endif

implicit none

! These module variables should be put into a derived type, and maybe
! not all of them are really necessary
!
integer ncid, norbs_id, nspin_id, nk_id, ncomplex_id, nvectors_id
integer wf_id    , eigval_id

public :: setup_wfs_netcdf_file, write_wfs_netcdf, get_wfs_netcdf
public :: open_wfs_netcdf_file, close_wfs_netcdf_file
public :: get_wfs_block_netcdf

private

CONTAINS

subroutine setup_wfs_netcdf_file( norbs, nk, nspin, nvectors)
                             
!  To be called only by the master node

      integer, intent(in)   ::    norbs  ! Number of atomic orbitals
      integer, intent(in)   ::    nspin  ! Number of spins 
      integer, intent(in)   ::    nk     ! Number of kpoints
      integer, intent(in)   ::    nvectors   ! Number of eigenvectors stored per spin and k-point

      call check( nf90_create('WFS.nc',NF90_CLOBBER,ncid))
!
!     Dimensions
!
      call check( nf90_def_dim(ncid,'norbs',norbs,norbs_id))  !"Number of basis orbitals"
      call check( nf90_def_dim(ncid,'nspin',nspin,nspin_id))   !"Number of spin components"
      call check( nf90_def_dim(ncid,'nk',nk,nk_id))     !"Number of k-points"
      call check( nf90_def_dim(ncid,'nvectors',nvectors,nvectors_id))     !"Number of k-points"
      call check( nf90_def_dim(ncid,'ncomplex',2,ncomplex_id))     !"Real/imaginary"

!
!     Variables
!
      call check( nf90_def_var(ncid,'eigval',nf90_float,(/nvectors_id,nk_id,nspin_id/),eigval_id))
      call check( nf90_put_att(ncid,eigval_id,'Description',"Eigenvalues in Ryd"))
      call check( nf90_def_var(ncid,'wf',nf90_float,    &
                            (/ncomplex_id,norbs_id,nvectors_id,nk_id,nspin_id/),wf_id))
      call check( nf90_put_att(ncid,wf_id,'Description',"Wavefunctions"))

      call check( nf90_enddef(ncid))

end subroutine setup_wfs_netcdf_file

subroutine write_wfs_netcdf( norbs, no_l, ik, ispin, psi, nvectors, eigvals )

use precision, only : dp

integer, intent(in)   ::    norbs  ! Total number of orbitals in unit cell (components of wf)
integer, intent(in)   ::    no_l   ! Number of eigenvectors in this node
integer, intent(in)   ::    ik     ! k-point index
integer, intent(in)   ::    ispin  ! spin index
integer, intent(in)   ::    nvectors ! total number of eigenvectors to be stored for this k-point and spin

real(dp), intent(in)  :: psi(2,norbs,no_l)    ! wave-functions stored in this node
real(dp), intent(in)  :: eigvals(nvectors)    ! wave-functions stored in this node

#ifdef MPI
      real(dp), dimension(:,:), pointer  ::  psi_buf  => null()
      integer  :: MPIerror, stat(MPI_STATUS_SIZE), count, BNode, tag, ivec_local
#endif

      integer  :: ivec

#ifdef MPI

if (Node == 0) then
   nullify( psi_buf )
   call re_alloc( psi_buf, 1, 2, 1, norbs, name='psi_buf', routine='iowfs_netcdf' )
   call check( nf90_put_var(ncid,eigval_id,eigvals(1:nvectors),         &
                            start=(/1, ik, ispin/),  count=(/nvectors, 1, 1 /) ) )
endif

    do ivec = 1, nvectors
       call WhichNodeOrb(ivec,Nodes,Bnode)
       call GlobalToLocalOrb(ivec,BNode,Nodes,ivec_local)
       tag = ivec

       if (Node == 0) then
          if (BNode == 0) then
             psi_buf(1:2,1:norbs) = psi(1:2,1:norbs,ivec_local)     ! Could do with pointer
          else
             call MPI_Recv(psi_buf(1,1),2*norbs,MPI_double_precision,BNode,tag,  &
                                                        MPI_Comm_World,stat,mpierror)
          endif
          call check( nf90_put_var(ncid,wf_id,psi_buf(1:2,1:norbs),         &
                                       start=(/1, 1, ivec, ik, ispin/),        &
                                       count=(/2, norbs, 1, 1, 1 /) ) )
       else      ! Non-master nodes

          if (Node == BNode) then
             call MPI_Send(psi(1,1,ivec_local),2*norbs, &
                           MPI_double_precision,0,tag,MPI_Comm_World,mpierror)
          endif

       endif

    enddo  ! loop over vectors

    call de_alloc(psi_buf,name="psi_buf", routine="iowfs_netcdf")

#else
    call check( nf90_put_var(ncid,eigval_id,eigvals(1:nvectors),         &
                             start=(/1, ik, ispin/),  count=(/nvectors, 1, 1 /) ) )
    do ivec = 1, nvectors
       call check( nf90_put_var(ncid, wf_id, psi(1:2,1:norbs,ivec),  & 
                              start = (/1, 1, ivec, ik, ispin /), &
                              count = (/2, norbs, 1, 1, 1 /) ))
    enddo
#endif

if (Node == 0) then
   call check( nf90_sync(ncid))
endif

end subroutine write_wfs_netcdf

subroutine close_wfs_netcdf_file()
if (Node == 0) then
   call check( nf90_close(ncid))
endif
end subroutine close_wfs_netcdf_file

subroutine open_wfs_netcdf_file()
if (Node == 0) then
   call check( nf90_open("WFS.nc",NF90_NOWRITE,ncid))
   call check( nf90_inq_varid(ncid,"wf",wf_id))
endif
end subroutine open_wfs_netcdf_file

subroutine get_wfs_netcdf(ivec,ik,ispin,norbs,psi_aux)

use precision, only : dp

integer, intent(in)   ::    ivec   ! Eigenvector index
integer, intent(in)   ::    ik     ! k-point index
integer, intent(in)   ::    ispin  ! spin index
integer, intent(in)   ::    norbs  ! Total number of orbitals in unit cell (components of wf)

real(dp), intent(out)  :: psi_aux(2,norbs)    ! wave-function

#ifdef MPI
      integer  :: MPIerror
#endif

if (Node == 0) then
       call check( nf90_get_var(ncid, wf_id, psi_aux(1:2,1:norbs),  & 
                              start = (/1, 1, ivec, ik, ispin /), &
                              count = (/2, norbs, 1, 1, 1 /) ))
endif
#ifdef MPI
call MPI_Bcast(psi_aux(1,1),2*norbs,MPI_Double_Precision,0,MPI_Comm_World, MPIerror)
#endif
end subroutine get_wfs_netcdf

subroutine get_wfs_block_netcdf(ivec_start,nvecs,ik,ispin,norbs,psi_block)

use precision, only : dp

integer, intent(in)   ::    ivec_start ! Index of first eigenvector to get
integer, intent(in)   ::    nvecs  ! Number of eigenvectors to get
integer, intent(in)   ::    ik     ! k-point index
integer, intent(in)   ::    ispin  ! spin index
integer, intent(in)   ::    norbs  ! Total number of orbitals in unit cell (components of wf)

real(dp), intent(out)  :: psi_block(2,norbs,nvecs)    ! wave-function block

#ifdef MPI
      integer  :: MPIerror
#endif

if (Node == 0) then
       call check( nf90_get_var(ncid, wf_id, psi_block(1:2,1:norbs,1:nvecs),  & 
                              start = (/1, 1, ivec_start, ik, ispin /), &
                              count = (/2, norbs, nvecs, 1, 1 /) ))
endif

#ifdef MPI
call MPI_Bcast(psi_block(1,1,1),2*norbs*nvecs,MPI_Double_Precision,0,MPI_Comm_World, MPIerror)
#endif

end subroutine get_wfs_block_netcdf

subroutine check(code)
use netcdf, only: nf90_noerr, nf90_strerror
integer, intent(in) :: code
if (code /= nf90_noerr) call die("netCDF error: " // NF90_STRERROR(code))
end subroutine check


#endif
end module iowfs_netcdf
