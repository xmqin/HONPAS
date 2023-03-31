!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been improved or fully created by:
! Nick Papior Andersen, 2014, nickpapior@gmail.com
!
module m_mesh_node

! Module for retaining information about the mesh.
! It is used by the Hartree module and the bias module.
!
! Created and copyrighted by: Nick Papior Andersen, 2014
! The use of this program is allowed for not-for-profit research only.
! Copy or disemination of all or part of this package is not
! permitted without prior and explicit authorization by the author.

  use precision, only : dp

  implicit none

  public
  save

  ! The offset for the current node
  real(dp) :: offset_r(3) = 0._dp
  ! the voxel-vectors for each sub-mesh element
  real(dp) :: dL(3,3) = 0._dp
  ! the voxel length along each direction
  real(dp) :: dMesh(3) = 0._dp

  ! offsets for the local node
  integer :: meshl(3) = 0
  integer :: offset_i(3) = 0

contains

  subroutine init_mesh_node( ucell , meshG , meshLim , nsm)

    use intrinsic_missing, only : VNORM
    use parallel, only : Node, Nodes, ProcessorY, IONode
    use units, only: Ang
#ifdef MPI
    use mpi_siesta
#endif

    ! The unit cell
    real(dp), intent(in) :: ucell(3,3)
    ! Number of mesh divisions of each lattice vector
    integer, intent(in) :: meshG(3), meshLim(2,3)
    ! Number of fine points per big point (see iogrid_netcdf)
    integer, intent(in) :: nsm 

    ! Number of big division points
    integer :: nm(3)
    ! Processor specifics
    integer :: ProcessorZ, blocY, blocZ, nremY, nremZ, idx(3)
    ! dimension tracking of the divisions
    integer :: iniX, iniY, iniZ, dimX, dimY, dimZ
    ! Local node dimensionality of the grid
    integer :: ldimX, ldimY, ldimZ
    ! Loop stuff
    integer :: node_Y, node_Z, cur_Node
    integer :: i
#ifdef MPI
    real(dp) :: ll(3)
    integer :: ix, iy, iz
    integer :: MPIerror, status(MPI_STATUS_SIZE)
#endif
    
    ! We calculate the spacing in each direction
    ! Notice that we now have:
    ! dL(1,1) = dX for stepping in the x-direction
    ! dL(1,2) = dX for stepping in the y-direction
    ! dL(1,3) = dX for stepping in the z-direction
    do i = 1 , 3 
       ldimX = max(meshG(i),1)
       ! The dimension stepping in each direction.
       dL(:,i) = ucell(:,i) / ldimX
    end do
    ! The voxel box-size in Cartesian coordinates 
    ! is calculated by adding all three vectors
    dMesh = matmul(dL,(/1._dp,1._dp,1._dp/))

    ! For nodes == 1 we have no offset
    ! (also some of the arrays are not initialized, which
    !  could lead to errors)
    if ( Nodes == 1 ) then
       meshl = meshG ! same as meshLim(2,:) * nsm
       return
    end if
       
    ! Now we need to calculate the offset of the local node

    ! Calculate the number of big-points
    nm(1:3) = meshG(1:3) / nsm

    ! In order to check that we do it correctly...
    dimX  = nm(1) * nsm
    ldimX = dimX

    ! Calculate number of z-processor divisions
    ProcessorZ = Nodes / ProcessorY
    
    ! Calculate the block-division in the Y-direction
    blocY = nm(2) / ProcessorY
    ! Calculate the block-division in the Z-direction
    blocZ = nm(3) / ProcessorZ
    ! If there are any remaining mesh-points, they will be added one 
    ! to each of the processors, in sequence
    nremY = nm(2) - blocY*ProcessorY
    nremZ = nm(3) - blocZ*ProcessorZ

    ! Initialize the loop construct variables
    cur_Node = 0
    ! Notice that we start from zero as we would like to calculate the
    ! offset
    iniX = 0
    iniY = 0
    do node_Y = 1, ProcessorY

       ! Initialize the dimension size in the Y-direction
       dimY = blocY
       ! The if there were too many points to be perfectly divisable
       ! we will add them.
       if ( node_Y <= nremY ) dimY = dimY + 1
       ! Extend into the mesh size
       dimY = dimY * nsm

       ! Initialize the Z-division (notice we index from zero
       ! as we would like to calculate the offset)
       iniZ = 0
       do node_Z = 1, ProcessorZ
          dimZ = blocZ
          if ( node_Z <= nremZ ) dimZ = dimZ + 1
          dimZ = dimZ * nsm ! For fine points
          
          if ( Node == cur_Node ) then
             ! Calculate the offset in the [YZ]-direction for the processor
             ! We know that iniX == 0, so no need, but we have it for
             ! consistency
             offset_r(:) = iniX * dL(:,1) + iniY * dL(:,2) + iniZ * dL(:,3)
             ldimY = dimY
             ldimZ = dimZ
          end if
          
          iniZ = iniZ + dimZ
          cur_Node = cur_Node + 1
       end do
       iniY = iniY + dimY
    end do

    ! Find quantities in mesh coordinates
    meshl(1) = (meshLim(2,1) - meshLim(1,1)+1)*nsm
    if ( ldimX /= meshl(1) ) &
         call die('Incorrect number of A divisions found.')
    meshl(2) = (meshLim(2,2) - meshLim(1,2)+1)*nsm
    if ( ldimY /= meshl(2) ) &
         call die('Incorrect number of B divisions found.')
    meshl(3) = (meshLim(2,3) - meshLim(1,3)+1)*nsm
    if ( ldimZ /= meshl(3) ) &
         call die('Incorrect number of C divisions found.')

    ! Calculate starting point for grid
    offset_i(:) = (meshLim(1,:) - 1)*nsm

#ifdef MESH_DEBUG
    ! We will print out information about the boxes that each processor
    ! has..
    if ( IONode ) then
       write(*,'(t2,a)') 'Printing offsets and corners in Ang'
       write(*,'(t3,a,3(tr1,f10.5))') &
               'dL-x: ',dL(:,1)/Ang
       write(*,'(t3,a,3(tr1,f10.5))') &
               'dL-y: ',dL(:,2)/Ang
       write(*,'(t3,a,3(tr1,f10.5))') &
               'dL-z: ',dL(:,3)/Ang
       write(*,'(t3,a,3(tr1,f10.5))') &
               'dMesh:',dMesh(:)/Ang
    end if
#ifdef MPI
    do i = 0 , Nodes - 1
       if ( Node == i ) then
          if ( IONode ) then
             ll = offset_r
             ix = ldimX
             iy = ldimY
             iz = ldimZ
             idx = offset_i
          else
             call MPI_Send(offset_r,3,MPI_Double_precision, &
                  0,0,MPI_Comm_World,MPIerror)
             call MPI_Send(ldimX,1,MPI_Integer, &
                  0,1,MPI_Comm_World,MPIerror)
             call MPI_Send(ldimY,1,MPI_Integer, &
                  0,2,MPI_Comm_World,MPIerror)
             call MPI_Send(ldimZ,1,MPI_Integer, &
                  0,3,MPI_Comm_World,MPIerror)
             call MPI_Send(offset_i,3,MPI_Integer, &
                  0,4,MPI_Comm_World,MPIerror)
          end if
       else if ( IONode ) then
          call MPI_Recv(ll,3,MPI_Double_precision, &
               i,0,MPI_Comm_World,status,MPIerror)
          call MPI_Recv(ix,1,MPI_Integer, &
               i,1,MPI_Comm_World,status,MPIerror)
          call MPI_Recv(iy,1,MPI_Integer, &
               i,2,MPI_Comm_World,status,MPIerror)
          call MPI_Recv(iz,1,MPI_Integer, &
               i,3,MPI_Comm_World,status,MPIerror)
          call MPI_Recv(idx,3,MPI_Integer, &
               i,4,MPI_Comm_World,status,MPIerror)
       end if
       if ( IONode ) then
          print *, i, ix,iy,iz
          write(*,'(t3,a,3(tr1,i8))') &
               'Indices:',idx(:)
          write(*,'(t3,a,3(tr1,f12.5))') &
               'Lower-left:',ll(:)/Ang
          write(*,'(t3,a,3(tr1,f12.5))') &
               'Upper-x   :',(ll(:)+ix * dL(:,1))/Ang
          write(*,'(t3,a,3(tr1,f12.5))') &
               'Upper-y   :',(ll(:)+iy * dL(:,2))/Ang
          write(*,'(t3,a,3(tr1,f12.5))') &
               'Upper-z   :',(ll(:)+iz * dL(:,3))/Ang
       end if
    end do
#else
    if ( IONode ) then
       write(*,'(t3,a,3(tr1,i8))') &
            'Indices:',idx(:)
       write(*,'(t3,a,3(tr1,f12.5))') &
            'Lower-left:',offset_r(:)/Ang
       write(*,'(t3,a,3(tr1,f12.5))') &
            'Upper-x   :',(offset_r(:)+ldimX * dL(:,1))/Ang
       write(*,'(t3,a,3(tr1,f12.5))') &
            'Upper-y   :',(offset_r(:)+ldimY * dL(:,2))/Ang
       write(*,'(t3,a,3(tr1,f12.5))') &
            'Upper-z   :',(offset_r(:)+ldimZ * dL(:,3))/Ang
    end if
#endif
#endif

  end subroutine init_mesh_node

  elemental subroutine mesh_correct_idx(mesh,idx)
    integer, intent(in) :: mesh
    integer, intent(inout) :: idx
    ! negative "supercell"
    do while ( idx <= 0 )
       idx = idx + mesh
    end do
    ! positive "supercell"
    do while ( mesh < idx )
       idx = idx - mesh
    end do
    
  end subroutine mesh_correct_idx
  
end module m_mesh_node

