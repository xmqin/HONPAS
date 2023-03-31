! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
subroutine setup_ordern_indexes(no_l,no_u,Nodes)

  use spatial, only:  nL2G, nG2L, nNode, nOrbPerNode
  use domain_decom, only :  use_dd_perm, dd_nuo, dd_perm, &
                            ulimit, llimit, dd_invp, dd_nnode
  use alloc, only: re_alloc
  use parallel, only: Node

#ifdef MPI
  use mpi_siesta
#endif

  implicit none

  integer, intent(in) :: no_l, no_u
  integer, intent(in) :: Nodes

#ifdef MPI
  integer, dimension(:), allocatable :: nl2gtmp, counts, displs
  integer :: i, MPIerr
#endif
  integer :: LOrb, GOrb, j, iu
  character(len=14) :: name

  call re_alloc( nL2G, 1, no_u, 1, Nodes, 'nL2G', 'setup_ordern_indexes' )
  call re_alloc( nG2L, 1, no_u, 'nG2L', 'setup_ordern_indexes' )
  call re_alloc( nNode, 1, no_u, 'nNode', 'setup_ordern_indexes' )
  call re_alloc( nOrbPerNode, 1, Nodes, 'nOrbPerNode','setup_ordern_indexes' )

! Global to local index and node mapping

! Note: nG2L does not need to be globalized, since it
!       is called by the local node, and all the information
!       has been already setup in domain_decom.F.
!       Julian Gale uses a globalized version, but checks
!       first whether the node asking is indeed the one
!       assigned to the orbital (see GlobalToLocalOrb in
!       parallelsubs.F)

!       nNode is already global (bcast in domain_decom.F)

  do Gorb = 1, no_u
     if (use_dd_perm) then
        Lorb = dd_perm(GOrb)
     else
        if (GOrb.ge.llimit .and. GOrb.lt.ulimit) then
           LOrb = GOrb - llimit + 1
        else
           LOrb = 0
        endif
     endif
     nG2L(Gorb) = LOrb
     nNode(Gorb) = dd_nnode(Gorb)
  enddo

! Local to global index
! (Local node only)

#ifdef MPI
  allocate(nl2gtmp(1:no_u))
  nl2gtmp(1:no_u) = 0
#else
  nL2G(:,1) = 0
#endif

  do LOrb = 1, no_l
     if (use_dd_perm) then
        GOrb = dd_invp(LOrb)
     else
        GOrb = LOrb + llimit - 1
     endif
     !     nL2G(LOrb,Node+1) = Gorb
#ifdef MPI
     nl2gtmp(LOrb) = Gorb
#else
     nL2G(LOrb,1) = Gorb
#endif
  enddo

#ifdef MPI
  ! Now all_gatherv
  allocate(counts(1:Nodes),displs(1:Nodes))
  do i = 1, Nodes
     counts(i) = no_u
     displs(i) = (i-1)*no_u
  enddo
  call mpi_allgatherv(nl2gtmp,no_u,MPI_Integer,nl2G(1,1),counts,displs, &
                      MPI_integer,MPI_Comm_World, MPIerr)
  deallocate(nl2gtmp,counts,displs)
#endif

! Number of orbitals per node

  ! nOrbPerNode(Node+1) = dd_nuo
#ifdef MPI
  ! Now all_gather...
  call mpi_allgather(dd_nuo,1,MPI_Integer,nOrbPerNode,1, &
                      MPI_integer,MPI_Comm_World, MPIerr)
#else
  nOrbPerNode(1) = dd_nuo
#endif

   write(name,"(a11,i3.3)")  "ON_INDEXES.", Node
   call io_assign(iu)
   open(unit=iu,file=name,form="formatted")
   write(iu,*) "nl2g"
   do Lorb = 1, no_u	
      write(iu, "(i6,8i9)") Lorb, (nL2G(Lorb,j),j=1,Nodes)
   enddo
   write(iu,*) "nNode"
   do Gorb = 1, no_u	
      write(iu, "(i6,i5)") Gorb, nNode(Gorb)
   enddo
   write(iu,*) "nG2L (node 0)"
   do Gorb = 1, no_u	
      write(iu, "(i6,i5)") Gorb, nG2L(Gorb)
   enddo
   write(iu,*) "nOrbPerNode"
   do j = 1, Nodes
      write(iu, "(i6,i8)") j-1, nOrbPerNode(j)
   enddo
   call io_close(iu)



end subroutine setup_ordern_indexes
