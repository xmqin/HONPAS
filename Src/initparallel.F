! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      subroutine initparallel(no_u,na_u,lasto,xa,ucell,rmaxh,rcoor,
     .                        isolve)
C
C  Initialises parallel parameters
C
C  Julian Gale, NRI, Curtin Uni, March 2004
C
      use precision
      use fdf
      use m_fdf_global, only: fdf_global_get
      use siesta_options, only: want_domain_decomposition
      use siesta_options, only: want_spatial_decomposition
      use parallel,     only : Node, Nodes, BlockSize, ProcessorY
      use parallelsubs, only : set_processorY, set_BlockSizeDefault
      use parallelsubs, only : set_processorYdefault
      use parallelsubs, only : WhichNodeOrb, GetNodeOrbs
      use parallelsubs, only : LocalToGlobalOrb, GlobalToLocalOrb
      use spatial
      use sys
      use alloc
      use domain_decom, only: use_dd, use_dd_perm, preSetOrbitLimits
#ifdef MPI
      use mpi_siesta, only: mpi_comm_world, mpi_integer
      use mpi_siesta, only: mpi_double_precision
#endif
      implicit none
C
C  Passed variables
C
      integer       :: no_u
      integer       :: na_u
      integer       :: isolve
      integer       :: lasto(0:na_u)
      real(dp)      :: xa(3,na_u)
      real(dp)      :: ucell(3,3)
      real(dp)      :: rcoor
      real(dp)      :: rmaxh
C
C  Local variables
C
#ifdef MPI
      integer       :: blocksizedefault
      integer       :: procYdefault
      integer       :: procYval
      integer       :: MPIerror
#endif
      integer       :: ii, j, LOrb, GOrb, iu
      integer       :: maxorb
      integer       :: ncell(3)
      integer       :: ncelltotal
      logical       :: lspatialok
      logical, save :: first = .true.
      real(dp)      :: acell
      real(dp)      :: bcell
      real(dp)      :: ccell
      real(dp)      :: cellmin
      real(dp)      :: cscale
      real(dp)      :: alpha
      real(dp)      :: beta
      real(dp)      :: gamma
      real(dp)      :: degtorad
      real(dp)      :: rcmax
      real(dp)      :: rcmaxopt
      real(dp), save:: tiny = 1.0d-6


C  Initialise on first call
      if (first) then
        nullify(natomsNode)
        nullify(natomsL2G)
        nullify(natomsG2L)
        nullify(ncellnodeptr)
        nullify(lbuffercell)
        nullify(nL2G)
        nullify(nG2L)
        nullify(nNode)
        nullify(nOrbPerNode)
        first = .false.
      endif

C Set spatial decomposition flag according to SCF method

      lspatial = .false.
#ifdef MPI
      if (isolve.eq.1) then
         if (want_domain_decomposition) then
            use_dd = .true.
            use_dd_perm = .false.
            lspatial = .false.
         else
            lspatial = want_spatial_decomposition
         endif
      else
         lspatial = .false.
      endif
#endif

C Spatial decomposition flag

      call fdf_global_get(rspatial, 'RcSpatial', 0.0_dp, 'Bohr' )

C Processor grid
      call fdf_global_get(npgrid(1), 'ProcessorGridX',1)
      call fdf_global_get(npgrid(2), 'ProcessorGridY',1)
      call fdf_global_get(npgrid(3), 'ProcessorGridZ',1)

C Set ProcessorY  
#ifdef MPI
      if (Node.eq.0) then
        call set_processorYdefault(Nodes,procYdefault)
        procYval = fdf_integer('processorY',procYdefault)
      endif
      call MPI_Bcast(procYval,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call set_processorY(procYval)
#else
      call set_processorY(1)
#endif

C Override spatial cut-off
      if (lspatial.and.rspatial.gt.0.0d0) then
        rcmax = rspatial
      else
        rcmax = max(rmaxh,rcoor)

C Check that this will lead to a reasonable number of domains relative to 
C number of processors and if not then modify
        ncell(1) = abs(ucell(1,1)/rcmax) + 1
        ncell(2) = abs(ucell(2,2)/rcmax) + 1
        ncell(3) = abs(ucell(3,3)/rcmax) + 1
        ncelltotal = ncell(1)*ncell(2)*ncell(3)
        if (ncelltotal.lt.Nodes) then
          cscale = dble(2*Nodes)/dble(ncelltotal)
          cscale = cscale**(1.0d0/3.0d0)
          ncell(1) = ncell(1)*cscale + 1
          ncell(2) = ncell(2)*cscale + 1
          ncell(3) = ncell(3)*cscale + 1
          rcmax = (ucell(1,1)-tiny)/dble(ncell(1))
          rcmax = min(rcmax,(ucell(2,2)-tiny)/dble(ncell(2)))
          rcmax = min(rcmax,(ucell(3,3)-tiny)/dble(ncell(3)))
        endif

        rcmaxopt = rcmax
      endif

C Find cell parameters
      degtorad = 4.0d0*atan(1.0d0)/180.0d0
      call uncell(ucell,acell,bcell,ccell,alpha,beta,gamma,degtorad)

C Check that cut-off is less than cell parameter
      cellmin = min(acell, bcell, ccell)
      if (cellmin.lt.2.0_dp*rcmax) then
        ! how big do we want the cells to be?
        !  Ideally, each node should have either
        ! one or zero atoms, which allows for the
        ! best division of cells to processors.
        !  Which rather suggests an alternative
        ! strategy, of course. But anyway.
        !  They should be small enough that there is
        ! a reasonable division of cells per node.
        !  However, they should not be so small that 
        ! we have so many cells that our accounting
        ! arrarys become stupidly long. Therefore
        ! I will arbitrarily decide (for the moment)
        ! that we should have 
        if (Nodes<16) then
          rcmax = 0.5_dp*cellmin/dble(Nodes)
        else
          rcmax = cellmin/32._dp
        endif
      endif

      if (.not.lspatial) then
C Set BlockSize
#ifdef MPI
        if (Node.eq.0) then
          call set_blocksizedefault(Nodes,no_u,blocksizedefault)
          BlockSize = fdf_integer('blocksize',blocksizedefault)
        endif
        call MPI_Bcast(BlockSize,1,MPI_integer,0,MPI_Comm_World,
     .    MPIerror)
#else
        call set_processorY(1)
        BlockSize = 8
#endif
      endif

#ifdef MPI
C Output indication of parallel parameters
      if (Node.eq.0) then
        if (lspatial) then
          write(6,'(/,a,f8.4,/)')
     .        '* Spatial decomposition: Cutoff = ',rcmax
        else
          write(6,'(/,a,2i4,/)')
     .        '* ProcessorY, Blocksize: ', procYval, Blocksize
        endif
      endif
#endif

C If spatial perform allocation of atoms to processors
      if (lspatial) then

        call setspatial(na_u,xa,ucell,rcmax,lspatialok)
        if (.not.lspatialok) then
          call die('Spatial decomposition failed')
        endif
        call setatomnodes(na_u,lasto,Node,Nodes)
        if (Nodes.gt.1) then
          ProcessorY = npgrid(2)
        else
          ProcessorY = 1
        endif

      elseif (.not. use_dd) then
C
C  Build dummy lists linking orbitals with parallel structure
C
        maxorb = lasto(na_u)
        call re_alloc(nL2G,1,maxorb,1,Nodes,name='nL2G')
        call re_alloc(nG2L,1,maxorb,name='nG2L')
        call re_alloc(nNode,1,maxorb,name='nNode')
        call re_alloc(nOrbPerNode,1,Nodes,name='nOrbPerNode')
C

        do j = 1, Nodes
           call GetNodeOrbs(maxorb,j-1,Nodes,nOrbPerNode(j))
        enddo
        nG2L(:) = 0
        do ii = 1,maxorb
          call GlobalToLocalOrb(ii,Node,Nodes,nG2L(ii))
          call WhichNodeOrb(ii, Nodes, nNode(ii))
          do j = 1, Nodes
             call LocalToGlobalOrb(ii,j-1,Nodes,nL2G(ii,j))
          enddo
        enddo

         if (node==0) then
            call io_assign(iu)
            open(unit=iu,file="BLOCK_INDEXES",form="formatted")
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
         endif

      endif

      if (use_dd) then
        call preSetOrbitLimits( no_u )
      endif

      return
      end
