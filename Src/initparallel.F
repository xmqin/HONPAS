! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine initparallel( no_u, na_u, lasto, xa,
     &                         ucell, rmaxh )
C
C  Initialises parallel parameters
C
C  Julian Gale, NRI, Curtin Uni, March 2004
C
      use precision
      use fdf
      use siesta_options, only: want_domain_decomposition
      use siesta_options, only: want_spatial_decomposition
      use siesta_options, only: isolve, SOLVE_ORDERN, SOLVE_PEXSI
      use siesta_options, only: rcoor
      use parallel,     only : Node, Nodes, BlockSize, ProcessorY
      use parallelsubs, only : WhichNodeOrb, GetNodeOrbs
      use parallelsubs, only : LocalToGlobalOrb, GlobalToLocalOrb
      use parallelsubs, only : pexsi_dist, pexsi_bs
      use spatial
      use sys
      use alloc,        only : re_alloc
      use domain_decom, only: use_dd, use_dd_perm, preSetOrbitLimits
      use parallelsubs, only : set_processorY, set_processorYdefault
#ifdef MPI
      use parallelsubs, only : set_BlockSizeDefault
      use mpi_siesta,   only: MPI_Comm_World, MPI_Comm_Self
#endif
      use class_OrbitalDistribution
      use sparse_matrices,     only : block_dist, single_dist

      implicit none

#ifdef MPI
      integer, external :: numroc
#endif
C
C  Passed variables
C
      integer       :: no_u
      integer       :: na_u
      integer       :: lasto(0:na_u)
      real(dp)      :: xa(3,na_u)
      real(dp)      :: ucell(3,3)
      real(dp)      :: rmaxh
C
C  Local variables
C
#ifdef MPI
      integer       :: blocksizedefault
      integer       :: procYdefault
      integer       :: procYval
      integer       :: n1, nn
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
      real(dp)      :: bbeta
      real(dp)      :: gamma
      real(dp)      :: degtorad
      real(dp)      :: rcmax
      real(dp)      :: rcmaxopt
      real(dp), save:: tiny = 1.0d-6

#ifdef DEBUG
      call write_debug( '  PRE initparallel' )
#endif

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
      if (isolve.eq.SOLVE_ORDERN) then
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

      rspatial = fdf_physical('RcSpatial',0.0d0,'Bohr')

C Processor grid
      npgrid(1) = fdf_integer('ProcessorGridX',1)
      npgrid(2) = fdf_integer('ProcessorGridY',1)
      npgrid(3) = fdf_integer('ProcessorGridZ',1)

C Set ProcessorY  
#ifdef MPI
      call set_processorYdefault(Nodes,procYdefault)
      procYval = fdf_integer('processorY',procYdefault)
      ! Sanity check
      if (procYval > Nodes) then
         if (Node.eq.0) then
            write(6,'(a)')
     .       '* Warning: ProcessorY > Nodes in fdf file. '//
     $       'Reset to estimated optimal value.'
         endif
         call set_processorY(procYdefault)
      else
         call set_processorY(procYval)
      endif
         
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
      call uncell(ucell,acell,bcell,ccell,alpha,bbeta,gamma,degtorad)

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

      if (.not. lspatial) then
C Set BlockSize
#ifdef MPI
        pexsi_dist = .false.
        call set_blocksizedefault(Nodes,no_u,blocksizedefault)
        BlockSize = fdf_integer('blocksize',blocksizedefault)
        call newDistribution(blocksize,MPI_Comm_World,block_dist)
#else
        BlockSize = 8
        call newDistribution(blocksize,-1,block_dist)
#endif
      endif

      ! Create single communicator
#ifdef MPI
      call newDistribution(no_u,MPI_Comm_Self,single_dist,
     &     name='Single-dist')
#else
      call newDistribution(no_u,-1           ,single_dist,
     &     name='Singledist')
#endif


#ifdef MPI
C Output indication of parallel parameters
      if (Node.eq.0) then
        if (lspatial) then
           write(6,'(/,a,f8.4,/)')
     .        '* Spatial decomposition: Cutoff = ',rcmax
        else
           write(6,'(/,a,2i4,/)')
     .        '* ProcessorY, Blocksize: ', processorY, Blocksize
        endif
        n1 = numroc(no_u,Blocksize,0,0,Nodes)
        nn = numroc(no_u,Blocksize,Nodes-1,0,Nodes)
        write(6,"(/,a,2i6,/)")
     $       "* Orbital distribution balance (max,min):", n1, nn
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
        if (Node.eq.0) then
           write(6,'(/,a,i4,/)')
     .        '* ProcessorY ', processorY
        endif

      elseif ((isolve.eq.SOLVE_ORDERN) .AND. (.not. use_dd)) then
C
C  Build dummy lists linking orbitals with parallel structure
C
        maxorb = lasto(na_u)
        call re_alloc( nL2G, 1, maxorb, 1, Nodes, 'nL2G',
     &                 'initparallel' )
        call re_alloc( nG2L, 1, maxorb, 'nG2L', 'initparallel' )
        call re_alloc( nNode, 1, maxorb, 'nNode', 'initparallel' )
        call re_alloc( nOrbPerNode, 1, Nodes, 'nOrbPerNode',
     &                 'initparallel' )
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

#ifdef DEBUG
      call write_debug( '  POS initparallel' )
#endif
      end subroutine initparallel
