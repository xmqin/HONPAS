! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!******************************************************************************
! MODULE m_fft
! 3-D fast complex Fourier transform
! Written by J.D.Gale (July 1999)
!******************************************************************************
!
!   PUBLIC procedures available from this module:
! fft   : 3-D complex FFT. Old version of J.D.Gale
!
!   PUBLIC parameters, types, and variables available from this module:
! none
!
!******************************************************************************
!
!   USED module procedures:
! use sys,          only: die           ! Termination routine
! use alloc,        only: de_alloc      ! De-allocation routine
! use alloc,        only: re_alloc      ! Re-allocation routine
! use parallelsubs, only: HowManyMeshPerNode
! use fft1d,        only: gpfa          ! 1-D FFT routine
! use fft1d,        only: setgpfa       ! Sets gpfa routine
! use m_timer,      only: timer_start   ! Start counting CPU time
! use m_timer,      only: timer_stop    ! Stop counting CPU time
!
!   USED MPI procedures
! use mpi_siesta
!
!   USED module parameters:
! use precision,    only: dp            ! Real double precision type
! use precision,    only: gp=>grid_p    ! Real precision type of mesh arrays
! use parallel,     only: Node, Nodes, ProcessorY
! use mesh,         only: nsm
!
!   EXTERNAL procedures used:
! gpfa    : 1D complex FFT
! setgpfa : Initializes gpfa
! timer   : CPU time counter
!
!******************************************************************************

      MODULE m_fft

      use precision,    only : dp, grid_p
      use parallel,     only : Node, Nodes, ProcessorY
      use parallelsubs, only : HowManyMeshPerNode
      use sys,          only : die
      use alloc,        only : re_alloc, de_alloc
      use mesh,         only : nsm
      use fft1d,        only : gpfa        ! 1-D FFT routine
      use fft1d,        only : setgpfa     ! Sets gpfa routine
      use m_timer,      only : timer_start ! Start counting CPU time
      use m_timer,      only : timer_stop  ! Stop counting CPU time
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

      PUBLIC:: 
     .  fft      ! 3D complex FFT

      PRIVATE ! Nothing is declared public beyond this point

      CONTAINS

!******************************************************************************

      subroutine fft( f, nMesh, isn )
C
C  Parallel FFT based on FFT routine of Templeton
C
C  On input :
C
C  real*8 f()       : contains data to be Fourier transformed
C  integer nMesh(3) : contains global dimensions of grid
C  integer isn      : indicates the direction of the transform
C
C  On exit :
C
C  real*8 f()      : contains data Fourier transformed
C
C  Julian Gale, July 1999
C
      implicit none
C
C  Passed arguments
C
      real(grid_p)
     .  f(*)
      integer
     .  nMesh(3), isn
C
C  Local variables
C 
      integer
     .  maxmaxtrigs, ntrigs,
     .  n1, n2, n3, n2l, n3l, CMeshG(3), CMeshL(3), i, n,
     .  nf, nft, ng, nmeshl, ProcessorZ, IOffSet, OldMesh(3)
      integer, save ::
     .  maxtrigs
      logical
     .  LocalPoints, lredimension
      real(dp)
     .  scale
      real(dp), dimension(:,:), pointer, save ::   trigs => null()
#ifdef MPI
      integer
     .  n1lf, nrem, Py, Pz
      real(grid_p), dimension(:), pointer  ::   ft
#endif

      save OldMesh
      data OldMesh/0,0,0/

      call timer( 'fft', 1 )
C
C  Set mesh size variables
C
      n1 = nMesh(1)
      n2 = nMesh(2)
      n3 = nMesh(3)
      CMeshG(1) = n1/nsm
      CMeshG(2) = n2/nsm
      CMeshG(3) = n3/nsm
      call HowManyMeshPerNode(CMeshG, Node, Nodes, nmeshl, CMeshL)
      n2l = CMeshL(2)*nsm
      n3l = CMeshL(3)*nsm

      n = n1*n2l*n3l
      ng = n1*n2*n3
      nf = 2*ng   ! size(f)

C
C  Set logical as to whether there are any locally stored points
C
      LocalPoints = (n.gt.0) 

C
C  Work out processor grid size
C
      ProcessorZ = Nodes/ProcessorY
      if (ProcessorY*ProcessorZ.ne.Nodes) 
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')

C
C  Allocate trigs array
C
      if (.not.associated(trigs)) then
        maxtrigs = 256
        nullify( trigs )
        call re_alloc( trigs, 1, maxtrigs, 1, 3, name='trigs',
     &                 routine='fft' )
      endif
C
C  Initialise the tables for the FFT if the mesh has changed
C
 10   lredimension = .false.
      maxmaxtrigs = maxtrigs
      do i=1,3
        if (OldMesh(i).ne.nMesh(i)) then
          call setgpfa(trigs(:,i), maxtrigs, ntrigs, nMesh(i))
          if (ntrigs.gt.maxmaxtrigs) then
            lredimension = .true.
            maxmaxtrigs = ntrigs
          else
            OldMesh(i) = nMesh(i)
          endif
        endif
      enddo
      if (lredimension) then
C
C  Resize FFT array for trig factors and set OldMesh to 0 to force recalculation
C
        maxtrigs = maxmaxtrigs
        call re_alloc( trigs, 1, maxtrigs, 1, 3, name='trigs',
     &                 routine='fft' )
        OldMesh(1:3) = 0
        goto 10
      endif

C
C  FFT in X direction
C
      if (LocalPoints) then
        call gpfa(f,f(2:nf),trigs(:,1),2,2*n1,n1,n2l*n3l,-isn)
      endif
#ifdef MPI
C***********************
C  2-D Processor Grid  *
C***********************
      n1lf = n1/ProcessorY
      nrem = n1 - n1lf*ProcessorY
      Py = (Node/ProcessorZ) + 1
      Pz = Node - (Py - 1)*ProcessorZ + 1
      if (Py.le.nrem) n1lf = n1lf + 1
C
C  Allocate local memory
C
      nullify( ft )
      call re_alloc( ft, 1, 2*n1lf*n2*n3l, name='ft', routine='fft' )
C
C  Redistribute data to be distributed by X and Z
C
      call redistribXZ(f,n1,n2l,n3l,ft,n1lf,n2,1,nsm,Node,Nodes)
C
C  FFT in Y direction
C
      if (LocalPoints) then
        nft = size(ft)
!$OMP parallel do default(shared), private(i,IOffSet)
        do i=0,n3l-1
          IOffSet=2*n1lf*n2*i
          call gpfa(ft(IOffSet+1:nft),ft(IOffSet+2:nft),trigs(:,2),
     .            2*n1lf,2,n2,n1lf,-isn)
        enddo
!$OMP end parallel do
      endif
C
C  Redistribute data back to original form
C
      call redistribXZ(f,n1,n2l,n3l,ft,n1lf,n2,-1,nsm,Node,Nodes)
C
C  Free local memory ready for re-use
C
      call de_alloc( ft, name='ft', routine='fft' )
C
C  Find new distributed x dimension
C
      n1lf = n1/ProcessorZ
      nrem = n1 - n1lf*ProcessorZ
      if (Pz.le.nrem) n1lf = n1lf + 1
C
C  Allocate local memory
C
      nullify( ft )   ! AG
      call re_alloc( ft, 1, 2*n1lf*n2l*n3, name='ft', routine='fft' )
C
C  Redistribute data to be distributed by X and Y
C
      call redistribXY(f,n1,n2l,n3l,ft,n1lf,n3,1,nsm,Node,Nodes)
C
C  FFT in Z direction
C
      if (LocalPoints) then
        nft = size(ft)
        call gpfa(ft,ft(2:nft),trigs(:,3),2*n1lf*n2l,2,n3,n1lf*n2l,-isn)
      endif
C
C  Redistribute data to be distributed by Z again
C
      call redistribXY(f,n1,n2l,n3l,ft,n1lf,n3,-1,nsm,Node,Nodes)
C
C  Free local memory
C
      call de_alloc( ft, name='ft', routine='fft' )
#else
C
C  FFT in Y direction
C
!$OMP parallel do default(shared), private(i,IOffSet)
      do i=0,n3-1
         IOffSet=2*n1*n2*i
         call gpfa(f(IOffSet+1:nf),f(IOffSet+2:nf),trigs(:,2),
     .             2*n1,2,n2,n1,-isn)
      enddo
!$OMP end parallel do
C
C  FFT in Z direction
C
      call gpfa(f,f(2:nf),trigs(:,3),2*n1*n2,2,n3,n1*n2,-isn)
#endif

C
C  Scale values
C
      if (LocalPoints.and.isn.gt.0) then
        scale=1.0_dp/dble(ng)
!$OMP parallel do default(shared), private(i)
        do i=1,2*n
          f(i)=f(i)*scale
        enddo
!$OMP end parallel do
      endif

      call timer( 'fft', 2 )

      return
      end subroutine fft


!******************************************************************************

#ifdef MPI
      subroutine redistribXZ(f,n1,n2l,n3l,ft,n1lf,n2,idir,nsm,
     .  Node,Nodes)
C
C  This routine redistributes the data over the Nodes as needed
C  for the FFT routines between the arrays f and ft
C
C  Array f is distributed in the Y/Z direction while ft is distributed
C  in the X/Z direction
C
C  idir = direction of redistribution : > 0 => f->ft
C                                       < 0 => ft->f
C
C  Use the processor grid to divide communicator according to Z
C
C  Julian Gale, July 1999
C

      implicit none

C
C  Passed arguments
C
      integer
     .  n1, n2, n2l, n1lf, n3l, Node, Nodes, idir, nsm
      real(grid_p)
     .  f(2,n1,n2l,n3l), ft(2,n1lf,n2,n3l)
C
C  Local variables
C
      integer
     .  BlockSizeY, BlockSizeYMax, jmin, jmax, jloc, n1lmax, NRemY,
     .  i, j, k, jl, kl, BNode, INode, SNode, jminS, jmaxS

      integer
     .  MPIerror, RequestR, RequestS, Status(MPI_Status_Size)

      integer, save ::
     .  ProcessorZ, Py, Pz, ZCommunicator, ZNode, ZNodes

      logical, save :: firsttimeZ = .true.

      real(grid_p), dimension(:,:,:,:), pointer  :: ftmp,ftmp2

      if (firsttimeZ) then
C
C  Determine processor grid coordinates
C
        ProcessorZ = Nodes/ProcessorY
        if (ProcessorY*ProcessorZ.ne.Nodes) 
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')
        Py = (Node/ProcessorZ) + 1
        Pz = Node - (Py - 1)*ProcessorZ + 1
C
C  Group processors into subsets by Z
C
        call MPI_Comm_Split(MPI_Comm_World, Pz, Py, ZCommunicator,
     .    MPIerror)
        call MPI_Comm_Rank(ZCommunicator,ZNode,MPIerror)
        call MPI_Comm_Size(ZCommunicator,ZNodes,MPIerror)
        firsttimeZ = .false.
      endif

      BlockSizeY = ((n2/nsm)/ProcessorY)*nsm
      NRemY = (n2 - ProcessorY*BlockSizeY)/nsm
      if (NRemY.gt.0) then
        BlockSizeYMax = BlockSizeY + nsm
      else
        BlockSizeYMax = BlockSizeY
      endif
      n1lmax = ((n1-1)/ProcessorY) + 1
C
C  Work out local dimensions
C
      jmin = (Py-1)*BlockSizeY + nsm*min(Py-1,NRemY) + 1
      jmax = jmin + BlockSizeY - 1
      if (Py-1.lt.NRemY) jmax = jmax + nsm
      jmax = min(jmax,n2)
      jloc = jmax - jmin + 1
C
C  Allocate local memory and initialise
C
      nullify( ftmp )
      call re_alloc( ftmp, 1, 2, 1, n1lmax, 1, BlockSizeYMax, 1, n3l,
     &               name='ftmp', routine='redistribXZ' )
      nullify( ftmp2 )
      call re_alloc( ftmp2, 1, 2, 1, n1lmax, 1, BlockSizeYMax, 1, n3l,
     &               name='ftmp2', routine='redistribXZ' )

!$OMP parallel do default(shared), private(i,j,k)
      do i = 1,n3l
        do j = 1,BlockSizeYMax
          do k = 1,n1lmax
            ftmp(1,k,j,i) = 0.0_grid_p
            ftmp(2,k,j,i) = 0.0_grid_p
            ftmp2(1,k,j,i) = 0.0_grid_p
            ftmp2(2,k,j,i) = 0.0_grid_p
          enddo
        enddo
      enddo
!$OMP end parallel do

      if (idir.gt.0) then
C***********
C  F -> FT *
C***********
C
C  Handle transfer of terms which are purely local
C
!$OMP parallel do default(shared), private(i,j,jl,kl,k)
        do i = 1,n3l
          do j = jmin,jmax
            jl = j - jmin + 1
            kl = 0
            do k = 1+ZNode, n1, ZNodes
              kl = kl + 1
              ft(1,kl,j,i) = f(1,k,jl,i)
              ft(2,kl,j,i) = f(2,k,jl,i)
            enddo
          enddo
        enddo
!$OMP end parallel do

C  Loop over all Node-Node vectors exchanging local data
C
        do INode = 1,ProcessorY-1
          BNode = (ZNode+INode)
          BNode = mod(BNode,ProcessorY)
          SNode = (ZNode-INode)
          SNode = mod(SNode+ProcessorY,ProcessorY)
C
C  Collect data to send
!$OMP parallel do default(shared), private(i,jl,kl,k)
          do i = 1,n3l
            do jl = 1,jloc
              kl = 0
              do k = 1+BNode, n1, ZNodes
                kl = kl + 1
                ftmp(1,kl,jl,i) = f(1,k,jl,i)
                ftmp(2,kl,jl,i) = f(2,k,jl,i)
              enddo
            enddo
          enddo
!$OMP end parallel do

C
C  Exchange data - send to right and receive from left
C
          call MPI_IRecv(ftmp2(1,1,1,1),2*n1lmax*BlockSizeYMax*n3l,
     .      MPI_grid_real,SNode,1,ZCommunicator,RequestR,MPIerror)
          call MPI_ISend(ftmp(1,1,1,1),2*n1lmax*BlockSizeYMax*n3l,
     .      MPI_grid_real,BNode,1,ZCommunicator,RequestS,MPIerror)
C
C  Wait for receive to complete
C
          call MPI_Wait(RequestR,Status,MPIerror)
C
C  Place received data into correct array
C
          jminS = SNode*BlockSizeY + nsm*min(SNode,NRemY) + 1
          jmaxS = jminS + BlockSizeY - 1
          if (SNode.lt.NRemY) jmaxS = jmaxS + nsm
          jmaxS = min(jmaxS,n2)
!$OMP parallel do default(shared), private(i,j,jl,kl,k)
          do i = 1,n3l
            do j = jminS,jmaxS
              jl = j - jminS + 1
              kl = 0
              do k = 1+ZNode, n1, ZNodes
                kl = kl + 1
                ft(1,kl,j,i) = ftmp2(1,kl,jl,i)
                ft(2,kl,j,i) = ftmp2(2,kl,jl,i)
              enddo
            enddo
          enddo
!$OMP end parallel do

C  Wait for send to complete
C
          call MPI_Wait(RequestS,Status,MPIerror)
        enddo
      elseif (idir.lt.0) then
C***********
C  FT -> F *
C***********
C
C  Handle transfer of terms which are purely local

!$OMP parallel do default(shared), private(i,j,jl,kl,k)
        do i = 1,n3l
          do j = jmin,jmax
            jl = j - jmin + 1
            kl = 0
            do k = 1+ZNode, n1, ZNodes
              kl = kl + 1
              f(1,k,jl,i) = ft(1,kl,j,i) 
              f(2,k,jl,i) = ft(2,kl,j,i) 
            enddo
          enddo
        enddo
!$OMP end parallel do

C  Loop over all Node-Node vectors exchanging local data
C
        do INode = 1,ProcessorY-1
          BNode = (ZNode+INode)
          BNode = mod(BNode,ProcessorY)
          SNode = (ZNode-INode)
          SNode = mod(SNode+ProcessorY,ProcessorY)
C
C  Collect data to send
C
          jminS = SNode*BlockSizeY + nsm*min(SNode,NRemY) + 1
          jmaxS = jminS + BlockSizeY - 1
          if (SNode.lt.NRemY) jmaxS = jmaxS + nsm
          jmaxS = min(jmaxS,n2)
!$OMP parallel do default(shared), private(i,j,jl,kl,k)
          do i = 1,n3l
            do j = jminS,jmaxS
              jl = j - jminS + 1
              kl = 0
              do k = 1+ZNode, n1, ZNodes
                kl = kl + 1
                ftmp(1,kl,jl,i) = ft(1,kl,j,i) 
                ftmp(2,kl,jl,i) = ft(2,kl,j,i) 
              enddo
            enddo
          enddo
!$OMP end parallel do

C  Exchange data - send to right and receive from left
C
          call MPI_IRecv(ftmp2(1,1,1,1),2*n1lmax*BlockSizeYMax*n3l,
     .      MPI_grid_real,BNode,1,ZCommunicator,RequestR,MPIerror)
          call MPI_ISend(ftmp(1,1,1,1),2*n1lmax*BlockSizeYMax*n3l,
     .      MPI_grid_real,SNode,1,ZCommunicator,RequestS,MPIerror)
C
C  Wait for receive to complete
C
          call MPI_Wait(RequestR,Status,MPIerror)
C
C  Place received data into correct array

!$OMP parallel do default(shared), private(i,jl,kl,k)
          do i = 1,n3l
            do jl = 1,jloc
              kl = 0
              do k = 1+BNode, n1, ZNodes
                kl = kl + 1
                f(1,k,jl,i) = ftmp2(1,kl,jl,i) 
                f(2,k,jl,i) = ftmp2(2,kl,jl,i) 
              enddo
            enddo
          enddo
!$OMP end parallel do

C  Wait for send to complete
C
          call MPI_Wait(RequestS,Status,MPIerror)
        enddo
      endif
C
C  Free local memory
C
      call de_alloc( ftmp2, name='ftmp2', routine='redistribXZ' )
      call de_alloc( ftmp, name='ftmp', routine='redistribXZ' )

      end subroutine redistribXZ

!******************************************************************************

      subroutine redistribXY(f,n1,n2l,n3l,ft,n1lf,n3,idir,nsm,
     .  Node,Nodes)
C
C  This routine redistributes the data over the Nodes as needed
C  for the FFT routines between the arrays f and ft
C
C  Array f is distributed in the Y/Z direction while ft is distributed
C  in the X/Y direction
C
C  idir = direction of redistribution : > 0 => f->ft
C                                       < 0 => ft->f
C
C  Use the processor grid to divide communicator according to Y
C
C  Currently written in not the most efficient but simple way! 
C  Need to improve communication later.
C
C  Julian Gale, July 1999
C

      implicit none

C
C  Passed arguments
C
      integer
     .  n1, n3, n2l, n1lf, n3l, nsm, Node, Nodes, idir
      real(grid_p)
     .  f(2,n1,n2l,n3l), ft(2,n1lf,n2l,n3)
C
C  Local variables
C
      integer
     .  i, j, k, il, kl, BNode, BlockSizeZ, imin, imax, iloc,
     .  INode, n1lmax, SNode, iminS, imaxS, NRemZ, BlockSizeZMax

      integer
     .  MPIerror, RequestR, RequestS, Status(MPI_Status_Size)

      integer, save ::
     .  ProcessorZ, Py, Pz, YCommunicator, YNode, YNodes

      logical, save :: firsttimeY = .true.

      real(grid_p), dimension(:,:,:,:), pointer :: ftmp,ftmp2

      if (firsttimeY) then
C
C  Determine processor grid coordinates
C
        ProcessorZ = Nodes/ProcessorY
        if (ProcessorY*ProcessorZ.ne.Nodes) 
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')
        Py = (Node/ProcessorZ) + 1
        Pz = Node - (Py - 1)*ProcessorZ + 1
C
C  Group processors into subsets by Y
C
        call MPI_Comm_Split(MPI_Comm_World, Py, Pz, YCommunicator,
     .    MPIerror)
        call MPI_Comm_Rank(YCommunicator,YNode,MPIerror)
        call MPI_Comm_Size(YCommunicator,YNodes,MPIerror)
        firsttimeY = .false.
      endif

      BlockSizeZ = ((n3/nsm)/ProcessorZ)*nsm
      NRemZ = (n3 - ProcessorZ*BlockSizeZ)/nsm
      if (NRemZ.gt.0) then
        BlockSizeZMax = BlockSizeZ + nsm
      else
        BlockSizeZMax = BlockSizeZ
      endif
      n1lmax = ((n1-1)/ProcessorZ) + 1
C
C  Allocate local memory and initialise
C
      nullify( ftmp )
      call re_alloc( ftmp, 1, 2, 1, n1lmax, 1, n2l, 1, BlockSizeZMax,
     &               name='ftmp', routine='redistribXY' )
      nullify( ftmp2 )
      call re_alloc( ftmp2, 1, 2, 1, n1lmax, 1, n2l, 1, BlockSizeZMax,
     &               name='ftmp2', routine='redistribXY' )

!$OMP parallel do default(shared), private(i,j,k)
      do i = 1,BlockSizeZMax
        do j = 1,n2l
          do k = 1,n1lmax
            ftmp(1,k,j,i) = 0.0_grid_p
            ftmp(2,k,j,i) = 0.0_grid_p
            ftmp2(1,k,j,i) = 0.0_grid_p
            ftmp2(2,k,j,i) = 0.0_grid_p
          enddo
        enddo
      enddo
!$OMP end parallel do

C  Work out local dimensions
C
      imin = (Pz-1)*BlockSizeZ + nsm*min(Pz-1,NRemZ) + 1
      imax = imin + BlockSizeZ - 1
      if (Pz-1.lt.NRemZ) imax = imax + nsm
      imax = min(imax,n3)
      iloc = imax - imin + 1

      if (idir.gt.0) then
C***********
C  F -> FT *
C***********
C
C  Handle transfer of terms which are purely local

!$OMP parallel do default(shared), private(i,il,j,kl,k)
        do i = imin,imax
          il = i - imin + 1
          do j = 1,n2l
            kl = 0
            do k = 1+YNode, n1, YNodes
              kl = kl + 1
              ft(1,kl,j,i) = f(1,k,j,il)
              ft(2,kl,j,i) = f(2,k,j,il)
            enddo
          enddo
        enddo
!$OMP end parallel do

C  Loop over all Node-Node vectors exchanging local data
C
        do INode = 1,ProcessorZ-1
          BNode = (YNode+INode)
          BNode = mod(BNode,ProcessorZ)
          SNode = (YNode-INode)
          SNode = mod(SNode+ProcessorZ,ProcessorZ)
C
C  Collect data to send

!$OMP parallel do default(shared), private(il,j,kl,k)
          do il = 1,iloc
            do j = 1,n2l
              kl = 0
              do k = 1+BNode, n1, YNodes
                kl = kl + 1
                ftmp(1,kl,j,il) = f(1,k,j,il)
                ftmp(2,kl,j,il) = f(2,k,j,il)
              enddo
            enddo
          enddo
!$OMP end parallel do

C  Exchange data - send to right and receive from left
C
          call MPI_IRecv(ftmp2(1,1,1,1),2*n1lmax*n2l*BlockSizeZMax,
     .      MPI_grid_real,SNode,1,YCommunicator,RequestR,MPIerror)
          call MPI_ISend(ftmp(1,1,1,1),2*n1lmax*n2l*BlockSizeZMax,
     .      MPI_grid_real,BNode,1,YCommunicator,RequestS,MPIerror)
C
C  Wait for receive to complete
C
          call MPI_Wait(RequestR,Status,MPIerror)
C
C  Place received data into correct array
C
          iminS = SNode*BlockSizeZ + nsm*min(SNode,NRemZ) + 1
          imaxS = iminS + BlockSizeZ - 1
          if (SNode.lt.NRemZ) imaxS = imaxS + nsm
          imaxS = min(imaxS,n3)
!$OMP parallel do default(shared), private(i,il,j,kl,k)
          do i = iminS,imaxS
            il = i - iminS + 1
            do j = 1,n2l
              kl = 0
              do k = 1+YNode, n1, YNodes
                kl = kl + 1
                ft(1,kl,j,i) = ftmp2(1,kl,j,il)
                ft(2,kl,j,i) = ftmp2(2,kl,j,il)
              enddo
            enddo
          enddo
!$OMP end parallel do

C  Wait for send to complete
C
          call MPI_Wait(RequestS,Status,MPIerror)
        enddo
      elseif (idir.lt.0) then
C***********
C  FT -> F *
C***********
C
C  Handle transfer of terms which are purely local

!$OMP parallel do default(shared), private(i,il,j,kl,k)
        do i = imin,imax
          il = i - imin + 1
          do j = 1,n2l
            kl = 0
            do k = 1+YNode, n1, YNodes
              kl = kl + 1
              f(1,k,j,il) = ft(1,kl,j,i)
              f(2,k,j,il) = ft(2,kl,j,i)
            enddo
          enddo
        enddo
!$OMP end parallel do

C  Loop over all Node-Node vectors exchanging local data
C
        do INode = 1,ProcessorZ-1
          BNode = (YNode+INode)
          BNode = mod(BNode,ProcessorZ)
          SNode = (YNode-INode)
          SNode = mod(SNode+ProcessorZ,ProcessorZ)
C
C  Collect data to send
C
          iminS = SNode*BlockSizeZ + nsm*min(SNode,NRemZ) + 1
          imaxS = iminS + BlockSizeZ - 1
          if (SNode.lt.NRemZ) imaxS = imaxS + nsm
          imaxS = min(imaxS,n3)
!$OMP parallel do default(shared), private(i,il,j,kl,k)
          do i = iminS,imaxS
            il = i - iminS + 1
            do j = 1,n2l
              kl = 0
              do k = 1+YNode, n1, YNodes
                kl = kl + 1
                ftmp(1,kl,j,il) = ft(1,kl,j,i)
                ftmp(2,kl,j,il) = ft(2,kl,j,i)
              enddo
            enddo
          enddo
!$OMP end parallel do

C
C  Exchange data - send to right and receive from left
C
          call MPI_IRecv(ftmp2(1,1,1,1),2*n1lmax*n2l*BlockSizeZMax,
     .      MPI_grid_real,BNode,1,YCommunicator,RequestR,MPIerror)
          call MPI_ISend(ftmp(1,1,1,1),2*n1lmax*n2l*BlockSizeZMax,
     .      MPI_grid_real,SNode,1,YCommunicator,RequestS,MPIerror)
C
C  Wait for receive to complete
C
          call MPI_Wait(RequestR,Status,MPIerror)
C
C  Place received data into correct array
C
!$OMP parallel do default(shared), private(il,j,kl,k)
          do il = 1,iloc
            do j = 1,n2l
              kl = 0
              do k = 1+BNode, n1, YNodes
                kl = kl + 1
                f(1,k,j,il) = ftmp2(1,kl,j,il)
                f(2,k,j,il) = ftmp2(2,kl,j,il)
              enddo
            enddo
          enddo
!$OMP end parallel do

C  Wait for send to complete
C
          call MPI_Wait(RequestS,Status,MPIerror)
        enddo
      endif
C
C  Free local memory
C
      call de_alloc( ftmp2, name='ftmp2', routine='redistribXY' )
      call de_alloc( ftmp, name='ftmp', routine='redistribXY' )

      return
      end subroutine redistribXY
#endif

!******************************************************************************

      END MODULE m_fft
