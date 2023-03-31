! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module globalise
#ifdef MPI

      use precision, only: dp
      use spatial, only: nL2G
      use spatial, only: nNode

      use on_core, only : numft, listft
      use on_core,   only : numhij, listhij
      use on_main,   only : ncL2G, ncG2L, nbandsloc, nbL2G
      use alloc,     only : re_alloc, de_alloc
      use mpi_siesta, only: mpi_comm_world
      use mpi_siesta, only: mpi_status_size, mpi_integer
      use mpi_siesta, only: mpi_sum, mpi_logical
      use mpi_siesta, only: mpi_double_precision

      implicit none

      public :: setglobalise, setglobaliseF, setglobaliseB,
     $          globalisestats, globalinitb, globalloadb2,
     $          globalcommb, globalreloadb2, globalrezerob2,
     $          globalisef, globalloadb3, globalreloadb3,
     $          globalrezerob3, globalloadb1, globalreloadb1,
     $          globalrezerob1, globalisec

C
C  NNodes                           = Number of neighbouring nodes
C  MaxBandsSendB                    = Maximum no. of basis functions to be transfered between nodes for bux1/2
C  MaxBasisSendC                    = Maximum no. of basis functions to be transfered between nodes for C
C  NNodeNo(NNodes)                  = Pointer to node numbers for each neighbouring node
C  Node2NNodePtr(Nodes)             = Pointer from Global to Local node number for communication
C  NNodeBandsSendB(NNodes)          = Number of bands to send to each neighbouring node for bux1/2
C  NNodeBasisSendC(NNodes)          = Number of basis functions to send to each neighbouring node for C
C  NNodeBandsSendListB(Ns,NNodes)   = List of Ns bands to send to each neighbouring node for bux1/2
C  NNodeBasisSendListC(Ns,NNodes)   = List of Ns basis functions to send to each neighbouring node for C
C
C  where i is a band index
C
        integer,          save :: NNodes
        integer,          save :: MaxBandsSendB
        integer,          save :: MaxBasisSendC
        integer, pointer, save :: NNodeNo(:) => null()
        integer, pointer, save :: NNodeInBux(:) => null()
        integer, pointer, save :: NNodeInBux2(:) => null()
        integer, pointer, save :: NNodeInBuxPreN(:) => null()
        integer, pointer, save :: NNodeInBux2PreN(:) => null()
        integer, pointer, save :: Node2NNodePtr(:) => null()
        integer, pointer, save :: NNodeBandsSendB(:) => null()
        integer, pointer, save :: NNodeBasisSendC(:) => null()
        integer, pointer, save :: NNodeBandsSendListB(:,:) => null()
        integer, pointer, save :: NNodeBandsSendListB2(:,:) => null()
        integer, pointer, save :: NNodeBasisSendListC(:,:) => null()
        integer, pointer, save, public :: numft2(:) => null()
        integer, pointer, save :: numftG(:,:) => null()
        integer, pointer, save :: numftL(:,:) => null()
        integer, pointer, save :: numhijG(:,:) => null()
        integer, pointer, save :: numhijL(:,:) => null()
        integer, pointer, save :: listhijG(:) => null()
        integer, pointer, save :: listhijL(:) => null()
        integer, pointer, save, public  :: listft2(:,:) => null()
        integer, pointer, save :: listftG(:,:) => null()
        integer, pointer, save :: listftL(:,:) => null()
        integer,          save :: MaxListHijG
        integer,          save :: MaxListHijL
        integer,          save :: MaxFtG
        integer,          save :: MaxFtL
        integer,          save, public  :: MaxFt2
C
C  Data store for communication of B
C
        integer, pointer, save :: ndataBR(:) => null()
        integer, pointer, save :: ndataBS(:) => null()
        integer, pointer, save :: ndataBR2(:) => null()
        integer, pointer, save :: ndataFR(:) => null()
        integer, pointer, save :: ndataFS(:) => null()
        real(dp),pointer, save :: buxstore(:) => null()
        real(dp),pointer, save :: buxstore2(:) => null()
C
C  Statistics for data transfer
C
        integer,          save :: nTransfer = 0
        integer,          save :: nTransferInteger = 0
        integer,          save :: nTransferDouble = 0

      PRIVATE

      contains

      subroutine setglobalise(no_u,no_l,no_cl,maxnh,numh,
     .                        listh,listhptr,Node,Nodes)
C
C  Initialise information required for the globalisation of arrays
C  by passing to neighbouring nodes.
C
C  Julian Gale, NRI, Curtin University of Technology, April 2004
C

      implicit none

C Passed variables
      integer ::
     .  no_u,no_l,no_cl,maxnh,Node,Nodes,
     .  numh(no_l),listhptr(no_l),listh(maxnh)

C Internal variables
      integer                 :: ind
      integer                 :: m
      integer                 :: mn
      integer                 :: MPIerror
      integer                 :: mu
      integer                 :: mug
      integer                 :: n
      integer                 :: nRequest
      integer, pointer, save :: Request(:) => null()
      integer, pointer, save :: Status(:,:) => null()
      logical, pointer, save :: lConnected(:,:) => null()
      logical, pointer, save :: lConnectedL(:) => null()
      logical, pointer, save :: lNodeNeeded(:) => null()

C Call timer
      call timer('setglobal',1)

C--------------------------
C Find neighbouring nodes -
C--------------------------
C Check orbitals to see which ones are connected to orbitals on other nodes
      call re_alloc( lConnected, 1, no_u, 1, Nodes, 'lConnected',
     &               'globalise' )

      lConnected(1:no_u,1:Nodes) = .false.
      do m = 1,no_l
        ind = listhptr(m)
        mug = nL2G(m,Node+1)
        do mn = 1,numh(m)
          n = listh(ind+mn)
          if (nNode(n).ne.Node.or.nNode(mug).ne.Node) then
            lConnected(mug,nNode(n)+1) = .true.
            lConnected(n,nNode(n)+1) = .true.
          endif
        enddo
      enddo

C Set flag as to whether a Node needs to be communicated with
      call re_alloc( lNodeNeeded, 1, Nodes, 'lNodeNeeded', 'globalise' )
      lNodeNeeded(1:Nodes) = .false.
      do n = 1,Nodes
        do mu = 1,no_u
          if (lConnected(mu,n)) lNodeNeeded(n) = .true.
        enddo
      enddo

C Set logical for current node to false since communication is not needed
      lNodeNeeded(Node+1) = .false.

C Count number of nodes that are neighbours
      call re_alloc(Node2NNodePtr,1,Nodes,name='Node2NNodePtr')
      NNodes = 0
      Node2NNodePtr(1:Nodes) = 0
      do n = 1,Nodes
        if (lNodeNeeded(n)) then
          NNodes = NNodes + 1
          Node2NNodePtr(n) = NNodes
        endif
      enddo
      call re_alloc(NNodeNo,1,NNodes,name='NNodeNo')

C Store node numbers
      NNodes = 0
      do n = 1,Nodes
        if (lNodeNeeded(n)) then
          NNodes = NNodes + 1
          NNodeNo(NNodes) = n - 1
        endif
      enddo
      call de_alloc(lNodeNeeded, 'lNodeNeeded', 'globalise' )

C Size data counter arrays
      call re_alloc(ndataBR,1,NNodes,name='ndataBR')
      call re_alloc(ndataBS,1,NNodes,name='ndataBS')
      call re_alloc(ndataBR2,1,NNodes,name='ndataBR2')
      call re_alloc(NNodeInBux,1,NNodes,name='NNodeInBux')
      call re_alloc(NNodeInBux2,1,NNodes,name='NNodeInBux2')
      call re_alloc(NNodeInBuxPreN,1,NNodes,name='NNodeInBuxPreN')
      call re_alloc(NNodeInBux2PreN,1,NNodes,name='NNodeInBux2PreN')

C----------------------------------------------
C Find elements of C/grad to be communicated  -
C----------------------------------------------
C Allocate array for local flags as to columns needed
      call re_alloc( lConnectedL, 1, no_u, 'lConnectedL', 'globalise')
      lConnectedL(1:no_u) = .false.

C Set local flags for C columns held locally
      do mu = 1,no_cl
        mug = ncL2G(mu)
        lConnectedL(mug) = .true.
      enddo

C Allocate communication arrays
      call re_alloc( Request, 1, 2*NNodes, 'Request', 'globalise')
      call re_alloc( Status, 1, MPI_Status_Size, 1, 2*NNodes, 'Request',
     &               'globalise')

C Communicate flags to neighbouring nodes
      nRequest = 0
      do n = 1,NNodes

C Post receives
        nRequest = nRequest + 1
        call MPI_IRecv(lConnected(1,n),no_u,MPI_logical,NNodeNo(n),
     .    NNodeNo(n),MPI_Comm_World,Request(nRequest),MPIerror)
      
C Post sends
        nRequest = nRequest + 1
        call MPI_ISend(lConnectedL,no_u,MPI_logical,NNodeNo(n),
     .    Node,MPI_Comm_World,Request(nRequest),MPIerror)
      
      enddo
     
C Statistics
      nTransfer = nTransfer + 2*NNodes
      nTransferInteger = nTransferInteger + 2*no_u*NNodes
     
C Wait for everything to finish
      call MPI_WaitAll(nRequest,Request,Status,MPIerror)

C Allocate arrays that depend on number of nodes to be communicated with
      call re_alloc(NNodeBasisSendC,1,NNodes,name='NNodeBasisSendC')

C Zero counters
      NNodeBasisSendC(1:NNodes) = 0

C Find number of basis functions that need to be transfered & size arrays
      MaxBasisSendC = 0
      do n = 1,NNodes
        do mu = 1,no_u
          if (lConnected(mu,n).and.lConnectedL(mu)) then
            NNodeBasisSendC(n) = NNodeBasisSendC(n) + 1
          endif
        enddo
        MaxBasisSendC = max(MaxBasisSendC,NNodeBasisSendC(n))
      enddo
      call re_alloc(NNodeBasisSendListC,1,MaxBasisSendC,
     .              1,NNodes,name='NNodeBasisSendListC')

C Find basis functions that need to be transfered
      NNodeBasisSendC(1:NNodes) = 0
      do n = 1,NNodes
        do mu = 1,no_u
          if (lConnected(mu,n).and.lConnectedL(mu)) then
            NNodeBasisSendC(n) = NNodeBasisSendC(n) + 1
            NNodeBasisSendListC(NNodeBasisSendC(n),n) = mu
          endif
        enddo
      enddo

C Free communication arrays
      call de_alloc( Status, 'Status', 'globalise' )
      call de_alloc( Request, 'Request', 'globalise' )
      call de_alloc( lConnectedL, 'lConnectedL', 'globalise' )
      call de_alloc( lConnected, 'lConnected', 'globalise' )

C Call timer
      call timer('setglobal',2)

      return

      end subroutine setglobalise
 
      subroutine setglobaliseB(nbands,Node)
C
C  Initialise information required for the globalisation of arrays
C  by passing to neighbouring nodes.
C
C  Julian Gale, NRI, Curtin University of Technology, April 2004
C

      implicit none

C Passed variables
      integer ::
     .  nbands,Node

C Internal variables
      integer                :: i
      integer                :: ig
      integer                :: il
      integer                :: jn
      integer                :: n
      integer                :: MPIerror
      integer                :: nRequest
      logical, pointer, save :: lBneededG(:,:) => null()
      logical, pointer, save :: lBneededL(:) => null()
      integer, pointer, save :: Request(:) => null()
      integer, pointer, save :: Status(:,:) => null()

C Call timer
      call timer('setglobalB',1)

C----------------------------------------------
C Find elements of bux1/2 to be communicated  -
C----------------------------------------------
C Allocate array for storing which bands are needed 
      call re_alloc( lBneededL, 1, nbands, 'lBneededL', 'globalise' )
      call re_alloc( lBneededG, 1, nbands, 1, NNodes, 'lBneededG',
     &               'globalise' )

C Allocate space for communication arrays
      call re_alloc( Request, 1, 2*NNodes, 'Request', 'globalise' )
      call re_alloc( Status, 1, MPI_Status_Size, 1, 2*NNodes, 'Status',
     &               'globalise' )

C Initialise local flag array as to whether band is needed
      lBneededL(1:nbands) = .false.

C Set local bands that are needed
      do il = 1,nbandsloc
        i = nbL2G(il)
        lBneededL(i) = .true.
      enddo

C Communicate flags to neighbouring nodes
      nRequest = 0
      do n = 1,NNodes

C Post receives
        nRequest = nRequest + 1
        call MPI_IRecv(lBneededG(1,n),nbands,MPI_logical,NNodeNo(n),
     .    NNodeNo(n),MPI_Comm_World,Request(nRequest),MPIerror)

C Post sends
        nRequest = nRequest + 1
        call MPI_ISend(lBneededL,nbands,MPI_logical,NNodeNo(n),
     .    Node,MPI_Comm_World,Request(nRequest),MPIerror)

      enddo

C Statistics
      nTransfer = nTransfer + 2*NNodes
      nTransferInteger = nTransferInteger + 2*nbands*NNodes

C Wait for everything to finish
      call MPI_WaitAll(nRequest,Request,Status,MPIerror)

C Initialise array holding number of bands to be sent/received
      call re_alloc(NNodeBandsSendB,1,NNodes,name='NNodeBandsSendB')
      NNodeBandsSendB(1:NNodes) = 0

      do il = 1,nbandsloc
        i = nbL2G(il)
        if (lBneededL(i)) then
          do n = 1,NNodes
            if (lBneededG(i,n)) then
              NNodeBandsSendB(n) = NNodeBandsSendB(n) + 1
            endif
          enddo
        endif
      enddo

C Make initial allocation of arrays for storing pointers to bux elements
      do n = 1,NNodes
        MaxBandsSendB = max(MaxBandsSendB,NNodeBandsSendB(n))
        NNodeBandsSendB(n) = 0
      enddo
      call re_alloc(NNodeBandsSendListB,1,MaxBandsSendB,1,NNodes,
     .              name='NNodeBandsSendListB')
      call re_alloc(NNodeBandsSendListB2,1,nbandsloc,1,NNodes,
     .              name='NNodeBandsSendListB2')

C Use flags to work out number of elements of bux that need to be sent locally
      NNodeBandsSendB(1:NNodes) = 0
      NNodeBandsSendListB2(1:nbandsloc,1:NNodes) = 0

      do il = 1,nbandsloc
        i = nbL2G(il)
        if (lBneededL(i)) then
          do n = 1,NNodes
            if (lBneededG(i,n)) then
              NNodeBandsSendB(n) = NNodeBandsSendB(n) + 1
              NNodeBandsSendListB(NNodeBandsSendB(n),n) = i
              NNodeBandsSendListB2(il,n) = NNodeBandsSendB(n)
            endif
          enddo
        endif
      enddo

C Allocate numhij arrays
      call re_alloc(numhijG,1,MaxBandsSendB,1,NNodes,name='numhijG')
      call re_alloc(numhijL,1,MaxBandsSendB,1,NNodes,name='numhijL')

C Initialise numhijL
      numhijL(1:MaxBandsSendB,1:NNodes) = 0

C Build array of numhij to be sent
      nRequest = 0
      do n = 1,NNodes

        do il = 1,nbandsloc
          ig = NNodeBandsSendListB2(il,n)
          if (ig.gt.0) then
            numhijL(ig,n) = numhij(il)
          endif
        enddo

C Post receives
        nRequest = nRequest + 1
        call MPI_IRecv(numhijG(1,n),NNodeBandsSendB(n),MPI_integer,
     .    NNodeNo(n),NNodeNo(n),MPI_Comm_World,Request(nRequest),
     .    MPIerror)

C Post sends
        nRequest = nRequest + 1
        call MPI_ISend(numhijL(1,n),NNodeBandsSendB(n),MPI_integer,
     .    NNodeNo(n),Node,MPI_Comm_World,Request(nRequest),MPIerror)

C Statistics
        nTransferInteger = nTransferInteger + NNodeBandsSendB(n)

      enddo

C Statistics
      nTransfer = nTransfer + 2*NNodes

C Wait for everything to finish
      call MPI_WaitAll(nRequest,Request,Status,MPIerror)

C Find total size of listhijG array
      MaxListHijG = 0
      do n = 1,NNodes
        NNodeInBux2PreN(n) = MaxListHijG
        do i = 1,NNodeBandsSendB(n)
          MaxListHijG = MaxListHijG + numhijG(i,n)
        enddo
      enddo

C Find total size of listhijL array
      MaxListHijL = 0
      do n = 1,NNodes
        NNodeInBuxPreN(n) = MaxListHijL
        do il = 1,nbandsloc
          ig = NNodeBandsSendListB2(il,n)
          if (ig.gt.0) then
            MaxListHijL = MaxListHijL + numhij(il)
          endif
        enddo
      enddo

C Allocate listhij arrays
      call re_alloc(listhijG,1,MaxListHijG,name='listhijG')
      call re_alloc(listhijL,1,MaxListHijL,name='listhijL')

C Build array of listhij to be sent
      nRequest = 0
      do n = 1,NNodes

        ndataBR(n) = 0
        ndataBS(n) = 0

        do il = 1,nbandsloc
          ig = NNodeBandsSendListB2(il,n)
          if (ig.gt.0) then
            do jn = 1,numhij(il)
              listhijL(NNodeInBuxPreN(n)+ndataBS(n)+jn) = listhij(jn,il)
            enddo
            ndataBS(n) = ndataBS(n) + numhij(il)
            ndataBR(n) = ndataBR(n) + numhijG(ig,n)
          endif
        enddo

C Post receives
        nRequest = nRequest + 1
        call MPI_IRecv(listhijG(NNodeInBux2PreN(n)+1),ndataBR(n),
     .    MPI_integer,NNodeNo(n),NNodeNo(n),MPI_Comm_World,
     .    Request(nRequest),MPIerror)

C Post sends
        nRequest = nRequest + 1
        call MPI_ISend(listhijL(NNodeInBuxPreN(n)+1),ndataBS(n),
     .    MPI_integer,NNodeNo(n),Node,MPI_Comm_World,Request(nRequest),
     .    MPIerror)

C Statistics
        nTransferInteger = nTransferInteger + ndataBS(n)

      enddo

C Statistics
      nTransfer = nTransfer + 2*NNodes

C Wait for everything to finish
      call MPI_WaitAll(nRequest,Request,Status,MPIerror)

C Free communication arrays
      call de_alloc( Status, 'Status', 'globalise' )
      call de_alloc( Request, 'Request', 'globalise' )

C Free workspace array
      call de_alloc( lBneededG, 'lBneededG', 'globalise' )
      call de_alloc( lBneededL, 'lBneededL', 'globalise' )

C Call timer
      call timer('setglobalB',2)

      return

      end subroutine setglobaliseB

C--------------------
C Globalise C/grad  -
C--------------------
      subroutine globaliseC(no_cl,ncmax,numc,c,nspin,Node)
C
C  Globalise the gradient array over adjacent nodes within spatial 
C  decomposition pattern.
C
C  Julian Gale, NRI, Curtin University of Technology, April 2004
C

      implicit none

C Passed variables
      integer
     .  ncmax,no_cl,nspin,Node

      integer
     .  numc(no_cl)

      real(dp) ::
     .  c(ncmax,no_cl,nspin)

C Internal variables
      integer                 :: i
      integer                 :: ik
      integer                 :: is
      integer                 :: mu
      integer                 :: mug
      integer                 :: n
      integer                 :: ndata
      integer                 :: MPIerror
      integer,  pointer, save :: Request(:,:) => null()
      integer,  pointer, save :: Status(:,:) => null()
      real(dp), pointer, save :: gtmp(:,:) => null()
      real(dp), pointer, save :: gtmp2(:,:) => null()

C Call timer
      call timer('globaliseC',1)

C
C Global reduction only involving neighbouring nodes
C

C Allocate space for communication arrays
      call re_alloc( Request, 1, 2, 1, NNodes, 'Request', 'globalise' )
      call re_alloc( Status, 1, MPI_Status_Size, 1, 2*NNodes, 'Status',
     &               'globalise' )
      call re_alloc( gtmp, 1, ncmax*MaxBasisSendC*nspin, 1, NNodes,
     &               'gtmp', 'globalise' )
      call re_alloc( gtmp2, 1, ncmax*MaxBasisSendC*nspin, 1, NNodes,
     &               'gtmp2', 'globalise' )

C Loop over neighbouring nodes
      do n = 1,NNodes

C Load data arrays
        ndata = 0
        do is = 1,nspin
          do i = 1,NNodeBasisSendC(n)
            mug = NNodeBasisSendListC(i,n)
            mu = ncG2L(mug)
            do ik = 1,numc(mu)
              ndata = ndata + 1
              gtmp(ndata,n) = c(ik,mu,is)
            enddo
          enddo
        enddo

C Post receives
        call MPI_IRecv(gtmp2(1,n),ndata,mpi_double_precision,NNodeNo(n),
     .    NNodeNo(n),MPI_Comm_World,Request(2,n),MPIerror)

C Post sends
        call MPI_ISend(gtmp(1,n),ndata,mpi_double_precision,NNodeNo(n),
     .    Node,MPI_Comm_World,Request(1,n),MPIerror)

C Statistics
        nTransferDouble = nTransferDouble + ndata

      enddo

C Statistics
      nTransfer = nTransfer + 2*NNodes

C Wait for everything to finish
      call MPI_WaitAll(2*NNodes,Request,Status,MPIerror)

C Loop over neighbouring nodes
      do n = 1,NNodes

C Load data arrays
        ndata = 0
        do is = 1,nspin
          do i = 1,NNodeBasisSendC(n)
            mug = NNodeBasisSendListC(i,n)
            mu = ncG2L(mug)
            do ik = 1,numc(mu)
              ndata = ndata + 1
              c(ik,mu,is) = c(ik,mu,is) + gtmp2(ndata,n)
            enddo
          enddo
        enddo
      enddo

C Deallocate space for communication arrays
      call de_alloc( gtmp2, 'gtmp2', 'globalise' )
      call de_alloc( gtmp, 'gtmp', 'globalise' )
      call de_alloc( Status, 'Status', 'globalise' )
      call de_alloc( Request, 'Request', 'globalise' )

C Call timer
      call timer('globaliseC',2)

      return

      end subroutine globaliseC

      subroutine globalisestats(Node)
C
C  Output statistics for globalisation process
C
C  Julian Gale, NRI, Curtin University of Technology, April 2004
C

C Passed variables
      integer
     .  Node

C Local variables
      integer          :: MPIerror
      integer          :: gtmp(3)
      integer          :: gtmp2(3)

C Globalise number of communications
      gtmp(1) = nTransfer
      gtmp(2) = nTransferInteger
      gtmp(3) = nTransferDouble
      call MPI_AllReduce(gtmp,gtmp2,3,MPI_integer,
     .                   MPI_sum,MPI_Comm_World,MPIerror)
      if (Node.eq.0) then
        write(6,'(/,''---------------------------------------'',
     .            ''---------------------------------------'')')
        write(6,'(''  Globalisation MPI statistics : '')')
        write(6,'(''---------------------------------------'',
     .            ''---------------------------------------'')')
        write(6,'(''  Total number of MPI communications = '',
     .            i20)') gtmp2(1)
        write(6,'(''  Number of integers communicated    = '',
     .            i20)') gtmp2(2)
        write(6,'(''  Number of doubles  communicated    = '',
     .            i20)') gtmp2(3)
        write(6,'(''---------------------------------------'',
     .            ''---------------------------------------'',/)')
      endif

      end subroutine globalisestats

C------------------------
C Globalise bux arrays  -
C------------------------
      subroutine globalinitB(narray)
C
C  Globalise the vectors from CtHC/CtSC over all nodes.
C  Initialisation of data.
C
C  Julian Gale, NRI, Curtin University of Technology, April 2004
C

C Passed variables
      integer
     .  narray

C Local variables
      integer
     .  n

C Copy pointers scaled by number of arrays
      do n = 1,NNodes
        NNodeInBux(n) = narray*NNodeInBuxPreN(n)
        NNodeInBux2(n) = narray*NNodeInBux2PreN(n)
      enddo

C Allocate space for storing array values to be communicated
      call re_alloc(buxstore,1,narray*MaxListHijL,
     .  name='buxstore',shrink=.false.)
      call re_alloc(buxstore2,1,narray*MaxListHijG,
     .  name='buxstore2',shrink=.false.)

C Initialise data counters
      ndataBS(1:NNodes) = 0
      ndataBR(1:NNodes) = 0
      ndataBR2(1:NNodes) = 0

      return

      end subroutine globalinitB

      subroutine globalloadB1(ibandloc,nbands,bux1)
C
C  Globalise the vectors from CtHC/CtSC over all nodes.
C  Load data into store prior to communication.
C
C  Julian Gale, NRI, Curtin University of Technology, April 2004
C

      implicit none

C Passed variables
      integer
     .  ibandloc,nbands
      real(dp) :: 
     .  bux1(nbands)

C Internal variables
      integer          :: i
      integer          :: iadd
      integer          :: ig
      integer          :: mu
      integer          :: n
      
C Loop over neighbouring nodes to distribute elements
      do n = 1,NNodes
      
C Set pointer to global band and NNodeSendB element
        ig = NNodeBandsSendListB2(ibandloc,n)
        if (ig.gt.0) then

C Load data arrays
          iadd = ndataBS(n) + NNodeInBux(n)
          do i = 1,numhij(ibandloc)
            mu = listhij(i,ibandloc)
            buxstore(iadd+i) = bux1(mu)
          enddo
          ndataBS(n) = ndataBS(n) + numhij(ibandloc)
          ndataBR(n) = ndataBR(n) + numhijG(ig,n)

        endif

      enddo
      
      return

      end subroutine globalloadB1

      subroutine globalloadB2(ibandloc,nbands,bux1,bux2)
C
C  Globalise the vectors from CtHC/CtSC over all nodes.
C  Load data into store prior to communication.
C
C  Julian Gale, NRI, Curtin University of Technology, April 2004
C

C Passed variables
      integer
     .  ibandloc,nbands
      real(dp) :: 
     .  bux1(nbands),bux2(nbands)

C Internal variables
      integer          :: i
      integer          :: iadd
      integer          :: ig
      integer          :: mu
      integer          :: n

C Loop over neighbouring nodes to distribute elements
      do n = 1,NNodes

C Set pointer to global band and NNodeSendB element
        ig = NNodeBandsSendListB2(ibandloc,n)
        if (ig.gt.0) then
  
C Load data arrays
          iadd = ndataBS(n) + NNodeInBux(n)
          do i = 1,numhij(ibandloc)
            mu = listhij(i,ibandloc)
            buxstore(iadd+i) = bux1(mu)
            buxstore(iadd+numhij(ibandloc)+i) = bux2(mu)
          enddo
          ndataBS(n) = ndataBS(n) + 2*numhij(ibandloc)
          ndataBR(n) = ndataBR(n) + 2*numhijG(ig,n)

        endif

      enddo

      return

      end subroutine globalloadB2

      subroutine globalcommB(Node)
C
C  Globalise the vectors from CtHC/CtSC over all nodes.
C  Communicate data for B matrix
C
C  Julian Gale, NRI, Curtin University of Technology, April 2004
C

C Passed variables
      integer
     .  Node

C Internal variables
      integer                 :: n
      integer                 :: MPIerror
      integer                 :: nRequest
      integer,  pointer, save :: Request(:) => null(),
     $                           Status(:,:) => null()

C Call timer
      call timer('globaliseB',1)

C Allocate space for communication arrays
      call re_alloc( Request, 1, 2*NNodes, 'Request', 'globalise' )
      call re_alloc( Status, 1, MPI_Status_Size, 1, 2*NNodes, 'Status',
     &               'globalise' )

C Loop over neighbouring nodes to distribute elements
      nRequest = 0
      do n = 1,NNodes

C Post receives
        if (ndataBR(n).gt.0) then
          nRequest = nRequest + 1
          call MPI_IRecv(buxstore2(NNodeInBux2(n)+1),ndataBR(n),
     .                   MPI_double_precision,NNodeNo(n),NNodeNo(n),
     .                   MPI_Comm_World,Request(nRequest),MPIerror)
        endif

C Post sends
        if (ndataBS(n).gt.0) then
          nRequest = nRequest + 1
          call MPI_ISend(buxstore(NNodeInBux(n)+1),ndataBS(n),
     .                   MPI_double_precision,NNodeNo(n),Node,
     .                   MPI_Comm_World,Request(nRequest),MPIerror)
        endif

C Statistics
        nTransferDouble = nTransferDouble + ndataBS(n)

      enddo

C Statistics
      nTransfer = nTransfer + 2*NNodes

C Wait for everything to finish
      call MPI_WaitAll(nRequest,Request,Status,MPIerror)

C Reinitialise data counters ready for unloading
      ndataBR(1:NNodes) = 0

C Deallocate space for communication arrays
      call de_alloc( Status, 'Status', 'globalise' )
      call de_alloc( Request, 'Request', 'globalise' )

C Call timer
      call timer('globaliseB',2)

      return

      end subroutine globalcommB

      subroutine globalreloadB1(ibandloc,nbands,bux1,buxg)
C
C  Globalise the vectors from CtHC/CtSC over all nodes.
C  Reload globalised data back out of arrays.
C
C  Julian Gale, NRI, Curtin University of Technology, April 2004
C     

C Passed variables
      integer
     .  ibandloc,nbands
      real(dp) :: 
     .  bux1(nbands),buxg(nbands)

C Internal variables
      integer          :: i
      integer          :: ig
      integer          :: j
      integer          :: jn
      integer          :: mu
      integer          :: n

C Initialise buxg
      do jn = 1,numhij(ibandloc)
        j = listhij(jn,ibandloc)
        buxg(j) = bux1(j)
      enddo

C Loop over neighbouring nodes
      do n = 1,NNodes

C Set pointer to global band and NNodeSendB element
        ig = NNodeBandsSendListB2(ibandloc,n)
        if (ig.gt.0) then

C Load globalised data arrays
          do i = 1,numhijG(ig,n)
            mu = listhijG(ndataBR2(n)+i+NNodeInBux2PreN(n))
            buxg(mu) = buxg(mu) + 
     .        buxstore2(ndataBR(n)+i+NNodeInBux2(n))
          enddo
          ndataBR(n) = ndataBR(n) + numhijG(ig,n)
          ndataBR2(n) = ndataBR2(n) + numhijG(ig,n)

        endif

      enddo

      return

      end subroutine globalreloadB1

      subroutine globalreloadB2(ibandloc,nbands,bux1,bux2,buxg)
C
C  Globalise the vectors from CtHC/CtSC over all nodes.
C  Reload globalised data back out of arrays.
C
C  Julian Gale, NRI, Curtin University of Technology, April 2004
C

C Passed variables
      integer
     .  ibandloc,nbands
      real(dp) :: 
     .  bux1(nbands),bux2(nbands),buxg(nbands,2)

C Internal variables
      integer          :: i
      integer          :: ig
      integer          :: j
      integer          :: jn
      integer          :: mu
      integer          :: n

C Initialise buxg
      do jn = 1,numhij(ibandloc)
        j = listhij(jn,ibandloc)
        buxg(j,1) = bux1(j)
        buxg(j,2) = bux2(j)
      enddo

C Loop over neighbouring nodes
      do n = 1,NNodes

C Set pointer to global band and NNodeSendB element
        ig = NNodeBandsSendListB2(ibandloc,n)
        if (ig.gt.0) then

C Load globalised data arrays
          do i = 1,numhijG(ig,n)
            mu = listhijG(ndataBR2(n)+i+NNodeInBux2PreN(n))
            buxg(mu,1) = buxg(mu,1) + 
     .        buxstore2(ndataBR(n)+i+NNodeInBux2(n))
            buxg(mu,2) = buxg(mu,2) + 
     .        buxstore2(ndataBR(n)+numhijG(ig,n)+i+NNodeInBux2(n))
          enddo
          ndataBR(n) = ndataBR(n) + 2*numhijG(ig,n)
          ndataBR2(n) = ndataBR2(n) + numhijG(ig,n)

        endif

      enddo

      return

      end subroutine globalreloadB2

      subroutine globalrezeroB1(ibandloc,nbands,buxg)
C
C  Re-zero globalised B vectors based on all values loaded
C
C  Julian Gale, NRI, Curtin University of Technology, May 2004
C

C Passed variables
      integer
     .  ibandloc,nbands
      real(dp) ::
     .  buxg(nbands)

C Internal variables
      integer          :: i
      integer          :: iadd
      integer          :: ig
      integer          :: j
      integer          :: jn
      integer          :: m
      integer          :: mu
      integer          :: n
          
C Rezero based on local elements
      do jn = 1,numhij(ibandloc)
        j = listhij(jn,ibandloc)
        buxg(j) = 0.0_dp
      enddo

C Loop over neighbouring nodes
      do n = 1,NNodes
            
C Set pointer to global band and NNodeSendB element
        ig = NNodeBandsSendListB2(ibandloc,n)
        if (ig.gt.0) then

C Find pointer to start of listhijG values
          m = ndataBR2(n) - numhijG(ig,n)
      
C Zero globalised data arrays
          iadd = m + NNodeInBux2PreN(n)
          do i = iadd+1,iadd+numhijG(ig,n)
            mu = listhijG(i)
            buxg(mu) = 0.0_dp
          enddo

        endif

      enddo

      return

      end subroutine globalrezeroB1

      subroutine globalrezeroB2(ibandloc,nbands,buxg)
C
C  Re-zero globalised B vectors based on all values loaded
C
C  Julian Gale, NRI, Curtin University of Technology, May 2004
C

C Passed variables
      integer
     .  ibandloc,nbands
      real(dp) :: 
     .  buxg(nbands,2)

C Internal variables
      integer          :: i
      integer          :: iadd
      integer          :: ig
      integer          :: j
      integer          :: jn
      integer          :: m
      integer          :: mu
      integer          :: n

C Rezero based on local elements
      do jn = 1,numhij(ibandloc)
        j = listhij(jn,ibandloc)
        buxg(j,1) = 0.0_dp
        buxg(j,2) = 0.0_dp
      enddo

C Loop over neighbouring nodes
      do n = 1,NNodes

C Set pointer to global band and NNodeSendB element
        ig = NNodeBandsSendListB2(ibandloc,n)
        if (ig.gt.0) then

C Find pointer to start of listhijG values
          m = ndataBR2(n) - numhijG(ig,n)

C Zero globalised data arrays
          iadd = m + NNodeInBux2PreN(n)
          do i = iadd+1,iadd+numhijG(ig,n)
            mu = listhijG(i)
            buxg(mu,1) = 0.0_dp
            buxg(mu,2) = 0.0_dp
          enddo

        endif

      enddo

      return

      end subroutine globalrezeroB2

      subroutine globalrezeroB3(ibandloc,nbands,buxg)
C
C  Re-zero globalised B vectors based on all values loaded
C
C  Julian Gale, NRI, Curtin University of Technology, May 2004
C

C Passed variables
      integer
     .  ibandloc,nbands
      real(dp) :: 
     .  buxg(nbands,3)

C Internal variables
      integer          :: i
      integer          :: iadd
      integer          :: ig
      integer          :: j
      integer          :: jn
      integer          :: m
      integer          :: mu
      integer          :: n

C Rezero based on local elements
      do jn = 1,numhij(ibandloc)
        j = listhij(jn,ibandloc)
        buxg(j,1) = 0.0_dp
        buxg(j,2) = 0.0_dp
        buxg(j,3) = 0.0_dp
      enddo

C Loop over neighbouring nodes
      do n = 1,NNodes

C Set pointer to global band and NNodeSendB element
        ig = NNodeBandsSendListB2(ibandloc,n)
        if (ig.gt.0) then

C Find pointer to start of listhijG values
          m = ndataBR2(n) - numhijG(ig,n)

C Zero globalised data arrays
          iadd = m + NNodeInBux2PreN(n)
          do i = iadd+1,iadd+numhijG(ig,n)
            mu = listhijG(i)
            buxg(mu,1) = 0.0_dp
            buxg(mu,2) = 0.0_dp
            buxg(mu,3) = 0.0_dp
          enddo

        endif

      enddo

      return

      end subroutine globalrezeroB3

C------------------------------------------------
C Globalise bux2/bux4/bux6  -  Blocked version  -
C------------------------------------------------
      subroutine globalloadB3(ibandloc,nbands,bux1,bux2,bux3)
C
C  Globalise the vectors from CtHC/CtSC over all nodes.
C  Load data into store prior to communication.
C
C  Julian Gale, NRI, Curtin University of Technology, April 2004
C

C Passed variables
      integer
     .  ibandloc,nbands
      real(dp) :: 
     .  bux1(nbands),bux2(nbands),bux3(nbands)

C Internal variables
      integer          :: i
      integer          :: iadd
      integer          :: ig
      integer          :: mu
      integer          :: n

C Loop over neighbouring nodes to distribute elements
      do n = 1,NNodes

C Set pointer to global band and NNodeSendB element
        ig = NNodeBandsSendListB2(ibandloc,n)
        if (ig.gt.0) then

C Load data arrays
          iadd = ndataBS(n) + NNodeInBux(n)
          do i = 1,numhij(ibandloc)
            mu = listhij(i,ibandloc)
            buxstore(iadd+i) = bux1(mu)
            buxstore(numhij(ibandloc)+iadd+i) = bux2(mu)
            buxstore(2*numhij(ibandloc)+iadd+i) = bux3(mu)
          enddo
          ndataBS(n) = ndataBS(n) + 3*numhij(ibandloc)
          ndataBR(n) = ndataBR(n) + 3*numhijG(ig,n)

        endif

      enddo

      return

      end subroutine globalloadB3

      subroutine globalreloadB3(ibandloc,nbands,bux1,bux2,bux3,buxg)
C
C  Globalise the vectors from CtHC/CtSC over all nodes.
C  Reload globalised data back out of arrays.
C
C  Julian Gale, NRI, Curtin University of Technology, April 2004
C


C Passed variables
      integer
     .  ibandloc,nbands
      real(dp) :: 
     .  bux1(nbands),bux2(nbands),bux3(nbands),buxg(nbands,3)

C Internal variables
      integer          :: i
      integer          :: ig
      integer          :: j
      integer          :: jn
      integer          :: mu
      integer          :: n

C Initialise buxg
      do jn = 1,numhij(ibandloc)
        j = listhij(jn,ibandloc)
        buxg(j,1) = bux1(j)
        buxg(j,2) = bux2(j)
        buxg(j,3) = bux3(j)
      enddo

C Loop over neighbouring nodes
      do n = 1,NNodes

C Set pointer to global band and NNodeSendB element
        ig = NNodeBandsSendListB2(ibandloc,n)
        if (ig.gt.0) then
            
C Load globalised data arrays
          do i = 1,numhijG(ig,n)
            mu = listhijG(ndataBR2(n)+i+NNodeInBux2PreN(n))
            buxg(mu,1) = buxg(mu,1) + 
     .        buxstore2(ndataBR(n)+i+NNodeInBux2(n))
            buxg(mu,2) = buxg(mu,2) + 
     .        buxstore2(ndataBR(n)+numhijG(ig,n)+i+NNodeInBux2(n))
            buxg(mu,3) = buxg(mu,3) + 
     .        buxstore2(ndataBR(n)+2*numhijG(ig,n)+i+NNodeInBux2(n))
          enddo
          ndataBR(n) = ndataBR(n) + 3*numhijG(ig,n)
          ndataBR2(n) = ndataBR2(n) + numhijG(ig,n)

        endif

      enddo

      return

      end subroutine globalreloadB3

C--------------------
C Globalise ft/fst  -
C--------------------
      subroutine globaliseF(no_cl,nbands,nftmax,numft,
     .                      listft,ft,ftG,Node)
C
C  Globalise the gradient array over adjacent nodes within spatial 
C  decomposition pattern.
C
C  Julian Gale, NRI, Curtin University of Technology, April 2004
C

C Passed variables
      integer
     .  nbands,nftmax,no_cl,Node

      integer
     .  listft(nftmax,no_cl),numft(no_cl)

      real(dp) ::
     .  ft(nftmax,no_cl),ftG(MaxFt2,no_cl)

C Internal variables
      integer                 :: i
      integer                 :: ik
      integer                 :: mu
      integer                 :: mug
      integer                 :: n
      integer                 :: MPIerror
      integer,  pointer, save :: Request(:,:) => null()
      integer,  pointer, save :: Status(:,:) => null()
      real(dp), pointer, save :: gtmp(:,:) => null()
      real(dp), pointer, save :: gtmp2(:,:) => null()
      real(dp), pointer, save :: gvec(:) => null()

C Call timer
      call timer('globaliseF',1)

C
C Global reduction only involving neighbouring nodes
C

C Allocate space for communication arrays
      call re_alloc( Request, 1, 2, 1, NNodes, 'Request', 'globalise' )
      call re_alloc( Status, 1, MPI_Status_Size, 1, 2*NNodes, 'Status',
     &               'globalise' )
      call re_alloc( gtmp, 1, MaxFtL, 1, NNodes, 'gtmp', 'globalise' )
      call re_alloc( gtmp2, 1, MaxFtG, 1, NNodes, 'gtmp2', 'globalise' )
      call re_alloc( gvec, 1, nbands, 'gvec', 'globalise' )

      ftG(1:MaxFt2,1:no_cl) = 0.0d0
      gvec(1:nbands) = 0.0d0
      do mu = 1,no_cl
        do ik = 1,numft(mu)
          gvec(listft(ik,mu)) = ft(ik,mu)
        enddo
        do ik = 1,numft2(mu)
          ftG(ik,mu) = gvec(listft2(ik,mu))
        enddo
        do ik = 1,numft(mu)
          gvec(listft(ik,mu)) = 0.0d0
        enddo
      enddo

C Loop over neighbouring nodes
      ndataFR(1:NNodes) = 0
      ndataFS(1:NNodes) = 0
      do n = 1,NNodes

C Load data arrays
        do i = 1,NNodeBasisSendC(n)
          mug = NNodeBasisSendListC(i,n)
          mu = ncG2L(mug)
          do ik = 1,numft(mu)
            gtmp(ndataFS(n)+ik,n) = ft(ik,mu)
          enddo
          ndataFS(n) = ndataFS(n) + numft(mu)
          ndataFR(n) = ndataFR(n) + numftG(i,n)
        enddo

C Post receives
        call MPI_IRecv(gtmp2(1,n),ndataFR(n),mpi_double_precision,
     $       NNodeNo(n),NNodeNo(n),MPI_Comm_World,Request(2,n),MPIerror)
C Post sends
        call MPI_ISend(gtmp(1,n),ndataFS(n),mpi_double_precision,
     $       NNodeNo(n),Node,MPI_Comm_World,Request(1,n),MPIerror)

C Statistics
        nTransferDouble = nTransferDouble + ndataFR(n) + ndataFS(n)

      enddo

C Statistics
      nTransfer = nTransfer + 2*NNodes

C Wait for everything to finish
      call MPI_WaitAll(2*NNodes,Request,Status,MPIerror)

C Loop over neighbouring nodes
      ndataFR(1:NNodes) = 0
      do n = 1,NNodes

C Load data arrays
        do i = 1,NNodeBasisSendC(n)
          mug = NNodeBasisSendListC(i,n)
          mu = ncG2L(mug)
          do ik = 1,numftG(i,n)
            gvec(listftG(ndataFR(n)+ik,n)) = gtmp2(ndataFR(n)+ik,n)
          enddo
          do ik = 1,numft2(mu)
            ftG(ik,mu) = ftG(ik,mu) + gvec(listft2(ik,mu))
          enddo
          do ik = 1,numftG(i,n)
            gvec(listftG(ndataFR(n)+ik,n)) = 0.0d0
          enddo
          ndataFR(n) = ndataFR(n) + numftG(i,n)
        enddo
      enddo

C Deallocate space for communication arrays
      call de_alloc( gvec, 'gvec', 'globalise' )
      call de_alloc( gtmp2, 'gtmp2', 'globalise' )
      call de_alloc( gtmp, 'gtmp', 'globalise' )
      call de_alloc( Status, 'Status', 'globalise' )
      call de_alloc( Request, 'Request', 'globalise' )

C Call timer
      call timer('globaliseF',2)

      return

      end subroutine globaliseF

      subroutine setglobaliseF(no_u,nbands,no_cl,maxnft,Node)
C
C  Initialise information required for the globalisation of arrays
C  relating to indexing of Ft by passing to neighbouring nodes.
C
C  Julian Gale, NRI, Curtin University of Technology, May 2004
C

C Passed variables
      integer ::
     .  no_u,nbands,no_cl,maxnft,Node

C Internal variables
      integer                :: i
      integer                :: m
      integer                :: MPIerror
      integer                :: mu
      integer                :: mug
      integer                :: muk
      integer                :: n
      integer                :: nRequest
      integer                :: nsizeL
      integer                :: nsizeG
      integer, pointer, save :: Request(:) => null(),
     $                          Status(:,:) => null()
      logical, pointer, save :: lNeeded(:) => null()

C Call timer
      call timer('setglobalF',1)

C----------------------------------------------
C Find elements of Ft/Fst to be communicated  -
C----------------------------------------------
C Size data arrays
      call re_alloc(ndataFR,1,NNodes,name='ndataFR')
      call re_alloc(ndataFS,1,NNodes,name='ndataFS')
      call re_alloc(numftL,1,MaxBasisSendC,1,NNodes,name='numftL')
      call re_alloc(numftG,1,MaxBasisSendC,1,NNodes,name='numftG')

C Allocate array for local flags as to columns needed
      call re_alloc( lNeeded, 1, no_u, 'Request', 'globalise' )

C Allocate communication arrays
      call re_alloc( Request, 1, 2*NNodes, 'Request', 'globalise' )
      call re_alloc( Status, 1, MPI_Status_Size, 1, 2*NNodes, 'Status',
     &               'globalise' )

C Loop over neighbouring nodes
      nRequest = 0
      do n = 1,NNodes

C Load data arrays
        do i = 1,NNodeBasisSendC(n)
          mug = NNodeBasisSendListC(i,n)
          mu = ncG2L(mug)
          numftL(i,n) = numft(mu)
        enddo

C Post receives
        nRequest = nRequest + 1
        call MPI_IRecv(numftG(1,n),NNodeBasisSendC(n),MPI_integer,
     .    NNodeNo(n),NNodeNo(n),MPI_Comm_World,Request(nRequest),
     .    MPIerror)

C Post sends
        nRequest = nRequest + 1
        call MPI_ISend(numftL(1,n),NNodeBasisSendC(n),MPI_integer,
     .    NNodeNo(n),Node,MPI_Comm_World,Request(nRequest),MPIerror)

C Statistics
        nTransferDouble = nTransferDouble + NNodeBasisSendC(n)

      enddo

C Statistics
      nTransfer = nTransfer + nRequest

C Wait for everything to finish
      call MPI_WaitAll(nRequest,Request,Status,MPIerror)

C Find maximum sizes of listft arrays to be transfered
      MaxFtL = 0
      MaxFtG = 0
      do n = 1,NNodes
        nsizeL = 0
        nsizeG = 0
        do i = 1,NNodeBasisSendC(n)
          nsizeL = nsizeL + numftL(i,n)
          nsizeG = nsizeG + numftG(i,n)
        enddo
        MaxFtL = max(MaxFtL,nsizeL)
        MaxFtG = max(MaxFtG,nsizeG)
      enddo

C Size data arrays
      call re_alloc(listftL,1,MaxFtL,1,NNodes,name='listftL')
      call re_alloc(listftG,1,MaxFtG,1,NNodes,name='listftG')

C Zero data counters
      ndataFR(1:NNodes) = 0
      ndataFS(1:NNodes) = 0

C Loop over neighbouring nodes
      nRequest = 0
      do n = 1,NNodes

C Load data arrays
        do i = 1,NNodeBasisSendC(n)
          mug = NNodeBasisSendListC(i,n)
          mu = ncG2L(mug)
          do muk = 1,numftL(i,n)
            listftL(ndataFS(n)+muk,n) = listft(muk,mu)
          enddo
          ndataFS(n) = ndataFS(n) + numftL(i,n)
          ndataFR(n) = ndataFR(n) + numftG(i,n)
        enddo

C Post receives
        nRequest = nRequest + 1
        call MPI_IRecv(listftG(1,n),ndataFR(n),MPI_integer,
     .    NNodeNo(n),NNodeNo(n),MPI_Comm_World,Request(nRequest),
     .    MPIerror)

C Post sends
        nRequest = nRequest + 1
        call MPI_ISend(listftL(1,n),ndataFS(n),MPI_integer,
     .    NNodeNo(n),Node,MPI_Comm_World,Request(nRequest),MPIerror)

C Statistics
        nTransferDouble = nTransferDouble + ndataFS(n)

      enddo

C Statistics
      nTransfer = nTransfer + nRequest

C Wait for everything to finish
      call MPI_WaitAll(nRequest,Request,Status,MPIerror)

C Create copy of numft/listft
      MaxFt2 = maxnft
      call re_alloc(numft2,1,no_cl,name='numft2')
      call re_alloc(listft2,1,MaxFt2,1,no_cl,name='listft2')
      do mu = 1,no_cl
        numft2(mu) = numft(mu)
        do muk = 1,numft(mu)
          listft2(muk,mu) = listft(muk,mu)
        enddo
      enddo

C Find true list of Ft elements
      ndataFR(1:NNodes) = 0
      lNeeded(1:no_u) = .false.
      do n = 1,NNodes
        do i = 1,NNodeBasisSendC(n)
          mug = NNodeBasisSendListC(i,n)
          mu = ncG2L(mug)
          do muk = 1,numft2(mu)
            lNeeded(listft2(muk,mu)) = .true.
          enddo
          do muk = 1,numftG(i,n)
            lNeeded(listftG(ndataFR(n)+muk,n)) = .true.
          enddo
          ndataFR(n) = ndataFR(n) + numftG(i,n)
          numft2(mu) = 0
          do m = 1,nbands
            if (lNeeded(m)) then
              numft2(mu) = numft2(mu) + 1
            endif
          enddo
          if (numft2(mu).gt.MaxFt2) then
            MaxFt2 = numft2(mu) + 20
            call re_alloc(listft2,1,MaxFt2,1,no_cl,name='listft2')
          endif
          numft2(mu) = 0
          do m = 1,nbands
            if (lNeeded(m)) then
              numft2(mu) = numft2(mu) + 1
              listft2(numft2(mu),mu) = m
              lNeeded(m) = .false.
            endif
          enddo
        enddo
      enddo

C Free communication arrays
      call de_alloc( Status, 'Status', 'globalise' )
      call de_alloc( Request, 'Request', 'globalise' )
      call de_alloc( lNeeded, 'lNeeded', 'globalise' )

C Call timer
      call timer('setglobalF',2)

      return

      end subroutine setglobaliseF

#endif
      end module globalise
