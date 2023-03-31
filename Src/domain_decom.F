! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

C This file contains module domain_decom.
C Creates a new distribution for the Hamiltonian and overlap
C matrices using domain decomposition. This data distribution is
C used in iterative methods, because it reduces the number of
C communications needed in a matrix per vector multiplication
C
C Written by Rogeli Grima (BSC) Oct.2008
C
C The following subroutines are public:
C   domainDecom        = Creates a Domain Decomposition of the
C                        matrix structure.
C   gatherGlobalMatrix = Process 0 gather the whole matrix.
C                        It is not mandatory for SIESTA (Only for debug)
C   preSetOrbitLimits  = Init the limits of the original distribution
C   resetDomainDecom   = Reset internal arrays of the module
C
C Includes the following internal variables:
C   logical use_dd      :: .true. if we want to use domain decomposition
C   logical use_dd_perm :: .true. if the domain decomposition permutation
C                          has been created.
C   integer llimit      :: Lower limit of the original distribution
C   integer ulimit      :: Upper limit of the original distribution
C   integer dd_nuo      :: Number of rows of current distribution
C   integer dd_ncolum   :: Number of columns of current distribution
C   integer dd_perm     :: Transforms a global index into a local index
C   integer dd_invp     :: Transforms a local index into a global index
C   integer dd_cperm    :: Transforms a global column index into a local index.
C   integer dd_nnode    :: Process number of every orbital.
C   integer dd_ncom     :: Number of communications
C   integer dd_comm     :: Communications (neighbour processes)
C   integer dd_bsiz     :: Size of communications (boundary size)
C   integer dd_lni      :: Local internal nodes. Orbitals that only appear
C                          in the current process.
C   integer dd_lnb      :: Local boundary nodes. Orbitals that belongs to
C                          current process but needed by other processes.
C   
      module domain_decom
      use alloc,         only : re_alloc, de_alloc
      use parallel,      only : Node, Nodes
      use printMat
      use scheComm
      implicit none

      PUBLIC :: domainDecom, gatherGlobalMatrix,
     &          resetDomainDecom, preSetOrbitLimits
!      gperm

      PRIVATE
      logical, public           :: use_dd
      logical, public           :: use_dd_perm
      integer, public           :: llimit
      integer, public           :: ulimit
      integer, public           :: dd_nuo
      integer, public           :: dd_ncolum
      integer, public           :: dd_ncom
      integer, public , pointer :: dd_perm(:)
      integer, public , pointer :: dd_invp(:)
      integer, public , pointer :: dd_cperm(:)
      integer, public , pointer :: dd_nnode(:)
      integer, public , pointer :: dd_comm(:)
      integer, public , pointer :: dd_bsiz(:)
      integer, public           :: dd_lni
      integer, public           :: dd_lnb
#ifdef _SLEPC_
      integer, public , pointer :: slepc_perm(:)
#endif

      integer          :: gni          ! Global Node internal
      integer          :: gnb          ! Global Node boundary
      integer, pointer :: listhPtrTot(:)
      integer, pointer :: listhTot(:)
      integer, pointer :: partition(:)
      integer, pointer :: ginde(:,:)
      integer, pointer :: gperm(:)
      integer, pointer :: ginvp(:)
      integer, pointer :: nnzInt(:)
      integer, pointer :: nnzBou(:)
      integer, pointer :: nNeigh(:)
      logical, pointer :: DomNeigh(:)
      type(COMM_T)     :: Gcomm

      contains

      subroutine domainDecom( nuotot, nuo, nnz )
C ==================================================================
C Redistribute a matrix using a domain decomposition. This method
C reduce the number and size of boundaries among the several domains.
C An scheduling of communications among processes is computed.
C ==================================================================
C SUBROUTINE domainDecom( nuotot, nuo, nnz )
C
C INPUT:
C integer  nuotot   : Total number of orbitals
C integer  nuo      : Local number of orbitals
C integer  nnz      : Total number of non-null elements in the
C                     matrix
C
C OUTPUT: matrix H structure is accesed via module.
C integer  numh     : number of non-null elements per row
C integer  listh    : column index of the non-null elements of the
C                     matrix.
C integer  listhptr : Start of every row of listh
C
C BEHAVIOR:
C Process 0 gather all the matrix structure from the other processes.
C This creates a data partition for each process, computes the
C boundary nodes and creates an scheduling of communications.
C Finally, process 0 send the new data distribution and the new
C scheduling to all the other processes.
C
C ==================================================================
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
C     Input variables
      integer          :: nuotot, nuo, nnz
      integer          :: MPIerror
      logical,    save :: firsttime = .true.

      if (firsttime) then
        firsttime = .false.
      else
        call resetDomainDecom( )
      endif

      call reduceGlobalMatrix( nuotot, nuo, nnz )
      if (node.eq.0) then
C       Set orbitals to partitions
        call makePartitions( nuotot, nodes )

        call boundaryNodes( nuotot )

        call checkPartitionSize( nuotot )

        call createPerm( nuotot )

        call compuNeigh( )

        call compuComm( )

      endif

      call distriGlobalMatrix( nuotot, nuo, nnz )

      if (node.eq.0) then
!       Free local memory
        if (associated(Gcomm%ind))
     &    call de_alloc( Gcomm%ind, 'comm%ind', 'scheComm' )
        call de_alloc( DomNeigh, 'DomNeigh', 'domainDecom' )
        call de_alloc( nNeigh, 'nNeigh', 'domainDecom' )
!        write(23,*) '  ******* ginvp & ginde & ginvp & ginde *******'
!        write(23,*) '  *******        Not deallocated        *******'
!        flush(23)
!        call de_alloc( nnzBou, 'nnzBou', 'domainDecom' )
!        call de_alloc( nnzInt, 'nnzInt', 'domainDecom' )
!        call de_alloc( ginvp, 'ginvp', 'domainDecom' )
!        call de_alloc( ginde, 'ginde', 'domainDecom' )
        call de_alloc( partition, 'partition', 'domainDecom' )
        call de_alloc( listhTot, 'listhTot', 'domainDecom' )
        call de_alloc( listhPtrTot, 'listhPtrTot', 'domainDecom' )

        call de_alloc( gperm, 'gperm', 'domainDecom' )

      endif

      end subroutine domainDecom

      subroutine reduceGlobalMatrix( nuotot, nuo, nnz )
C ==================================================================
C Process 0 gathers the whole matrix structure from all other
C processes
C ==================================================================
C SUBROUTINE reduceGlobalMatrix( nuotot, nuo, nnz )
C
C INPUT:
C integer  nuotot   : Total number of orbitals
C integer  nuo      : Local number of orbitals
C integer  nnz      : Total number of non-null elements in the
C                     matrix
C matrix H structure is accesed via module.
C integer  numh     : number of non-null elements per row
C integer  listh    : column index of the non-null elements of the
C                     matrix.
C
C OUTPUT: Whole H matrix structure is stored in current module.
C integer  numhTot     : number of non-null elements per row
C integer  listhTot    : column index of the non-null elements of
C                        the matrix.
C integer  listhPtrTot : Start of every row of listh
C
C BEHAVIOR:
C Process 0 computes the number of rows that should receive from
C every process.
C Process 0 receive numh from all the processes and compute the
C number of non-zero elements in other processes.
C Process 0 receive the non-zero elements 
C
C ==================================================================
      use sparse_matrices, only : numh, listh
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
C     Input variables
      integer,           intent(in) :: nuotot, nuo, nnz
C     Local variables
      integer                       :: blocks, remain, PP, io, ind, ii,
     &                                 nnzG, nnzL
      integer,              pointer :: nuos(:)
#ifdef MPI
      integer                       :: MPIerror, Status(MPI_Status_Size)
#endif

      if (node.eq.0) then
C       Compute the number of rows that has every process
        nullify(nuos)
        call re_alloc( nuos, 0, Nodes-1, 'nuos', 'domainDecom' )
        blocks = nuotot/Nodes
        remain = nuotot - blocks*Nodes
        do PP= 0, Nodes-1
          if (PP.lt.remain) then
            nuos(PP) = blocks + 1
          else
            nuos(PP) = blocks
          endif
        enddo

C       Receive the number of non-null elements of every row
        nullify(listhPtrTot)
        call re_alloc( listhPtrTot, 0, nuotot, 'listhPtrTot',
     &                 'domainDecom' )

        do ii= 1, nuo
          listhPtrTot(ii) = numh(ii)
        enddo
        ind = nuo + 1
#ifdef MPI
        do PP= 1, Nodes-1
          call MPI_recv( listhPtrTot(ind), nuos(PP), MPI_INTEGER,
     &                   PP, 1, MPI_Comm_world, Status, MPIerror )
          ind = ind + nuos(PP)
        enddo
#endif
        listhPtrTot(0) = 1
        do io= 1, nuotot
          listhPtrTot(io) = listhPtrTot(io) + listhPtrTot(io-1)
        enddo
        nnzG = listhPtrTot(nuotot) - 1

        nullify(listhTot)
        call re_alloc( listhTot, 1, nnzG, 'listhTot',
     &                 'domainDecom' )

C       Receive the non-null elements
        do ii= 1, nnz
          listhTot(ii) = listh(ii)
        enddo
        ind = nuo
#ifdef MPI
        do PP= 1, Nodes-1
          nnzL = listhPtrTot(ind+nuos(PP)) - listhPtrTot(ind)
          ii = listhPtrTot(ind)
          call MPI_recv( listhTot(ii), nnzL, MPI_INTEGER,
     &                   PP, 1, MPI_Comm_world, Status, MPIerror )
          ind = ind + nuos(PP)
        enddo
#endif
!        call printMatrix( nuotot, nuotot, listhTot, listhPtrTot,
!     &                    'globmat1.ps', 'Global_matrix' )
        call de_alloc( nuos, 'nuos', 'domainDecom' )
      else
#ifdef MPI
C       Send the number of non-null elements of every row
        call MPI_send( numh, nuo, MPI_INTEGER, 0, 1,
     &                 MPI_Comm_world, MPIerror )
C       Send the non-null elements
        call MPI_send( listh, nnz, MPI_INTEGER, 0, 1,
     &                 MPI_Comm_world, MPIerror )
#endif
      endif
      end subroutine reduceGlobalMatrix


      subroutine makePartitions( no, nparts )
C ==================================================================
C Create a domain decomposition of a matrix graph.
C ==================================================================
C SUBROUTINE makePartitions( no, nparts )
C
C INPUT:
C integer  no        : Number of nodes of the matrix graph
C integer  nparts    : Number of domains
C
C OUTPUT: 
C integer  partition : Domain of every node of the graph
C
C BEHAVIOR:
C Creates a domain decomposition of a graph. Every node is assigned
C to a domain.
C
C ==================================================================
      use sys,           only : die
      implicit none
!     Input/Output variables
      integer,           intent(in) :: no, nparts
!     Local variables
      integer                       :: ii, jj, ind, nnz, dummy(1),
     &                                 metis_options(1), edgecut
      integer,              pointer :: iia(:), jja(:)
      logical                       :: notfound

      nullify(partition)
      call re_alloc( partition, 1, no, 'partition', 'domainDecom' )

      if (nodes.eq.1) then
        do ii= 1, no
          partition(ii) = 1
        enddo
      else
!       Drop diagonal nodes (autoconnection)
        nnz = listhPtrTot(no) - listhPtrTot(0) - no
        nullify(jja,iia)
        call re_alloc( jja, 0, no, 'jja', 'domainDecom' )
        call re_alloc( iia, 1, nnz, 'iia', 'domainDecom' )

        jja(0) = 1
        ind    = 1
        do jj= 1, no
          notfound = .true.
          DO ii= listhPtrTot(jj-1), listhPtrTot(jj)-1
            if (listhTot(ii).eq.jj) then
              notfound = .false.
            else
              iia(ind) = listhTot(ii)
              ind = ind + 1
            endif
          ENDDO
          if (notfound) call die( 'not diagonal element in this row' )
          jja(jj) = ind
        enddo

!        call printMatrix( no, no, iia, jja,
!       &                  'globmat2.ps', 'Global_matrix' )

!       Use metis library to create a Domain decomposition
        metis_options(1) = 0
#if defined(ON_DOMAIN_DECOMP) || defined(SIESTA__METIS)
        if (nparts .gt. 8 ) then
          call METIS_PartGraphKway( no, jja, iia, dummy, dummy, 0,
     &                              1, nparts, metis_options,
     &                              edgecut, partition )
        else
          call METIS_PartGraphRecursive( no, jja, iia, dummy, dummy, 0,
     &                                   1, nparts, metis_options,
     &                                   edgecut, partition )
        endif
#else
        ! This statement should not be reached if we caught
        ! the initial attempt to request domain decomposition without
        ! the proper symbol (in siesta_options)
        !
        call die("Need to compile with ON_DOMAIN_DECOMP support")
#endif

        call de_alloc( jja, 'jja', 'domainDecom' )
        call de_alloc( iia, 'iia', 'domainDecom' )
      endif
      end subroutine makePartitions


      subroutine boundaryNodes( no )
C ==================================================================
C Mark nodes with conections with different domains as boundary nodes.
C ==================================================================
C SUBROUTINE boundaryNodes( no )
C
C INPUT:
C integer  no        : Number of nodes of the matrix graph
C
C OUTPUT: 
C integer  partition : Domain of every node. Positive if it is an
C                      internal node and Negative for boundary nodes
C integer  ginde     : Contains the number of internal and boundary
C                      nodes for every process.
C integer  gni       : Total number of internal nodes.
C integer  gnb       : Total number of boundary nodes.
C
C BEHAVIOR:
C For every node, check if it is a boundary node. Update the number
C of internal or boundary nodes.
C
C ==================================================================
      implicit none
!     Input/Output variables
      integer,           intent(in) :: no
!     Local variables
      logical                       :: notfound
      integer                       :: ii, jj

      nullify(ginde)
      call re_alloc( ginde, 1, 4, 1, Nodes, 'ginde', 'domainDecom' )
      ginde = 0

      do jj = 1, no
!       look if all niegh. are in the same partition
        notfound = .true.
        ii = listhPtrTot(jj-1)
        do while (ii.lt.listhPtrTot(jj) .and. notfound)
          if (partition(jj).ne.abs(partition(listhTot(ii))))
     &      notfound = .false.
          ii = ii + 1
        enddo
        if ( .not.notfound ) then
          ginde(4,partition(jj)) = ginde(4,partition(jj)) + 1
          partition(jj) = -partition(jj)
        else
          ginde(2,partition(jj)) = ginde(2,partition(jj)) + 1
        endif
      enddo

      ginde(1,1) = 1
      do jj= 2, Nodes
        ginde(1,jj) = ginde(2,jj-1) + ginde(1,jj-1)
      enddo
      gni = ginde(1,Nodes) + ginde(2,Nodes) - 1

      ginde(3,1) = gni + 1
      do jj= 2, Nodes
        ginde(3,jj) = ginde(4,jj-1) + ginde(3,jj-1)
      enddo
      gnb = no - gni
      end subroutine boundaryNodes 


      subroutine checkPartitionSize( no )
C ==================================================================
C Compute the number of nonzeros of every domain for internal and
C boundary nodes.
C ==================================================================
C SUBROUTINE checkPartitionSize( no )
C
C INPUT:
C integer  no     : Number of nodes of the matrix graph
C
C OUTPUT: 
C integer  nnzInt : number of non-null elements for internal nodes
C integer  nnzBou : number of non-null elements for boundary nodes
C
C BEHAVIOR:
C For every node, check if it is an internal or a boundary node and
C count the number of nonzero elements.
C
C ==================================================================
      implicit none
!     Input/Output variables
      integer,           intent(in) :: no
!     Local variables
      integer                       :: io, PP

      nullify(nnzInt,nnzBou)
      call re_alloc( nnzInt, 1, Nodes, 'nnzInt', 'domainDecom' )
      call re_alloc( nnzBou, 1, Nodes, 'nnzBou', 'domainDecom' )
      nnzInt = 0
      nnzBou = 0
      do io= 1, no
        PP = partition(io)
        if (PP.gt.0) then
          nnzInt(PP) = nnzInt(PP)+listhPtrTot(io)-listhPtrTot(io-1)
        else
          PP = -PP
          nnzBou(PP) = nnzBou(PP)+listhPtrTot(io)-listhPtrTot(io-1)
        endif
      enddo
      end subroutine checkPartitionSize


      subroutine createPerm( no )
      implicit none
!     Input/Output variables
      integer,           intent(in) :: no
!     Local variables
      integer               :: dom, maxni, maxnb, maxnz, noloc, boloc,
     &                         metisOpt(1), ii, jj, kk, offsetI,
     &                         offsetB
      integer,      pointer :: adj(:), xadj(:), perI(:), invI(:),
     &                         perB(:), invB(:), perR(:), invR(:)

      nullify(gperm,ginvp)
      call re_alloc( gperm, 1, no, 'gperm', 'domainDecom' )
      call re_alloc( ginvp, 1, no, 'ginvp', 'domainDecom' )

      maxni = 1
      do dom= 1, Nodes
        maxni = max(maxni,ginde(2,dom))
      enddo
      maxnb = 1
      do dom= 1, Nodes
        maxnb = max(maxnb,ginde(4,dom))
      enddo
      maxnz = 1
      do dom= 1, Nodes
        maxnz = max(maxnz,nnzInt(dom))
      enddo

      nullify(adj,xadj,perI,invI,perR,invR,perB,invB)
      call re_alloc(  adj, 1, maxnz,  'adj', 'domainDecom' )
      call re_alloc( xadj, 0, maxni, 'xadj', 'domainDecom' )
      call re_alloc( perI, 1,    no, 'perI', 'domainDecom' )
      call re_alloc( invI, 1, maxni, 'invI', 'domainDecom' )
      call re_alloc( perR, 1, maxni, 'perR', 'domainDecom' )
      call re_alloc( invR, 1, maxni, 'invR', 'domainDecom' )
      call re_alloc( perB, 1,    no, 'perB', 'domainDecom' )
      call re_alloc( invB, 1, maxnb, 'invB', 'domainDecom' )

      metisOpt = 0
      offsetI  = 0
      offsetB  = gni

      do dom= 1, Nodes
!       Loop over each subdomain: separate interior and boundary nodes
        call createSubgraph( no, listhTot, listhPtrTot, dom,
     &                       adj, xadj, perI, invI, perB, invB )
!       Number interior nodes
        noloc = ginde(2,dom)
        if (noloc.gt.0) then
#if defined(ON_DOMAIN_DECOMP) || defined(SIESTA__METIS)
          call METIS_NodeND( noloc, xadj, adj, 1, metisOpt, perR, invR )
#endif

          do ii = 1, noloc
            kk        = invI(ii)
            jj        = perR(ii) + offsetI
            gperm(kk) = jj
            ginvp(jj) = kk
          end do
          offsetI = offsetI + noloc
        endif
!       Number boundary nodes
        boloc = ginde(4,dom)
        do ii = 1, boloc
          kk       = invB(ii)
          jj       = ii + offsetB
          gperm(kk) = jj
          ginvp(jj) = kk
        end do
        offsetB = offsetB + boloc
      enddo
      call de_alloc(  adj,  'adj', 'domainDecom' )
      call de_alloc( xadj, 'xadj', 'domainDecom' )
      call de_alloc( perI, 'perI', 'domainDecom' )
      call de_alloc( invI, 'invI', 'domainDecom' )
      call de_alloc( perR, 'perR', 'domainDecom' )
      call de_alloc( invR, 'invR', 'domainDecom' )
      call de_alloc( perB, 'perB', 'domainDecom' )
      call de_alloc( invB, 'invB', 'domainDecom' )

!      call printPermMatrix( no, listhTot, listhPtrTot, gperm, ginvp,
!     &                      'globmat3.ps', 'Global matrix' )

      end subroutine createPerm

      subroutine compuNeigh( )
      implicit none
!     Local variables
      integer                       :: tsize, io, jo, ko, lo, PP, QQ,
     &                                 ini, fin, ind

      tsize = (Nodes-1)*Nodes/2 
      nullify(DomNeigh)
      call re_alloc( DomNeigh, 1, tsize, 'DomNeigh', 'domainDecom' )

      DomNeigh = .false.

      do PP=1, Nodes
        ini = ginde(3,PP)
        fin = ginde(3,PP) + ginde(4,PP) - 1

        do io= ini, fin
          jo = ginvp(io)
          do ko= listhPtrTot(jo-1), listhPtrTot(jo)-1
            lo = listhTot(ko)
            QQ = abs(partition(lo))
            if (PP.lt.QQ) then
              ind = PP + ((QQ-1)*(QQ-2))/2
              DomNeigh(ind) = .true.
            else if (PP.gt.QQ) then
              ind = QQ + ((PP-1)*(PP-2))/2
              DomNeigh(ind) = .true.
            endif
          enddo
        enddo
      enddo
      end subroutine compuNeigh

      subroutine compuComm( )
      implicit none
!     Local variables
      integer          :: ind, ncom, PP, QQ, KK
      integer, pointer :: src(:), dst(:)

      ncom   = 0
      ind    = 0
      do PP=2, Nodes
        do QQ=1, PP-1
          ind = ind + 1
          if (DomNeigh(ind)) ncom = ncom + 1
        enddo
      enddo

      nullify(nNeigh)
      call re_alloc( nNeigh, 1, Nodes, 'nNeigh', 'domainDecom' )

      if (ncom.gt.0) then
        nullify(src,dst)
        call re_alloc( src, 1, ncom, 'src', 'domainDecom' )
        call re_alloc( dst, 1, ncom, 'dst', 'domainDecom' )
        ncom = 0
        ind  = 0
        do PP=2, Nodes
          do QQ=1, PP-1
            ind = ind + 1
            if (DomNeigh(ind)) then
              ncom = ncom + 1
              src(ncom) = QQ
              dst(ncom) = PP
            endif
          enddo
        enddo

        Gcomm%np = Nodes
C       reschedule the communications in order to minimize the time
        call scheduleComm( ncom, src, dst, Gcomm )

        do PP=1, Nodes
          KK = 0
          do QQ= 1, Gcomm%ncol
            ind = Gcomm%ind(QQ,PP)
            if (ind.ne.0) then
              KK = KK + 1
              if (src(ind).eq.PP) then
                Gcomm%ind(KK,PP) = dst(ind)
              else
                Gcomm%ind(KK,PP) = src(ind)
              endif
            endif
          enddo
          nNeigh(PP) = KK
        enddo

        call de_alloc( dst, 'src', 'domainDecom' )
        call de_alloc( src, 'dst', 'domainDecom' )
      else
        nNeigh     = 0
        Gcomm%np   = Nodes
        Gcomm%ncol = 0
        nullify(Gcomm%ind)
      endif

      end subroutine compuComm


      subroutine createSubgraph( no, ia, ja, dom, adj, xadj,
     &                           permI, invI, permB, invB )
      implicit none
!     Input Variables
      integer, intent(in)  :: no, dom
      integer, intent(in)  :: ia(*), ja(0:no)
!     Output Variables
      integer, intent(out) :: adj(*), xadj(0:*), permI(*), invI(*),
     &                        permB(*), invB(*)
!     Local Variables
      integer              :: next, vv, ww, ii, nodeI, nodeB

      nodeI   = 0
      nodeB   = 0
      xadj(0) = 1
      next    = 1
      do vv = 1, no
        if (partition(vv) == dom) then
          nodeI       = nodeI + 1
          permI(vv)   = nodeI
          invI(nodeI) = vv
          do ii = ja(vv-1), ja(vv)-1
            ww = ia(ii)
            if (ww/=vv) then
              if (partition(ww) == dom) then
                adj(next) = ww
                next = next + 1
              endif
            endif
          enddo
          xadj(nodeI) = next

        else if (partition(vv) == -dom) then
          nodeB       = nodeB + 1
          permB(vv)   = nodeB
          invB(nodeB) = vv
        endif
      enddo

      do ii = 1, next-1
        adj(ii) = permI(adj(ii))
      enddo

      end subroutine createSubgraph


      subroutine distriGlobalMatrix( notot, no, nnz )
      use sparse_matrices, only : numh, listhptr, listh
#ifdef MPI
      use mpi_siesta
#endif
      implicit none
!     Input Variables
      integer,  intent(in) :: notot
      integer, intent(out) :: no, nnz
!     Local Variables
      integer              :: PP, QQ, RR, l_no, l_lni, l_lnb, l_nnz,
     &                        l_nei, l_bou, io, jo, ll, ind, BUFF(5)
      integer,     pointer :: linvp(:), lbinvp(:), lbsize(:), iia(:),
     &                        jja(:)
#ifdef MPI
      integer              :: MPIerror, Status(MPI_Status_Size)
#endif

      nullify(dd_nnode)
      call re_alloc( dd_nnode, 1, notot, 'dd_nnode', 'domainDecom' )
      if (Node.eq.0) dd_nnode = abs(partition) - 1
#ifdef MPI
      call MPI_BCast( dd_nnode, notot, MPI_INTEGER, 0,
     &                MPI_Comm_world, MPIerror )
#endif
      if (Node.eq.0) then
        do PP=2, Nodes
          l_lni = ginde(2,PP)
          l_lnb = ginde(4,PP)
          l_no  = l_lni + l_lnb
          l_nnz = nnzInt(PP) + nnzBou(PP)
          l_nei = nNeigh(PP)
          l_bou = 0
          do QQ= 1, l_nei
            l_bou = l_bou + ginde(4,Gcomm%ind(QQ,PP))
          enddo

          nullify(linvp)
          call re_alloc( linvp, 1, l_no, 'linvp', 'domainDecom' )
          linvp(1:ginde(2,PP)) =
     &      ginvp(ginde(1,PP):ginde(1,PP)+ginde(2,PP)-1)
          linvp(ginde(2,PP)+1:l_no) =
     &      ginvp(ginde(3,PP):ginde(3,PP)+ginde(4,PP)-1)

          nullify(lbinvp)
          call re_alloc( lbinvp, 1, l_bou, 'lbinvp', 'domainDecom' )
          ind = 1
          do QQ= 1, l_nei
            RR = Gcomm%ind(QQ,PP)
            ll = ginde(4,RR)
            lbinvp(ind:ind+ll-1) =
     &        ginvp(ginde(3,RR):ginde(3,RR)+ll-1)
            ind = ind + ll
          enddo

          nullify(lbsize)
          if (l_nei.gt.0) then
            call re_alloc( lbsize, 1, l_nei, 'lbsize', 'domainDecom' )
            do QQ= 1, l_nei
              lbsize(QQ) = ginde(4,Gcomm%ind(QQ,PP))
            enddo
          endif

          nullify(jja,iia)
          call re_alloc( jja, 1, l_no, 'jja', 'domainDecom' )
          call re_alloc( iia, 1, l_nnz, 'iia', 'domainDecom' )

          ind = 1
          DO io= 1, l_no
            jo                = linvp(io)
            ll                = listhPtrTot(jo) - listhPtrTot(jo-1)
            jja(io)           = ll
            iia(ind:ind+ll-1) = listhTot(listhPtrTot(jo-1):
     &                                   listhPtrTot(jo)-1)
            ind               = ind + ll
          ENDDO

          BUFF(1) = l_lni
          BUFF(2) = l_lnb
          BUFF(3) = l_nnz
          BUFF(4) = l_nei
          BUFF(5) = l_bou

#ifdef MPI
          call MPI_Send( BUFF, 5, MPI_INTEGER, PP-1, 0,
     &                   MPI_Comm_world, MPIerror )

          call MPI_Send( linvp, l_no, MPI_INTEGER, PP-1, 0,
     &                   MPI_Comm_world, MPIerror )

          call MPI_Send( lbinvp, l_bou, MPI_INTEGER, PP-1, 0,
     &                   MPI_Comm_world, MPIerror )

          if (l_nei.gt.0) then
            call MPI_Send( Gcomm%ind(1,PP), l_nei, MPI_INTEGER, PP-1,
     &                     0, MPI_Comm_world, MPIerror )

            call MPI_Send( lbsize, l_nei, MPI_INTEGER, PP-1, 0,
     &                     MPI_Comm_world, MPIerror )
         endif

          call MPI_Send( jja, l_no, MPI_INTEGER, PP-1, 0,
     &                   MPI_Comm_world, MPIerror )

          call MPI_Send( iia, l_nnz, MPI_INTEGER, PP-1, 0,
     &                   MPI_Comm_world, MPIerror )
#endif

          call de_alloc( jja, 'jja', 'domainDecom' )
          call de_alloc( iia, 'iia', 'domainDecom' )
          if (l_nei.gt.0)
     &      call de_alloc( lbsize, 'lbsize', 'domainDecom' )
          call de_alloc( lbinvp, 'lbinvp', 'domainDecom' )
          call de_alloc( linvp, 'linvp', 'domainDecom' )
        enddo
C       Compute data for local Process
        l_lni = ginde(2,1)
        l_lnb = ginde(4,1)
        l_no  = l_lni + l_lnb
        l_nnz = nnzInt(1) + nnzBou(1)
        l_nei = nNeigh(1)
        l_bou = 0
        do QQ= 1, l_nei
          l_bou = l_bou + ginde(4,Gcomm%ind(QQ,1))
        enddo

        nullify(dd_invp)
        call re_alloc( dd_invp, 1, l_no, 'dd_invp', 'domainDecom' )
        dd_invp(1:ginde(2,1)) =
     &    ginvp(ginde(1,1):ginde(1,1)+ginde(2,1)-1)
        dd_invp(ginde(2,1)+1:l_no) =
     &    ginvp(ginde(3,1):ginde(3,1)+ginde(4,1)-1)

        nullify(lbinvp)
        call re_alloc( lbinvp, 1, l_bou, 'lbinvp', 'domainDecom' )
        ind = 1
        do QQ= 1, l_nei
          RR = Gcomm%ind(QQ,1)
          ll = ginde(4,RR)
          lbinvp(ind:ind+ll-1) =
     &      ginvp(ginde(3,RR):ginde(3,RR)+ll-1)
          ind = ind + ll
        enddo

        nullify(dd_comm,dd_bsiz)
        if (l_nei.gt.0) then
          call re_alloc( dd_comm, 1, l_nei, 'dd_comm', 'domainDecom' )
          dd_comm = Gcomm%ind(1:l_nei,1)

          call re_alloc( dd_bsiz, 1, l_nei, 'dd_bsiz', 'domainDecom' )
          do QQ= 1, l_nei
            dd_bsiz(QQ) = ginde(4,Gcomm%ind(QQ,1))
          enddo
        endif

        call re_alloc( numh, 1, l_no, 'numh', 'sparseMat' )
        call re_alloc( listh, 1, l_nnz, 'listh', 'sparseMat' )

        ind = 1
        DO io= 1, l_no
          jo                  = dd_invp(io)
          ll                  = listhPtrTot(jo) - listhPtrTot(jo-1)
          numh(io)            = ll
          listh(ind:ind+ll-1) = listhTot(listhPtrTot(jo-1):
     &                                 listhPtrTot(jo)-1)
          ind                 = ind + ll
        ENDDO
      else
C       Receive Data from process 0
#ifdef MPI
        call MPI_recv( BUFF, 5, MPI_INTEGER, 0, 0, MPI_Comm_world,
     &                 Status, MPIerror )
#endif

        l_lni = BUFF(1)
        l_lnb = BUFF(2)
        l_nnz = BUFF(3)
        l_nei = BUFF(4)
        l_bou = BUFF(5)
        l_no  = l_lni + l_lnb

        nullify(dd_invp)
        call re_alloc( dd_invp, 1, l_no, 'dd_invp', 'domainDecom' )
#ifdef MPI
        call MPI_Recv( dd_invp, l_no, MPI_INTEGER, 0, 0,
     &                 MPI_Comm_world, Status, MPIerror )
#endif

        nullify(lbinvp)
        call re_alloc( lbinvp, 1, l_bou, 'lbinvp', 'domainDecom' )
#ifdef MPI
        call MPI_Recv( lbinvp, l_bou, MPI_INTEGER, 0, 0,
     &                 MPI_Comm_world, Status, MPIerror )
#endif

        nullify(dd_comm,dd_bsiz)
        if (l_nei.gt.0) then

          call re_alloc( dd_comm, 1, l_nei, 'dd_comm', 'domainDecom' )
#ifdef MPI
          call MPI_Recv( dd_comm, l_nei, MPI_INTEGER, 0, 0,
     &                   MPI_Comm_world, Status, MPIerror )
#endif

          call re_alloc( dd_bsiz, 1, l_nei, 'dd_bsiz', 'domainDecom' )
#ifdef MPI
          call MPI_Recv( dd_bsiz, l_nei, MPI_INTEGER, 0, 0,
     &                  MPI_Comm_world, Status, MPIerror )
#endif
        endif

        call re_alloc( numh, 1, l_no, 'numh', 'sparseMat' )
        call re_alloc( listh, 1, l_nnz, 'listh', 'sparseMat' )

#ifdef MPI
        call MPI_Recv( numh, l_no, MPI_INTEGER, 0, 0,
     &                 MPI_Comm_world, Status, MPIerror )

        call MPI_Recv( listh, l_nnz, MPI_INTEGER, 0, 0,
     &                 MPI_Comm_world, Status, MPIerror )
#endif
      endif

      no        = l_no
      nnz       = l_nnz
      dd_lni    = l_lni
      dd_lnb    = l_lnb
      dd_ncom   = l_nei
      dd_ncolum = l_no + l_bou
      call re_alloc( listhptr, 1, no, 'listhptr', 'sparseMat' )
      ind = 0
      do io= 1, no
        listhptr(io) = ind
        ind          = ind + numh(io)
      enddo

      nullify(dd_perm)
      call re_alloc( dd_perm, 1, notot, 'dd_perm','domainDecom' )
      dd_perm = 0
      do io= 1, no
        ind          = dd_invp(io)
        dd_perm(ind) = io
      enddo


      nullify(dd_cperm)
      call re_alloc( dd_cperm, 1, notot, 'dd_cperm',
     &               'domainDecom' )
      dd_cperm = dd_perm
      do io= 1, l_bou
        ind           = lbinvp(io)
        dd_cperm(ind) = io + l_no
      enddo

      call de_alloc( lbinvp, 'lbinvp', 'domainDecom' )

#ifdef _SLEPC_
      if (node.ne.0)
     &  call re_alloc(slepc_perm, 1, notot,'slepc_perm','domainDecom')
#ifdef MPI
      call MPI_BCast( slepc_perm, notot, MPI_INTEGER, 0,
     &                MPI_Comm_world, MPIerror )
#endif
#endif

      use_dd_perm = .true.
      dd_nuo      = no
      end subroutine distriGlobalMatrix

!
!     This is just a debug function to gather the global matrix into
!     process 0 when the files are distributed using DomainDecomposition
!
      subroutine gatherGlobalMatrix( nuotot, nuo, nnz, ML, filename )
#ifdef MPI
      use mpi_siesta
#endif
      use precision
      use sparse_matrices, only : numh, listh
      implicit none
!     Input Variables
      integer,  intent(in) :: nuotot, nuo, nnz
      real(dp),     target :: ML(:)
      character(*)         :: filename
!     Local Variables
      integer              :: PP, nuo_l, io, jo, ko, nnzG, nnz_l,
     &                        ll, ind, fd
      integer,     pointer :: ja(:), ia(:), jndx(:)
      real(dp),    pointer :: M(:), an(:)
#ifdef MPI
      integer              :: MPIerror, Status(MPI_Status_Size)
#endif
      integer              :: ival(1)
      complex              :: cval(1)
!
!     Need arrays: ginde, ginvp
!
      if (node.eq.0) then
        nullify(listhPtrTot)
        call re_alloc( listhPtrTot, 0, nuotot, 'listhPtrTot',
     &                 'domainDecom' )

        do PP=1, Nodes
          nuo_l = ginde(2,PP) + ginde(4,PP)
          if (PP.eq.1) then
            ja => numh
          else
            nullify(ja)
            call re_alloc( ja, 1, nuo_l, 'ja', 'domainDecom' )
#ifdef MPI
            call MPI_recv( ja, nuo_l, MPI_INTEGER,
     &                     PP-1, 1, MPI_Comm_world, Status, MPIerror )
#endif
          endif
          jo = ginde(1,PP)
          do io= 1, ginde(2,PP)
            ko              = ginvp(jo)
            listhPtrTot(ko) = ja(io)
            jo              = jo + 1
          enddo
          jo = ginde(3,PP)
          do io= io, nuo_l
            ko              = ginvp(jo)
            listhPtrTot(ko) = ja(io)
            jo              = jo + 1
          enddo
          if (PP.eq.1) then
            ja => null()
          else
            call de_alloc( ja, 'ja', 'domainDecom' )
          endif
        enddo
        listhPtrTot(0) = 1
        do io= 1, nuotot
          listhPtrTot(io) = listhPtrTot(io) + listhPtrTot(io-1)
        enddo
        nnzG = listhPtrTot(nuotot)-1

        nullify(listhTot,M)
        call re_alloc( listhTot, 1, nnzG, 'listhTot',
     &                 'domainDecom' )
        call re_alloc( M, 1, nnzG, 'M', 'domainDecom' )
        do PP=1, Nodes
          nuo_l = ginde(2,PP) + ginde(4,PP)
          nnz_l = nnzInt(PP) + nnzBou(PP)
          if (PP.eq.1) then
            ia => listh
            an => ML
          else
            nullify(ia,an)
            call re_alloc( ia, 1, nnz_l, 'ia', 'domainDecom' )
            call re_alloc( an, 1, nnz_l, 'an', 'domainDecom' )
#ifdef MPI
            call MPI_recv( ia, nnz_l, MPI_INTEGER,
     &                     PP-1, 1, MPI_Comm_world, Status, MPIerror )
            call MPI_recv( an, nnz_l, MPI_DOUBLE_PRECISION,
     &                     PP-1, 1, MPI_Comm_world, Status, MPIerror )
#endif
          endif
          jo   = ginde(1,PP)
          ind = 1
          do io= 1, ginde(2,PP)
            ko  = ginvp(jo)
            ll  = listhPtrTot(ko)-listhPtrTot(ko-1)
            listhTot(listhPtrTot(ko-1):listhPtrTot(ko)-1) =
     &        ia(ind:ind+ll-1)
            M(listhPtrTot(ko-1):listhPtrTot(ko)-1) =
     &        an(ind:ind+ll-1)
            jo  = jo + 1
            ind = ind + ll
          enddo
          jo = ginde(3,PP)
          do io= io, nuo_l
            ko  = ginvp(jo)
            ll  = listhPtrTot(ko)-listhPtrTot(ko-1)
            listhTot(listhPtrTot(ko-1):listhPtrTot(ko)-1) =
     &        ia(ind:ind+ll-1)
            M(listhPtrTot(ko-1):listhPtrTot(ko)-1) =
     &        an(ind:ind+ll-1)
            jo  = jo + 1
            ind = ind + ll
          enddo
          if (PP.eq.1) then
            ia => null()
            an => null()
          else
            call de_alloc( ia, 'ia', 'domainDecom' )
            call de_alloc( an, 'an', 'domainDecom' )
          endif
        enddo

        call io_assign( fd )
        open( unit=fd, file=filename, status='unknown' )

        call re_alloc( jndx, 1, nnzG, 'jndx', 'domainDecom' )
        do io= 1, nuotot
          jndx(listhPtrTot(io-1):listhPtrTot(io)-1) = io
        enddo

        call mmwrite( fd, 'coordinate', 'real', 'general', nuotot,
     &                nuotot, nnzG, listhTot, jndx, ival, M, cval )
        call de_alloc( jndx, 'jndx', 'domainDecom' )


!        write(fd,*) nuotot, nuotot, nnzG
!        write(fd,*) listhPtrTot
!        write(fd,*) listhTot
!        write(fd,*) M

!        do io= 1, nuotot
!          do jo= listhPtrTot(io-1), listhPtrTot(io)-1
!            write(fd,*) io, listhTot(jo), M(jo)
!          enddo
!        enddo
!        call pxfflush(fd)

        call io_close( fd )

        call de_alloc( ja, 'ja', 'domainDecom' )
        call de_alloc( an, 'an', 'domainDecom' )
        call de_alloc( M, 'M', 'domainDecom' )
        call de_alloc( listhPtrTot, 'listhPtrTot','domainDecom' )
        call de_alloc( listhTot, 'listhTot', 'domainDecom' )
      else
#ifdef MPI
        call MPI_send( numh, nuo, MPI_INTEGER, 0, 1,
     &                 MPI_Comm_world, MPIerror )
        call MPI_send( listh, nnz, MPI_INTEGER, 0, 1,
     &                 MPI_Comm_world, MPIerror )
        call MPI_send( ML, nnz, MPI_DOUBLE_PRECISION, 0, 1,
     &                 MPI_Comm_world, MPIerror )
#endif
      endif
      end subroutine gatherGlobalMatrix

      subroutine resetDomainDecom(  )
      implicit none

      if (associated(dd_perm))
     &  call de_alloc( dd_perm, 'dd_perm','domainDecom' )

      if (associated(dd_invp))
     &  call de_alloc( dd_invp, 'dd_invp','domainDecom' )

      if (associated(dd_cperm))
     &  call de_alloc( dd_cperm, 'dd_cperm','domainDecom' )

      if (associated(dd_comm))
     &  call de_alloc( dd_comm, 'dd_comm','domainDecom' )

      if (associated(dd_bsiz))
     &  call de_alloc( dd_bsiz, 'dd_bsiz','domainDecom' )

      if (associated(dd_nnode))
     &  call de_alloc( dd_nnode, 'dd_nnode','domainDecom' )

#ifdef _SLEPC_
      if (associated(slepc_perm))
     &  call de_alloc( slepc_perm, 'slepc_perm','domainDecom' )
#endif

      end subroutine resetDomainDecom


      subroutine preSetOrbitLimits( no )
      implicit none
C     Input variables
      integer      :: no
C     Local variables
      integer      :: blocks, remain
      
      blocks = no/Nodes
      remain = no - blocks*Nodes
      
      llimit = blocks*Node + 1 + min(remain,Node)
      ulimit = blocks*(Node+1) + 1 + min(remain,Node+1)
      dd_nuo = ulimit - llimit
      end subroutine preSetOrbitLimits

      end module domain_decom
