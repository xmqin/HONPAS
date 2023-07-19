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
      module meshdscf
C
C Stores quantities that are connected with Dscf in mesh local
C form when data is distributed for parallel execution
C
      use precision, only: dp
      use alloc
      implicit none

C ----------------------------------------------------------------------
C Dscf related variables for parallel distributed form
C ----------------------------------------------------------------------
C integer listdl(sum(numdl))           : List of non-zero elements in
C                                      : a row of DscfL
C integer listdlptr(nrowsDscfL)        : Pointer to row in listdl
C integer NeedDscfL(nuotot)            : Pointer as to whether a row of
C                                      : Dscf is needed in DscfL
C integer nrowsDscfL                   : Number of rows of DscfL
C integer numdl(nrowsDscfL)            : Number of non-zero elements in
C                                      : a row of DscfL
C real(dp)  DscfL(maxndl,nrowsDscfL)     : Local copy of Dscf elements
C                                      : needed for the local mesh
C ----------------------------------------------------------------------

      integer, save :: nrowsDscfL

      integer, dimension(:), pointer, save :: listdl
      integer, dimension(:), pointer, save :: listdlptr
      integer, dimension(:), pointer, save :: NeedDscfL
      integer, dimension(:), pointer, save :: numdl

      real(dp), dimension(:,:), pointer, save :: DscfL

      logical, public, save :: nullified_pointers = .false.

      end module meshdscf

      subroutine CreateLocalDscfPointers( nmpl, nuotot, numd,
     .                                    listdptr, listd )
C
C Calculates the values of the orbitals at the mesh points
C
C ----------------------------------------------------------------------
C Input :
C ----------------------------------------------------------------------
C integer nmpl          : Number of mesh points in unit cell locally
C integer nuotot        : Total number of basis orbitals in unit cell
C integer numd(nuo)     : Number of nonzero density-matrix
C                       : elements for each matrix row
C integer listdptr(nuo) : Pointer to start of rows of density-matrix
C integer listd(maxnh)  : Nonzero-density-matrix-element column
C                       : indexes for each matrix row
C ----------------------------------------------------------------------
C Output :
C ----------------------------------------------------------------------
C All output quantities are in the module meshdscf
C ----------------------------------------------------------------------

C
C Modules
C
      use atomlist,     only : indxuo
      use meshdscf,     only : listdl, listdlptr, NeedDscfL, 
     .                         nrowsDscfL, numdl, nullified_pointers
      use meshphi,      only : endpht, lstpht
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb
      use precision
      use alloc
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

C
C Passed arguments
C
      integer, intent(in) :: nmpl
      integer, intent(in) :: nuotot
      integer, intent(in) :: numd(*)
      integer, intent(in) :: listdptr(*)
      integer, intent(in) :: listd(*)

C
C Local variables
C
      integer
     .  BNode, i, ii, io, iio, ip, imp, iu, numdele, maxndmax,
     .  nsize
#ifdef MPI
      integer
     .  MPIerror
#endif

      integer, dimension(:), pointer :: ibuffer

      if (.not. nullified_pointers) then
        nullify( NeedDscfL )
        nullify( listdl )
        nullify( listdlptr )
        nullify( numdl )
        nullified_pointers = .true.
      endif

C Create pointer as to whether a given row of DscfL is needed in NeedDscfL
C This pointer is never deallocated...

      call re_alloc( NeedDscfL, 1, nuotot, name='NeedDscfL',
     &                 routine='CreateLocalDscfPointers' )

      NeedDscfL(1:nuotot) = 0
      do ip = 1,nmpl
        do imp = 1+endpht(ip-1), endpht(ip)
          i = lstpht(imp)
          iu = indxuo(i)
          NeedDscfL(iu) = 1
        enddo
      enddo
      nrowsDscfL = 0
      do i = 1,nuotot
        if (NeedDscfL(i).eq.1) then
          nrowsDscfL = nrowsDscfL + 1
          NeedDscfL(i) = nrowsDscfL
        endif
      enddo

C Allocate/reallocate memory for numdl and listdlptr

      call re_alloc( numdl, 1, max(1,nrowsDscfL), name='numdl',
     &                 routine='CreateLocalDscfPointers' )
      call re_alloc( listdlptr, 1, max(1,nrowsDscfL), 
     &                 name='listdlptr',
     &                 routine='CreateLocalDscfPointers' )

C Distribute information about numd globally
      maxndmax = 0
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          numdele = numd(iio)
        endif
#ifdef MPI
        call MPI_Bcast(numdele,1,MPI_integer,BNode,
     .    MPI_Comm_World,MPIerror)
#endif
        if (NeedDscfL(io).gt.0) numdl(NeedDscfL(io)) = numdele
        if (numdele.gt.maxndmax) maxndmax = numdele
      enddo

C Create listdlptr using numdl
      listdlptr(1) = 0
      do io = 2,nrowsDscfL
        listdlptr(io) = listdlptr(io-1) + numdl(io-1)
      enddo

C Allocate/reallocate listdl
      if (nrowsDscfL.gt.0) then
        nsize = listdlptr(nrowsDscfL)+numdl(nrowsDscfL)
      else
        nsize = 1
      endif
      call re_alloc( listdl, 1, nsize, name='listdl',
     &               routine='CreateLocalDscfPointers' )

C Distribute information about listd globally
      nullify( ibuffer )
      call re_alloc( ibuffer, 1, maxndmax, name='ibuffer',
     &               routine='CreateLocalDscfPointers' )
      ibuffer(1:maxndmax) = 0      ! AG: superfluous
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do ii = 1,numd(iio)
            ibuffer(ii) = listd(listdptr(iio)+ii)
          enddo
        endif
#ifdef MPI
        call MPI_Bcast(ibuffer,maxndmax,MPI_integer,BNode,
     .    MPI_Comm_World,MPIerror)
#endif
        iio = NeedDscfL(io)
        if (iio.gt.0) then
          do ii = 1,numdl(iio)
            listdl(listdlptr(iio)+ii) = ibuffer(ii)
          enddo
        endif
      enddo
      call de_alloc( ibuffer, name='ibuffer' )

      end subroutine CreateLocalDscfPointers

      subroutine matrixOtoM( maxnd, numd, listdptr, maxndl, nuo, 
     .                       nuotot, nspin, Dscf, DscfL )
C ********************************************************************
C Transforms a matrix which is distributed by block cyclic 
C distribution of orbitals to a matrix that contains all
C the orbital rows needed for a mesh point distribution 
C over the nodes.
C Created by J.D.Gale, February 2000
C *********************** INPUT **************************************
C integer maxnd         : First dimension of Dscf
C integer numd(nuo)     : Number of non-zero elements in row of Dscf
C integer listdptr(nuo) : Pointer to start of rows in Dscf
C integer maxndl        : First dimension of DscfL
C integer nuo           : Local no. of orbitals in unit cell
C integer nuotot        : Total no. of orbitals in unit cell
C integer nspin         : Number of spin components
C real*8  Dscf(maxnd,nspin) : Matrix in orbital distributed form
C *********************** OUTPUT *************************************
C real*8  DscfL(maxndl,nspin) : Matrix in mesh distributed form
C ********************************************************************

C  Modules
      use precision
      use meshdscf, only : nrowsDscfL, numdl, listdlptr, NeedDscfL
      use alloc,    only : re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
      use parallel,     only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
#endif

      implicit none

C Argument types and dimensions
      integer
     .   maxnd, maxndl, nspin, nuo, nuotot, numd(nuo), 
     .   listdptr(nuo)

      real(dp)
     .   Dscf(maxnd,nspin), DscfL(maxndl,nspin)

      external
     .   memory

C Internal variables and arrays
      integer
     .  io, il, ispin

#ifdef MPI
      integer 
     .  BNode, ii, iio, maxno, maxnog, MPIerror
      real(dp), dimension(:), pointer ::  buffer
#endif

C***********************
C  Parallel execution  *
C***********************
#ifdef MPI
C Find local and global size of maxno
      maxno = 0
      do io = 1,nrowsDscfL
        maxno = max(maxno,numdl(io))
      enddo
      call MPI_AllReduce(maxno,maxnog,1,MPI_integer,MPI_max,
     .  MPI_Comm_World,MPIerror)

C Allocate local Dscf storage array
      nullify(buffer)
      call re_alloc( buffer, 1, maxnog*nspin, name='buffer',
     &               routine='matrixOtoM' )

C Zero buffer as we will be passing the complete array
C as the number of explicitly non-zero elements will
C not be known on all nodes
      buffer(1:maxnog*nspin) = 0.0d0  ! AG: superfluous

C Loop over rows of Dscf broadcasting to all other nodes
      do io = 1,nuotot
        call WhichNodeOrb(io,Nodes,BNode)

C If this row is local to this node then copy into buffer
        if (Node.eq.BNode) then
          call GlobalToLocalOrb(io,Node,Nodes,iio)
          do ispin = 1,nspin
            ii = (ispin-1)*maxnog
            buffer(ii+1:ii+numd(iio)) = 
     .        Dscf(listdptr(iio)+1:listdptr(iio)+numd(iio),ispin)
          enddo
        endif

C Broadcast buffer array
        call MPI_Bcast(buffer,maxnog*nspin,MPI_double_precision,BNode,
     .                 MPI_Comm_World,MPIerror)

C Get pointer for this row of Dscf and see if it is needed for DscfL
        il = NeedDscfL(io)
        if (il.gt.0) then
          do ispin = 1,nspin
            ii = (ispin-1)*maxnog
            DscfL(listdlptr(il)+1:listdlptr(il)+numdl(il),ispin) = 
     .        buffer(ii+1:ii+numdl(il))
          enddo
        endif
      enddo

C Deallocate buffer memory
      call de_alloc(buffer, name="buffer")
#else
C*********************
C  Serial execution  *
C*********************
C Loop over rows of Dscf checking to see if they are in DscfL
      do ispin = 1,nspin
        do io = 1,nuotot

C Get pointer for this row of Dscf and see if it is needed for DscfL
          il = NeedDscfL(io)
          if (il.gt.0) then
            DscfL(listdlptr(il)+1:listdlptr(il)+numdl(il),ispin) = 
     .        Dscf(listdptr(io)+1:listdptr(io)+numdl(il),ispin)
          endif

        enddo
      enddo
#endif

      end subroutine matrixOtoM

      subroutine matrixMtoO( maxnvl, maxnv, numVs, listVsptr, nuo, 
     .                       nuotot, nspin, VsL, Vs )

C ********************************************************************
C Transforms a matrix which is distributed by mesh points to a matrix
C that is distributed by a block cyclic distribution over the orbitals
C and adds it to an existing array of this form.
C Created by J.D.Gale, February 2000
C *********************** INPUT **************************************
C integer maxnvl          : First dimension of VsL and maximum number
C                           of nonzero elements in VsL
C integer maxnv           : First dimension of Vs and maximum number
C                           of nonzero elements in Vs
C integer numVs(nuo)      : Number of non-zero elements in row of Vs
C integer listVsptr(nuo)  : Pointer to start of rows in Vs
C integer nuo             : Local no. of orbitals in unit cell
C integer nuotot          : Total no. of orbitals in unit cell
C integer nspin           : Number of spin components
C real*8  VsL(maxnvl,nspin) : Mesh contribution to be added to Vs
C ******************** INPUT AND OUTPUT *******************************
C real*8  Vs(maxnv,nspin) : Value of nonzero elements in each row of Vs
C                           to which the potential matrix elements are
C                           summed up
C *********************************************************************

C  Modules
      use precision
      use meshdscf
      use alloc, only: re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
      use parallel,     only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb,
     .                         LocalToGlobalOrb
#endif

      implicit none

C Argument types and dimensions
      integer
     .   maxnv, maxnvl, nspin, nuo, nuotot, numVs(nuo), 
     .   listVsptr(nuo)
      real(dp)
     .   Vs(maxnv,nspin), VsL(maxnvl,nspin)

C Internal variables and arrays
      integer
     .  i, iu, ispin

#ifdef MPI
      integer :: ii, in, iul, maxnvg, MPIerror, nVsL, nVsLmax
      real(dp), dimension(:), pointer ::  Vi
      integer, dimension(:), pointer  ::  nVsLPtr, listViptr
#endif
      
C***********************
C  Parallel execution  *
C***********************
#ifdef MPI
C Find the maximum number of rows on anyone Node
      call MPI_AllReduce(nrowsDscfL,nVsLmax,1,MPI_integer,
     .  MPI_max,MPI_Comm_World,MPIerror)

C Find the maximum value of maxno on any Node and globally
      call MPI_AllReduce(maxnvl,maxnvg,1,MPI_integer,
     .  MPI_max,MPI_Comm_World,MPIerror)

C Allocate buffer memory for transfer of matrix and pointer
      nullify(Vi)
      call re_alloc( Vi, 1, maxnvg*nspin, name='Vi',
     &               routine='matrixMtoO' )
      nullify(nVsLPtr)
      call re_alloc( nVsLPtr, 1, nuotot, name='nVsLPtr',
     &               routine='matrixMtoO' )
      nullify(listViptr)
      call re_alloc( listViptr, 1, nuotot, name='listViptr',
     &               routine='matrixMtoO' )

C Initialise buffers to avoid number trapping errors during transfer
      Vi(1:maxnvg) = 0.0d0  ! AG : superfluous
      listViptr(1:nuotot) = 0

C Loop over Nodes for broadcasting of local data
      do in = 0,Nodes-1

C Broadcast the number of rows stored locally
        if (in.eq.Node) nVsL = nrowsDscfL
        call MPI_Bcast(nVsL,1,MPI_integer,in,MPI_Comm_World,MPIerror)

C If this is the broadcasting Node then copy the pointer to the buffer
        if (in.eq.Node) then
          nVsLPtr(1:nuotot) = NeedDscfL(1:nuotot)
        endif

C Broadcast the pointer information
        call MPI_Bcast(nVsLPtr,nuotot,MPI_integer,in,MPI_Comm_World,
     .                 MPIerror)

C Need to broadcast the row pointer info for Vi (listdlptr)
        if (in.eq.Node.and.nrowsDscfL.gt.0) then
          listViptr(1:nrowsDscfL) = listdlptr(1:nrowsDscfL)
        endif
        call MPI_Bcast(listViptr,nuotot,MPI_integer,in,MPI_Comm_World,
     .                 MPIerror)

C If this is the broadcasting Node then copy the data to the buffer
        if (in.eq.Node) then
          do ispin = 1,nspin
            ii = (ispin-1)*maxnvg
            Vi(ii+1:ii+maxnvl) = VsL(1:maxnvl,ispin)
          enddo
        endif

C Broadcast the buffer information
        call MPI_Bcast(Vi,maxnvg*nspin,MPI_double_precision,in,
     .                 MPI_Comm_World,MPIerror)

C Add those elements that are needed locally to the values already
C stored in the orbital oriented array
        do ispin = 1,nspin
          ii = (ispin-1)*maxnvg
          do i = 1,nuo
            call LocalToGlobalOrb(i,Node,Nodes,iu)
            iul = nVsLPtr(iu)
            if (iul.gt.0) then
              Vs(listVsptr(i)+1:listVsptr(i)+numVs(i),ispin) = 
     .          Vs(listVsptr(i)+1:listVsptr(i)+numVs(i),ispin) + 
     .          Vi(ii+listViptr(iul)+1:ii+listViptr(iul)+numVs(i))
            endif
          enddo
        enddo

      enddo

C Deallocate buffer and pointer memory
      call de_alloc(listViptr, name="listViptr")
      call de_alloc(nVsLPtr, name="nVsLPtr")
      call de_alloc(Vi, name="Vi")
#else
C*********************
C  Serial execution  *
C*********************
C Add those elements that are needed locally to the values already
C stored in the orbital oriented array
      do ispin = 1,nspin
        do i = 1,nuotot
          iu = NeedDscfL(i)
          if (iu.gt.0) then
            Vs(listVsptr(i)+1:listVsptr(i)+numVs(i),ispin) = 
     .        Vs(listVsptr(i)+1:listVsptr(i)+numVs(i),ispin) + 
     .        VsL(listdlptr(iu)+1:listdlptr(iu)+numVs(i),ispin)
          endif
        enddo
      enddo
#endif

      end
