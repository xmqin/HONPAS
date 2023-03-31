! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module meshdscf
C
C Stores quantities that are connected with Dscf in mesh local
C form when data is distributed for parallel execution
C
      use precision, only: dp
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
      integer, public            :: nrowsDscfL
      integer,  pointer, public  :: listdl(:), listdlptr(:),
     $                              NeedDscfL(:), numdl(:)
      real(dp), pointer, public  :: DscfL(:,:)
      logical           :: first_time = .true.

      public ::  matrixOtoM, matrixMtoO, matrixMtoOC
      public ::  resetDscfPointers
      public ::  CreateLocalDscfPointers

      private

      CONTAINS

      subroutine resetDscfPointers( )
      use alloc, only : de_alloc
      use m_dscfcomm, only: resetdscfComm
      implicit none
      call resetdscfComm( )
      call de_alloc( listdl,    'listdl',    'meshdscf' )
      call de_alloc( listdlptr, 'listdlptr', 'meshdscf' )
      call de_alloc( NeedDscfL, 'NeedDscfL', 'meshdscf' )
      call de_alloc( numdl,     'numdl',     'meshdscf' )
      call de_alloc( DscfL,     'DscfL',     'meshdscf' )
      end subroutine resetDscfPointers

      subroutine CreateLocalDscfPointers( nmpl, nuotot, numd,
     &                                    listdptr, listd )
C
C Calculates the values of the orbitals at the mesh points
C
C Update: Computes the communications needed to move data ordered
C by orbitals to data ordered by mesh (function dscfComm). All-to-all
C communications has been substituted by point-to-point communications.
C Written by Rogeli Grima (BSC) Dec.2007
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
      use meshphi,      only : endpht, lstpht
      use parallel,     only : Node, Nodes
      use parallelsubs, only : GlobalToLocalOrb, WhichNodeOrb
      use precision
      use alloc
      use m_dscfComm,   only : dscfcomm
      use m_dscfComm,   only : DCsize, DCself, DCtotal, DCmax
      use m_dscfComm,   only : DCmaxnd
      use m_dscfComm,   only : DCsrc, DCdst, DCsic, DCinvp


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
      integer          :: i, ii, io, iio, ip, imp, iu, numdele,
     &                    maxndmax, nsize, nn, mm 
      integer, pointer :: ibuffer(:)
#ifdef MPI
      integer          :: MPIerror, Status(MPI_Status_Size)
#endif

#ifdef DEBUG
      call write_debug( '      PRE CreateLocalDscfPointers' )
#endif

      if (first_time) then
        first_time = .false.
      else
        call resetDscfPointers( )
      endif

      nullify( NeedDscfL, listdl, listdlptr, numdl )

C     Create pointer as to whether a given row of DscfL is needed in NeedDscfL
C     This pointer is never deallocated...
      call re_alloc( NeedDscfL, 1, nuotot, 'NeedDscfL', 'meshdscf' )

      do io= 1, nuotot
        NeedDscfL(io) = 0
      enddo

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

C     Computes the communications needed to move data ordered
C     by orbitals to data ordered by mesh.
      call dscfComm( nuotot, nrowsDscfL, NeedDscfL )

C     Allocate/reallocate memory for numdl and listdlptr
      call re_alloc( numdl, 1, max(1,nrowsDscfL), 'numdl', 'meshdscf' )
      call re_alloc( listdlptr, 1, max(1,nrowsDscfL), 
     &               'listdlptr', 'meshdscf' )

C     Distribute information about numd globally
C     Use communications scheduling precomputed in dscfComm
      nullify( ibuffer )
      call re_alloc( ibuffer, 1, DCmax, 'ibuffer', 'meshdscf' )

C     This data is in the current process. No communications needed
      do ii= 1, DCself
        io = DCinvp(ii)
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        numdele              = numd(iio)
        numdl(NeedDscfL(io)) = numdele
      enddo

      nn       = DCself + 1
      maxndmax = 0
#ifdef MPI
      do i= 1, DCsize
C       For all the needed communications
        numdele = 0
        if (DCsrc(i).eq.Node) then
C         If this is the source node, copy data into a buffer
C         and send it to the receiver.
          do ii= 1, DCsic(i)
            io = DCinvp(nn)
            call GlobalToLocalOrb(io,Node,Nodes,iio)
            ibuffer(ii) = numd(iio)
            numdele = numdele + ibuffer(ii)
            nn = nn + 1
          enddo
          call MPI_Send( ibuffer, DCsic(i), MPI_integer,
     &                   DCdst(i), 1, MPI_Comm_World, MPIerror )
        else
C         If this is the destination node, receive data from source node.
          call mpi_recv( ibuffer, DCsic(i), MPI_integer,
     &                   DCsrc(i),  1, MPI_Comm_world, Status,
     &                   MPIerror )
          do ii= 1, DCsic(i)
            io                   = DCinvp(nn)
            numdl(NeedDscfL(io)) = ibuffer(ii)
            numdele              = numdele + ibuffer(ii)
            nn                   = nn + 1
          enddo
        endif
        if (numdele.gt.maxndmax) maxndmax = numdele
      enddo
#endif
      call de_alloc( ibuffer, 'ibuffer', 'meshdscf' )
      DCmaxnd = maxndmax

C     Create listdlptr using numdl
      listdlptr(1) = 0
      do io = 2,nrowsDscfL
        listdlptr(io) = listdlptr(io-1) + numdl(io-1)
      enddo

C     Allocate/reallocate listdl
      if (nrowsDscfL.gt.0) then
        nsize = listdlptr(nrowsDscfL)+numdl(nrowsDscfL)
      else
        nsize = 1
      endif
      call re_alloc( listdl, 1, nsize, 'listdl', 'meshdscf' )

C     Distribute information about listd globally
C     Use communications scheduling precomputed in dscfComm
      nullify( ibuffer )
      call re_alloc( ibuffer, 1, DCmaxnd, 'ibuffer', 'meshdscf' )

C     This data is in the current process. No communications needed
      do ii= 1, DCself
        io = DCinvp(ii)
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        iu = NeedDscfL(io)
        do i = 1,numd(iio)
          listdl(listdlptr(iu)+i) = listd(listdptr(iio)+i)
        enddo
      enddo

      nn = DCself + 1
#ifdef MPI
      do i= 1, DCsize
C       For all the needed communications
        if (DCsrc(i).eq.Node) then
C         If this is the source node, copy data into a buffer
C         and send it to the receiver.
          mm = 0
          do ii= 1, DCsic(i)
            io = DCinvp(nn)
            call GlobalToLocalOrb(io,Node,Nodes,iio)
            ibuffer(mm+1:mm+numd(iio)) = listd(listdptr(iio)+1:
     &                                         listdptr(iio)+numd(iio))
            mm = mm + numd(iio)
            nn = nn + 1
          enddo
          call MPI_Send( ibuffer, mm, MPI_integer, DCdst(i), 1, 
     &                   MPI_Comm_World, MPIerror )
        else
C         If this is the destination node, receive data from source node.
          mm = 0
          do ii= 0, DCsic(i)-1
            io = DCinvp(nn+ii)
            mm = mm + numdl(NeedDscfL(io))
          enddo
          call mpi_recv( ibuffer, mm, MPI_integer, DCsrc(i), 1,
     &                   MPI_Comm_world, Status, MPIerror )
          mm = 0
          do ii= 1, DCsic(i)
            io = DCinvp(nn)
            iu = NeedDscfL(io)
            listdl(listdlptr(iu)+1:listdlptr(iu)+numdl(iu)) =
     &              ibuffer(mm+1:mm+numdl(iu))
            mm = mm + numdl(iu)
            nn = nn + 1
          enddo
        endif
      enddo
#endif
      call de_alloc( ibuffer, 'ibuffer', 'meshdscf' )

#ifdef DEBUG
      call write_debug( '      POS CreateLocalDscfPointers' )
#endif
      end subroutine CreateLocalDscfPointers

      subroutine matrixOtoM( maxnd, numd, listdptr, maxndl, nuo, 
     &                       nspin, Dscf, DscfL )
C ********************************************************************
C Transforms a matrix which is distributed by block cyclic 
C distribution of orbitals to a matrix that contains all
C the orbital rows needed for a mesh point distribution 
C over the nodes.
C Created by J.D.Gale, February 2000
C
C Update: All-to-all communications has been substituted by
C point-to-point communications. These have been precomputed in 
C dscfComm.
C Written by Rogeli Grima (BSC) Dec.2007
C
C *********************** INPUT **************************************
C integer maxnd         : First dimension of Dscf
C integer numd(nuo)     : Number of non-zero elements in row of Dscf
C integer listdptr(nuo) : Pointer to start of rows in Dscf
C integer maxndl        : First dimension of DscfL
C integer nuo           : Local no. of orbitals in unit cell
C integer nspin         : Number of spin components
C real*8  Dscf(maxnd,nspin) : Matrix in orbital distributed form
C *********************** OUTPUT *************************************
C real*8  DscfL(maxndl,nspin) : Matrix in mesh distributed form
C ********************************************************************

C  Modules
      use precision
      use alloc,        only : re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
      use parallel,     only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
      use m_dscfComm
#endif
      implicit none

C Argument types and dimensions
      integer           :: maxnd, maxndl, nspin, nuo, numd(nuo),
     &                     listdptr(nuo)
      real(dp)          :: Dscf(maxnd,nspin), DscfL(maxndl,nspin)
      external             memory

C Internal variables and arrays
      integer           :: io, ispin
#ifdef MPI
      integer           :: i, ii, iu, mm, nn, iio, MPIerror,
     &                     Status(MPI_Status_Size)
      real(dp), pointer :: buffer(:)
#else
      integer           :: il
#endif

C***********************
C  Parallel execution  *
C***********************
#ifdef MPI
C Allocate local Dscf storage array
      nullify(buffer)
      call re_alloc( buffer, 1, DCmaxnd*nspin, 'buffer', 'meshdscf' )

C     Copy the data that is in the current process
      do ii= 1, DCself
        io = DCinvp(ii)
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        iu = NeedDscfL(io)
        do ispin = 1,nspin
          DscfL(listdlptr(iu)+1:listdlptr(iu)+numdl(iu),ispin) =
     &    Dscf(listdptr(iio)+1:listdptr(iio)+numd(iio),ispin)
        enddo
      enddo
      nn = DCself

C     Send or receive the data from/to other processes
      do i= 1, DCsize
        if (DCsrc(i).eq.Node) then
          mm = 0
          do ispin = 1,nspin
            do ii= 1, DCsic(i)
              io = DCinvp(nn+ii)
              call GlobalToLocalOrb(io,Node,Nodes,iio)
              buffer(mm+1:mm+numd(iio)) = Dscf(listdptr(iio)+1:
     &                           listdptr(iio)+numd(iio),ispin)
              mm = mm + numd(iio)
            enddo
          enddo
          call MPI_Send( buffer, mm, MPI_double_precision, DCdst(i), 1,
     &                   MPI_Comm_World, MPIerror )
        else
          mm = 0
          do ii= 1, DCsic(i)
            io = DCinvp(nn+ii)
            iu = NeedDscfL(io)
            mm = mm + numdl(iu)
          enddo
          mm = mm*nspin
          call mpi_recv( buffer, mm, MPI_double_precision, DCsrc(i), 1,
     &                   MPI_Comm_world, Status, MPIerror )
          mm = 0
          do ispin = 1,nspin
            do ii= 1, DCsic(i)
              io = DCinvp(nn+ii)
              iu = NeedDscfL(io)
              DscfL(listdlptr(iu)+1:listdlptr(iu)+numdl(iu),ispin) =
     &              buffer(mm+1:mm+numdl(iu))
              mm = mm + numdl(iu)
            enddo
          enddo
        endif
        nn = nn + DCsic(i)
      enddo
      call de_alloc( buffer, 'buffer', 'meshdscf' )
#else
C*********************
C  Serial execution  *
C*********************
C Loop over rows of Dscf checking to see if they are in DscfL
      do ispin = 1,nspin
        do io = 1,nuo

C Get pointer for this row of Dscf and see if it is needed for DscfL
          il = NeedDscfL(io)
          if (il.gt.0) then
            DscfL(listdlptr(il)+1:listdlptr(il)+numdl(il),ispin) = 
     &        Dscf(listdptr(io)+1:listdptr(io)+numdl(il),ispin)
          endif

        enddo
      enddo
#endif
      end subroutine matrixOtoM

      subroutine matrixMtoO( maxnvl, maxnv, numVs, listVsptr, nuo, 
     &                       nspin, VsL, Vs )
C ********************************************************************
C Transforms a matrix which is distributed by mesh points to a matrix
C that is distributed by a block cyclic distribution over the orbitals
C and adds it to an existing array of this form.
C Created by J.D.Gale, February 2000
C
C Update: All-to-all communications has been substituted by
C point-to-point communications. These have been precomputed in 
C dscfComm.
C Written by Rogeli Grima (BSC) Dec.2007
C
C *********************** INPUT **************************************
C integer maxnvl          : First dimension of VsL and maximum number
C                           of nonzero elements in VsL
C integer maxnv           : First dimension of Vs and maximum number
C                           of nonzero elements in Vs
C integer numVs(nuo)      : Number of non-zero elements in row of Vs
C integer listVsptr(nuo)  : Pointer to start of rows in Vs
C integer nuo             : Local no. of orbitals in unit cell
C integer nspin           : Number of spin components
C real*8  VsL(maxnvl,nspin) : Mesh contribution to be added to Vs
C ******************** INPUT AND OUTPUT *******************************
C real*8  Vs(maxnv,nspin) : Value of nonzero elements in each row of Vs
C                           to which the potential matrix elements are
C                           summed up
C *********************************************************************

C  Modules
      use precision
      use alloc, only: re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
      use parallel,     only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb,
     &                         LocalToGlobalOrb
      use m_dscfComm
#endif

      implicit none

C Argument types and dimensions
      integer
     &   maxnv, maxnvl, nspin, nuo, numVs(nuo), 
     &   listVsptr(nuo)
      real(dp)
     &   Vs(maxnv,nspin), VsL(maxnvl,nspin)

C Internal variables and arrays
      integer
     &  i, iu, ispin

#ifdef MPI
      integer           :: ii, MPIerror, Status(MPI_Status_Size), nn,
     &                     mm, io, iio
      real(dp), pointer :: buffer(:)
#endif

C***********************
C  Parallel execution  *
C***********************
#ifdef MPI
C     Allocate a buffer to Send/Receive data
      nullify(buffer)
      call re_alloc( buffer, 1, DCmaxnd*nspin, 'buffer', 'meshdscf' )

C     Copy the data that is in the current process from VsL to Vs
      do ii= 1, DCself
        io = DCinvp(ii)
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        iu = NeedDscfL(io)
        do ispin = 1,nspin
          Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio),ispin) = 
     &    Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio),ispin) + 
     &    VsL(listdlptr(iu)+1:listdlptr(iu)+numdl(iu),ispin)
        enddo
      enddo
      nn = DCself

C     Send or receive the data from/to other processes
      do i= 1, DCsize
C       The communication is from dst to src
        if (DCsrc(i).eq.Node) then
          mm = 0
          do ii= 1, DCsic(i)
            io = DCinvp(nn+ii)
            call GlobalToLocalOrb(io,Node,Nodes,iio)
            mm = mm + numVs(iio)
          enddo
          mm = mm*nspin

          call mpi_recv( buffer, mm, MPI_double_precision, DCdst(i), 1,
     &                   MPI_Comm_world, Status, MPIerror )
          mm = 0
          do ispin = 1,nspin
            do ii= 1, DCsic(i)
              io = DCinvp(nn+ii)
              call GlobalToLocalOrb(io,Node,Nodes,iio)
              Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio),ispin) = 
     &        Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio),ispin) + 
     &        buffer(mm+1:mm+numVs(iio))
              mm = mm + numVs(iio)
            enddo
          enddo
        else
          mm = 0
          do ispin = 1,nspin
            do ii= 1, DCsic(i)
              io = DCinvp(nn+ii)
              iio = NeedDscfL(io)
              buffer(mm+1:mm+numdl(iio)) =
     &        VsL(listdlptr(iio)+1:listdlptr(iio)+numdl(iio),ispin)
              mm = mm + numdl(iio)
            enddo
          enddo
          call MPI_Send( buffer, mm, MPI_double_precision, DCsrc(i), 1,
     &                   MPI_Comm_World, MPIerror )
        endif
        nn = nn + DCsic(i)
      enddo
      call de_alloc( buffer, 'buffer', 'meshdscf' )
#else
C*********************
C  Serial execution  *
C*********************
C Add those elements that are needed locally to the values already
C stored in the orbital oriented array
      do ispin = 1,nspin
        do i = 1,nuo
          iu = NeedDscfL(i)
          if (iu.gt.0) then
            Vs(listVsptr(i)+1:listVsptr(i)+numVs(i),ispin) = 
     &        Vs(listVsptr(i)+1:listVsptr(i)+numVs(i),ispin) + 
     &        VsL(listdlptr(iu)+1:listdlptr(iu)+numVs(i),ispin)
          endif
        enddo
      enddo
#endif
      end subroutine matrixMtoO

      subroutine matrixMtoOC( maxnvl, maxnv, numVs, listVsptr, nuo, 
     &                        VsL, Vs )
C ********************************************************************
C Transforms a matrix which is distributed by mesh points to a matrix
C that is distributed by a block cyclic distribution over the orbitals
C and adds it to an existing array of this form.
C Created by J.D.Gale, February 2000
C
C Update: All-to-all communications has been substituted by
C point-to-point communications. These have been precomputed in 
C dscfComm.
C Written by Rogeli Grima (BSC) Dec.2007
C Generalized by Javier Junquera in May. 2012 to the case where the 
C potential is a complex variable.
C
C *********************** INPUT **************************************
C integer maxnvl          : First dimension of VsL and maximum number
C                           of nonzero elements in VsL
C integer maxnv           : First dimension of Vs and maximum number
C                           of nonzero elements in Vs
C integer numVs(nuo)      : Number of non-zero elements in row of Vs
C integer listVsptr(nuo)  : Pointer to start of rows in Vs
C integer nuo             : Local no. of orbitals in unit cell
C real*8  VsL(maxnvl,2)   : Mesh contribution to be added to Vs
C ******************** INPUT AND OUTPUT *******************************
C complex*8  Vs(maxnv)    : Value of nonzero elements in each row of Vs
C                           to which the matrix elements are
C                           summed up
C *********************************************************************

C  Modules
      use precision
      use alloc, only: re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
      use parallel,     only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb,
     &                         LocalToGlobalOrb
      use m_dscfComm
#endif

      implicit none

C Argument types and dimensions
      integer
     &   maxnv, maxnvl, nuo, numVs(nuo), 
     &   listVsptr(nuo)
      real(dp)
     &   VsL(maxnvl,2)
      complex(dp)
     &   Vs(maxnv)

C Internal variables and arrays
      integer
     &  i, iu, irelim

#ifdef MPI
      integer           :: ii, MPIerror, Status(MPI_Status_Size), nn,
     &                     mm, io, iio
      real(dp), pointer :: buffer(:)
#endif

C***********************
C  Parallel execution  *
C***********************
#ifdef MPI
C     Allocate a buffer to Send/Receive data
      nullify(buffer)
      call re_alloc( buffer, 1, DCmaxnd*2, 'buffer', 'meshdscf' )

C     Copy the data that is in the current process from VsL to Vs
      do ii= 1, DCself
        io = DCinvp(ii)
        call GlobalToLocalOrb(io,Node,Nodes,iio)
        iu = NeedDscfL(io)
        Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio)) = 
     &    Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio)) + 
     &    cmplx(VsL(listdlptr(iu)+1:listdlptr(iu)+numdl(iu),1),
     &          VsL(listdlptr(iu)+1:listdlptr(iu)+numdl(iu),2), kind=dp)
      enddo
      nn = DCself

C     Send or receive the data from/to other processes
      do i= 1, DCsize
C       The communication is from dst to src
        if (DCsrc(i).eq.Node) then
          mm = 0
          do ii= 1, DCsic(i)
            io = DCinvp(nn+ii)
            call GlobalToLocalOrb(io,Node,Nodes,iio)
            mm = mm + numVs(iio)
          enddo
          mm = mm*2

          call mpi_recv( buffer, mm, MPI_double_precision, DCdst(i), 1,
     &                   MPI_Comm_world, Status, MPIerror )
          mm = 0
          do irelim = 1, 2
            do ii= 1, DCsic(i)
              io = DCinvp(nn+ii)
              call GlobalToLocalOrb(io,Node,Nodes,iio)
              if( irelim .eq. 1) then
                Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio)) = 
     &          Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio)) + 
     &          buffer(mm+1:mm+numVs(iio)) * cmplx(1.0,0.0,kind=dp)
              elseif( irelim .eq. 2) then
                Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio)) = 
     &          Vs(listVsptr(iio)+1:listVsptr(iio)+numVs(iio)) + 
     &          buffer(mm+1:mm+numVs(iio)) * cmplx(0.0,1.0,kind=dp)
              endif
              mm = mm + numVs(iio)
            enddo
          enddo
        else
          mm = 0
          do irelim = 1, 2
            do ii= 1, DCsic(i)
              io = DCinvp(nn+ii)
              iio = NeedDscfL(io)
              buffer(mm+1:mm+numdl(iio)) =
     &        VsL(listdlptr(iio)+1:listdlptr(iio)+numdl(iio),irelim)
              mm = mm + numdl(iio)
            enddo
          enddo
          call MPI_Send( buffer, mm, MPI_double_precision, DCsrc(i), 1,
     &                   MPI_Comm_World, MPIerror )
        endif
        nn = nn + DCsic(i)
      enddo
      call de_alloc( buffer, 'buffer', 'meshdscf' )
#else
C*********************
C  Serial execution  *
C*********************
C Add those elements that are needed locally to the values already
C stored in the orbital oriented array
      do i = 1,nuo
        iu = NeedDscfL(i)
        if (iu.gt.0) then
          Vs(listVsptr(i)+1:listVsptr(i)+numVs(i)) = 
     &      Vs(listVsptr(i)+1:listVsptr(i)+numVs(i)) + 
     &      cmplx(VsL(listdlptr(iu)+1:listdlptr(iu)+numVs(i),1),
     &            VsL(listdlptr(iu)+1:listdlptr(iu)+numVs(i),2),kind=dp)
        endif
      enddo
#endif
      end subroutine matrixMtoOC

      end module meshdscf

