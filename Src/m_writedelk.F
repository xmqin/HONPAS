! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_writedelk
      public :: write_delk
      CONTAINS
      subroutine write_delk( gamma, no_u, no_s, nspin, indxuo, maxnh, 
     .             numh, listhptr, listh, delkmat, S, qtot, temp, xij)
C *********************************************************************
C Saves the hamiltonian and overlap matrices, and other data required
C to obtain the bands and density of states
C Writen by J.Soler July 1997.
C Note because of the new more compact method of storing H and S
C this routine is NOT backwards compatible
C *************************** INPUT **********************************
C logical       gamma         : Is only gamma point used?
C ******************** INPUT or OUTPUT (depending on task) ***********
C integer no_u                : Number of basis orbitals per unit cell
C integer no_s                : Number of basis orbitals per supercell
C integer nspin               : Spin polarization (1 or 2)
C integer indxuo(no_s)        : Index of orbitals in supercell
C integer maxnh               : First dimension of listh, delkmat, S and
C                               second of xij
C integer numh(nuo)           : Number of nonzero elements of each row
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to the start of each row (-1)
C                               of hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element column
C                               indexes for each matrix row
C real*8  delkmat(maxnh)      : Hamiltonian in sparse form
C real*8  S(maxnh)            : Overlap in sparse form
C real*8  qtot                : Total number of electrons
C real*8  temp                : Electronic temperature for Fermi smearing
C real*8  xij(3,maxnh)        : Vectors between orbital centers (sparse)
C                               (not read/written if only gamma point)

C
C  Modules
C
      use precision, only: dp, sp
      use parallel,     only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, LocalToGlobalOrb,
     .                         GlobalToLocalOrb, GetNodeOrbs
      use atm_types,    only : nspecies
      use atomlist,     only : iphorb, iaorb
      use siesta_geom,  only : na_u, isa
      use atmfuncs,     only : nofis, labelfis, zvalfis
      use atmfuncs,     only : cnfigfio, lofio, zetafio
      use fdf
      use files,        only : slabel, label_length
      use sys,          only : die
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

      logical           gamma
      integer           maxnh, no_u, no_s, nspin
      integer           indxuo(no_s), listh(maxnh), numh(*), listhptr(*)
      real(dp)          S(maxnh),
     .                  qtot, temp, xij(3,maxnh)
      complex(dp)       delkmat(maxnh)

      external          io_assign, io_close

C Internal variables and arrays
      character(len=label_length+4) :: fname
      integer    im, is, iu, ju, k, mnh, ns, ia, io
      integer    ih,hl,nuo,maxnhtot,maxhg
      integer, dimension(:), allocatable :: numhg
#ifdef MPI
      integer    MPIerror, Request, Status(MPI_Status_Size), BNode
      integer,  dimension(:),   allocatable :: ibuffer
      real(dp), dimension(:),   allocatable :: buffer
      complex(dp), dimension(:),   allocatable :: buffer3
      real(dp), dimension(:,:), allocatable :: buffer2
#endif
      logical    baddim, gammaonfile, write_xijk

C Find name of file
      fname = trim(slabel) // '.dek'

C Find total numbers over all Nodes
#ifdef MPI
      call MPI_AllReduce(maxnh,maxnhtot,1,MPI_integer,MPI_sum,
     .  MPI_Comm_World,MPIerror)
#else
      maxnhtot = maxnh
#endif

        if (Node.eq.0) then
C Open file
          call io_assign( iu )
          open( iu, file=fname, form='formatted', status='unknown' )      

C Write overall data
          write(iu,*) " no_u, no_s, nspin, maxnhtot = ",
     .      no_u, no_s, nspin, maxnhtot

C Read logical
          write(iu,*) " gamma = ", gamma

C Allocate local array for global numh
          allocate(numhg(no_u))
          call memory('A','I',no_u,'iohs')

C Write out indxuo
          if (.not.gamma) then
            do ih = 1, no_s
              write(iu,*)" ih, indxuo(ih) = ", ih, indxuo(ih)
            enddo
          endif

        endif

C Create globalised numh
        do ih = 1,no_u
#ifdef MPI
          call WhichNodeOrb(ih,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
            hl = ih
#endif
            numhg(ih) = numh(hl)
#ifdef MPI
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            call MPI_ISend(numh(hl),1,MPI_integer,
     .        0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.0) then
            call MPI_IRecv(numhg(ih),1,MPI_integer,
     .        BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
          endif
#endif
        enddo

        if (Node.eq.0) then
C Write numh
          maxhg = 0
          do ih = 1,no_u
            maxhg = max(maxhg,numhg(ih))
          enddo
          do ih = 1, no_u
            write(iu,*)" ih, numhg(ih) = ", ih, numhg(ih)
          enddo 
#ifdef MPI
          allocate(buffer(maxhg))
          call memory('A','D',maxhg,'iohs')
          allocate(buffer3(maxhg))
          call memory('A','D',maxhg,'iohs')
          allocate(ibuffer(maxhg))
          call memory('A','I',maxhg,'iohs')
#endif
        endif

C Write listh
        do ih = 1,no_u
#ifdef MPI
          call WhichNodeOrb(ih,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
            hl = ih
#endif
            do im = 1, numh(hl)
              write(iu,*) ' im, listh(listhptr(hl)+im) = ',
     .          im, listh(listhptr(hl)+im)
            enddo
#ifdef MPI
          elseif (Node.eq.0) then
            call MPI_IRecv(ibuffer,numhg(ih),MPI_integer,BNode,1,
     .        MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            call MPI_ISend(listh(listhptr(hl)+1),numh(hl),MPI_integer,
     .        0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
            if (Node.eq.0) then
               do im = 1, numhg(ih)
                 write(iu,*) " im, ibuffer(im) = ",
     .             im, ibuffer(im)
               enddo
            endif
          endif
#endif
        enddo

#ifdef MPI
        if (Node.eq.0) then
          call memory('D','I',size(ibuffer),'iohs')
          deallocate(ibuffer)
        endif
#endif


C Write delkmat
        do ih=1,no_u
#ifdef MPI
          call WhichNodeOrb(ih,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
            hl = ih
#endif
            do im = 1, numh(hl)
!              write(iu,'(a,2i5,2f12.5)')
!     .          ' im, delkmat(listhptr(hl)+im) = ',
!     .          im, listhptr(hl)+im, delkmat(listhptr(hl)+im)
              write(iu,'(2f12.5)')
     .           delkmat(listhptr(hl)+im)
            enddo
#ifdef MPI
          elseif (Node.eq.0) then
            call MPI_IRecv(buffer3,numhg(ih),MPI_double_complex,
     .        BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            call MPI_ISend(delkmat(listhptr(hl)+1),numh(hl),
     .        MPI_double_complex,0,1,MPI_Comm_World,
     .        Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
            if (Node.eq.0) then
               do im = 1, numhg(ih)
!                 write(iu,'(a,i5,2f12.5)') ' im, buffer3(im) = ',
!     .             im, buffer3(im)
                 write(iu,'(2f12.5)') buffer3(im)
               enddo 
            endif
          endif
#endif
        enddo

C Write Overlap matrix
        do ih = 1,no_u
#ifdef MPI
          call WhichNodeOrb(ih,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
            hl = ih
#endif
            do im = 1, numh(hl)
              write(iu,'(a,2i5,2f12.5)') ' im, S(listhptr(hl)+im) = ',
     .          im, listhptr(hl)+im, (real(S(listhptr(hl)+im),kind=sp))
            enddo 
#ifdef MPI
          elseif (Node.eq.0) then
            call MPI_IRecv(buffer,numhg(ih),MPI_double_precision,
     .        BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
            call MPI_ISend(S(listhptr(hl)+1),numh(hl),
     .        MPI_double_precision,0,1,MPI_Comm_World,
     .        Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
            if (Node.eq.0) then
               do im = 1, numhg(ih)
                 write(iu,'(a,i5,2f12.5)')
     .             ' im, (real(buffer(im),kind=sp) = ',
     .             im, (real(buffer(im),kind=sp))
               enddo 
            endif
          endif
#endif
        enddo

#ifdef MPI
          if (Node .eq. 0) then
C Free buffer array
             call memory('D','D',size(buffer),'iohs')
             deallocate(buffer)
             call memory('D','D',size(buffer3),'iohs')
             deallocate(buffer3)
          endif
#endif

        if (Node.eq.0) then
          write(iu,*) ' qtot,temp = ', qtot,temp
        endif

        write_xijk = .TRUE.

        if (write_xijk) then
#ifdef MPI
C Allocate buffer array
          if (Node .eq. 0) then
             allocate(buffer2(3,maxhg))
             call memory('A','D',3*maxhg,'iohs')
          endif
#endif
          do ih = 1,no_u
#ifdef MPI
            call WhichNodeOrb(ih,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
              hl = ih
#endif
              do im = 1, numh(hl)
                write(iu,*) im, (real(xij(k,listhptr(hl)+im),kind=sp),
     $                   k=1,3)
              enddo 
#ifdef MPI
            elseif (Node.eq.0) then
              call MPI_IRecv(buffer2(1,1),3*numhg(ih),
     .          MPI_double_precision,BNode,1,MPI_Comm_World,
     .          Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
              call MPI_ISend(xij(1,listhptr(hl)+1),3*numh(hl),
     .          MPI_double_precision,0,1,MPI_Comm_World,
     .          Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
              if (Node.eq.0) then
                 do im = 1, numhg(ih)
                 write(iu,*) im, (real(buffer2(k,im),kind=sp),
     $                k=1,3)
                 enddo 
              endif
            endif
#endif
          enddo
#ifdef MPI
          if (Node .eq. 0) then
C Free buffer array
             call memory('D','D',size(buffer2),'iohs')
             deallocate(buffer2)
          endif
#endif
        endif   ! write_xijk

        if (Node.eq.0) then

!
!       Write other useful info
!
          write(iu,*)' nspecies = ', nspecies
          write(iu,*)(labelfis(is), zvalfis(is),nofis(is),is=1,nspecies)
           do is = 1, nspecies
              do io=1,nofis(is)
                 write(iu,*)cnfigfio(is,io),lofio(is,io),zetafio(is,io)
              enddo
           enddo
           write(iu,*) na_u
           write(iu,*) (isa(ia),ia=1,na_u)
           write(iu,*) (iaorb(io), iphorb(io), io=1,no_u)

C Deallocate local array for global numh
          call memory('D','I',size(numhg),'iohs')
          deallocate(numhg)
C Close file
          call io_close( iu )
        endif

      end subroutine write_delk
      end module m_writedelk
