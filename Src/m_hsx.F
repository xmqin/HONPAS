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
      module m_hsx
      public :: write_hsx
      CONTAINS
      subroutine write_hsx( gamma, no_u, no_s, nspin, indxuo, maxnh, 
     .                 numh, listhptr, listh, H, S, qtot, temp, xij)
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
C integer maxnh               : First dimension of listh, H, S and
C                               second of xij
C integer numh(nuo)           : Number of nonzero elements of each row
C                               of hamiltonian matrix
C integer listhptr(nuo)       : Pointer to the start of each row (-1)
C                               of hamiltonian matrix
C integer listh(maxnh)        : Nonzero hamiltonian-matrix element column
C                               indexes for each matrix row
C real*8  H(maxnh,nspin)      : Hamiltonian in sparse form
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

      character(len=label_length+3) :: paste
      logical           gamma
      integer           maxnh, no_u, no_s, nspin
      integer           indxuo(no_s), listh(maxnh), numh(*), listhptr(*)
      real(dp)          H(maxnh,nspin), S(maxnh),
     .                  qtot, temp, xij(3,maxnh)

      external          io_assign, io_close, paste

C Internal variables and arrays
      character(len=label_length+3) :: fname
      integer    im, is, iu, ju, k, mnh, ns, ia, io
      integer    ih,hl,nuo,maxnhtot,maxhg
      integer, dimension(:), allocatable :: numhg
#ifdef MPI
      integer    MPIerror, Request, Status(MPI_Status_Size), BNode
      integer,  dimension(:),   allocatable :: ibuffer
      real(dp), dimension(:),   allocatable :: buffer
      real(dp), dimension(:,:), allocatable :: buffer2
#endif
      logical    baddim, frstme, gammaonfile, write_xijk
      save       frstme, fname
      data frstme /.true./

C Find name of file
      if (frstme) then
        fname = paste( slabel, '.HSX' )
        frstme = .false.
      endif

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
          open( iu, file=fname, form='unformatted', status='unknown' )      

C Write overall data
          write(iu) no_u, no_s, nspin, maxnhtot

C Read logical
          write(iu) gamma

C Allocate local array for global numh
          allocate(numhg(no_u))
          call memory('A','I',no_u,'iohs')

C Write out indxuo
          if (.not.gamma) then
            write(iu) (indxuo(ih),ih=1,no_s)
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
          write(iu) (numhg(ih),ih=1,no_u)
#ifdef MPI
          allocate(buffer(maxhg))
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
            write(iu) (listh(listhptr(hl)+im),im = 1,numh(hl))
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
               write(iu) (ibuffer(im),im = 1,numhg(ih))
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


C Write Hamiltonian
        do is=1,nspin
          do ih=1,no_u
#ifdef MPI
            call WhichNodeOrb(ih,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
              hl = ih
#endif
              write(iu) (real(H(listhptr(hl)+im,is),kind=sp),
     $             im=1,numh(hl))
#ifdef MPI
            elseif (Node.eq.0) then
              call MPI_IRecv(buffer,numhg(ih),MPI_double_precision,
     .          BNode,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(ih,Node,Nodes,hl)
              call MPI_ISend(H(listhptr(hl)+1,is),numh(hl),
     .          MPI_double_precision,0,1,MPI_Comm_World,
     .          Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
              if (Node.eq.0) then
                 write(iu) (real(buffer(im),kind=sp),im=1,numhg(ih))
              endif
            endif
#endif
          enddo
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
            write(iu) (real(S(listhptr(hl)+im),kind=sp),im = 1,numh(hl))
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
               write(iu) (real(buffer(im),kind=sp),im=1,numhg(ih))
            endif
          endif
#endif
        enddo

#ifdef MPI
          if (Node .eq. 0) then
C Free buffer array
             call memory('D','D',size(buffer),'iohs')
             deallocate(buffer)
          endif
#endif

        if (Node.eq.0) then
          write(iu) qtot,temp
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
              write(iu) ((real(xij(k,listhptr(hl)+im),kind=sp),
     $                   k=1,3),im =1,numh(hl))
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
                 write(iu) ((real(buffer2(k,im),kind=sp),
     $                k=1,3),im=1,numhg(ih))
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
           write(iu) nspecies
           write(iu) (labelfis(is), zvalfis(is),nofis(is),is=1,nspecies)
           do is = 1, nspecies
              do io=1,nofis(is)
                 write(iu) cnfigfio(is,io), lofio(is,io), zetafio(is,io)
              enddo
           enddo
           write(iu) na_u
           write(iu) (isa(ia),ia=1,na_u)
           write(iu) (iaorb(io), iphorb(io), io=1,no_u)

C Deallocate local array for global numh
          call memory('D','I',size(numhg),'iohs')
          deallocate(numhg)
C Close file
          call io_close( iu )
        endif

      end subroutine write_hsx
      end module m_hsx
