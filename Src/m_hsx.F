! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_hsx
      public :: write_hsx
      CONTAINS
      subroutine write_hsx( no_u, no_s, nspin, indxuo, maxnh, 
     .                 numh, listhptr, listh, H, S, qtot, temp, xij)
C *********************************************************************
C Saves the hamiltonian and overlap matrices, and other data required
C to obtain the bands and density of states
C Writen by J.Soler July 1997.
C Note because of the new more compact method of storing H and S
C this routine is NOT backwards compatible
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

      ! This really means that supercell info is supported by default
      logical, parameter ::  gamma = .false. 
      logical, parameter ::  write_xijk = .true.

      integer           maxnh, no_u, no_s, nspin
      integer           indxuo(no_s), listh(maxnh), numh(*), listhptr(*)
      real(dp)          H(maxnh,nspin), S(maxnh),
     .                  qtot, temp, xij(3,maxnh)

      external          io_assign, io_close

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

      call timer("writeHSX",1)
      
C Find name of file
      fname = trim(slabel) // '.HSX'

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

C Write logical (always .false.)
          write(iu) gamma

C Allocate local array for global numh
          allocate(numhg(no_u))
          call memory('A','I',no_u,'iohs')

C Write out indxuo always
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

        if (write_xijk) then           ! Always
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

      call timer("writeHSX",2)

      end subroutine write_hsx

!-----------------------------------------------------------------
      subroutine write_hs_formatted(no_u, nspin,
     $                 maxnh, numh, listhptr, listh, H, S)
C *********************************************************************
C Saves the hamiltonian and overlap matrices in formatted form
C ONLY for gamma case
C ******************** INPUT
C integer no_u                : Number of basis orbitals per unit cell
C integer nspin               : Spin polarization (1 or 2)
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

C
C  Modules
C
      use precision, only: dp, sp
      use parallel,     only : Node, Nodes
      use parallelsubs, only : WhichNodeOrb, LocalToGlobalOrb,
     .                         GlobalToLocalOrb, GetNodeOrbs
      use atm_types,    only : nspecies
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

      integer           maxnh, no_u, nspin
      integer           listh(maxnh), numh(*), listhptr(*)
      real(dp)          H(maxnh,nspin), S(maxnh)
      external          io_assign, io_close


      integer    im, is, iu, ju, k, mnh, ns, ia, io
      integer    ih,hl,nuo,maxnhtot,maxhg
      integer, dimension(:), allocatable :: numhg, hg_ptr
#ifdef MPI
      integer    MPIerror, Request, Status(MPI_Status_Size), BNode
      integer,  dimension(:),   allocatable :: ibuffer
      real(dp), dimension(:),   allocatable :: buffer
#endif


      call timer("write_HS_fmt",1)

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
          open( iu, file="S.matrix", form='formatted', status='unknown',
     $          position="rewind")      

C Write overall data
          write(iu,*) no_u, no_u, maxnhtot

C Allocate local array for global numh
          allocate(numhg(no_u))
          allocate(hg_ptr(no_u+1))
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
C Write row pointers
         maxhg = 0
         hg_ptr(1) = 1
         do ih = 1,no_u
           maxhg = max(maxhg,numhg(ih))
           hg_ptr(ih+1) = hg_ptr(ih) + numhg(ih)
         enddo
         write(iu,*) (hg_ptr(ih),ih=1,no_u+1)
         deallocate(hg_ptr)
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
            write(iu,*) (listh(listhptr(hl)+im),im = 1,numh(hl))
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
               write(iu,*) (ibuffer(im),im = 1,numhg(ih))
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
        if (Node.eq.0) then
           call io_close(iu)
           call system("cp S.matrix H.matrix")
           call io_assign( iu )
           open( iu, file="H.matrix", form='formatted', status='old',
     $          action="write", position="append")      
        endif


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
              write(iu,*) (real(H(listhptr(hl)+im,is),kind=sp),
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
                 write(iu,*) (real(buffer(im),kind=sp),im=1,numhg(ih))
              endif
            endif
#endif
          enddo
        enddo


        if (node == 0) then
           call io_close(iu)
           call io_assign( iu )
           open( iu, file="S.matrix", form='formatted', status='old',
     $          action="write", position="append")      
        endif

C Write Overlap matrix
        do ih = 1,no_u
#ifdef MPI
          call WhichNodeOrb(ih,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(ih,Node,Nodes,hl)
#else
            hl = ih
#endif
            write(iu,*) (real(S(listhptr(hl)+im),kind=sp),
     $                   im = 1,numh(hl))
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
               write(iu,*) (real(buffer(im),kind=sp),im=1,numhg(ih))
            endif
          endif
#endif
        enddo

#ifdef MPI
        if (Node .eq. 0) then
           call memory('D','D',size(buffer),'iohs')
           deallocate(buffer)
           deallocate(numhg)   
        endif
#endif

        if (node == 0) then
           call io_close( iu )
        endif

      call timer("write_HS_fmt",2)

      end subroutine write_hs_formatted

      end module m_hsx
