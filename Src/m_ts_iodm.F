!     
!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996- .
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!    

! NOTE: Module adapted from m_iodm.F to handle the .TSDE file. The difference with
!       respect to the SIESTA .DM file is that, at the end of the file are added
!       the EDM sparse matrix and the Fermi Energy. In this way, the .TSDE file can
!       be used for example by Denchar in order to plot the TS density matrix.
! F.D. Novaes
      module m_ts_iodm

      use precision,     only: dp
      use parallel,     only : ionode,  Node, Nodes
      use parallelsubs, only : LocalToGlobalOrb
      use sys,          only : die
      use fdf,          only : fdf_boolean, fdf_string
      use files,        only : slabel, label_length
      use alloc,        only : re_alloc, de_alloc
      use m_mpi_utils,  only : broadcast
#ifdef MPI
      use mpi_siesta
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
#endif

      implicit none

      private
      PUBLIC :: write_ts_dm, read_ts_dm, read_dynamic_ts_dm, ts_init_dm

C     Saved internal variables:
      logical,           save :: frstme = .true., scndme = .false.
      character(len=label_length+3), save :: fnameu
      character(len=label_length+4), save :: fnamef
      character(len=label_length+4), save :: fnamei, fnameo
      logical,                       save :: fmti, fmto
      character(len=11),             save :: formin, formout

!     Character formats for formatted I/O
!     We assume that we have no integers greater than (1e10-1)
!     (which is bigger than the largest 32-bit integer)
!     and that floating point accuracy is sufficiently represented
!     by 16 decimal digits of precision (which suffices for 64-bit
!     IEEE floats)

      character(len=*), parameter :: intfmt = '(I11)'
      character(len=*), parameter :: floatfmt = '(ES22.14)'

      CONTAINS

      subroutine read_dynamic_ts_dm ( maxnd, no_l, nspin, numd, 
     .                     listdptr, listd, dm, edm, ef, found )

      ! The only "in" parameter is no_l, to determine the distribution
      integer, intent(in) :: no_l

      ! These are all "output" sparsity parameters
      ! We read the file info into a temporary
      ! set of arrays and change the structure if needed
      ! All pointers should be nullified on entry

      integer, intent(out)  :: maxnd
      integer, intent(out)  :: nspin

      integer, pointer      :: numd(:)
      integer, pointer      :: listdptr(:)
      integer, pointer      :: listd(:)
      real(dp), pointer     :: dm(:,:)
      real(dp), pointer     :: edm(:,:)
      real(dp), intent(out) :: ef

      logical, intent(out) :: found

C     Internal variables
      logical   exist3
      integer   im, is, unit1, m, nb
      integer   no_u, ml, ndmaxg
      integer, dimension(:), pointer  :: numdg
#ifdef MPI
      integer   MPIerror, Request, Status(MPI_Status_Size)
      integer   BNode

      real(dp), dimension(:), pointer :: buffer
      integer,  dimension(:), pointer :: ibuffer
#endif
      
      external          chkdim

      call ts_setup_file_modes()

C     Find file name

      if (Node.eq.0) then
         inquire (file=fnamei,        exist=exist3)
      endif
      call broadcast(exist3)

      found = .false.
      if ( .not. exist3) RETURN

C     Find total number of basis functions over all Nodes
#ifdef MPI
      call MPI_AllReduce(no_l,no_u,1,MPI_integer,MPI_sum,
     .     MPI_Comm_World,MPIerror)
#else
      no_u = no_l
#endif

      if (Node.eq.0) then
         write(6,'(/,a)') 'ts_iodm: Reading Density Matrix from files'
         call io_assign(unit1)
         open( unit1, file=fnamei, form=formin, status='old' )
         rewind(unit1)
         if (fmti) then
            read(unit1, intfmt) nb, nspin
         else
            read(unit1) nb, nspin
         endif
      endif

C     Communicate the values to all Nodes and adjust to allow for
C     distributed memory before checking the dimensions
#ifdef MPI
      call MPI_Bcast(nb,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(nspin,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

      ! If the total number of orbitals does not match, bail out
      if (no_u /= nb) then
         if (Node.eq.0) then
           write(6,"(a,i6,/,a)")
     $      "WARNING: Wrong number of orbitals in TSDE file: ",
     $       nb,
     $      "WARNING: Falling back to atomic initialization of DM."
           write(0,"(a,i6,/,a)")
     $      "WARNING: Wrong number of orbitals in TSDE file: ",
     $       nb,
     $      "WARNING: Falling back to atomic initialization of DM."
            call io_close(unit1)
         endif
         found = .false.
         RETURN
      endif

C     Allocate local buffer array for globalised numd
      nullify(numdg)
      call re_alloc(numdg, 1, no_u, name="numdg", routine="read_ts_dm")
      if (Node.eq.0) then
         if (fmti) then
            read(unit1, intfmt) (numdg(m),m=1,no_u)
         else
            read(unit1) (numdg(m),m=1,no_u)
         endif
      endif
      call broadcast(numdg(1:no_u))

C     Convert global numd pointer to local form and generate listdptr

      maxnd = 0
      do m = 1,no_l
         call LocalToGlobalOrb(m,Node,Nodes,ml)
         numd(m) = numdg(ml)
         maxnd = maxnd + numdg(ml)
         if (m .eq. 1) then
            listdptr(1) = 0
         else
            listdptr(m) = listdptr(m-1) + numd(m-1)
         endif
      enddo
      ndmaxg = maxval(numdg(1:no_u))

      ! Allocate density matrix, energy density matrix and listd
      call re_alloc(listd,1,maxnd,name="listd",routine="read_ts_dm")
      call re_alloc(dm,1,maxnd,1,nspin,name="dm",routine="read_ts_dm")
      call re_alloc(edm,1,maxnd,1,nspin,name="edm",routine="read_ts_dm")

#ifdef MPI
C     Create buffer arrays for transfering density matrix between nodes and lists
      nullify(buffer,ibuffer)
      call re_alloc(buffer,1,ndmaxg,name="buffer",routine="read_ts_dm")
      call re_alloc(ibuffer,1,ndmaxg,name="ibuffer",
     .                               routine="read_ts_dm")
#endif

      do m = 1,no_u
#ifdef MPI
         call WhichNodeOrb(m,Nodes,BNode)
         if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
            ml = m
#endif
            if (fmti) then
               read(unit1, intfmt) 
     .              (listd(listdptr(ml)+im),im=1,numd(ml))
            else
               read(unit1) (listd(listdptr(ml)+im),im=1,numd(ml))
            endif
#ifdef MPI
         elseif (Node.eq.0) then
            if (fmti) then
               read(unit1, intfmt) (ibuffer(im),im=1,numdg(m))
            else
               read(unit1) (ibuffer(im),im=1,numdg(m))
            endif
            call MPI_ISend(ibuffer,numdg(m),MPI_integer,
     .           BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            call MPI_IRecv(listd(listdptr(ml)+1),numd(ml),
     .           MPI_integer,0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         endif
         if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
         endif
#endif
      enddo

#ifdef MPI
      call de_alloc(ibuffer,name="ibuffer",routine="read_ts_dm")
#endif

! Read DM
      do is = 1,nspin
         do m = 1,no_u
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
               ml = m
#endif
               if (fmti) then
                  read(unit1, floatfmt)
     .                 (dm(listdptr(ml)+im,is),im=1,numd(ml))
               else
                  read(unit1) (dm(listdptr(ml)+im,is),im=1,numd(ml))
               endif
#ifdef MPI
            elseif (Node.eq.0) then
               if (fmti) then
                  read(unit1, floatfmt) (buffer(im),im=1,numdg(m))
               else
                  read(unit1) (buffer(im),im=1,numdg(m))
               endif
               call MPI_ISend(buffer,numdg(m),MPI_double_precision,
     .              BNode,1,MPI_Comm_World,Request,MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
               call MPI_IRecv(dm(listdptr(ml)+1,is),numd(ml),
     .              MPI_double_precision,0,1,MPI_Comm_World,Request,
     .              MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
               call MPI_Barrier(MPI_Comm_World,MPIerror)
            endif
#endif
         enddo
      enddo

! Read EDM
      do is = 1,nspin
         do m = 1,no_u
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
               ml = m
#endif
               if (fmti) then
                  read(unit1, floatfmt)
     .                 (edm(listdptr(ml)+im,is),im=1,numd(ml))
               else
                  read(unit1) (edm(listdptr(ml)+im,is),im=1,numd(ml))
               endif
#ifdef MPI
            elseif (Node.eq.0) then
               if (fmti) then
                  read(unit1, floatfmt) (buffer(im),im=1,numdg(m))
               else
                  read(unit1) (buffer(im),im=1,numdg(m))
               endif
               call MPI_ISend(buffer,numdg(m),MPI_double_precision,
     .              BNode,1,MPI_Comm_World,Request,MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
               call MPI_IRecv(edm(listdptr(ml)+1,is),numd(ml),
     .              MPI_double_precision,0,1,MPI_Comm_World,Request,
     .              MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
               call MPI_Barrier(MPI_Comm_World,MPIerror)
            endif
#endif
         enddo
      enddo

#ifdef MPI
      call de_alloc(buffer,name="buffer",routine="read_ts_dm")
#endif
      call de_alloc(numdg,name="numdg",routine="read_ts_dm")

! Read EF and close unit1
      if (Node.eq.0) then
         read(unit1) ef
         call io_close(unit1)
      endif

#ifdef MPI
      call MPI_BCast(ef,1,MPI_double_precision,0,MPI_Comm_World,
     .              MPIerror)
#endif

      found = .true.

      end subroutine read_dynamic_ts_dm
!-----------------------      
      subroutine read_ts_dm ( maxnd, no_l, nspin, numd, 
     .                     listdptr, listd, dm, edm, ef, found )

C     Reads density matrix from file
C     Written by P.Ordejon and J.M.Soler. May 1997.
C     Made into a separate routine by Alberto Garcia, August 2007
C     ********* INPUT ***************************************************
C     integer   maxnd    : First dimension of listd and dm
C     integer   no_l   : Number of atomic orbitals
C     integer   nspin    : Number of spins (1 or 2)
C     ********* OUTPUT ************
C     integer numd(no_l)     : Control vector of DM matrix
C     (number of nonzero elements of each row)
C     integer listdptr(no_l) : Control vector of DM matrix
C     (pointer to the start of each row)
C     integer listd(maxnd)     : Control vector of DM matrix
C     (list of nonzero elements of each row)
C     real*8  dm(maxnd,nspin)  : Density matrix
C     logical found : Has DM been found in disk? 
C     *******************************************

      ! These are the "input" sparsity parameters
      ! We might want to read the file info into a temporary
      ! set of arrays and change the structure if needed

      integer, intent(in) :: maxnd
      integer, intent(in) :: no_l
      integer, intent(in) :: nspin

      integer, intent(out) :: numd(no_l)
      integer, intent(out) :: listdptr(no_l)
      integer, intent(out) :: listd(maxnd)
      real(dp), intent(out) :: dm(maxnd, nspin)
      real(dp), intent(out) :: edm(maxnd, nspin)
      real(dp), intent(out) :: ef

      logical, intent(out) :: found

C     Internal variables
      logical   exist3
      integer   im, is, unit1, unit2, m, nb, ndmax, ns
      integer   no_u, ml, ndmaxg
      integer, dimension(:), pointer  :: numdg
#ifdef MPI
      integer   MPIerror, Request, Status(MPI_Status_Size)
      integer   BNode

      real(dp), dimension(:), pointer :: buffer
      integer,  dimension(:), pointer :: ibuffer
#endif
      
      external          chkdim

      call ts_setup_file_modes()

C     Find file name

      if (Node.eq.0) then
         inquire (file=fnamei,        exist=exist3)
      endif
      call broadcast(exist3)

      found = .false.
      if ( .not. exist3) RETURN

C     Find total number of basis functions over all Nodes
#ifdef MPI
      call MPI_AllReduce(no_l,no_u,1,MPI_integer,MPI_sum,
     .     MPI_Comm_World,MPIerror)
#else
      no_u = no_l
#endif

      if (Node.eq.0) then
         write(6,'(/,a)') 'ts_iodm: Reading Density Matrix from files'
         call io_assign(unit1)
         open( unit1, file=fnamei, form=formin, status='old' )
         rewind(unit1)
         if (fmti) then
            read(unit1, intfmt) nb, ns
         else
            read(unit1) nb, ns
         endif
      endif

C     Communicate the values to all Nodes and adjust to allow for
C     distributed memory before checking the dimensions
#ifdef MPI
      call MPI_Bcast(nb,1,MPI_integer,0,MPI_Comm_World,MPIerror)
      call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

      ! This checks for equality of spin and number of rows.
      call chkdim( 'iodm', 'no_u', no_u, nb, 0 )
      call chkdim( 'iodm', 'nspin',  nspin,  ns, 0 )

C     Allocate local buffer array for globalised numd
      nullify(numdg)
      call re_alloc(numdg, 1, no_u, name="numdg", routine="read_ts_dm")
      if (Node.eq.0) then
         if (fmti) then
            read(unit1, intfmt) (numdg(m),m=1,no_u)
         else
            read(unit1) (numdg(m),m=1,no_u)
         endif
      endif
      call broadcast(numdg(1:no_u))

C     Convert global numd pointer to local form and generate listdptr
!      nullify(numd_tmp,listhptr_tmp)
!      call re_alloc(numd_tmp,1,no_l,name="numd_tmp")
!      call re_alloc(listhptr_tmp,1,no_l,name="listhptr_tmp")
      ndmax = 0
      do m = 1,no_l
         call LocalToGlobalOrb(m,Node,Nodes,ml)
         numd(m) = numdg(ml)
!         numd_tmp(m) = numdg(ml)
         ndmax = ndmax + numdg(ml)
         if (m .eq. 1) then
!            listdptr_tmp(1) = 0
            listdptr(1) = 0
         else
!            listdptr_tmp(m) = listdptr_tmp(m-1) + numd_tmp(m-1)
            listdptr(m) = listdptr(m-1) + numd(m-1)
         endif
      enddo
      ndmaxg = maxval(numdg(1:no_u))

C     Check size of first dimension of dm
!      nullify(dm_tmp)
!      call re_alloc(dm_tmp,1,ndmax,name="dm_local")
      call chkdim( 'iodm', 'maxnd', maxnd, ndmax, 1 )

#ifdef MPI
C     Create buffer arrays for transfering density matrix between nodes and lists
      nullify(buffer,ibuffer)
      call re_alloc(buffer,1,ndmaxg,name="buffer",routine="read_ts_dm")
      call re_alloc(ibuffer,1,ndmaxg,name="ibuffer",
     .                               routine="read_ts_dm")
#endif

      do m = 1,no_u
#ifdef MPI
         call WhichNodeOrb(m,Nodes,BNode)
         if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
            ml = m
#endif
            if (fmti) then
               read(unit1, intfmt) 
     .              (listd(listdptr(ml)+im),im=1,numd(ml))
            else
               read(unit1) (listd(listdptr(ml)+im),im=1,numd(ml))
            endif
#ifdef MPI
         elseif (Node.eq.0) then
            if (fmti) then
               read(unit1, intfmt) (ibuffer(im),im=1,numdg(m))
            else
               read(unit1) (ibuffer(im),im=1,numdg(m))
            endif
            call MPI_ISend(ibuffer,numdg(m),MPI_integer,
     .           BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            call MPI_IRecv(listd(listdptr(ml)+1),numd(ml),
     .           MPI_integer,0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         endif
         if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
         endif
#endif
      enddo

#ifdef MPI
      call de_alloc(ibuffer,name="ibuffer",routine="read_ts_dm")
#endif

! Read DM
      do is = 1,nspin
         do m = 1,no_u
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
               ml = m
#endif
               if (fmti) then
                  read(unit1, floatfmt)
     .                 (dm(listdptr(ml)+im,is),im=1,numd(ml))
               else
                  read(unit1) (dm(listdptr(ml)+im,is),im=1,numd(ml))
               endif
#ifdef MPI
            elseif (Node.eq.0) then
               if (fmti) then
                  read(unit1, floatfmt) (buffer(im),im=1,numdg(m))
               else
                  read(unit1) (buffer(im),im=1,numdg(m))
               endif
               call MPI_ISend(buffer,numdg(m),MPI_double_precision,
     .              BNode,1,MPI_Comm_World,Request,MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
               call MPI_IRecv(dm(listdptr(ml)+1,is),numd(ml),
     .              MPI_double_precision,0,1,MPI_Comm_World,Request,
     .              MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
               call MPI_Barrier(MPI_Comm_World,MPIerror)
            endif
#endif
         enddo
      enddo

! Read EDM
      do is = 1,nspin
         do m = 1,no_u
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
               ml = m
#endif
               if (fmti) then
                  read(unit1, floatfmt)
     .                 (edm(listdptr(ml)+im,is),im=1,numd(ml))
               else
                  read(unit1) (edm(listdptr(ml)+im,is),im=1,numd(ml))
               endif
#ifdef MPI
            elseif (Node.eq.0) then
               if (fmti) then
                  read(unit1, floatfmt) (buffer(im),im=1,numdg(m))
               else
                  read(unit1) (buffer(im),im=1,numdg(m))
               endif
               call MPI_ISend(buffer,numdg(m),MPI_double_precision,
     .              BNode,1,MPI_Comm_World,Request,MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
               call MPI_IRecv(edm(listdptr(ml)+1,is),numd(ml),
     .              MPI_double_precision,0,1,MPI_Comm_World,Request,
     .              MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
               call MPI_Barrier(MPI_Comm_World,MPIerror)
            endif
#endif
         enddo
      enddo

#ifdef MPI
      call de_alloc(buffer,name="buffer",routine="read_ts_dm")
#endif
      call de_alloc(numdg,name="numdg",routine="read_ts_dm")

! Read EF and close unit1
      if (Node.eq.0) then
         read(unit1) ef
         call io_close(unit1)
      endif

#ifdef MPI
      call MPI_BCast(ef,1,MPI_double_precision,0,MPI_Comm_World,
     .              MPIerror)
#endif

      found = .true.

      end subroutine read_ts_dm
!-----------------------      
      subroutine write_ts_dm (maxnd, no_l, nspin,
     $     numd, listdptr, listd, dm, edm, ef)

      integer, intent(in) :: maxnd
      integer, intent(in) :: no_l
      integer, intent(in) :: nspin
      integer, intent(in) :: numd(1:no_l)
      integer, intent(in) :: listdptr(1:no_l)
      integer, intent(in) :: listd(maxnd)
      real(dp), intent(in) :: dm(maxnd,nspin)
      real(dp), intent(in) :: edm(maxnd,nspin)
      real(dp), intent(in) :: ef

      integer :: no_u, m, ml, im, ndmaxg, unit1, is
#ifdef MPI
      integer   MPIerror, Request, Status(MPI_Status_Size)
      integer   BNode
      real(dp), dimension(:), pointer :: buffer
      integer,  dimension(:), pointer :: ibuffer
#endif
      integer, dimension(:), pointer  :: numdg

      call ts_setup_file_modes()

C     Find total number of basis functions over all Nodes
#ifdef MPI
      call MPI_AllReduce(no_l,no_u,1,MPI_integer,MPI_sum,
     .     MPI_Comm_World,MPIerror)
#else
      no_u = no_l
#endif

      if (Node.eq.0) then
         call io_assign(unit1)
         open( unit1, file=fnameo, form=formout, status='unknown' )
         rewind(unit1)
         if (fmto) then
            write(unit1, intfmt) no_u, nspin
         else
            write(unit1) no_u, nspin
         endif
      endif

      nullify(numdg)
      call re_alloc(numdg, 1, no_u, name="numdg", routine="read_ts_dm")

C     Create globalised numd
      do m = 1,no_u
#ifdef MPI
         call WhichNodeOrb(m,Nodes,BNode)
         if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
            ml = m
#endif
            numdg(m) = numd(ml)
#ifdef MPI
         elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            call MPI_ISend(numd(ml),1,MPI_integer,
     .           0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         elseif (Node.eq.0) then
            call MPI_IRecv(numdg(m),1,MPI_integer,
     .           BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         endif
         if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
         endif
#endif
      enddo

C     Write out numd array
      if (Node.eq.0) then
         ndmaxg = maxval(numdg(1:no_u))
         if (fmto) then
            write(unit1, intfmt) (numdg(m),m=1,no_u)
         else
            write(unit1) (numdg(m),m=1,no_u)
         endif
#ifdef MPI
         nullify(buffer,ibuffer)
         call re_alloc(buffer,1,ndmaxg,name="buffer",
     $                          routine="write_ts_dm")
         call re_alloc(ibuffer,1,ndmaxg,name="ibuffer",
     $                          routine="write_ts_dm")
#endif
      endif

C     Write out listd array
      do m = 1,no_u
#ifdef MPI
         call WhichNodeOrb(m,Nodes,BNode)
         if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
            ml = m
#endif
            if (fmto) then
               write(unit1, intfmt)
     .              (listd(listdptr(ml)+im),im=1,numd(ml))
            else
               write(unit1) (listd(listdptr(ml)+im),im=1,numd(ml))
            endif
#ifdef MPI
         elseif (Node.eq.0) then
            call MPI_IRecv(ibuffer,numdg(m),MPI_integer,BNode,1,
     .           MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            call MPI_ISend(listd(listdptr(ml)+1),numd(ml),MPI_integer,
     .           0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
         endif
         if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
            if (Node.eq.0) then
               if (fmto) then
                  write(unit1, intfmt) (ibuffer(im),im=1,numdg(m))
               else
                  write(unit1) (ibuffer(im),im=1,numdg(m))
               endif
            endif
         endif
#endif
      enddo

#ifdef MPI
      if (Node.eq.0) then
         call de_alloc(ibuffer,name="ibuffer",routine="write_ts_dm")
      endif
#endif

C     Write density matrix
      do is=1,nspin
         do m=1,no_u
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
               ml = m
#endif
               if (fmto) then
                  write(unit1, floatfmt) 
     .                 (dm(listdptr(ml)+im,is),im=1,numd(ml))
               else
                  write(unit1) (dm(listdptr(ml)+im,is),im=1,numd(ml))
               endif
#ifdef MPI
            elseif (Node.eq.0) then
               call MPI_IRecv(buffer,numdg(m),MPI_double_precision,
     .              BNode,1,MPI_Comm_World,Request,MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
               call MPI_ISend(dm(listdptr(ml)+1,is),numd(ml),
     .              MPI_double_precision,0,1,MPI_Comm_World,Request,
     .              MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
               call MPI_Barrier(MPI_Comm_World,MPIerror)
               if (Node.eq.0) then
                  if (fmto) then
                     write(unit1, floatfmt) (buffer(im),im=1,numdg(m))
                  else
                     write(unit1) (buffer(im),im=1,numdg(m))
                  endif
               endif
            endif
#endif
         enddo
      enddo

C     Write energy density matrix
      do is=1,nspin
         do m=1,no_u
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
               ml = m
#endif
               if (fmto) then
                  write(unit1, floatfmt) 
     .                 (edm(listdptr(ml)+im,is),im=1,numd(ml))
               else
                  write(unit1) (edm(listdptr(ml)+im,is),im=1,numd(ml))
               endif
#ifdef MPI
            elseif (Node.eq.0) then
               call MPI_IRecv(buffer,numdg(m),MPI_double_precision,
     .              BNode,1,MPI_Comm_World,Request,MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
               call GlobalToLocalOrb(m,Node,Nodes,ml)
               call MPI_ISend(edm(listdptr(ml)+1,is),numd(ml),
     .              MPI_double_precision,0,1,MPI_Comm_World,Request,
     .              MPIerror)
               call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
               call MPI_Barrier(MPI_Comm_World,MPIerror)
               if (Node.eq.0) then
                  if (fmto) then
                     write(unit1, floatfmt) (buffer(im),im=1,numdg(m))
                  else
                     write(unit1) (buffer(im),im=1,numdg(m))
                  endif
               endif
            endif
#endif
         enddo
      enddo

! Write EF and close
      if (Node.eq.0) then
#ifdef MPI
         call de_alloc(buffer,name="buffer",routine="write_ts_dm")
#endif
         write(unit1) ef
         call io_close(unit1)
      endif

      call de_alloc(numdg,name="numdg",routine="write_ts_dm")

      end subroutine write_ts_dm

!-------------------
    
      subroutine ts_init_dm(found)

      use parallel, only: IOnode 
      use m_ts_global_vars, only : TSinit, TSrun

      logical, intent(in) :: found

         if(.not. found )  then ! not a TS continuation run
           TSinit = .true.  ! start to converge a diagon run
           TSrun = .false.

           if(IONode) then
              write(*,'(a)')
     .          'TRANSIESTA: No TS-DensityMatrix file found'
              write(*,'(a)')
     .           'TRANSIESTA: Initialization runs using diagon'
           endif
         else
           TSinit = .false.
           TSrun = .true.

           if(IONode) then
             write(*,'(a)') 'TRANSIESTA: Continuation run'
             write(*,'(/a)')
     .                 '                     ************************  '
             write(*,'(a)')
     .                 '                     *   TRANSIESTA BEGIN   *  '
             write(*,'(a)')
     .                 '                     ************************  '
           endif
         end if              !found TSDM-file

 
      end subroutine ts_init_dm

!-------------------
      subroutine ts_setup_file_modes()
!
!     Utility routine to deal with different options
!     (AG, after first implementation by Toby White)

!     We might want to use formatted DM files in order to
!     transfer them between computers. 
!     However, whenever the settings are different for input/
!     output, this is only the case for the first step.
!     Thereafter, they must be the same; otherwise we will
!     end up reading from the wrong file.

      if (ionode) then
         if (frstme) then
            fmto = fdf_boolean('DM.FormattedFiles', .false.)
            fmti = fdf_boolean('DM.FormattedInput', fmto)
            fmto = fdf_boolean('DM.FormattedOutput', fmto)
            frstme = .false.
            scndme = .true.
         elseif (scndme) then
            fmti = fmto
            scndme = .false.
         endif
         fnameu = trim(slabel) // '.TSDE'
         fnamef = trim(slabel) // '.TSDEF'
         if (fmto) then
            formout = 'formatted'
            fnameo = fnamef
         else
            formout = 'unformatted'
            fnameo = fnameu
         endif
         if (fmti) then
            formin = 'formatted'
            fnamei = fnamef
         else
            formin = 'unformatted'
            fnamei = fnameu
         endif
      endif

      end subroutine ts_setup_file_modes

      end module m_ts_iodm
