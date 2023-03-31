!     
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!     
      module m_iodm_old

      use precision,     only: dp
      use parallel,     only : ionode,  Node, Nodes
      use parallelsubs, only : LocalToGlobalOrb
      use sys,          only : die
      use fdf,          only : fdf_boolean, fdf_string, fdf_get
      use files,        only : slabel, label_length
      use alloc,        only : re_alloc, de_alloc
      use m_mpi_utils,  only : broadcast
#ifdef MPI
      use mpi_siesta
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
#endif

      implicit none

      private

      ! aliases for old routines

      interface read_spmatrix
         module procedure read_dm
      end interface
      interface write_spmatrix
         module procedure write_spmatrix1, write_spmatrix2
      end interface
      interface read_dynamic_spmatrix
         module procedure read_dynamic_dm
      end interface

      PUBLIC :: write_dm, read_dm, read_dynamic_dm
      PUBLIC :: write_spmatrix, read_spmatrix, read_dynamic_spmatrix

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

      subroutine read_dynamic_dm ( maxnd, no_l, nspin, numd, 
     .                     listdptr, listd, dm, found, userfile )

      ! The only "in" parameter is no_l, to determine the distribution
      integer, intent(in) :: no_l

      ! These are all "output" sparsity parameters
      ! We read the file info into a temporary
      ! set of arrays and change the structure if needed
      ! All pointers should be nullified on entry

      integer, intent(out) :: maxnd
      integer, intent(out) :: nspin

      integer, pointer     :: numd(:)
      integer, pointer     :: listdptr(:)
      integer, pointer     :: listd(:)
      real(dp), pointer    :: dm(:,:)

      logical, intent(out) :: found
      character(len=*), intent(in), optional :: userfile


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

      call setup_file_modes(userfile)

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
         write(6,'(/,a)')
     $        'Reading matrix from file ' // trim(fnamei)
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
     $      "WARNING: Wrong number of orbitals in DM file: ",
     $       nb,
     $      "WARNING: Falling back to atomic initialization of DM."
           write(0,"(a,i6,/,a)")
     $      "WARNING: Wrong number of orbitals in DM file: ",
     $       nb,
     $      "WARNING: Falling back to atomic initialization of DM."
            call io_close(unit1)
         endif
         found = .false.
         RETURN
      endif

C     Allocate local buffer array for globalised numd
      nullify(numdg)
      call re_alloc( numdg, 1, no_u, 'numdg', 'read_dm' )
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

      ! Allocate density matrix and listd
      call re_alloc( listd, 1, maxnd, 'listd', 'read_dm' )
      call re_alloc( dm, 1, maxnd, 1, nspin, 'dm', 'read_dm' )

#ifdef MPI
C     Create buffer arrays for transfering density matrix between nodes and lists
      nullify(buffer,ibuffer)
      call re_alloc( buffer,  1, ndmaxg, 'buffer'  , 'read_dm' )
      call re_alloc( ibuffer, 1, ndmaxg, 'ibuffer' , 'read_dm' )
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
      call de_alloc( ibuffer, 'ibuffer', 'read_dm' )
#endif

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

#ifdef MPI
      call de_alloc( buffer, 'buffer', 'read_dm' )
#endif
      call de_alloc( numdg, 'numdg', 'read_dm' )

      if (Node.eq.0) then
         call io_close(unit1)
      endif

      found = .true.

      end subroutine read_dynamic_dm
!-----------------------      
      subroutine read_dm ( maxnd, no_l, nspin, numd, 
     .                     listdptr, listd, dm, found, userfile )

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

      logical, intent(out) :: found
      character(len=*), intent(in), optional :: userfile

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

      call timer("ReadM",1)

      call setup_file_modes(userfile)

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
         write(6,'(/,a)') 'Reading Matrix from file ' // trim(fnamei)
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
      call re_alloc( numdg, 1, no_u, 'numdg', 'read_dm' )
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
      call re_alloc( buffer,  1, ndmaxg, 'buffer',  'read_dm' )
      call re_alloc( ibuffer, 1, ndmaxg, 'ibuffer', 'read_dm' )
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
      call de_alloc( ibuffer, 'ibuffer', 'read_dm' )
#endif

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

#ifdef MPI
      call de_alloc( buffer, 'buffer', 'read_dm' )
#endif
      call de_alloc( numdg, 'numdg', 'read_dm' )

      if (Node.eq.0) then
         call io_close(unit1)
      endif

      found = .true.
      call timer("ReadM",2)

      end subroutine read_dm
!-----------------------      
      subroutine write_dm (maxnd, no_l, nspin,
     $     numd, listdptr, listd, dm, userfile)

      integer, intent(in) :: maxnd
      integer, intent(in) :: no_l
      integer, intent(in) :: nspin
      integer, intent(in) :: numd(1:no_l)
      integer, intent(in) :: listdptr(1:no_l)
      integer, intent(in) :: listd(maxnd)
      real(dp), intent(in) :: dm(maxnd,nspin)

      character(len=*), intent(in), optional :: userfile

      integer :: no_u, m, ml, im, ndmaxg, unit1, is
#ifdef MPI
      integer   MPIerror, Request, Status(MPI_Status_Size)
      integer   BNode
      real(dp), dimension(:), pointer :: buffer
      integer,  dimension(:), pointer :: ibuffer
#endif
      integer, dimension(:), pointer  :: numdg

      call timer("WriteM",1)

      call setup_file_modes(userfile)

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
      call re_alloc( numdg, 1, no_u, 'numdg', 'write_dm' )

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
         call re_alloc( buffer,  1, ndmaxg, 'buffer',  'write_dm' )
         call re_alloc( ibuffer, 1, ndmaxg, 'ibuffer', 'write_dm' )
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
         call de_alloc( ibuffer, 'ibuffer', 'write_dm' )
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

      if (Node.eq.0) then
#ifdef MPI
         call de_alloc( buffer, 'buffer', 'write_dm' )
#endif
         call io_close(unit1)
      endif

      call de_alloc( numdg, 'numdg', 'write_dm' )

      call timer("WriteM",2)

      end subroutine write_dm

!-----------------------------------------------------
      subroutine write_spmatrix2(M,file,when)
      !
      ! Convenience wrapper
      !
      use sparse_matrices, only: numh, listhptr, listh
      use m_matio,         only: write_mat

      real(dp), intent(in) :: M(:,:)
      character(len=*), intent(in), optional :: file
      logical, intent(in), optional          :: when

      integer :: no_l, maxnh, nspin

      if (present(when)) then
         if (.not. when) return
      endif

      no_l = size(numh)
      maxnh = size(M,dim=1)
      nspin = size(M,dim=2)

      if (fdf_get("Use.Blocked.WriteMat",.false.)) then
         call write_mat (maxnh, no_l, nspin,
     $     numh, listhptr, listh, M, userfile=trim(file)//'.blocked',
     $                               compatible=.false.)
      else
         call write_dm (maxnh, no_l, nspin,
     $     numh, listhptr, listh, M, userfile=file)
      endif
      end subroutine write_spmatrix2

!-----------------------------------------------------
      subroutine write_spmatrix1(M,file,when)
      !
      ! Convenience wrapper for 1-dim sparse data
      ! (i.e., no spin, as in S matrix)
      !
      use sparse_matrices, only: numh, listhptr, listh
      use m_matio,         only: write_mat

      real(dp), intent(in) :: M(:)
      character(len=*), intent(in), optional :: file
      logical, intent(in), optional          :: when

      integer :: no_l, maxnh, nspin

      if (present(when)) then
         if (.not. when) return
      endif

      no_l = size(numh)
      maxnh = size(M,dim=1)
      nspin = 1
      
      if (fdf_get("Use.Blocked.WriteMat",.false.)) then
         call write_mat (maxnh, no_l, nspin,
     $     numh, listhptr, listh, M, userfile=trim(file)//'.blocked',
     $                               compatible=.false.)
      else
         call write_dm (maxnh, no_l, nspin,
     $     numh, listhptr, listh, M, userfile=file)
      endif
      end subroutine write_spmatrix1

!-------------------
      subroutine setup_file_modes(userfile)

      character(len=*), intent(in), optional :: userfile
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
         fnameu = trim(slabel) // '.DM'
         fnamef = trim(slabel) // '.DMF'
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

         if (present(userfile)) then
            fnamei = userfile
            fnameo = userfile
            formout = 'unformatted'
            formin  = 'unformatted'
         endif
      endif

      end subroutine setup_file_modes

      end module m_iodm_old
