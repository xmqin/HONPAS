! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine iodm( task, maxnd, nbasis, nspin, numd, 
     .                 listdptr, listd, dm, found )
C *******************************************************************
C Reads/writes density matrix from/to file
C Written by P.Ordejon and J.M.Soler. May 1997.
C ********* INPUT ***************************************************
C character task*(*) : 'read' or 'write'
C integer   maxnd    : First dimension of listd and dm
C integer   nbasis   : Number of atomic orbitals
C integer   nspin    : Number of spins (1 or 2)
C ********* INPUT OR OUTPUT (depending on task) *********************
C integer numd(nbasis)     : Control vector of DM matrix
C                            (number of nonzero elements of each row)
C integer listdptr(nbasis) : Control vector of DM matrix
C                            (pointer to the start of each row)
C integer listd(maxnd)     : Control vector of DM matrix
C                            (list of nonzero elements of each row)
C real*8  dm(maxnd,nspin)  : Density matrix
C ********* OUTPUT *************************************************
C logical found : Has DM been found in disk? (Only when task='read')
C ******************************************************************
      
C
C  Modules
C
      use precision,    only : dp
      use parallel,     only : ionode,  Node, Nodes
      use sys,          only : die
      use fdf,          only : fdf_boolean, fdf_string
      use files,        only : slabel, label_length
      use alloc,        only : re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
      use parallelsubs, only : WhichNodeOrb, GlobalToLocalOrb
      use parallelsubs, only : LocalToGlobalOrb
#endif

      implicit  none

      character(len=*), intent(in) :: task
      integer, intent(in) :: maxnd
      integer, intent(in) :: nbasis
      integer, intent(in) :: nspin

      integer, intent(inout) :: numd(nbasis)
      integer, intent(inout) :: listdptr(nbasis)
      integer, intent(inout) :: listd(maxnd)
      real(dp), intent(inout) :: dm(maxnd, nspin)

      logical, intent(out) :: found

C Internal variables
      logical   file_exists
      integer   im, is, unit1, m, nb, ndmax, ns
      integer   nbasistot, ml, ndmaxg
      integer, pointer :: numdg(:) => null()
#ifdef MPI
      integer   MPIerror, Request, Status(MPI_Status_Size)
      integer   BNode
      real(dp), pointer :: buffer(:) => null()
      integer,  pointer :: ibuffer(:) => null()
#endif
      

C Saved internal variables:
      logical,           save :: frstme = .true., scndme = .false.
      character(len=label_length+3), save :: fnameu
      character(len=label_length+4), save :: fnamef
      character(len=label_length+4), save :: fnamei, fnameo
      logical,                       save :: fmti, fmto
      character(len=11),             save :: formin, formout

! Character formats for formatted I/O
! We assume that we have no integers greater than (1e10-1)
! (which is bigger than the largest 32-bit integer)
!
! and that floating point accuracy is sufficiently represented
! by 16 decimal digits of precision (which suffices for 64-bit
! IEEE floats)

      character(len=*), parameter :: intfmt = '(I11)'
      character(len=*), parameter :: floatfmt = '(ES22.14)'

C External procedures
      external          chkdim, memory

! We might want to use formatted DM files in order to
! transfer them between computers. 

! However, whenever the settings are different for input/
! output, this is only the case for the first step.

! Thereafter, they must be the same; otherwise we will
! end up reading from the wrong file.

! Therefore, on all steps after the first, the 

C Find file name
#ifdef IO_TIMING
      call timer( 'iodm', 1 )
#endif
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
      endif

C Find total number of basis functions over all Nodes
#ifdef MPI
      call MPI_AllReduce(nbasis,nbasistot,1,MPI_integer,MPI_sum,
     .     MPI_Comm_World,MPIerror)
#else
      nbasistot = nbasis
#endif

C Allocate local buffer array for globalised numd
      call re_alloc( numdg, 1, nbasistot, 'numdg', 'iodm' )

      if (task.eq.'read' .or. task.eq.'READ') then
        if (Node.eq.0) then
          inquire (file=fnamei,  exist=file_exists)
        endif

#ifdef MPI
C Broadcast logicals so that all processors take the same route
        call MPI_Bcast(file_exists,1,MPI_logical,0,
     $                             MPI_Comm_World,MPIerror)
#endif

        if (file_exists) then

          if (Node.eq.0) then
            write(6,'(/,a)') 'iodm: Reading Density Matrix from file'
            call io_assign(unit1)
            open( unit1, file=fnamei, form=formin, status='old' )
            rewind(unit1)
            if (fmti) then
              read(unit1, intfmt) nb, ns
            else
              read(unit1) nb, ns
            endif
          endif

C Communicate the values to all Nodes and adjust to allow for
C distributed memory before checking the dimensions
#ifdef MPI
          call MPI_Bcast(nb,1,MPI_integer,0,MPI_Comm_World,MPIerror)
          call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif

C Check dimensions
          call chkdim( 'iodm', 'nbasis', nbasistot, nb, 0 )
          call chkdim( 'iodm', 'nspin',  nspin,  ns, 0 )

          if (Node.eq.0) then
            if (fmti) then
              read(unit1, intfmt) (numdg(m),m=1,nbasistot)
            else
              read(unit1) (numdg(m),m=1,nbasistot)
            endif
          endif
#ifdef MPI
          call MPI_Bcast(numdg,nbasistot,MPI_integer,0,MPI_Comm_World,
     .      MPIerror)
#endif

C Convert global numd pointer to local form and generate listdptr
          ndmax = 0
          do m = 1,nbasis
#ifdef MPI
            call LocalToGlobalOrb(m,Node,Nodes,ml)
#else
            ml = m
#endif
            numd(m) = numdg(ml)
            ndmax = ndmax + numd(m)
            if (m .eq. 1) then
              listdptr(1) = 0
            else
              listdptr(m) = listdptr(m-1) + numd(m-1)
            endif
          enddo
          ndmaxg = 0
          do m = 1,nbasistot
            ndmaxg = max(ndmaxg,numdg(m))
          enddo

C Check size of first dimension of dm
          call chkdim( 'iodm', 'maxnd', maxnd, ndmax, 1 )

#ifdef MPI
C Create buffer arrays for transfering density matrix between nodes and lists
          call re_alloc( buffer, 1, ndmaxg, 'buffer', 'iodm' )
          call re_alloc( ibuffer, 1, ndmaxg, 'ibuffer', 'iodm' )
#endif

          do m = 1,nbasistot
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
              ml = m
#endif
              if (fmti) then
                read(unit1, intfmt) 
     .               (listd(listdptr(ml)+im),im=1,numd(ml))
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
     .          BNode,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
              call MPI_IRecv(listd(listdptr(ml)+1),numd(ml),
     .          MPI_integer,0,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            endif
            if (BNode.ne.0) then
              call MPI_Barrier(MPI_Comm_World,MPIerror)
            endif
#endif
          enddo

#ifdef MPI
          call de_alloc( ibuffer, 'ibuffer', 'iodm' )
#endif

          do is = 1,nspin
            do m = 1,nbasistot
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
     .            BNode,1,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)
              elseif (Node.eq.BNode) then
                call GlobalToLocalOrb(m,Node,Nodes,ml)
                call MPI_IRecv(dm(listdptr(ml)+1,is),numd(ml),
     .            MPI_double_precision,0,1,MPI_Comm_World,Request,
     .            MPIerror)
                call MPI_Wait(Request,Status,MPIerror)
              endif
              if (BNode.ne.0) then
                call MPI_Barrier(MPI_Comm_World,MPIerror)
              endif
#endif
            enddo
          enddo

#ifdef MPI
C Free buffer array
          call de_alloc( buffer, 'buffer', 'iodm' )
#endif
          if (Node.eq.0) then
            call io_close(unit1)
          endif

          found = .true.

        else

          found = .false.

        endif

      elseif (task.eq.'write' .or. task.eq.'WRITE') then

        if (Node.eq.0) then
          call io_assign(unit1)

          open( unit1, file=fnameo, form=formout, status='unknown' )
          rewind(unit1)
          if (fmto) then
            write(unit1, intfmt) nbasistot, nspin
          else
            write(unit1) nbasistot, nspin
          endif
        endif

C Create globalised numd
        do m = 1,nbasistot
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
     .        0,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.0) then
            call MPI_IRecv(numdg(m),1,MPI_integer,
     .        BNode,1,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          endif
          if (BNode.ne.0) then
            call MPI_Barrier(MPI_Comm_World,MPIerror)
          endif
#endif
        enddo

C Write out numd array
        if (Node.eq.0) then
          ndmaxg = 0
          do m = 1,nbasistot
            ndmaxg = max(ndmaxg,numdg(m))
          enddo
          if (fmto) then
            write(unit1, intfmt) (numdg(m),m=1,nbasistot)
          else
            write(unit1) (numdg(m),m=1,nbasistot)
          endif
#ifdef MPI
          call re_alloc( buffer, 1, ndmaxg, 'buffer', 'iodm' )
          call re_alloc( ibuffer, 1, ndmaxg, 'ibuffer', 'iodm' )
#endif
        endif

C Write out listd array
        do m = 1,nbasistot
#ifdef MPI
          call WhichNodeOrb(m,Nodes,BNode)
          if (BNode.eq.0.and.Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
            ml = m
#endif
            if (fmto) then
              write(unit1, intfmt)
     .             (listd(listdptr(ml)+im),im=1,numd(ml))
            else
              write(unit1) (listd(listdptr(ml)+im),im=1,numd(ml))
            endif
#ifdef MPI
          elseif (Node.eq.0) then
            call MPI_IRecv(ibuffer,numdg(m),MPI_integer,BNode,1,
     .        MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.BNode) then
            call GlobalToLocalOrb(m,Node,Nodes,ml)
            call MPI_ISend(listd(listdptr(ml)+1),numd(ml),MPI_integer,
     .        0,1,MPI_Comm_World,Request,MPIerror)
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
          call de_alloc( ibuffer, 'ibuffer', 'iodm' )
        endif
#endif

C Write density matrix
        do is=1,nspin
          do m=1,nbasistot
#ifdef MPI
            call WhichNodeOrb(m,Nodes,BNode)
            if (BNode.eq.0.and.Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
#else
              ml = m
#endif
              if (fmto) then
                write(unit1, floatfmt) 
     .               (dm(listdptr(ml)+im,is),im=1,numd(ml))
              else
                write(unit1) (dm(listdptr(ml)+im,is),im=1,numd(ml))
              endif
#ifdef MPI
            elseif (Node.eq.0) then
              call MPI_IRecv(buffer,numdg(m),MPI_double_precision,
     .          BNode,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.BNode) then
              call GlobalToLocalOrb(m,Node,Nodes,ml)
              call MPI_ISend(dm(listdptr(ml)+1,is),numd(ml),
     .          MPI_double_precision,0,1,MPI_Comm_World,Request,
     .          MPIerror)
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
          call de_alloc( buffer, 'buffer', 'iodm' )
#endif
          call io_close(unit1)
        endif

      else
        if (Node.eq.0) then
          call die('iodm: incorrect task')
        endif
      endif

C Deallocate local buffer array for globalised numd
      call de_alloc( numdg, 'numdg', 'iodm' )

#ifdef IO_TIMING
      call timer( 'iodm', 2 )
#endif
      end subroutine iodm
