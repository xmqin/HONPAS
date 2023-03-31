! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine iolwf( task, nbasis, nbasisCloc, nbasisloc, maxnc, 
     .                  found, nspin)
C *******************************************************************
C Reads/writes localized wave functions from/to file
C Written by P.Ordejon and J.M.Soler. May 1997.
C ********* INPUT ***************************************************
C character task*3     : 'read' or 'write'
C integer   nbasis     : Number of atomic orbitals globally
C integer   nbasisCloc : Number of rows stored locally in C
C integer   nbasisloc  : Number of orbitals for which the current node
C                      : is responsible
C integer   nspin      : Number of spins (1 or 2)
C ********* INPUT OR OUTPUT (depending on task) *********************
C integer numc(nbasis)       : Control vector of c matrix
C                              (number of nonzero elements of each row)
C integer listc(maxc,maxoloc): Control vector of c matrix
C                              (list of nonzero elements of each row)
C integer maxnc              : Dimension of C matrix
C real*8  c(maxc,maxoloc,nspin) : Density matrix
C ********* OUTPUT *************************************************
C logical found : Have LWF's been found in disk? (Only when task='read')
C ******************************************************************

C
C  Modules
C
      use precision,     only : dp
      use parallel,      only : Node, Nodes
      use parallelsubs
      use files,         only : slabel, label_length
      use on_main,       only : c, cold, numc, numcold, listc, listcold
      use on_main,       only : xi, g, hg
      use sys,           only : die
#ifdef MPI
      use mpi_siesta,    only : MPI_logical, MPI_Status_Size
      use mpi_siesta,    only : MPI_Comm_World, MPI_integer
      use mpi_siesta,    only : MPI_double_precision
#endif
      use alloc

      implicit  none

      character task*(*)
      logical   found
      integer   nbasis, nbasisCloc, nbasisloc, nspin, maxnc
      external  chkdim, io_assign, io_close, timer

C Internal variables and arrays
      character(len=label_length+4) :: fname
      logical                       :: exist3
      integer, pointer :: numcg(:) => null()
      integer :: im, is, unit1, m, nb, ncmax, ns
#ifdef MPI
      integer :: MPIerror, n, ndata, Request, Status(MPI_Status_Size)
      real(dp), pointer :: bdens(:) => null()
      integer, pointer :: idens(:) => null()
#endif

      call timer( 'iolwf', 1 )

#ifdef MPI
C Allocate local buffer array
      call re_alloc( idens, 1, nbasis, 'idens', 'iolwf' )
      call re_alloc( bdens, 1, nbasis, 'bdens', 'iolwf' )
#endif
      
C     Find file name
      fname = trim(slabel) // '.LWF'

      if (task.eq.'read' .or. task.eq.'READ') then
        if (Node.eq.0) then
          inquire(file=fname,exist=exist3)
        endif

#ifdef MPI
C Broadcast logicals
        call MPI_Bcast(exist3,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif

C Look now for new-format files
        if (exist3) then

C Allocate local workspace
          call re_alloc( numcg, 1, nbasis, 'numcg', 'iolwf' )

          if (Node.eq.0) then
            call io_assign(unit1)
            open( unit1, file=fname,
     .          form='unformatted', status='unknown' )
            rewind(unit1)
            read(unit1) nb, ns
          endif

#ifdef MPI        
C Broadcast values for basis set and spin so that cross-checks can be performed
          call MPI_Bcast(nb,1,MPI_integer,0,MPI_Comm_World,MPIerror)
          call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
          call chkdim( 'iolwf', 'nbasis', nbasis, nb, 0 )
          call chkdim( 'iolwf', 'nspin',  nspin,  ns, 0 )

C Read in global numc array
          if (Node.eq.0) then
            read(unit1) (numcg(m),m=1,nbasis)
          endif

C Collect elements of numc array on IO node from self
          if (Node.eq.0) then
            do m = 1,nOrbPerNode(1)
              numcold(m) = numcg(nL2G(m,1))
            enddo
          endif

#ifdef MPI
C Collect elements of numc array on IO node from other nodes
          do n = 1,Nodes-1
            if (Node.eq.0) then
              ndata = 0
              do m = 1,nbasis
                if (nNode(m).eq.n) then
                  ndata = ndata + 1
                  idens(ndata) = numcg(m)
                endif
              enddo
              call MPI_ISend(idens,ndata,MPI_integer,
     .          n,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.n) then
              ndata = 0
              do m = 1,nbasis
                if (nNode(m).eq.n) then
                  ndata = ndata + 1
                endif
              enddo
              call MPI_IRecv(idens,ndata,MPI_integer,
     .          0,1,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
              ndata = 0
              do m = 1,nbasis
                if (nNode(m).eq.n) then
                  ndata = ndata + 1
                  numcold(nG2L(m)) = idens(ndata)
                endif
              enddo
            endif
          enddo
#endif

C Check size of arrays
          ncmax = 0
          do m = 1,nOrbPerNode(Node+1)
            ncmax = max(ncmax,numc(m))
          enddo
          if (ncmax.gt.maxnc) then
            maxnc = ncmax 
            call re_alloc(listc,1,maxnc,1,nbasisCloc,name='listc')
            call re_alloc(listcold,1,maxnc,1,nbasisloc,name='listcold')
            call re_alloc(c,1,maxnc,1,nbasisCloc,1,nspin,name='c')
            call re_alloc(cold,1,maxnc,1,nbasisloc,1,nspin,name='cold')
            call re_alloc(xi,1,maxnc,1,nbasisCloc,1,nspin,name='xi')
            call re_alloc(g,1,maxnc,1,nbasisCloc,1,nspin,name='g')
            call re_alloc(hg,1,maxnc,1,nbasisCloc,1,nspin,name='hg')
          endif

C Initialise C before reading
          cold(1:maxnc,1:nbasisloc,1:nspin) = 0.0d0

C Loop over basis functions reading in rows of listc and distributing
          do m = 1,nbasis 
#ifdef MPI
            if (nNode(m).eq.0.and.Node.eq.0) then
#endif
C Data is on IO node so just read it in
              read(unit1)(listcold(im,nG2L(m)),im=1,numcg(m))
#ifdef MPI
            elseif (Node.eq.0) then
C Send data to IO node
              read(unit1)(idens(im),im=1,numcg(m))
              call MPI_ISend(idens,numcg(m),MPI_integer,
     .          nNode(m),m,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.nNode(m)) then
C Receive data on IO node
              ndata = numcold(nG2L(m))
              call MPI_IRecv(idens,ndata,MPI_integer,
     .          0,m,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
              do im = 1,ndata
                listcold(im,nG2L(m)) = idens(im)
              enddo
            endif
#endif
          enddo

C Loop over basis functions reading in rows of c and distributing
          do is = 1,nspin
            do m = 1,nbasis 
#ifdef MPI 
              if (nNode(m).eq.0.and.Node.eq.0) then
#endif    
C Data is on IO node so just read it in
                read(unit1)(cold(im,nG2L(m),is),im=1,numcg(m))
#ifdef MPI  
              elseif (Node.eq.0) then
C Send data to IO node
                read(unit1)(bdens(im),im=1,numcg(m))
                call MPI_ISend(bdens,numcg(m),MPI_double_precision,
     .            nNode(m),m,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)
              elseif (Node.eq.nNode(m)) then
C Receive data on IO node
                ndata = numcold(nG2L(m))
                call MPI_IRecv(bdens,ndata,MPI_double_precision,
     .            0,m,MPI_Comm_World,Request,MPIerror)
                call MPI_Wait(Request,Status,MPIerror)
                do im = 1,ndata
                  cold(im,nG2L(m),is) = bdens(im)
                enddo
              endif
#endif    
            enddo
          enddo

          if (Node.eq.0) then
            call io_close(unit1)
          endif

C Free local workspace
          call de_alloc( numcg, 'numcg', 'iolwf' )
          found = .true.
        else
          found = .false.
        endif

      elseif (task.eq.'write' .or. task.eq.'WRITE') then

C Allocate local workspace
        call re_alloc( numcg, 1, nbasis, 'numcg', 'iolwf' )

        if (Node.eq.0) then
          call io_assign(unit1)
          open( unit1, file=fname,
     .        form='unformatted', status='unknown' )
          rewind(unit1)
          write(unit1) nbasis, nspin
        endif

C Collect elements of numc array on IO node from self
        if (Node.eq.0) then
          do m = 1,nOrbPerNode(1)
            numcg(nL2G(m,1)) = numc(m)
          enddo
        endif

#ifdef MPI
C Collect elements of numc array on IO node from other nodes
        do n = 1,Nodes-1
          if (Node.eq.n) then
            ndata = 0
            do m = 1,nbasis
              if (nNode(m).eq.n) then
                ndata = ndata + 1
                idens(ndata) = numc(nG2L(m))
              endif
            enddo
            call MPI_ISend(idens,ndata,MPI_integer,
     .        0,n,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.0) then
            ndata = 0
            do m = 1,nbasis
              if (nNode(m).eq.n) then
                ndata = ndata + 1
              endif
            enddo
            call MPI_IRecv(idens,ndata,MPI_integer,
     .        n,n,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
            ndata = 0
            do m = 1,nbasis
              if (nNode(m).eq.n) then
                ndata = ndata + 1
                numcg(m) = idens(ndata)
              endif
            enddo
          endif
        enddo
#endif

C Write out numc array
        if (Node.eq.0) then
          write(unit1) (numcg(m),m=1,nbasis)
        endif

C Loop over basis functions collecting rows of listc and writing out
        do m = 1,nbasis 
#ifdef MPI
          if (nNode(m).eq.0.and.Node.eq.0) then
#endif
C Data is on IO node so just write it out
            write(unit1)(listc(im,nG2L(m)),im=1,numcg(m))
#ifdef MPI
          elseif (Node.eq.nNode(m)) then
C Send data to IO node
            ndata = numc(nG2L(m))
            do im = 1,ndata
              idens(im) = listc(im,nG2L(m))
            enddo
            call MPI_ISend(idens,ndata,MPI_integer,
     .        0,m,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
          elseif (Node.eq.0) then
C Receive data on IO node
            call MPI_IRecv(idens,numcg(m),MPI_integer,
     .        nNode(m),m,MPI_Comm_World,Request,MPIerror)
            call MPI_Wait(Request,Status,MPIerror)
            write(unit1)(idens(im),im=1,numcg(m))
          endif
#endif
        enddo

C Loop over basis functions collecting rows of c and writing out
        do is = 1,nspin
          do m = 1,nbasis 
#ifdef MPI 
            if (nNode(m).eq.0.and.Node.eq.0) then
#endif    
C Data is on IO node so just write it out
              write(unit1)(c(im,nG2L(m),is),im=1,numcg(m))
#ifdef MPI  
            elseif (Node.eq.nNode(m)) then
C Send data to IO node
              ndata = numc(nG2L(m))
              do im = 1,ndata
                bdens(im) = c(im,nG2L(m),is)
              enddo
              call MPI_ISend(bdens,ndata,MPI_double_precision,
     .          0,m,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
            elseif (Node.eq.0) then
C Receive data on IO node
              call MPI_IRecv(bdens,numcg(m),MPI_double_precision,
     .          nNode(m),m,MPI_Comm_World,Request,MPIerror)
              call MPI_Wait(Request,Status,MPIerror)
              write(unit1)(bdens(im),im=1,numcg(m))
            endif
#endif    
          enddo
        enddo

        if (Node.eq.0) then
          call io_close(unit1)
        endif

C Free local workspace
        call de_alloc( numcg, 'numcg', 'iolwf' )

      else
         call die('iolwf: incorrect task')
      endif

#ifdef MPI
C Free local buffer memory
      call de_alloc( bdens, 'bdens', 'iolwf' )
      call de_alloc( idens, 'idens', 'iolwf' )
#endif

      call timer( 'iolwf', 2 )
      end
