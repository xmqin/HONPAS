! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_iorho
!
!     To have both reading and writing, allowing a change in the 
!     values of variables used as array dimensions, is not very
!     safe. Thus, the old iorho has been split in three routines:

!      write_rho: Writes grid magnitudes to file
!      read_rho: Reads grid magnitudes from file
!      check_rho: Checks the appropriate dimensions for reading.
!
!     The magnitudes are saved in SINGLE PRECISION (sp), regardless
!     of the internal precision used in the program. Historically,
!     the internal precision has been "single". 
!

      use precision,    only : dp, grid_p, sp
      use parallel,     only : Node, ProcessorY
      use sys,          only : die
      use alloc,        only : re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
      use parallel,     only : Nodes
      use parallelsubs, only : HowManyMeshPerNode
#endif

      implicit          none

      character(len=*), parameter  :: fform = "unformatted"

      public :: write_rho, read_rho, check_rho
         
      private

      CONTAINS

      subroutine write_rho( fname, cell, mesh, nsm, maxp, nspin, rho)

! Writes the electron density or potential at the mesh points.
!
C Writen by J.Soler July 1997.
C Parallel modifications added, while maintaining independence
C of density matrix on disk from parallel distribution. Uses a
C block distribution for density matrix. It is important to
C remember that the density matrix is divided so that all
C sub-points reside on the same node. Modified by J.D.Gale March 1999.
C NOTE : in order to achieve a consistent format of the disk file
C each record in the unformatted file corresponds to one pair of
C Y and Z values. Hence there will be a total of mesh(2) x mesh(3)
C records.

      character(len=*), intent(in) ::     fname       ! File name
      integer, intent(in)          ::     nsm         
          !Number of sub-mesh points per mesh point
          !  (not used in this version) 
      integer, intent(in)          ::     mesh(3)
          ! Number of mesh divisions of each lattice vector
      integer, intent(in)          ::     maxp
      integer, intent(in)          ::     nspin
      real(grid_p), intent(in)     ::     rho(maxp,nspin)
      real(dp), intent(in)         ::     cell(3,3)

      external          io_assign, io_close, memory

C Internal variables and arrays
      integer    i, ip, iu, is, j, np,
     .           BlockSizeY, BlockSizeZ, ProcessorZ,
     .           meshnsm(3), NRemY, NRemZ,
     .           iy, iz, izm, Ind, ir

#ifdef MPI
      integer    Ind2, MPIerror, Request,
     .           Status(MPI_Status_Size), BNode, NBlock
      real(grid_p), pointer :: bdens(:) => null()
#endif
      real(sp),     pointer :: temp(:) => null()

#ifdef MPI
#ifdef DEBUG
        call write_debug( 'ERROR: write_rho not ready yet' )
        call write_debug( fname )
#endif

C Work out density block dimensions
      if (mod(Nodes,ProcessorY).gt.0)
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')
      ProcessorZ = Nodes/ProcessorY
      BlockSizeY = ((((mesh(2)/nsm)-1)/ProcessorY) + 1)*nsm
      call re_alloc(bdens, 1, BlockSizeY*mesh(1), 'bdens', 'write_rho')
#else
      ProcessorZ = 1
#endif

      call re_alloc( temp, 1, mesh(1), 'temp', 'write_rho' )
C Open file
      if (Node.eq.0) then
        call io_assign( iu )
        open( iu, file=fname, form=fform, status='unknown' )      
      endif

      np = mesh(1) * mesh(2) * mesh(3)

      meshnsm(1) = mesh(1)/nsm
      meshnsm(2) = mesh(2)/nsm
      meshnsm(3) = mesh(3)/nsm

C Write data
      if (Node.eq.0) then
        if (fform .eq. 'formatted') then
          do i = 1,3
            write(iu,*) (cell(j,i),j=1,3)
          enddo
          write(iu,*) mesh, nspin
        else
          write(iu) cell
          write(iu) mesh, nspin
        endif
      endif

C Outer loop over spins
        do is = 1,nspin

          Ind = 0

C Loop over Z dimension of processor grid
          do iz = 1,ProcessorZ

            BlockSizeZ = (meshnsm(3)/ProcessorZ)
            NRemZ = meshnsm(3) - BlockSizeZ*ProcessorZ
            if (iz-1.lt.NRemZ) BlockSizeZ = BlockSizeZ + 1
            BlockSizeZ = BlockSizeZ*nsm

C Loop over local Z mesh points
            do izm = 1,BlockSizeZ

C Loop over blocks in Y mesh direction
              do iy = 1,ProcessorY

C Work out size of density sub-matrix to be transfered
                BlockSizeY = (meshnsm(2)/ProcessorY)
                NRemY = meshnsm(2) - BlockSizeY*ProcessorY
                if (iy-1.lt.NRemY) BlockSizeY = BlockSizeY + 1
                BlockSizeY = BlockSizeY*nsm

#ifdef MPI
                NBlock = BlockSizeY*mesh(1)
C Work out which node block is stored on
                BNode = (iy-1)*ProcessorZ + iz - 1

                if (BNode.eq.0.and.Node.eq.BNode) then
#endif
C If density sub-matrix is local Node 0 then just write it out
                  if (fform .eq. 'formatted') then
                    do ir = 1,BlockSizeY
                      write(iu,'(e15.6)') (rho(Ind+ip,is),
     .                  ip=1,mesh(1))
                      Ind = Ind + mesh(1)
                    enddo
                  else
!AG**
                    do ir = 1,BlockSizeY
                      temp(1:mesh(1)) =
     $                      real(rho(Ind+1:Ind+mesh(1),is),kind=sp)
                      write(iu) temp(1:mesh(1))
                      Ind = Ind + mesh(1)
                    enddo
                  endif

#ifdef MPI
                elseif (Node.eq.0) then
C If this is Node 0 then recv and write density sub-matrix
                  call MPI_IRecv(bdens,NBlock,MPI_grid_real,BNode,1,
     .              MPI_Comm_World,Request,MPIerror)
                  call MPI_Wait(Request,Status,MPIerror)

                elseif (Node.eq.BNode) then
C If this is the Node where the density sub-matrix is, then send
                  call MPI_ISend(rho(Ind+1,is),NBlock,MPI_grid_real,0,1,
     .              MPI_Comm_World,Request,MPIerror)
                  call MPI_Wait(Request,Status,MPIerror)
                  Ind = Ind + NBlock

                endif

                if (BNode.ne.0) then
                  call MPI_Barrier(MPI_Comm_World,MPIerror)
                  if (Node.eq.0) then
                    Ind2 = 0
                    if (fform .eq. 'formatted') then
                      do ir = 1,BlockSizeY
                        write(iu,'(e15.6)') (bdens(Ind2+ip),ip=1,
     .                    mesh(1))
                        Ind2 = Ind2 + mesh(1)
                      enddo
                    else
!AG**
                      do ir = 1,BlockSizeY
                      temp(1:mesh(1)) =
     $                      real(bdens(Ind2+1:Ind2+mesh(1)),kind=sp)
                        write(iu) temp(1:mesh(1))
                        Ind2 = Ind2 + mesh(1)
                      enddo

                    endif
                  endif
                endif
#endif

              enddo

            enddo

          enddo

        enddo

        if (Node.eq.0) then
C Close file
          call io_close( iu )
        endif

#ifdef MPI
C Deallocate density buffer memory
      call de_alloc( bdens, 'bdens', 'write_rho' )
#endif
      call de_alloc( temp, 'temp', 'write_rho' )
      end subroutine write_rho

!------------------------------------------------------------------------

      subroutine read_rho( fname, cell, mesh, nsm, maxp, nspin, rho)

C *********************************************************************
C Reads the electron density or potential at the mesh points.
!
! AG: The magnitudes are saved in SINGLE PRECISION (sp), regardless
!     of the internal precision used in the program. Historically,
!     the internal precision has been "single". 
!
C Writen by J.Soler July 1997.
C Parallel modifications added, while maintaining independence
C of density matrix on disk from parallel distribution. Uses a
C block distribution for density matrix. It is important to
C remember that the density matrix is divided so that all
C sub-points reside on the same node. Modified by J.D.Gale March 1999.
C NOTE : in order to achieve a consistent format of the disk file
C each record in the unformatted file corresponds to one pair of
C Y and Z values. Hence there will be a total of mesh(2) x mesh(3)
C records.
C *************************** INPUT **********************************
C character*(*) fname     : File name for input or output
C integer nsm             : Number of sub-mesh points per mesh point
C                           (not used in this version)
C integer maxp            : First dimension of array rho
C integer nspin           : Second dimension of array rho
C ************************** OUTPUT **********************************
C real*8  cell(3,3)       : Lattice vectors
C integer mesh(3)         : Number of mesh divisions of each
C                           lattice vector
C real    rho(maxp,nspin) : Electron density
C ******************** BEHAVIOUR **************************************
C If the values of maxp or nspin on input are less than
C those required to copy the array rho from the file, the subroutine
C stops. Use subroutine check_rho to find the right maxp and nspin.
C *********************************************************************

C Arguments
      character(len=*), intent(in) ::     fname
      integer, intent(in)          ::     nsm, maxp, nspin
      integer, intent(out)         ::     mesh(3)
      real(grid_p), intent(out)    ::     rho(maxp,nspin)
      real(dp), intent(out)        ::     cell(3,3)

      external          io_assign, io_close, memory

C Internal variables and arrays
      integer    ip, iu, is, np, ns, 
     .           BlockSizeY, BlockSizeZ, ProcessorZ,
     .           meshnsm(3), npl, NRemY, NRemZ,
     .           iy, iz, izm, Ind, ir

#ifdef MPI
      integer    Ind2, MPIerror, Request, meshl(3),
     .           Status(MPI_Status_Size), BNode, NBlock
      logical    ltmp
      real(grid_p), pointer :: bdens(:) => null()
#endif
      real(sp),     pointer :: temp(:) => null()

      logical    baddim, found


#ifdef MPI
C Work out density block dimensions
      if (mod(Nodes,ProcessorY).gt.0)
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')
      ProcessorZ = Nodes/ProcessorY
      BlockSizeY = ((((mesh(2)/nsm)-1)/ProcessorY) + 1)*nsm
      call re_alloc( bdens, 1, BlockSizeY*mesh(1), 'bdens', 'read_rho' )
#else
      ProcessorZ = 1
#endif

C Check if input file exists
        if (Node.eq.0) then
          inquire( file=fname, exist=found )
          if (.not. found) call die("Cannot find file " // trim(fname))
        endif

C Open file
          if (Node.eq.0) then
            call io_assign( iu )
            open( iu, file=fname, form=fform, status='old' )      

C Read cell vectors and number of points
            if (fform .eq. 'formatted') then
              read(iu,*) cell
              read(iu,*) mesh, ns
            else
              read(iu) cell
              read(iu) mesh, ns
            endif
          endif
          call re_alloc( temp, 1, mesh(1), 'temp', 'read_rho' )
#ifdef MPI
          call MPI_Bcast(cell(1,1),9,MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)
          call MPI_Bcast(mesh,3,MPI_integer,0,MPI_Comm_World,MPIerror)
          call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
          np = mesh(1) * mesh(2) * mesh(3)

C Get local dimensions
          meshnsm(1) = mesh(1)/nsm
          meshnsm(2) = mesh(2)/nsm
          meshnsm(3) = mesh(3)/nsm
#ifdef MPI
          call HowManyMeshPerNode(meshnsm,Node,Nodes,npl,meshl)
#else
          npl = np
#endif

C Check dimensions
          baddim = .false.
          if (ns .gt. nspin) baddim = .true.
          if (npl .gt. maxp) baddim = .true.

#ifdef MPI
C Globalise dimension check
          call MPI_AllReduce(baddim,ltmp,1,MPI_logical,MPI_Lor,
     .      MPI_Comm_World,MPIerror)
          baddim = ltmp
#endif

          if (baddim) then
             call die("Dimensions of array rho too small in read_rho")
          endif

C Outer loop over spins
          do is = 1,ns

          Ind = 0

C Loop over Z mesh direction
          do iz = 1,ProcessorZ

C Work out number of mesh points in Z direction
            BlockSizeZ = (meshnsm(3)/ProcessorZ)
            NRemZ = meshnsm(3) - BlockSizeZ*ProcessorZ
            if (iz-1.lt.NRemZ) BlockSizeZ = BlockSizeZ + 1
            BlockSizeZ = BlockSizeZ*nsm

C Loop over local Z mesh points
            do izm = 1,BlockSizeZ

C Loop over blocks in Y mesh direction
              do iy = 1,ProcessorY

C Work out size of density sub-matrix to be transfered
                BlockSizeY = (meshnsm(2)/ProcessorY)
                NRemY = meshnsm(2) - BlockSizeY*ProcessorY
                if (iy-1.lt.NRemY) BlockSizeY = BlockSizeY + 1
                BlockSizeY = BlockSizeY*nsm

#ifdef MPI
                NBlock = BlockSizeY*mesh(1)
C Work out which node block is stored on
                BNode = (iy-1)*ProcessorZ + iz - 1

                if (BNode.eq.0.and.Node.eq.BNode) then
#endif
C If density sub-matrix is local Node 0 then just read it in
                  if (fform .eq. 'formatted') then
                    do ir = 1,BlockSizeY
                      read(iu,*) (rho(Ind+ip,is),ip=1,mesh(1))
                      Ind = Ind + mesh(1)
                    enddo
                  else
!
! AG**: Check sizes of records here -- grid-precision issues
!       Either agree on a "single" format for file storage,
!       or put in some intelligence in all post-processors

                    do ir = 1,BlockSizeY
                      read(iu) temp(1:mesh(1))
                      rho(Ind+1:Ind+mesh(1),is) =
     $                              real(temp(1:mesh(1)), kind=grid_p)
                      Ind = Ind + mesh(1)
                    enddo
                  endif

#ifdef MPI
                elseif (Node.eq.0) then
C If this is Node 0 then read and send density sub-matrix
                  Ind2 = 0
                  if (fform .eq. 'formatted') then
                    do ir = 1,BlockSizeY
                      read(iu,*) (bdens(Ind2+ip),ip=1,mesh(1))
                      Ind2 = Ind2 + mesh(1)
                    enddo
                  else
!AG**
                    do ir = 1,BlockSizeY
                      read(iu) temp(1:mesh(1))
                      bdens(Ind2+1:Ind2+mesh(1)) =
     $                              real(temp(1:mesh(1)), kind=grid_p)
                      Ind2 = Ind2 + mesh(1)
                    enddo
                  endif
                  call MPI_ISend(bdens,NBlock,MPI_grid_real,BNode,1,
     .              MPI_Comm_World,Request,MPIerror)
                  call MPI_Wait(Request,Status,MPIerror)

                elseif (Node.eq.BNode) then
C If this is the Node where the density sub-matrix is, then receive
                  call MPI_IRecv(rho(Ind+1,is),NBlock,MPI_grid_real,
     .              0,1,MPI_Comm_World,Request,MPIerror)
                  call MPI_Wait(Request,Status,MPIerror)
                  Ind = Ind + NBlock

                endif

                if (BNode.ne.0) then
                  call MPI_Barrier(MPI_Comm_World,MPIerror)
                endif
#endif

              enddo

            enddo

          enddo

          enddo

C Close file
        if (Node.eq.0) then
          call io_close( iu )
        endif

#ifdef MPI
C Deallocate density buffer memory
      call de_alloc( bdens, 'bdens', 'read_rho' )
#endif
      call de_alloc( temp, 'temp', 'read_rho' )
      end subroutine read_rho

!------------------------------------------------------------------
      subroutine check_rho(fname,maxp,nspin,nsm,found,overflow)

C integer maxp            : Required first dimension of array rho,
C                           equal to mesh(1)*mesh(2)*mesh(3)
C integer nspin           : Number of spin polarizations (1 or 2)
C integer nsm             : Number of sub-mesh points per mesh point
C                           (not used in this version)  ??
C logical found           : Were data found? 
C logical overflow        : True if maxp or nspin were changed

      character(len=*), intent(in) ::     fname
      integer, intent(in)          ::     nsm
      integer, intent(inout)       ::     maxp, nspin
      logical, intent(out)         ::     found, overflow

      real(dp) :: cell(3,3)
      integer  :: mesh(3)
      integer  :: meshnsm(3), npl, ns, iu, np, npmax
#ifdef MPI
      integer  :: meshl(3)
      logical  :: ltmp
      integer  :: MPIerror
#endif

C Check if input file exists
        if (Node.eq.0) then
          inquire( file=fname, exist=found )
        endif
#ifdef MPI
        call MPI_Bcast(found,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
        if (.not. found) return

C Open file
          if (Node.eq.0) then
            call io_assign( iu )
            open( iu, file=fname, form=fform, status='old' )      

C Read cell vectors and number of points
            if (fform .eq. 'formatted') then
              read(iu,*) cell
              read(iu,*) mesh, ns
            else
              read(iu) cell
              read(iu) mesh, ns
            endif
          endif
#ifdef MPI
          call MPI_Bcast(cell(1,1),9,MPI_double_precision,0,
     .      MPI_Comm_World,MPIerror)
          call MPI_Bcast(mesh,3,MPI_integer,0,MPI_Comm_World,MPIerror)
          call MPI_Bcast(ns,1,MPI_integer,0,MPI_Comm_World,MPIerror)
#endif
          np = mesh(1) * mesh(2) * mesh(3)

C Get local dimensions
          meshnsm(1) = mesh(1)/nsm
          meshnsm(2) = mesh(2)/nsm
          meshnsm(3) = mesh(3)/nsm
#ifdef MPI
          call HowManyMeshPerNode(meshnsm,Node,Nodes,npl,meshl)
#else
          npl = np
#endif

C Check dimensions
          overflow = .false.
          if (ns .gt. nspin) overflow = .true.
          if (npl .gt. maxp) overflow = .true.

#ifdef MPI
C Globalise dimension check
          call MPI_AllReduce(overflow,ltmp,1,MPI_logical,MPI_Lor,
     .      MPI_Comm_World,MPIerror)
          overflow = ltmp
#endif

          if (overflow) then    ! Some processor has npl > maxp ...
#ifdef MPI
C Find largest value of npl
            call MPI_AllReduce(npl,npmax,1,MPI_integer,MPI_Max,
     .        MPI_Comm_World,MPIerror)
#else
            npmax = np
#endif
            maxp = npmax
            nspin = ns
          endif

          if (Node.eq.0) call io_close( iu )

          end subroutine check_rho

      end module m_iorho
