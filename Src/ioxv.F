! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!!**    There should be two separate routines for reading and writing

      module m_ioxv

      logical, public, save   :: xv_file_read = .false.

      public :: ioxv

      CONTAINS

      subroutine ioxv( task, cell, vcell,
     .                 na, isa, iza, xa, va, found )
c *******************************************************************
c Saves positions and velocities.
c J.M.Soler. July 1997.
c ********** INPUT **************************************************
c character task*(*) : 'read' or 'write' (or 'READ' or 'WRITE')
c ********** INPUT OR OUTPUT (depending of task) ********************
c real*8  cell(3,3)  : Unit cell vectors
c real*8  vcell(3,3) : Unit cell vector velocities (Parrinello-Rahman)
c integer na         : Number of atoms
c integer isa(na)    : Atomic species index
c integer iza(na)    : Atomic numbers
c real*8  xa(3,na)   : Atomic positions
c real*8  va(3,na)   : Atomic velocities
c ********** OUTPUT *************************************************
c logical found      : Has input file been found
c                      (only for task='read')
c ********** UNITS **************************************************
c Units are arbitrary, but the use with task='write' and task='read'
c must be consistent
c *******************************************************************

C
C  Modules
C
      use files,      only : slabel, label_length
      use parallel,   only : Node
      use precision,  only : dp
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

      logical           found
      integer           na, isa(na), iza(na)
      real(dp)          cell(3,3), va(3,na), vcell(3,3), xa(3,na)
      character(len=*) ::         task

      external          io_assign, io_close

c Internal variables and arrays
      character(len=label_length+3) ::   fname

      integer    ia, iu, iv, ix
#ifdef MPI
      integer    MPIerror
#endif

C     Find name of file
      fname = trim(slabel) // '.XV'
      
C Only do reading and writing for IOnode
      if (Node.eq.0) then

C Choose between read or write
        if (task.eq.'read' .or. task.eq.'READ') then

C Check if input file exists
          inquire( file=fname, exist=found )
          if (found) then

C Open file
            call io_assign( iu )
            open( iu, file=fname, status='old' )      

C Read data
            write(6,'(/,a)') 
     .       'ioxv: Reading coordinates and velocities from file'
            do iv = 1,3
              read(iu,*) (cell(ix,iv),ix=1,3),(vcell(ix,iv),ix=1,3)
            enddo
            read(iu,*) na
            do ia = 1,na
              read(iu,*)
     .          isa(ia),iza(ia),(xa(ix,ia),ix=1,3),(va(ix,ia),ix=1,3)
            enddo

C Close file
            call io_close( iu )

          else
C If input file not found, go to exit point
            goto 999
          endif

        elseif (task.eq.'write' .or. task.eq.'WRITE') then

C Open file
          call io_assign( iu )
          open( iu, file=fname, form='formatted', status='unknown' )

C Write data on file
          write(iu,'(2(3x,3f18.9))')
     .      ((cell(ix,iv),ix=1,3),(vcell(ix,iv),ix=1,3),iv=1,3)
          write(iu,*) na
          do ia = 1,na
            write(iu,'(i3,i6,3f18.9,3x,3f18.9)')
     .        isa(ia),iza(ia),(xa(ix,ia),ix=1,3),(va(ix,ia),ix=1,3)
          enddo

          call io_close( iu )

        endif
      endif

  999 continue

      if (task.eq.'read' .or. task.eq.'READ') then
         
C If data has been read in then broadcast the values to all Nodes
#ifdef MPI
         call MPI_Bcast(found,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
         if (found) then
            xv_file_read = .true.
#ifdef MPI
            call MPI_Bcast(na,1,MPI_integer,0,MPI_Comm_World,MPIerror)
            call MPI_Bcast(cell(1,1),9,MPI_double_precision,0,
     .           MPI_Comm_World,MPIerror)
            call MPI_Bcast(vcell(1,1),9,MPI_double_precision,0,
     .           MPI_Comm_World,MPIerror)
            call MPI_Bcast(isa,na,MPI_integer,0,MPI_Comm_World,MPIerror)
            call MPI_Bcast(iza,na,MPI_integer,0,MPI_Comm_World,MPIerror)
            call MPI_Bcast(xa(1,1),3*na,MPI_double_precision,0,
     .           MPI_Comm_World,MPIerror)
            call MPI_Bcast(va(1,1),3*na,MPI_double_precision,0,
     .           MPI_Comm_World,MPIerror)
#endif
         endif
      endif

      end subroutine ioxv

      end module m_ioxv

