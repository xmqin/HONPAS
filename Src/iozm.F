! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      subroutine iozm( task, cell, vcell, xa, found )
c *******************************************************************
c Saves Z-matrix coordinates and forces
c R. C. Hoft    January 2006
c ********** INPUT **************************************************
c character task*(*) : 'read' or 'write' (or 'READ' or 'WRITE')
c ********** INPUT OR OUTPUT (depending of task) ********************
c real(dp) cell(3,3)  : Unit cell vectors
c real(dp) vcell(3,3)  : Unit cell vector velocities (Parrinello-Rahman)
c real(dp) Zmat(3*na) : Z-matrix coordinates
c real(dp) ZmatForce(3*na) : Z-matrix forces
c real(dp) ZmatForceVar(3*na) : Forces on additional constrained coordinates
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
      use precision
      use parallel,  only : Node
      use files,     only : slabel, label_length
      use zmatrix,   only : Zmat, ZmatForce, ZmatForceVar
      use zmatrix,   only : Zmat_to_Cartesian
      use siesta_geom,  only : na_u
#ifdef MPI
      use mpi_siesta
#endif

      implicit          none

      character(len=*)                    :: task
      logical                             :: found
      real(dp)                            :: cell(3,3)
      real(dp)                            :: vcell(3,3)
      real(dp)                            :: xa(3,na_u)
      external          io_assign, io_close

c Internal variables and arrays
      character(len=label_length+3) :: fname
      integer    ia, iu, iv, ix, k
#ifdef MPI
      integer    MPIerror
#endif

C Find name of file
      fname = trim(slabel) // '.ZM'

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
     .       'iozm: Reading Z-matrix coordinates and forces from file'
            do iv = 1,3
              read(iu,*) (cell(ix,iv),ix=1,3),(vcell(ix,iv),ix=1,3)
            enddo
c           read(iu,*) na_u
            do ia = 1,na_u
              read(iu,*) (Zmat(3*(ia-1)+k),k=1,3),
     .                   (ZmatForce(3*(ia-1)+k),k=1,3)
            enddo
            do ia = 1,3*na_u
              read(iu,*) ZmatForceVar(ia)
            enddo
            call Zmat_to_Cartesian(xa)

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
c         write(iu,*) na
          do ia = 1,na_u
            write(iu,'(6f18.9)') (Zmat(3*(ia-1)+k),k=1,3),
     .                     (ZmatForce(3*(ia-1)+k),k=1,3)
          enddo
          do ia=1,3*na_u
            write(iu,'(f18.9)') ZmatForceVar(ia)
          enddo

C Close file
          call io_close( iu )

        endif
      endif

  999 continue

C If data has been read in then broadcast the values to all Nodes
#ifdef MPI
      call MPI_Bcast(found,1,MPI_logical,0,MPI_Comm_World,MPIerror)
      if (found.and.(task.eq.'read' .or. task.eq.'READ')) then
        call MPI_Bcast(cell(1,1),9,MPI_double_precision,0,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(vcell(1,1),9,MPI_double_precision,0,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(xa(1,1),3*na_u,MPI_double_precision,0,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(Zmat(1),3*na_u,MPI_double_precision,0,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(ZmatForce(1),3*na_u,MPI_double_precision,0,
     .    MPI_Comm_World,MPIerror)
        call MPI_Bcast(ZmatForceVar(1),3*na_u,MPI_double_precision,0,
     .    MPI_Comm_World,MPIerror)
      endif
#endif

      return
      end
