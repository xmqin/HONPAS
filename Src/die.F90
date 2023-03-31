! Stand-alone 'die' routine
!--------------------------------------------------
      subroutine die(str)

      use parallel, only : Node

#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in)  :: str

#ifdef MPI
      integer MPIerror
#endif

! Even though formally (in MPI 1.X), only the master node
! can do I/O, in those systems that allow it having each
! node state its complaint can be useful.

      write(6,'(a)') trim(str)
      write(0,'(a)') trim(str)
      write(6,'(a,i4)') 'Stopping Program from Node: ', Node
      write(0,'(a,i4)') 'Stopping Program from Node: ', Node

      if (Node .eq. 0) then
         call pxfflush(6)
         call pxfflush(0)
      endif

#ifdef MPI
      call MPI_Abort(MPI_Comm_World,1,MPIerror)
      stop
#else
      call pxfabort()
#endif
      end subroutine die

