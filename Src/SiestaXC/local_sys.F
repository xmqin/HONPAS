      module sys
      public :: die
      contains
       subroutine die(str)
       character(len=*), intent(in), optional :: str
       if (present(str)) then
          call local_die(str)
       else
          call local_die('')
       endif
       end subroutine die
      end module sys
     
         
! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!--------------------------------------------------
! Stand-alone 'die' routine for use by libraries and
! low-level modules.
!
! Each program using the module or library needs to
! provide a routine with the proper interface, but
! accomodating the needs and conventions of the program.
! For example, in Siesta:
!
!   - The use of a Siesta-specific 'mpi_siesta' module.
!   - The need to have the pxf functionality.
!   - The use of 'unit 6' as output.
!
! Routines using this functionality should include
! the following
!
!     interface
!      subroutine die(str)
!      character(len=*), intent(in)  :: str
!      end subroutine die
!     end interface
!
!------------------------------------------------------

      subroutine local_die(str)

      use parallel, only : Node
#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in)  :: str
      
      external ::  io_assign, io_close
      integer  ::  lun
#ifdef MPI
      integer MPIerror
#endif

! Even though formally (in MPI 1.X), only the master node
! can do I/O, in those systems that allow it having each
! node state its complaint can be useful.

!!                                       if (Node.eq.0) then
      write(6,'(a)') trim(str)
      write(0,'(a)') trim(str)
      write(6,'(a,i4)') 'Stopping Program from Node: ', Node
      write(0,'(a,i4)') 'Stopping Program from Node: ', Node
!!                                       endif
      if (Node .eq. 0) then
         call io_assign( lun )
         open(lun,file="MESSAGES",status="unknown",
     $        position="append",action="write")
         write(lun,"(a)") 'FATAL: ' // trim(str)
         call io_close(lun)
         call pxfflush(6)
         call pxfflush(0)
      endif

#ifdef MPI
      call MPI_Abort(MPI_Comm_World,1,MPIerror)
      stop
#else
      call pxfabort()
#endif
      end subroutine local_die
      
      subroutine message(level,str)

      use parallel, only : Node

      ! One of INFO, WARNING, FATAL
      character(len=*), intent(in)  :: level

      character(len=*), intent(in)  :: str

      external ::  io_assign, io_close
      integer  ::  lun
      
      if (Node .eq. 0) then
         write(6,'(a)') trim(str)
         write(0,'(a)') trim(str)
         call io_assign(lun)
         open(lun,file="MESSAGES",status="unknown",
     $        position="append",action="write")
         write(lun,"(a)") trim(level) // ": " // trim(str)
         call io_close(lun)
         call pxfflush(6)
         call pxfflush(0)
      endif

      end subroutine message

      subroutine die(str)
      character(len=*), intent(in)  :: str
      call local_die(str)
      end subroutine die
