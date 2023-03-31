! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!!@LICENSE

      module sys
!
!     Termination and messaging routines, MPI aware
!
      implicit none

      public :: die      ! Prints an error message and calls MPI_Abort
      public :: bye      ! Prints an error message and calls MPI_Finalize
      public :: message  ! Prints a message string if node==0

      private

      CONTAINS

      subroutine message(level,str)

      use parallel, only : Node

      character(len=*), intent(in)  :: level
      ! One of INFO, WARNING, FATAL
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
!
!--------------------------------------------------
      subroutine die(str)

      use parallel, only : Node
      use siesta_cml
#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in), optional   :: str

      external ::  io_assign, io_close
      integer  ::  lun
#ifdef MPI
      integer MPIerror
#endif

! Even though formally (in MPI 1.X), only the master node
! can do I/O, in those systems that allow it having each
! node state its complaint can be useful.

!!                                       if (Node.eq.0) then
      if (present(str)) then
         write(6,'(a)') trim(str)
         write(0,'(a)') trim(str)
      endif
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
         If (cml_p) Then
            Call cmlFinishFile(mainXML)
         Endif                  !cml_p
      endif

#ifdef MPI
      call MPI_Abort(MPI_Comm_World,1,MPIerror)
      stop
#else
      call pxfabort()
#endif
      end subroutine die

!---------------------------------------------------------
      subroutine bye(str)

      use parallel, only : Node
      use siesta_cml
#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in), optional   :: str

#ifdef MPI
      integer rc
#endif

      if (Node.eq.0) then
         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
         write(6,'(a)') 'Requested End of Run. Bye!!'
         call pxfflush(6)
         If (cml_p) Then
            Call cmlFinishFile(mainXML)
         Endif                  !cml_p
      endif

#ifdef MPI
      call MPI_Finalize(rc)
      stop
#else
      stop
#endif
      end subroutine bye

      end module sys

