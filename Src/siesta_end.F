      MODULE m_siesta_end
      private
      public :: siesta_end

      CONTAINS

      subroutine siesta_end()
      use alloc, only: alloc_report
      use m_timestamp, only : timestamp
      use siesta_cmlsubs, only: siesta_cml_exit
      use siesta_cml,     only: cml_p, cmlEndModule, mainXML
      use m_wallclock, only : wallclock
      use parallel, only: IOnode
#ifdef MPI
      use mpi_siesta, only: MPI_Finalize
#endif
      implicit none

#ifdef MPI
      integer MPIerror
#endif

!------------------------------------------------------------------------- BEGIN
! Output memory use up to the end of the program
      call printmemory( 6, 1 )

! Print allocation report
      call alloc_report( printNow=.true. )

! Stop time counter
      call timer( 'siesta', 2 )
      call timer( 'all', 3 )

! Print final date and time
      if (IOnode) then
        call timestamp('End of run')
        call wallclock('End of run')
      endif

! Finalize MPI
#ifdef MPI
      call MPI_Finalize( MPIerror )
#endif

      if (cml_p) then
        call cmlEndModule(mainXML)
        call siesta_cml_exit()
      endif      
!--------------------------------------------------------------------------- END
      END subroutine siesta_end
      END MODULE m_siesta_end
