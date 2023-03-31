! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_siesta_end
      private
      public :: siesta_end

      CONTAINS

      subroutine siesta_end( )
      use fdf,             only : fdf_shutdown
      use precision,       only : dp
      use alloc,           only : alloc_report
      use siesta_cml,      only : cml_p, mainXML, cmlEndModule
      use siesta_cmlsubs,  only : siesta_cml_exit
      use siesta_master,   only : siesta_subroutine  ! Is siesta a subroutine?
#ifdef DEBUG_XC
      USE siestaXC,        only : closeDebugOutputFile
#endif /* DEBUG_XC */
      use m_timestamp,     only : timestamp
      use m_wallclock,     only : wallclock, wallclock_shutdown
      use parallel,        only : IOnode, SIESTA_worker, Node
#ifdef MPI
      use mpi_siesta
#endif
      use moreMeshSubs,    only : resetMeshDistr
      use densematrix,     only : resetDenseMatrix
      use sparse_matrices, only : resetSparseMatrices
      use m_fixed,         only : resetFixedPointers
      use band,            only : reset_bands
      use meshphi,         only : resetMeshPhi
      use meshdscf,        only : resetDscfPointers
      use m_dscfcomm,      only : resetdscfComm
      use m_rhog,          only : resetRhog
      ! NP, new mixing
      use m_mixing_scf,    only : mixers_scf_reset

#ifdef SIESTA__FLOOK
      use flook_siesta,    only : slua_call, LUA_FINALIZE
      use flook_siesta,    only : slua_close
      use siesta_options,  only : LUA
#endif

      use m_diag, only: diag_exit

#ifdef TRACING
      use extrae_eventllist, only : writeList, deleteList, eventlist
#endif
      use m_io,            only : io_assign, io_close
      use m_cite, only: announce_citations
      use atomlist,        only : reset_atomlists
      use m_new_matel,     only : new_MATEL
      use listsc_module,   only : LISTSC_RESET

#ifdef HAVE_LIBINT
      use siesta_options, only: hfx_wanted
      use nao2gto_data, only: hfx_gtos, hfx_libderiv, hfx_libint
      use nao2gto_wrappers, only: nao2gto_libderiv_free,
     &                            nao2gto_libint_free
#endif

      implicit none

      real(dp) :: R12(3), S12, DSDR(3)
      integer :: iu
#ifdef MPI
      integer MPIerror
#endif
!------------------------------------------------------------------------- BEGIN
      if (SIESTA_worker) then

      if ( ionode ) call announce_citations()

#ifdef SIESTA__FLOOK
      ! Final call for the end-user
      call slua_call(LUA,LUA_FINALIZE)
      ! Clean up the lua environment
      call slua_close(LUA)
#endif

!     Free allocated memory
      ! Reset matrix element tables
      R12 = 0._dp
      call new_MATEL('S', 0, 0, R12, S12, DSDR)
      call LISTSC_reset()
      call reset_atomlists()
      call resetSparseMatrices( )
      call resetMeshDistr( )
      call resetdscfComm( )
      call resetDenseMatrix( )
      call resetFixedPointers( )
      call reset_bands( )
      call resetMeshPhi( )
      call resetDscfPointers( )
      call resetRhog()
      call mixers_scf_reset()

      ! Clean diagonalization
      call diag_exit()
     
#ifdef HAVE_LIBINT
      if ( hfx_wanted ) then
        call nao2gto_libint_free(hfx_libint)
        call nao2gto_libderiv_free(hfx_libderiv)
        deallocate(hfx_gtos)
      end if
#endif
!     Output memory use up to the end of the program
!!      call printmemory( 6, 1 )

!     Print allocation report
      call alloc_report( printNow=.true. , shutdown=.true.)

#ifdef DEBUG_XC
! Close debug outup file
      call closeDebugOutputFile()
#endif /* DEBUG_XC */

!     Stop time counter
      call timer( 'siesta', 2 )
      call timer( 'all', 3 )

#ifdef TRACING
!     Write extrae user function list
      if (IOnode) then
        call writeList(eventlist)
      end if
      call deleteList(eventlist)
#endif /* TRACING */


!     Print final date and time
      if (IOnode) then
        call timestamp('End of run')
        call wallclock('End of run')
        ! Ensure clock file is closed
        call wallclock_shutdown()
      endif

      if (cml_p) then
        call cmlEndModule(mainXML)
        call siesta_cml_exit()
      endif      
      endif  ! SIESTA_worker

! Finalize MPI, unless siesta is running as a subroutine
      !
      ! Add stamp to MESSAGES file
      call message('INFO','Job completed')

      ! Clean the fdf-files
      call fdf_shutdown()

#ifdef MPI
      call MPI_Barrier(MPI_Comm_World, MPIerror)
      MPIerror = MPI_SUCCESS
      if (.not.siesta_subroutine) call MPI_Finalize( MPIerror )
#endif

      if (IOnode) then
         call io_assign(iu)
         open(iu, file='0_NORMAL_EXIT', form='FORMATTED')
#ifdef MPI
         write(iu,'(a,i0)') 'SIESTA completed, MPI exit: ',
     &        MPIerror
#else
         write(iu,'(a)') 'SIESTA completed successfully'
#endif
         call io_close(iu)
      endif
!--------------------------------------------------------------------------- END
      END subroutine siesta_end
      END MODULE m_siesta_end
