! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code has been fully implemented by:
! Nick Papior, 2014
!
! Please attribute the original author in case of dublication
subroutine tbt_end()

#ifdef MPI
  use mpi_siesta, only : MPI_Finalize
#endif
  use parallel, only : Node
  use alloc, only   : alloc_report
  use m_timestamp, only : timestamp
  use m_wallclock, only : wallclock
#ifdef NCDF_4
  use m_tbt_delta, only: delete_delta
  use m_tbt_dH, only: dH
  use m_tbt_dSE, only: dSE
#endif

#ifdef MPI
  integer :: MPIerror
#endif

#ifdef NCDF_4
  call delete_delta(dH)
  call delete_delta(dSE)
#endif
  
  ! Stop time counter
  call timer( 'tbtrans', 2 )
  call timer( 'all', 3 )

  ! Print allocation report
  call alloc_report( printNow=.true. )

  if ( Node == 0 ) then
     call timestamp("End of run")
     call wallclock("End of run")
  end if

#ifdef MPI
  call MPI_Finalize(MPIerror)
#endif

end subroutine tbt_end
