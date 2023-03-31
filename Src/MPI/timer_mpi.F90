! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
MODULE timer_mpi_m

! Disconnectable interface to siesta's timer routine. 
! J.M.Soler. May.2009

public :: timer_mpi

CONTAINS

  SUBROUTINE timer_mpi( name, opt )
    character(len=*), intent(in):: name
    integer,          intent(in):: opt

#ifdef MPI_TIMING
    external timer
    call timer( name, opt )
#endif

  END SUBROUTINE timer_mpi

END MODULE timer_mpi_m
