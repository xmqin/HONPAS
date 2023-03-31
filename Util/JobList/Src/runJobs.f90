! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
program runJobs

! Runs siesta jobs contained in a list, read from standard input
! J.M.Soler. Nov.2012

  use jobList, only: runFile => runJobs

  implicit none
  integer :: unit

  unit = 5
  call runFile( unit )

end program runJobs


