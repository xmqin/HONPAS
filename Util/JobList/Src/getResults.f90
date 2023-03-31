! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
program getResults

! Collects results of siesta jobs contained in a list, read from standard input
! J.M.Soler. Nov.2012

  use jobList, only: collectResults => getResults

  implicit none
  integer :: unit

  unit = 5
  call collectResults( unit )

end program getResults


