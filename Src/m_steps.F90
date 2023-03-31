! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_steps
  implicit none
  public
  integer:: inicoor           ! Initial step in geommetry iteration
  integer:: fincoor           ! Final step in geommetry iteration
  integer:: istp              ! Geommetry iteration step starting in istp=1
  logical:: final=.false.     ! Last geometry step?

end module m_steps
