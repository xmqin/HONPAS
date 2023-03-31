! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_eo

  use precision, only: dp
	implicit none

  public

  real(dp), pointer, save :: eo(:,:,:)  ! Hamiltonian eigenvalues
  real(dp), pointer, save :: qo(:,:,:)  ! Occupations of eigenstates

end module m_eo





