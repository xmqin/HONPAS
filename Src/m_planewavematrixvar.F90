! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! ---------------------------------------------------------------------
! MODULE m_planewavematrixvar
! In this module we define all the variables required
! to compute the matrix elements of a plane wave in a basis of
! numerical atomic orbitals
! Implemented by J. Junquera, June-July 2013
! ---------------------------------------------------------------------

module m_planewavematrixvar

  use precision, only: dp

  implicit none

  complex(dp), pointer, save ::   delkmatgen(:,:)

  complex(dp), pointer, save ::   delkmat(:)

  real(dp)                   ::   wavevector(3)

! ---------------------------------------------------------------------
! complex(dp) delkmat(:)     : matrix elements of a plane wave
! <\phi_{\mu}|exp^( isigneikr * i * \vec{kptpw} \cdot \vec{r} )|\phi_{\nu}>
!                              This pointer has to be allocated in 
!                              the calling routine.
!                              This is a sparse matrix, whose structure
!                              follows the same scheme as the
!                              hamiltonian or overlap matrix elements.
! complex(dp) delkmatgen(:,:): The matrix elements of a planewave are not
!                              self-consistent. They can be computed
!                              once for a given k-point and stored in
!                              delkmatgen.
!                              The first index refers to the number of 
!                              different k-points for which the planewave
!                              will be computed.
!                              The second index is the index of the
!                              sparse matrix.
!                              This pointer has to be allocated in 
!                              the calling routine.
! real(dp)    wavevector(3)  : Wave vector of the plane wave.
! ---------------------------------------------------------------------

endmodule m_planewavematrixvar



