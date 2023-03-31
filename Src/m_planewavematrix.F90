! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! ---------------------------------------------------------------------
! MODULE m_planewavematrix
! In this module we define all the variables and subroutines required
! to compute the matrix elements of a plane wave in a basis of
! numerical atomic orbitals
! Implemented by J. Junquera, June-July 2013
! ---------------------------------------------------------------------

module m_planewavematrix

  CONTAINS

  subroutine planewavematrix( isigneikr, kptpw )

! ------------------------------ INPUT ---------------------------------------
!   integer     : isigneikr    ! Sign of the exponent of the planewave
!   real(dp)    : kptpw(3)     ! Wavevector of the planewave
!                                (Units in Bohr^-1) 
! ------------------------------ OUTPUT --------------------------------------
!   complex(dp) : delkmat(:)   ! Matrix element of the plane wave in sparse
!                                format (same as in the Hamiltonian or Overlap
!                                matrices). 
!                                Note that the output is not part of the 
!                                variables of the subroutine, but is defined
!                                as a variable of this module.
! ------------------------------ BEHAVIOUR -----------------------------------
!   This subroutine calls delk_wrapper to compute the matrix elements of a
!   plane wave, with wave vector kptpw and sign isigneikr.
!   The matrix elements are computed in a subroutine called delk,
!   essentially a copy of vmat.
!   delk is only called if isigneikr = +1 or -1. 
! ----------------------------------------------------------------------------

!   Used module variables

    use precision,        only: dp       ! Real double precision type
    use siesta_geom,      only: isa      ! Species index of each atom
    use atomlist,         only: no_s     ! Number of orbitals in the supercell
    use atomlist,         only: iaorb    ! Atom to which each orbital belongs
    use atomlist,         only: iphorb   ! Orbital index within atom
    use atomlist,         only: no_l     ! Number of orbitals in the local node
    use atomlist,         only: no_u     ! Number of orbitals in the unit cell

    use sparse_matrices,  only: maxnh    ! Maximum number of orbitals
                                         !   interacting
    use sparse_matrices,  only: numh     ! Number of nonzero elements of each
                                         !   row of the hamiltonian matrix 
    use sparse_matrices,  only: listh    ! Nonzero hamiltonian-matrix elements
    use sparse_matrices,  only: listhptr ! Pointer to start of each row of H

!   Used module procedures
    use m_dhscf,          only: delk_wrapper    ! Mesh subroutine
    use m_planewavematrixvar, only: wavevector ! Wave vector of the plane wave
 
    implicit none

    integer,  intent(in) :: isigneikr   ! Sign of the exponent of the planewave
    real(dp), intent(in) :: kptpw(3)    ! Wavevector of the planewave
                                        !   (Units in Bohr^-1) 

!   Initialize the wave vector
    wavevector = kptpw

    call delk_wrapper(isigneikr, no_s, maxnh,  &
                            numh, listhptr, listh,  &
                            no_l, no_u, iaorb, iphorb, isa )
 
  endsubroutine planewavematrix

endmodule m_planewavematrix



