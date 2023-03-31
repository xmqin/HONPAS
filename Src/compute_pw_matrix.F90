! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
subroutine compute_pw_matrix( nncount, bvectorsfrac )
!
! In this subroutine we compute the matrix elements of the plane wave,
! for all the wave vectors that connect a given k-point to its nearest
! neighbours
!
! ------------------------------ INPUT -----------------------------------------
! integer nncount            : The number of nearest neighbours belonging to 
!                                each k-point of the Monkhorst-Pack mesh
! real(dp) bvectorsfrac(:,:) : Vectors that connect each mesh k-point 
!                              to its nearest neighbours.
! ------------------------------ OUTPUT ----------------------------------------
! Matrix elements of the plane wave for all the required wave vectors
! connecting a given k-point to its nearest neighbours
! This array (delkmatgen) is defined in the module m_planewavematrixvar,
! file m_planewavematrixvar.F90
! ------------------------------ BEHAVIOUR -------------------------------------
! First of all, we allocate the arrays where the matrix elements of the
! plane waves will be stored
!
! Loop on all the neighbour k-points
!
!   For each neighbour, compute the vector that connects each mesh point
!   to its nearest neighbour in bohr^-1.
!   Remember that bvectorsfrac are read in reduced units, so we have
!   to multiply then by the reciprocal lattice vector.
!   This is done in the function getkvector.
!
!   Call to the subroutine that computes the matrix elements of the plane
!   wave in the grid (planewavematrix)
! 
!   Store the matrix elements of the plane wave for a given wave vector
!   (array delkmat)
!   in a permanent array,
!   (array delkmatgen)
!   So delkmat can be rewritten for the next wave vector.
!
! ------------------------------ UNITS -----------------------------------------
!   bvectorsfrac are given in reduced units.
! ------------------------------------------------------------------------------

! Used module procedures
  use alloc,             only: re_alloc        ! Re-allocation routine
  use m_planewavematrix, only: planewavematrix ! To compute the plane wave
                                               !    matrix elements in the grid

! Used module variables
  use precision,            only: dp         ! Real douple precision type     
  use m_planewavematrixvar, only: delkmat    ! Matrix elements of a plane wave
  use m_planewavematrixvar, only: delkmatgen ! Matrix elements of a plane wave
                                             !   (saved)
  use sparse_matrices,      only: maxnh      ! Maximum number of orbitals
                                             !   overlapping
  use parallel,             only: IOnode     ! Input/Output node

!! For debugging
!  use m_writedelk,          only: write_delk
!  use atomlist,             only: no_u, no_s, iaorb, iphorb
!  use atomlist,             only: indxuo, qtot
!  use sparse_matrices,      only: maxnh, listh, listhptr, numh, xijo, S
!  use m_spin,               only: nspin
!  use siesta_options,       only: temp
!  use sys,                  only: die
!! End debugging

  implicit none

! The number of nearest neighbours belonging to 
!   each k-point of the Monkhorst-Pack mesh
  integer,  intent(in)          :: nncount           
! The vectors b that connect each mesh-point k to its nearest neighbours
  real(dp), intent(in)          :: bvectorsfrac(3,nncount) 

!
! Internal variables
!
  integer  :: inn
  real(dp) :: bvector(3)

!! For debugging
!  logical  :: onlygamma
!! End debugging

  if( IOnode )              &
 &        write(6,'(/,a)')  &
 &       'compute_pw_matrix: Computing the matrix elements of a plane wave'

! First of all, we allocate the arrays where the matrix elements of the
! plane waves will be stored
  if ( .not. associated(delkmatgen) ) then
    nullify(delkmatgen)
    call re_alloc(delkmatgen, 1, nncount, 1, maxnh,                  &
 &                name='delkmatgen', routine='compute_pw_matrix')
  endif
  if ( .not. associated(delkmat) ) then
    nullify(delkmat)
    call re_alloc(delkmat, 1, maxnh,  &
 &                name='delkmat',routine='compute_pw_matrix')
  endif

! Loop on neighbour k-points
  do inn = 1, nncount

!   For each neighbour, compute the vector that connects each mesh point
!   to its nearest neighbour in bohr^-1 (done in the subroutine getkvector).
!   Remember that bvectorsfrac are read in reduced units, so we have
!   to multiply then by the reciprocal lattice vector.
    call getkvector( bvectorsfrac(:,inn), bvector )

!!$    if( IOnode ) then  
!!$      write(6,'(a)')         &
!!$ &      'compute_pw_matrix: Vector connecting a k-point with its neighbours'
!!$      write(6,'(a,3f12.5)')  &
!!$ &      'compute_pw_matrix:', bvector(:)
!!$    endif

!   Call to the subroutine that computes the matrix elements of the plane
!   wave in the grid (planewavematrix)
    call planewavematrix( -1, bvector )

!   Store the matrix elements of the plane wave for a given wave vector
!   (array delkmat)
!   in a permanent array,
!   (array delkmatgen)
!   So delkmat can be rewritten for the next wave vector.
    delkmatgen(inn,:) = delkmat(:)

!! For debugging
!      onlygamma = .false.
!      call write_delk( onlygamma, no_u, no_s, nspin, indxuo, maxnh, &
! &                     numh, listhptr, listh, delkmat, S, qtot, temp, xijo)
!      call die()
!! End debugging

  enddo

end subroutine compute_pw_matrix

subroutine getkvector(frac, vector)
!
! Converts from fractional to Bohr^-1 coordinates in reciprocal space.
!

use m_siesta2wannier90, only : reclatvec ! Reciprocal lattice vectors
use precision,          only : dp        ! Real double precision type

integer              :: i,j
real(dp),intent(in)  :: frac(3)
real(dp),intent(out) :: vector(3)


  do i = 1, 3
    vector(i) = 0.0_dp
    do j = 1, 3
      vector(i) = frac(j)*reclatvec(i,j) + vector(i)
    enddo
  enddo

end subroutine getkvector

