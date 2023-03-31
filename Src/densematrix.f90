! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module densematrix
!
!  Contains the dense matrix arrays used within SIESTA
!
  use precision
  
  implicit none

  private
  
  ! Ensure they are initialially nullified
  real(dp), public, pointer :: Haux(:) => null()
  real(dp), public, pointer :: Saux(:) => null()
  real(dp), public, pointer :: psi(:) => null()

  public :: allocDenseMatrix
  public :: resetDenseMatrix

contains

  subroutine allocDenseMatrix(nHaux, nSaux, npsi)
    use alloc, only : re_alloc
    integer, intent(in) :: nHaux, nSaux, npsi

    !> If the arrays are already allocated with the same
    !> bounds nothing will be done
    call re_alloc(Haux, 1, nHaux, 'Haux', 'densematrix', copy=.false., shrink=.false.)
    call re_alloc(Saux, 1, nSaux, 'Saux', 'densematrix', copy=.false., shrink=.false.)
    call re_alloc(psi, 1, npsi, 'psi', 'densematrix', copy=.false., shrink=.false.)

  end subroutine allocDenseMatrix

  !> Deallocates auxiliary arrays.
  !> Note that it is safe to call the routine even if
  !> (some) arrays are not associated. Nothing will be
  !> done in that case.
  subroutine resetDenseMatrix(dealloc_Haux, dealloc_Saux, dealloc_psi)
    use alloc, only : de_alloc

    !> This flag is used in connection with the OMM
    !> module: it needs diagon-computed eigenvectors
    !> as seeds for the first few iterations.
    !> [[diagon]] will not deallocate psi in that case
    logical, intent(in), optional :: dealloc_Haux, dealloc_Saux, dealloc_psi

    logical :: ldealloc

    ldealloc = .true.
    if ( present(dealloc_Haux) ) ldealloc = dealloc_Haux
    if ( ldealloc ) then
      call de_alloc(Haux, 'Haux', 'densematrix')
      nullify(Haux)
    end if

    ldealloc = .true.
    if ( present(dealloc_Saux) ) ldealloc = dealloc_Saux
    if ( ldealloc ) then
      call de_alloc(Saux, 'Saux', 'densematrix')
      nullify(Saux)
    end if

    ldealloc = .true.
    if ( present(dealloc_psi) ) ldealloc = dealloc_psi
    if ( ldealloc ) then
      call de_alloc(psi, 'psi', 'densematrix')
      nullify(psi)
    end if
    
  end subroutine resetDenseMatrix

end module densematrix
