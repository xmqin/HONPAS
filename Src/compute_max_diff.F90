! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_compute_max_diff
  use precision,         only: dp

  !> Temporary for storing the old maximum change
  real(dp), public, save :: dDmax_current

  interface compute_max_diff
     module procedure compute_max_diff_1d
     module procedure compute_max_diff_2d
  end interface compute_max_diff
  public   :: compute_max_diff
  
contains
  
  subroutine compute_max_diff_2d(X1,X2, max_diff)
    
#ifdef MPI
    use m_mpi_utils, only: globalize_max
#endif
    
    real(dp), intent(in)  :: X1(:,:), X2(:,:)
    real(dp), intent(out) :: max_diff
    
    integer :: n1, n2
    integer :: i1, i2
#ifdef MPI
    real(dp) :: buffer1
#endif
    
    n1 = size(X1, 1)
    n2 = size(X1, 2)
    if ( size(X2, 1) /= n1 ) then
       call die('compute_max_diff: Sizes of the arrays are not &
            &conforming (1-D)')
    end if
    if ( size(X2, 2) /= n2 ) then
       call die('compute_max_diff: Sizes of the arrays are not &
            &conforming (2-D)')
    end if
    
    max_diff = 0.0_dp
!$OMP parallel do default(shared), private(i2,i1), &
!$OMP& reduction(max:max_diff)
    do i2 = 1 , n2
       do i1 = 1 , n1
          max_diff = max(max_diff, abs(X1(i1,i2) - X2(i1,i2)) )
       end do
    end do
!$OMP end parallel do
    
#ifdef MPI
    ! Ensure that max_diff is the same on all nodes 
    call globalize_max(max_diff, buffer1)
    max_diff = buffer1
#endif

    dDmax_current = max_diff

  end subroutine compute_max_diff_2d

  subroutine compute_max_diff_1d(X1, X2, max_diff)
    
#ifdef MPI
    use m_mpi_utils, only: globalize_max
#endif
    
    real(dp), intent(in)  :: X1(:), X2(:)
    real(dp), intent(out) :: max_diff
    
    integer :: n1
    integer :: i1
#ifdef MPI
    real(dp) :: buffer1
#endif
    
    n1 = size(X1, 1)
    if ( size(X2, 1) /= n1 ) then
       call die('compute_max_diff: Sizes of the arrays are not &
            &conforming (1-D)')
    end if
    
    max_diff = 0.0_dp
!$OMP parallel do default(shared), private(i1), reduction(max:max_diff)
    do i1 = 1 , n1
       max_diff = max(max_diff, abs(X1(i1) - X2(i1)) )
    end do
!$OMP end parallel do
    
#ifdef MPI
    ! Ensure that max_diff is the same on all nodes 
    call globalize_max(max_diff, buffer1)
    max_diff = buffer1
#endif

    dDmax_current = max_diff

  end subroutine compute_max_diff_1d
  
end module m_compute_max_diff
