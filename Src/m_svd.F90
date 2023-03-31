! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_svd
public :: solve_with_svd
CONTAINS
subroutine solve_with_svd(ain,rhs,c,info,rcond,rank_out,sigma)
  use precision, only: dp

  real(dp), intent(in)  :: ain(:,:)
  real(dp), intent(in)  :: rhs(:)
  real(dp), intent(out) :: c(:)
  integer, intent(out)  :: info
  real(dp), intent(in), optional   :: rcond
  integer, intent(out), optional  :: rank_out
  real(dp), intent(out), optional :: sigma(:)

  real(dp), allocatable :: a(:,:), s(:), work(:), b(:)

  integer, parameter :: nb = 64
  integer :: n, lwork

  REAL(DP) RCOND_local, RNORM
  INTEGER  ::    I, J, M, RANK, LDA

  REAL(DP), external ::  DNRM2
  EXTERNAL ::      DGELSS

!
  n = size(ain,dim=1)
  lda = n
  m = n
  allocate(a(n,n))
  a = ain
  allocate(b(n),s(n))

  lwork = 3*n+nb*(n+m)
  allocate(work(lwork))

  b(:) = rhs(:)

!        Choose RCOND to reflect the relative accuracy of the input data

  if (present(rcond)) then
     rcond_local = rcond
  else
     RCOND_local = 1.0e-6_dp
  endif
  ! Singular values s_i < rcond*s_1 will be neglected for
  ! the estimation of the rank.

!        Solve the least squares problem min( norm2(b - Ax) ) for the x
!        of minimum norm.

         CALL DGELSS(M,N,1,A,LDA,B,M,S,RCOND_local,RANK,WORK,LWORK,INFO)
!
         IF (INFO.EQ.0) THEN

            c(1:n) = b(1:n)

            if (present(sigma)) then
               sigma(1:n) = s(1:n)
            endif
            if (present(rank_out)) then
               rank_out = rank
            endif

         END IF

         deallocate(a,s,work,b)

END subroutine solve_with_svd
end module m_svd
