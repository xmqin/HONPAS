! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_broyddj

!     Accelerated and thrifty Broyden scheme.
!     Reference: D.D. Johnson, PRB 38, 12807 (1988)
!     Full update as in Kresse and Furthmuller, Comp. Mat. Sci 6, 15 (1996)

!     Alberto Garcia (wdpgaara@lg.ehu.es), May 2005

use precision, only: dp
use precision, only: wp=>broyden_p       ! Precision of work arrays

use m_mpi_utils, only: Globalize_sum
use parallel, only: ionode
use alloc, only: re_alloc, de_alloc

use sys, only: message, die

implicit none

type, public :: broyden_t
private
   integer  :: n
   integer  :: maxit
   integer  :: it
   real(dp) :: wprime
   real(dp) :: jinv0
   real(dp), dimension(:), pointer                 :: w
   real(dp), dimension(:,:), pointer               :: dFdF
   real(wp), dimension(:,:), pointer               :: u
   real(wp), dimension(:,:), pointer               :: dF
   logical  :: setup
   logical  :: cycle_on_maxit
   logical  :: debug
end type broyden_t

public :: broyden_reset, broyden_step , broyden_destroy, broyden_init
public :: broyden_is_setup

private :: dlinds
private 

CONTAINS

function broyden_is_setup(br) result(is_setup)
type(broyden_t), intent(inout) :: br
logical                        :: is_setup

is_setup = br%setup
end function broyden_is_setup

!--------------------------------------------------------------
subroutine broyden_reset(br,n,maxit,cycle_on_maxit, &
                         jinv0,wprime)
!     Jinit is the starting INVERSE "spring constant".
!     wprime is the (small) initial weight
!     Typical value: 0.01

type(broyden_t), intent(inout) :: br
integer, intent(in)  ::  n
integer, intent(in)  ::  maxit
logical, intent(in), optional       :: cycle_on_maxit
real(dp), intent(in) ::  jinv0, wprime

br%cycle_on_maxit = .false.
if (present(cycle_on_maxit)) then
   br%cycle_on_maxit = cycle_on_maxit
endif
   
br%n = n
br%maxit = maxit
br%it = -1             ! Mark as just reset
br%jinv0 = jinv0
br%wprime = wprime
br%setup = .true.

! if (associated(br%dF)) deallocate(br%dF)
! allocate(br%dF(1:n,0:maxit))
! if (associated(br%u)) deallocate(br%u)
! allocate(br%u(1:n,0:maxit))
! if (associated(br%w)) deallocate(br%w)
! allocate(br%w(0:maxit))
! if (associated(br%dFdF)) deallocate(br%dFdF)
! allocate(br%dFdF(0:maxit,0:maxit))

call re_alloc(br%dF,i1min=1,i1max=n, &
              i2min=0, i2max=maxit, copy=.false., shrink=.false.)
call re_alloc(br%u ,i1min=1,i1max=n, &
              i2min=0, i2max=maxit, copy=.false., shrink=.false.)
call re_alloc(br%w ,i1min=0, i1max=maxit, copy=.false., shrink=.false.)
call re_alloc(br%dFdF ,i1min=0,i1max=maxit, &
              i2min=0, i2max=maxit, copy=.false., shrink=.false.)

end subroutine broyden_reset

!----------------------------------------------------------

subroutine broyden_step(br,x,F,w,newx)
type(broyden_t), intent(inout)      :: br
real(dp), dimension(:), intent(in)  :: x, F
real(dp), dimension(:), intent(out) :: newx
real(dp), intent(in)                :: w

integer :: i, j, n, nn, m, k, l, maxit
real(dp) local_norm, norm, gamma

real(dp), dimension(:), allocatable    :: local_dFF, dFF, local_dFdF
real(dp), dimension(:,:), allocatable  :: a, beta, betabar
real(dp), dimension(:,:), allocatable  :: alpha, alpha_old

n = br%n

!
!  If just reset, or if we went over the maximum number
!  of history slots (i.e., if we did not recycle the oldest
!  positions), restart
!
if (br%it == -1 .or. br%it > br%maxit) then
    if (br%debug .and. ionode) print *,"(Re)starting the Broyden process."
    br%it = 0
    br%dF(1:n,0) = F(1:n)
    newx(1:n) = x(1:n) + br%jinv0*F(1:n)
    br%u(1:n,0) = newx(1:n) - x(1:n)
    RETURN
endif 

m = br%it                ! still previous iteration , ie: 0 for first call
maxit = br%maxit
br%w(m) = w

br%dF(1:n,m) = F(1:n) - br%dF(1:n,m)
local_norm = dot_product(br%dF(1:n,m),br%dF(1:n,m))
!
! Add up norms from all processors and distribute
!
call Globalize_sum(local_norm,norm)
!
norm = sqrt(norm)
!
br%dF(1:n,m) = br%dF(1:n,m)/norm
br%u(1:n,m) = br%u(1:n,m)/norm           ! now contains dx(:,m)
br%u(1:n,m) = br%jinv0 * br%dF(1:n,m) + br%u(1:n,m)
! 

allocate (local_dFF(0:br%it))
allocate(dFF(0:br%it))
allocate(local_dFdF(0:br%it))
allocate(a(0:br%it,0:br%it))
allocate(beta(0:br%it,0:br%it))
allocate(betabar(0:br%it,0:br%it))
allocate(alpha(0:br%it,0:br%it))
allocate(alpha_old(0:br%it,0:br%it))

! Compute new dot products with newly available term
!
do j = 0, m
   local_dFdF(j) = dot_product(br%dF(1:n,j),br%dF(1:n,m))
enddo
! Sum contributions from all processors and distribute

call Globalize_sum(local_dFdF(0:m),br%dFdF(0:m,m))

!
! Symmetry
!
do i = 0, m-1
   br%dFdF(m,i) = br%dFdF(i,m)
enddo
! 
! add wprime**2 in the diagonal to get beta by inversion
!
do i = 0, m
   do j= 0, m
      a(i,j) = br%w(i)*br%w(j) * br%dFdF(i,j)
   enddo
enddo
do i = 0, m
   a(i,i) = a(i,i) + br%wprime**2
enddo
!
! Invert matrix (small...)
!
call dlinds(m+1,a,m+1,beta,m+1)
!
! Partial update with original inverse jacobian
!
newx(1:n) = x(1:n) + br%jinv0*F(1:n) 
!
do i = 0, m
   local_dFF(i) = dot_product(br%dF(1:n,i),F(1:n))
enddo
!
!! Sum contributions from all processors and distribute

call Globalize_sum(local_dFF(0:m),dFF(0:m))

!-------------------------------------------------------
do k = 0, m
   do nn = 0, m
      betabar(k,nn) = 0.0_dp
      do j = 0, m
         betabar(k,nn) = betabar(k,nn)  - &
                  br%w(k) * br%w(j) * beta(k,j) * br%dFdF(nn,j)
      enddo
   enddo
enddo
do k = 0, m
   betabar(k,k) = betabar(k,k) + 1.0_dp
enddo
!-------------------------------------------------------
! Compute alpha (i.e., the coefficients of the Zk in the ui)
! recursively starting from 0
!
l = 0          
do
   do k=0,l
      do i=0,l
         alpha(k,i) = beta(k,i)*br%w(k)*br%w(i)
      enddo
      do i=0,l-1   ! previous alpha only reaches up to m-1...
         do nn=0, l-1
            alpha(k,i) = alpha(k,i) + betabar(k,nn)*alpha_old(nn,i)
         enddo
      enddo
   enddo
   if (l == m) exit
   alpha_old(0:l,0:l) = alpha(0:l,0:l)
   l = l + 1
enddo  
   
if (br%debug) then
   if (ionode) then
      print *, "Alpha:"
      do k = 0, m
         print *, (alpha(k,i), i=0,m)
      enddo
   endif
endif

do l = 0, m
   gamma = 0.0_dp
   do k=0,m
      gamma = gamma + alpha(k,l)*dFF(k)
   enddo
   if (br%debug) then
      if (ionode) print *, "gamma ", l, ": ", gamma
   endif
   newx(1:n) = newx(1:n) - gamma * br%u(1:n,l)
enddo

!
deallocate (local_dFF)
deallocate(dFF)
deallocate(local_dFdF)
deallocate(a)
deallocate(beta)
deallocate(betabar)
deallocate(alpha)
deallocate(alpha_old)

!
! Increase iteration counter and store magnitudes for next one
!
br%it = br%it + 1

if (br%it > br%maxit) then
   
   if (br%cycle_on_maxit) then
         !call message("Cycling the Broyden process...")
         br%dFdF(0:maxit-1,0:maxit-1) = br%dFdF(1:maxit,1:maxit)
         br%w(0:maxit-1) = br%w(1:maxit)
         br%u(1:n,0:maxit-1) = br%u(1:n,1:maxit)
         br%dF(1:n,0:maxit-1) = br%dF(1:n,1:maxit)
         br%it = br%maxit  ! note fixing of marker, to cycle over and over
   else
      RETURN   ! will restart on next entry, as there is no space
               ! to allocate more history records
   endif
endif      

br%dF(1:n,br%it) = F(1:n)
br%u(1:n,br%it) = newx(1:n) - x(1:n)

end subroutine broyden_step

!------------------------------------------------------------
subroutine broyden_destroy(br)
type(broyden_t), intent(inout) :: br

call de_alloc(br%dF)
call de_alloc(br%u)
call de_alloc(br%w)
call de_alloc(br%dFdF)

end subroutine broyden_destroy

!------------------------------------------------------------
subroutine broyden_init(br,debug)
type(broyden_t), intent(inout) :: br
logical, intent(in), optional :: debug

   nullify(br%dF, br%u)
   nullify(br%w, br%dFdF)

   br%debug = .false.
   if (present(debug)) then
      br%debug = debug
   endif
   br%setup = .false.

end subroutine broyden_init

!------------------------------------------------------------
subroutine dlinds(n,a,np1,ainv,np)

!     Inverse of a matrix. LAPACK version

    integer, intent(in)   ::  n,np1,np
    real(dp), intent(in)  ::  a(np1,*)
    real(dp), intent(out) ::  ainv(np,*)


      real(dp), dimension(n) ::  work   ! Automatic arrays
      integer, dimension(n)  ::  ipiv

      integer info, lwork

      ainv(1:n,1:n) = a(1:n,1:n)

!     USE LAPACK

      call dgetrf(n,n,ainv,np,ipiv,info)
      if (info.ne.0) then
         write(6,*) 'Error in DGETRF. INFO:',info
         call die()
      endif

      lwork = n
      call dgetri(n,ainv,np,ipiv,work,lwork,info)
      if (info.ne.0) then
         write(6,*) 'Error in DGETRI. INFO:',info
         call die()
      endif

      end subroutine dlinds


end module m_broyddj
