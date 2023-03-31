! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!!@LICENSE
!
! ==============================================================================
! MODULE interpolation
! ==============================================================================
! Cubic spline interpolation utilities
! ==============================================================================
! Public types in this module:
!   spline_t     : derived type to hold spline info of a function
! Public parameters, variables and arrays in this module:
!   none
! Public procedures provided by this module:
!   generate_spline : generate spline interpolation info of a function
!   evaluate_spline : evaluate a function at given point(s) using splines
! ==============================================================================
! SUBROUTINE generate_spline( dat, x, y, n, dydx1, dydxn, d2ydx2, stat )
!   Generate data required to interpolate a 1D function with cubic splines
! Required input:
!   real(dp) x(n)  : independent variable at each point
!   real(dp) y(n)  : function values at each point
!   integer  n     : number of points
! Optional input:
!   real(dp) dydx1 : dy/dx at x(1) (if not present, assumes d2y/dx2=0 at x(1))
!   real(dp) dydxn : dy/dx at x(n) (if not present, assumes d2y/dx2=0 at x(n))
! Output:
!   type(spline_t) dat : spline info of function y(x)
! Optional output:
!   real(dp) d2ydx2(n) : d2y/dx2 at mesh points
!   integer  stat      : error status
! Behaviour:
!   If dydx1 and/or dydxn are not present, assumes d2y/dx2=0 at x(1) and/or x(n)
!   Returns with stat=-1 if mesh is not monotonic, without any other output
!   Returns with stat=-2 if n<2, without any other output
! Algorithm:
!   Computes d2y/dx2 at each point, to make y(x) and dy/dx continuous.
!   The x values are analyzed to check if the mesh is linear, logarithmic,
!   or general. This is used by evaluate_spline to find the mesh interval at
!   which the evaluation points lie.
! Ref: "Numerical Recipes", W.H. Press et al, Cambridge U.P.
! ==============================================================================
! SUBROUTINE evaluate_spline( dat, x, y )
!   Evaluate function at given point(s) using spline data info
! Input:
!   type(spline_t) dat  : spline info of function y(x), from generate_spline
!   real(dp)       x(:) : point(s) at which function must be evaluated
! Output:
!   real(dp)       y(:) : value of function at given point(s)
! Behaviour:
!   Single- and multiple-point evaluators are overloaded with the same name,
!   so that x and y may be also scalar values
!   Stops with an error message if any point x is out of the interp. interval
! ==============================================================================
! subroutine clean_spline(dat)
!   type(spline_t), intent(inout) :: dat
! Behavior: 
!     Deallocates the array components and resets %n to zero
! end subroutine clean_spline
! ==============================================================================
! SUBROUTINE spline( x, y, n, dydx1, dydxn, d2ydx2 )
!   Included for compatibility with Numerical Recipes interface
! Input:
!   real(dp) x(n)      ! mesh points
!   real(dp) y(n)      ! function value at mesh points
!   integer  n         ! number of mesh points
!   real(dp) dydx1     ! dy/dx at x(1). If >1.e30, assumes d2y/dx2=0 at x(1)
!   real(dp) dydxn     ! dy/dx at x(n). If >1.e30, assumes d2y/dx2=0 at x(n)
! Output:
!   real(dp) d2ydx2(n) ! d2y/dx2 at mesh points
! ==============================================================================
! SUBROUTINE spline( dx, y, n, dydx1, dydxn, d2ydx2 )
!   Included for compatibility with interface used in siesta
! Input:
!   real(dp) dx        ! mesh interval of a uniform mesh starting at x=0
!   real(dp) y(n)      ! function value at mesh points
!   integer  n         ! number of mesh points
!   real(dp) dydx1     ! dy/dx at x(1). If >1.e30, assumes d2y/dx2=0 at x(1)
!   real(dp) dydxn     ! dy/dx at x(n). If >1.e30, assumes d2y/dx2=0 at x(n)
! Output:
!   real(dp) d2ydx2(n) ! d2y/dx2 at mesh points
! ==============================================================================
! SUBROUTINE splint( xi, yi, d2ydx2, n, x, y, dydx )
!   Included for compatibility with Numerical Recipes interface
! Input:
!   real(dp) xi(n)     ! mesh points
!   real(dp) yi(n)     ! function value at mesh points
!   real(dp) d2ydx2(n) ! second derivative at mesh points
!   integer  n         ! number of mesh points
!   real(dp) x         ! point at which function is needed
! Output:
!   real(dp) y         ! function value at point x
!   real(dp) dydx      ! function derivative at point x
! ==============================================================================
! SUBROUTINE splint( dx, yi, d2ydx2, n, x, y, dydx )
!   Included for compatibility with interface used in siesta
! Input:
!   real(dp) dx        ! mesh interval of a uniform mesh starting at x=0
!   real(dp) yi(n)     ! function value at mesh points
!   real(dp) d2ydx2(n) ! second derivative at mesh points
!   integer  n         ! number of mesh points
!   real(dp) x         ! point at which function is needed
! Output:
!   real(dp) y         ! function value at point x
!   real(dp) dydx      ! function derivative at point x
! ==============================================================================
! SUBROUTINE polint(XA,YA,N,X,Y,DYDX) 
!   Lagrange interpolation
! Input:
!   real*8  XA(N) : x values of the function y(x) to interpolate
!   real*8  YA(N) : y values of the function y(x) to interpolate
!   integer N     : Number of data points
!   real*8  X     : x value at which the interpolation is desired
! Output:
!   real*8  Y     : interpolated value of y(x) at X
!   real*8  DYDX  : interpolated derivative dy/dx at X
!                   Notice: this argument has a different meaning
!                   in the Numerical Recipes' polint subroutine
! Ref:
!   W.H.Press et al, Numerical Recipes, Cambridge U.P.
! ==============================================================================
! Written by J.M.Soler and A.Garcia. May.2015
! ==============================================================================
!
! This module calls an external subroutine 'die', with interface
!
!     interface
!      subroutine die(str)
!      character(len=*), intent(in)  :: str
!      end subroutine die
!     end interface

MODULE interpolation

  implicit none

  integer, parameter :: dp = selected_real_kind(10,100)


PUBLIC :: &
  spline_t,        &! derived type to hold spline info of a function
  generate_spline, &! generate spline info
  evaluate_spline, &! evaluate a function at given point(s)
  clean_spline,    &! deallocate the components of a spline_t object
  spline,          &! compatible with Numerical Recipes interface
  splint,          &! compatible with Numerical Recipes interface
  polint            ! Lagrange polynomial interpolation

PRIVATE ! nothing is declared public beyond this point

! Derived type to hold spline info of a function
type spline_t
  private
  character(len=3):: mesh          ! mesh type ('lin'|'log'|'gen')
  real(dp)        :: xmin          ! x(1)-xtol
  real(dp)        :: xmax          ! x(n)+xtol
  real(dp)        :: x1            ! x(1)
  real(dp)        :: dx            ! linear mesh interval
  real(dp)        :: a,b           ! log-mesh parameters x(k)=b*[exp(a*(k-1))-1]
  integer         :: n             ! number of points
  real(dp),allocatable:: x(:)      ! mesh points
  real(dp),allocatable:: y(:)      ! function value at mesh points
  real(dp),allocatable:: d2ydx2(:) ! 2nd derivative at mesh points
end type

! Overloaded spline generators and evaluators
interface generate_spline
  module procedure generate_spline_master    ! new interface
  module procedure generate_spline_x  ! Numerical Recipes interface
  module procedure generate_spline_dx ! older interface
end interface generate_spline
interface evaluate_spline
  module procedure evaluate_spline    ! evaluate function at single point
  module procedure evaluate_spline_n  ! evaluate function at multiple points
  module procedure evaluate_spline_x  ! Numerical Recipes interface
  module procedure evaluate_spline_dx ! siesta interface
end interface evaluate_spline
interface spline
  module procedure generate_spline_x  ! Numerical Recipes interface
  module procedure generate_spline_dx ! siesta interface
end interface spline
interface splint
  module procedure evaluate_spline_x  ! Numerical Recipes interface
  module procedure evaluate_spline_dx ! siesta interface
end interface splint

! Internal module parameters
real(dp),parameter:: dydxMax = 0.99e30_dp  ! max. dy/dx at x(1) and x(n)
real(dp),parameter:: tol = 1.e-6_dp   ! tolerance for 'within range' condition

CONTAINS  

!-------------------------------------------------------------------------------

SUBROUTINE generate_spline_master( dat, x, y, n, dydx1, dydxn, d2ydx2, store, stat )

! Generate data for cubic spline interpolation of a function

implicit none
type(spline_t),   intent(out):: dat    ! spline interpolation info
integer,          intent(in) :: n      ! number of mesh points
real(dp),         intent(in) :: x(n)   ! independent variable at mesh points
real(dp),         intent(in) :: y(n)   ! function value at mesh points
real(dp),optional,intent(in) :: dydx1  ! dy/dx at x(1)
real(dp),optional,intent(in) :: dydxn  ! dy/dx at x(n)
real(dp),optional,intent(out):: d2ydx2(n) ! d2y/dx2 at mesh points
logical, optional, intent(in):: store  ! Store data in the spline object
integer, optional,intent(out):: stat   ! error status:
                                       !  (-1 if nonmonotonic mesh)
                                       !  (-2 if n < 2)

! Internal variables and arrays
character(len=*),parameter:: myName = 'generate_spline_master '
real(dp):: a, b, dx, dx1, dx2, dxn, dxm, dxp, dy1, dyn, dym, dyp, &
           p, s, u(n), v(n), xtol, ypp(n)
integer :: flag, k
character(len=3):: meshType

! Check that mesh is monotonic
if (n<2) then
  flag = -2   ! single-point
elseif (n>2 .and. any( (x(2:n-1)-x(1:n-2)) * (x(3:n)-x(2:n-1)) <= 0.0_dp )) then
  flag = -1   ! non-monotonic mesh
else
  flag = 0
endif
if (present(stat)) stat = flag
if (flag<0) return

! Find if mesh is linear (a=0) or logarithmic
!   x(k)=x(1)+b*[exp(a*(k-1))-1] => (x(k+1)-x(k))/(x(k)-x(k-1))=exp(a)
if (n<3) then
  meshType = 'lin'
  dx = x(2)-x(1)
else
  a = log( (x(3)-x(2)) / (x(2)-x(1)) )
  if (all(abs( (x(3:n)-x(2:n-1))/(x(2:n-1)-x(1:n-2))-exp(a) ) < 1.e-12_dp)) then
    if (abs(a)<1.e-12_dp) then
      meshType = 'lin'
      dx = x(2)-x(1)
    else  ! try x(k) = x(1) + b*[exp(a*(k-1))-1]
      meshType = 'log'
      b = (x(2)-x(1)) / (exp(a)-1._dp)
    endif
  else
    meshType = 'gen'
  endif
endif

! Set boundary conditions at end points
if (present(dydx1)) then
    dx1 = x(2)-x(1)
    dy1 = y(2)-y(1)
    v(1) = -0.5_dp
    u(1) = (3.0_dp/dx1)*(dy1/dx1-dydx1)
else
    v(1) = 0._dp
    u(1) = 0._dp
endif
if (present(dydxn)) then
    dxn = x(n)-x(n-1)
    dyn = y(n)-y(n-1)
    v(n) = 0.5_dp
    u(n) = (3.0_dp/dxn)*(dydxn-dyn/dxn)
else
    v(n) = 0.0_dp
    u(n) = 0.0_dp
endif

! Decomposition loop of the tridiagonal equations
do k = 2,n-1
    dxm = x(k)-x(k-1)
    dym = y(k)-y(k-1)
    dxp = x(k+1)-x(k)
    dyp = y(k+1)-y(k)
    s = dxm / (dxm+dxp)
    p = s*v(k-1) + 2.0_dp
    v(k) = (s-1.0_dp)/p
    u(k) = ( 6.0_dp*(dyp/dxp-dym/dxm)/(dxm+dxp) - s*u(k-1) )/p
enddo

! Backsubstitution loop of the tridiagonal equations
ypp(n) = (u(n)-v(n)*u(n-1)) / (v(n)*v(n-1)+1.0_dp)
do k = n-1,1,-1
  ypp(k) = v(k)*ypp(k+1) + u(k)
enddo

! Store output, if requested
if (present(d2ydx2)) d2ydx2 = ypp

if ( present(store) ) then
  if ( .not. store ) return
end if

! Store spline info in output variable
allocate(dat%x(n), dat%y(n), dat%d2ydx2(n))

dx = (x(n)-x(1))/(n-1)
xtol = abs(dx)*tol
dat%mesh = meshType
dat%xmin = min(x(1),x(n)) - xtol
dat%xmax = max(x(1),x(n)) + xtol
dat%x1 = x(1)
dat%dx = dx
dat%a = a
if (dat%mesh == 'log') dat%b = b
dat%n = n
dat%x(:) = x(:)
dat%y(:) = y(:)
dat%d2ydx2(:) = ypp(:)

end subroutine generate_spline_master

!-------------------------------------------------------------------------------

SUBROUTINE evaluate_spline( dat, x, y, dydx ) ! Evaluate function at one point

implicit none
type(spline_t),   intent(in) :: dat  ! info for spline interpolation of y(x)
real(dp),         intent(in) :: x    ! point at which function is needed
real(dp),         intent(out):: y    ! function value at point x
real(dp),optional,intent(out):: dydx ! function derivative at point x

! Internal variables
integer :: kh, kl
real(dp):: xh, xl

! Find mesh interval of point x
call find_interval( dat, dat%x, x, kl, kh, xl, xh )

! Find interpolation within given interval
call interpolate_interval( xl, xh, dat%y(kl), dat%y(kh), &
                           dat%d2ydx2(kl), dat%d2ydx2(kh), x, y, dydx ) 

END SUBROUTINE evaluate_spline

!-------------------------------------------------------------------------------

SUBROUTINE find_interval( dat, xi, x, kl, kh, xl, xh )

! Find interpolation interval

implicit none
type(spline_t),   intent(in) :: dat   ! spline info
real(dp),         intent(in) :: xi(:) ! mesh points
real(dp),         intent(in) :: x     ! point at which function is needed
integer,          intent(out):: kl    ! lower interval index
integer,          intent(out):: kh    ! higher interval index
real(dp),         intent(out):: xl    ! lower interval mesh point
real(dp),         intent(out):: xh    ! higher interval mesh point

! Internal variables
character(len=*),parameter:: myName = 'evaluate_spline/find_interval '
real(dp)::  a, b, dx, h, x1, xmin, xmax, xtol
integer ::  k, n

! Get some parameters
xmin = dat%xmin
xmax = dat%xmax
n = size(xi)

! Check that point is within interpolation range
if (x<dat%xmin .or. x>dat%xmax) &
  call die(myName//'ERROR: x out of range')

! Find mesh interval of point x
select case(dat%mesh)
case('lin')      ! linear mesh
  x1 = dat%x1
  dx = dat%dx
  kl = 1+floor((x-x1)/dx)
  kl = min(n-1,max(1,kl))
  kh = kl+1
  xl = x1+(kl-1)*dx
  xh = x1+(kh-1)*dx
case ('log')     ! logarithmic mesh
  x1 = dat%x1
  a  = dat%a
  b  = dat%b
  kl = 1+floor(log(1+(x-x1)/b)/a)
  kl = min(n-1,max(1,kl))
  kh = kl+1
  xl = x1+b*(exp((kl-1)*a)-1)
  xh = x1+b*(exp((kh-1)*a)-1)
case('gen')      ! general mesh => use bisection
  kl = 1
  kh = n
  do while(kh-kl>1)
    k = (kh+kl)/2
    if (xi(k)>x) then
      kh = k
    else
      kl = k
    endif
  enddo
  xl = xi(kl)
  xh = xi(kh)
case default
  call die(myName//'ERROR: unknown mesh type')
end select

END SUBROUTINE find_interval

!-------------------------------------------------------------------------------

SUBROUTINE interpolate_interval( xl, xh, yl, yh, d2ydx2l, d2ydx2h, x, y, dydx )

! Evaluate function at point x, within a given interval

implicit none
real(dp),         intent(in) :: xl      ! lower interval point
real(dp),         intent(in) :: xh      ! higher interval point
real(dp),         intent(in) :: yl      ! function value at xl
real(dp),         intent(in) :: yh      ! function value at xh
real(dp),         intent(in) :: d2ydx2l ! d2y/dx2 at at xl
real(dp),         intent(in) :: d2ydx2h ! d2y/dx2 at at xh
real(dp),         intent(in) :: x       ! point at which function is needed
real(dp),         intent(out):: y       ! function value at point x
real(dp),optional,intent(out):: dydx    ! function derivative at point x

! Internal variables
real(dp)::  A, B, dAdx, dBdx, dx

! Find interpolated value
dx = xh-xl
A = (xh-x)/dx
B = (x-xl)/dx
y = A*yl + B*yh + ( (A**3-A)*d2ydx2l + (B**3-B)*d2ydx2h ) * (dx**2)/6.0_dp

! Find interpolated derivative
if (present(dydx)) then
  dAdx = -1._dp/dx
  dBdx = +1._dp/dx
  dydx = dAdx*yl + dBdx*yh + ( dAdx*(3*A**2-1._dp)*d2ydx2l + & 
                               dBdx*(3*B**2-1._dp)*d2ydx2h ) * (dx**2)/6.0_dp
endif

END SUBROUTINE interpolate_interval

!-------------------------------------------------------------------------------

SUBROUTINE evaluate_spline_n( dat, x, y, dydx )

! Evaluate a function at several points

implicit none
type(spline_t),   intent(in) :: dat     ! info for spline interpolation of y(x)
real(dp),         intent(in) :: x(:)    ! point at which function is needed
real(dp),         intent(out):: y(:)    ! function value at point x
real(dp),optional,intent(out):: dydx(:) ! function derivative at point x

integer:: ix, nx

nx = size(x)
do ix = 1,nx
  call evaluate_spline( dat, x(ix), y(ix), dydx(ix) )
enddo

END SUBROUTINE evaluate_spline_n

!-------------------------------------------------------------------------------

SUBROUTINE generate_spline_x( x, y, n, dydx1, dydxn, d2ydx2 )

! Included for compatibility with an older interface

implicit none
integer, intent(in) :: n         ! number of mesh points
real(dp),intent(in) :: x(n)      ! mesh of independent variable
real(dp),intent(in) :: y(n)      ! function value at mesh points
real(dp),intent(in) :: dydx1     ! dy/dx at x(1)
real(dp),intent(in) :: dydxn     ! dy/dx at x(n)
real(dp),intent(out):: d2ydx2(n) ! d2y/dx2 at mesh points

! Internal variables
integer :: stat
type(spline_t):: dat

! Generate spline dat
if (dydx1>dydxMax .and. dydxn>dydxMax) then
  call generate_spline_master( dat, x, y, n, d2ydx2=d2ydx2, store=.false., stat=stat)
elseif (dydxn>dydxMax) then
  call generate_spline_master( dat, x, y, n, dydx1=dydx1, d2ydx2=d2ydx2, store=.false., stat=stat)
elseif (dydx1>dydxMax) then
  call generate_spline_master( dat, x, y, n, dydxn=dydxn, d2ydx2=d2ydx2, store=.false., stat=stat)
else
  call generate_spline_master( dat, x, y, n, dydx1, dydxn, d2ydx2=d2ydx2, store=.false., stat=stat)
endif

! If the status is faulty the resulting d2y/dx2 will be
! forcefully set to 0
! Then the user has to do something differently
if ( stat /= 0 ) d2ydx2(:) = 0._dp

! Deallocate arrays in the spline data-structure
call clean_spline(dat)

END SUBROUTINE generate_spline_x

!-------------------------------------------------------------------------------

SUBROUTINE generate_spline_dx( dx, y, n, dydx1, dydxn, d2ydx2 )

! Included for compatibility with an older interface

implicit none
integer, intent(in) :: n         ! number of mesh points
real(dp),intent(in) :: dx        ! mesh interval
real(dp),intent(in) :: y(n)      ! function value at mesh points
real(dp),intent(in) :: dydx1     ! dy/dx at x(1)
real(dp),intent(in) :: dydxn     ! dy/dx at x(n)
real(dp),intent(out):: d2ydx2(n) ! d2y/dx2 at mesh points

! Internal variables
integer :: ix
real(dp):: x(n)

! Generate spline data
x = (/( (ix-1)*dx, ix=1,n )/)
call generate_spline_x( x, y, n, dydx1, dydxn, d2ydx2 )

END SUBROUTINE generate_spline_dx

!-------------------------------------------------------------------------------

SUBROUTINE evaluate_spline_x( xi, yi, d2ydx2, n, x, y, dydx )

! Included for compatibility with an older interface

implicit none
integer,          intent(in) :: n         ! number of mesh points
real(dp),         intent(in) :: xi(n)     ! mesh of independent variable
real(dp),         intent(in) :: yi(n)     ! function value at mesh points
real(dp),         intent(in) :: d2ydx2(n) ! function value at point x
real(dp),         intent(in) :: x         ! point at which function is needed
real(dp),         intent(out):: y         ! function value at point x
real(dp),optional,intent(out):: dydx      ! function derivative at point x

integer :: kh, kl
real(dp):: xh, xl, xtol
type(spline_t):: dat

! Find mesh interval of point x (warning: no check that x is within range)
dat%mesh = 'gen'
dat%x1 = xi(1)
dat%dx = (xi(n)-xi(1))/(n-1)
xtol = abs(dat%dx)*tol
dat%xmin = min(xi(1),xi(n)) - xtol
dat%xmax = max(xi(1),xi(n)) + xtol
call find_interval( dat, xi, x, kl, kh, xl, xh )

! Find interpolation within given interval
call interpolate_interval( xl, xh, yi(kl), yi(kh), &
                           d2ydx2(kl), d2ydx2(kh), x, y, dydx ) 

END SUBROUTINE evaluate_spline_x

!-------------------------------------------------------------------------------

SUBROUTINE evaluate_spline_dx( dx, yi, d2ydx2, n, x, y, dydx )

! Included for compatibility with an older interface

implicit none
integer,          intent(in) :: n         ! number of mesh points
real(dp),         intent(in) :: dx        ! mesh interval
real(dp),         intent(in) :: yi(n)     ! function value at mesh points
real(dp),         intent(in) :: d2ydx2(n) ! function value at point x
real(dp),         intent(in) :: x         ! point at which function is needed
real(dp),         intent(out):: y         ! function value at point x
real(dp),optional,intent(out):: dydx      ! function derivative at point x

integer :: kh, kl
real(dp):: xh, xl, xmax, xtol

! Check that point is within interpolation range
xmax = (n-1)*dx
xtol = dx*tol
if (x<-xtol .or. x>xmax+xtol) &
  call die('evaluate_spline ERROR: x out of range')

! Find mesh interval of point x (warning: no check that x is within range)
kl = 1+floor(x/dx)
kl = min(n-1,max(1,kl))
kh = kl+1
xl = (kl-1)*dx
xh = (kh-1)*dx

! Find interpolation within given interval
call interpolate_interval( xl, xh, yi(kl), yi(kh), &
                           d2ydx2(kl), d2ydx2(kh), x, y, dydx ) 

END SUBROUTINE evaluate_spline_dx

!-------------------------------------------------------------------------------

SUBROUTINE polint( xi, yi, n, x, y, dydx )  ! Lagrange interpolation

  implicit none
  integer, intent(in) :: n     ! number of mesh points
  real(dp),intent(in) :: xi(*) ! mesh points
  real(dp),intent(in) :: yi(*) ! function values at mesh points
  real(dp),intent(in) :: x     ! point at which function is required
  real(dp),intent(out):: y     ! function value at point x
  real(dp),intent(out):: dydx  ! function derivative at point x

  integer :: i0, m, nm
  real(dp):: c(n), d(n), e(n), dcdx(n), dddx(n), dedx(n), dx1(n), dx2(n), dxi(n)

! Neville's algorithm, as described in NR
  c = yi(1:n)
  d = yi(1:n)
  dcdx = 0._dp
  dddx = 0._dp
  dydx = 0._dp
  i0 = n/2
  y = yi(i0)
  do m = 1,n-1               ! loop on increasing polynomial order
    nm = n-m
    dxi(1:nm) = xi(1:nm)-xi(m+1:n)
    if (any(dxi(1:nm)==0.0_dp)) &
      call die('polint: ERROR: two mesh points are equal')
    dx1(1:nm) = xi(1:nm)-x
    dx2(1:nm) = xi(m+1:n)-x
    e(1:nm) = (c(2:nm+1)-d(1:nm))/dxi(1:nm)
    c(1:nm) = dx1(1:nm)*e(1:nm)
    d(1:nm) = dx2(1:nm)*e(1:nm)
    dedx(1:nm) = (dcdx(2:nm+1) - dddx(1:nm)) / dxi(1:nm)
    dcdx(1:nm) = - e(1:nm) + dx1(1:nm)*dedx(1:nm)
    dddx(1:nm) = - e(1:nm) + dx2(1:nm)*dedx(1:nm)
    i0 = i0-1
    if (2*i0 < nm) then
      i0 = i0+1
      y = y+c(i0)
      dydx = dydx+dcdx(i0)
    else
      y = y+d(i0)
      dydx = dydx+dddx(i0)
    endif
  end do

END SUBROUTINE polint

subroutine clean_spline(dat)
type(spline_t), intent(inout) :: dat

integer :: n

n = dat%n
if (allocated(dat%x)) then
   deallocate(dat%x, dat%y, dat%d2ydx2)
endif
dat%n = 0
end subroutine clean_spline

END MODULE interpolation

