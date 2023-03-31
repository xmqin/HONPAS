! ---
!  This file is distributed under the terms of the
!  Modified BSD License: see the BSD_LICENSE file
! ---
!
! *******************************************************************
! module mesh1D
!-----------------------------------------------------------------
! Contains utility routines to manipulate one-dimensional functions defined 
! in a mesh, and to solve some differential equations by Numerov's method
!-----------------------------------------------------------------
! Used module procedures:
!   interpolation, only: polint  ! Polynomial (Lagrange) interpolation
!   interpolation, only: spline  ! Sets spline interpolation
!   interpolation, only: splint  ! Finds spline interpolation
! Used module types and variables:
!   precision, only: dp        ! Real double precision type
!-----------------------------------------------------------------
! Public module procedures:
!   derivative,        &! Returns function derivatives the same mesh points
!   get_mesh,          &! Returns a previously-set 1D mesh
!   get_n,             &! Returns the number of mesh points
!   integral,          &! Returns the integral of a function defined in a mesh
!   interpolation_local,     &! Returns interpolated values at arbitrary points
!   locate,            &! Given x0, it returns real value i0 such that x(i0)=x0
!   numerov,           &! Solves d2y/dx2 = f(x,y) = f0(x) + f1(x)*y
!   set_interpolation, &! Sets interpolation method (lagrange|spline)
!   set_mesh            ! Sets a uniform, logarithmic, or arbitrary 1D mesh
! Public module types and variables: none
!-----------------------------------------------------------------
! get_n: returns the number of mesh points required to meet some conditions
! Usage:
!   n = get_n( xmin, xmax, dxmin, dxmax )
!     integer  n     ! Number of logarithmic-mesh points
!     real(dp) xmin  ! First mesh point
!     real(dp) xmax  ! Last mesh point
!     real(dp) dxmin ! First mesh interval: x(2)-x(1)
!     real(dp) dxmax ! Next to last mesh interval: x(n+1)-x(n)
! Notes:
! - All arguments are input and required
! - No mesh needs to have been previously set. Instead, this
!   function is to help to find the argument n needed to set it.
!-----------------------------------------------------------------
! set_mesh: sets a uniform, logarithmic, or arbitrary 1D mesh
! Usage:
!   call set_mesh( n, x, xmin, xmax, a, dxndx1 )
!     integer  n     : Number of mesh points
!     real(dp) x(n)  : Mesh point values
!     real(dp) xmin  : Value of first mesh point
!     real(dp) xmax  : Value of last mesh point
!     real(dp) a     : Fractional increment of succesive mesh 
!                      intervals of a logarithmic mesh:
!                        x(i) = b * (exp(a*(i-1)) - 1)
!     real(dp) dxndx1: Ratio between last and first intervals
!                      in a logarithmic mesh
! Notes:
! - All arguments are input
! - All the arguments, except n, are optional.
! - If x is present, no other optional arguments may be present
! - If x is not present, xmax must be present, and only one of
!   a or dxndx1 may be present
! - If xmin is not present, xmin=0 is assumed
! - If only xmax (and optionally xmin) is present, a regular mesh 
!   is generated
! - Stops with an error message if x is present but x(i) is not monotonic
! - If x is present and it is a uniform or a logarithmic mesh, it will be
!   identified as such, and analytic formulas will be used accordingly.
!   However, this takes extra time, so it is preferable to use the
!   other arguments to set a uniform or a logarithmic mesh.
! Examples:
! - To initialize a previously generated mesh:
!      call set_mesh( n, x=x_mesh(1:n) )
! - To generate a regular mesh:
!      call set_mesh( n, xmax=x_max )
! - To generate a logarithmic mesh:
!      call set_mesh( n, dxndx1=delta_x_max/delta_x_min )
! Warning:
! - One should not fix parameter a while increasing nx to 
!   convergence, because this may easily lead to an extremely  
!   small dx1 and to roundoff errors. Generally, it is safer
!   and easier to fix dxndx1, specially for this purpose.
!-----------------------------------------------------------------
! get_mesh: returns a previously-set 1D mesh
! Usage:
!   call get_mesh( n, nx, x, dxdi, d2xdi2, d3xdi3, d4xdi4 )
!     integer  n         : Size of array arguments (input)
!     integer  nx        : Number of mesh points (ouput)
!     real(dp) x(n)      : Mesh point values
!     real(dp) dxdi(n)   : dx/di
!     real(dp) d2xdi2(n) : d2x/di2
!     real(dp) d3xdi3(n) : d3x/di3
!     real(dp) d4xdi4(n) : d4x/di4
! Notes:
! - All arguments are output, except n
! - All the arguments, except n and nx, are optional.
! - If any optional array is present, nx is limited by its size, 
!   i.e. nx.le.n
! - If no optional array is present,  present, nx is the number 
!   of mesh points stored
! Examples:
! - To get the first derivatives:
!      call get_mesh( n, nx, dxdi=dxdi(1:n) )
!----------------------------------------------------------------
! set_interpolation: sets interpolation method and end-point derivatives
! Usage:
!   call set_interpolation( method, ypleft, ypright )
!     character(len=*) method  : Interpolation method (lagrange|spline)
!     real(dp) ypleft, ypright : dy/dx values at end points
! Notes:
! - Stops with an error message if method /= ('lagrange'|'spline')
! - Arguments ypleft and ypright are optional and they are used only 
!   when method='spline'. If absent, natural spline conditions (zero  
!   second derivative at corresponding end points) are used.
! Examples:
! - To set natural spline at left and zero derivative at right point:
!      call set_interpolation( 'spline', ypright=0._dp )
!----------------------------------------------------------------
! locate: given x0, returns real value i0 such that x(i0)=x0
! Usage:
!   i0 = locate( x0 )
!     real(dp) x0 : point to locate, within range x(1):x(n)
!     real(dp) i0 : real value such taht interpolated x(i0)=x0
! Notes:
! - Stops with an error message if mesh is not defined
! - For numeriacal meshes stops with an error message if x0 is out of 
!   range x(1):x(n). For uniform and logarithmic meshes it returns an
!   extrapolated value.
!----------------------------------------------------------------
! interpolation: returns the interpolation of a function at one 
! or several points
! Usage:
!   ynew = interpolation( nnew, xnew, n, y, x, dx )
!     real(dp) interpolation(nnew) : Interpolated function values
!     integer  nnew          : Number of interpolation points
!     real(dp) xnew(nnew)    : Interpolation points
!     integer  n             : Number of mesh points
!     real(dp) y(n)          : Values of function to interpolate
!     real(dp) x(n)          : Mesh values (optional)
!     real(dp) dx            : Mesh interval (optional)
! Notes:
! - All arguments are input but the function is not pure because,
!   if x or dx are present, they change the default mesh.
! - The two optional arguments x and dx are incompatible.
!   If none of them is present, the last defined mesh is used.
! Examples:
! - To change the interpoltion mesh of y(x):
!      y(1:n) = interpolation( n, xnew(1:n), n, y(1:n), xold(1:n) )
!----------------------------------------------------------------
! derivative: returns the derivative of a function at the same
! mesh points
! Usage:
!   dydx = derivative( n, y, x, dx, order )
!     real(dp) derivative(n) : Return function derivative
!     integer  n             : Number of points
!     real(dp) y(n)          : Values of function to derivate
!     real(dp) x(n)          : Mesh values (optional)
!     real(dp) dx            : Mesh interval (optional)
!     integer  order         : order of the derivative (optional)
! Notes:
! - All arguments are input but the function is not pure because,
!   if x or dx are present, they change the default mesh.
! - The two optional arguments x and dx are incompatible. They
!   are used as in routine numerov. If none is present, the last
!   defined mesh is used.
! - If order is not present, order=1 is assumed.
! - If n is smaller than the number of mesh points stored, the
!   first n of them are used.
! Examples:
! - To find the Laplacian of y(x), using a previously defined mesh:
!      d2ydx2(1:n) = derivative( n, y(1:n), order=2 )
!----------------------------------------------------------------
! integral: returns the integral of a function defined in a mesh
! Usage:
!   the_integral = integral( n, y, x, dx )
!     real(dp) integral : Return value
!     integer  n        : Number of points
!     real(dp) y(n)     : Values of function to integrate
!     real(dp) x(n)     : Mesh values (optional)
!     real(dp) dx       : Mesh interval (optional)
! Notes:
! - All arguments are input but the function is not pure because,
!   if x or dx are present, they change the default mesh.
! - The two optional arguments x and dx are incompatible. They
!   are used as in routine numerov. If none is present, the last
!   defined mesh is used.
! - If n is smaller than the number of mesh points stored, the
!   first n of them are used.
! Examples:
! - To find the norm of psi, using a previously defined mesh:
!      norm = integral( n, y=psi(1:n)**2 )
!-----------------------------------------------------------------
! numerov: solves an ordinary differential equation of the form
!    d2y/dx2 = f(x,y) = f0(x) + f1(x)*y
! by the Numerov method:
!    (y(x+dx)-2*y(x)+y(x-dx))/dx2 = (f(x+dx)+10*f(x)+f(x-dx))/12
! from which y(x+dx) can be solved from y(x) and y(x-dx).
! Typical cases are f0=-4*pi*rho(x), f1=0 (Poisson eq.) 
!               and f0=0, f1(x)=2*(V(x)-E) (Schroedinger equation)
! Notice that f must not depend on y' (first derivative)
! Usage:
!   call numerov( n, y, yp, ypp, f0, f1, x, dx, norm )
!     integer  n     : Number of mesh points
!     real(dp) y(n)  : Solution
!     real(dp) yp(n) : First derivative y' = dy/dx
!     real(dp) ypp(n): Second derivative y'' = d2y/dx2 = f(x,y)
!     real(dp) f0(n) : Component of f(x,y) independent of y
!     real(dp) f1(n) : Componnet of f(x,y) proportional to y
!     real(dp) x(n)  : Mesh points
!     real(dp) dx    : Mesh interval (only x OR dx allowed)
!     real(dp) norm  : Norm for solution of homogeneous eqs.
! Notes:
! - All arguments, except y and ypp, are input.
! - All the arguments, except n, are optional.
! - If y is not present, only the mesh is initialized. If the
!   differential equation is solved repeatedly with the same mesh,
!   it is recommended to initialize it only once, making further
!   calls without x or dx present. In that case, the mesh used is
!   defined by the last call to any routine (set_mesh, numerov,
!   derivative, or integral) with a mesh-setting argument
! - If n is smaller than the number of mesh points stored, the
!   first n of them are used.
! - Only one of x OR dx are allowed to define the mesh
! - If f0 or f1 are not present, they are assumed to be zero.
! - The first two values of y must be given on input, to specify 
!   the boundary condition. If f0=0 (homogeneous equation), at
!   least one of them must be nonzero.
! - For f0 not present (homogeneous equation) the normalization
!   of the solution (integral of y**2) is fixed by argument norm.
! - If norm is not present, the norm of y is fixed by the input 
!   values of y(1) and y(2)
! - Output array yp is useful to normalize unbound states
! - Output array ypp is useful for a spline interpolation of y(x)
! Examples:
! - To initialize a variable radial mesh for inwards integration
!     call numerov( n+1, x=r_mesh(n:0:-1) )
! - To solve Poisson's equation with the previously defined mesh:
!     call numerov( n+1, y=rV(n:0:-1), f0=-4*pi*rrho(n:0:-1) )
!   with rV=r*V, rrho=r*rho. rV(n) and rV(n-1) may be initialized
!   to Q (integral of rho) to set the zero of potential at infty.
! - To integrate Schroedinger's equation with a regular mesh
!   (psi(1) and psi(2) required on input):
!     call numerov( n, y=psi(1:n), f1=V(1:n)-E, &
!                   dx=delta_x, norm=1.d0 )
!----------------------------------------------------------------
! Units: 
! Units of x and y in all routines are arbitrary, but they must be
! consistent for all arguments within the same call. Units of x 
! must be also consistent in different calls with the same mesh.
!----------------------------------------------------------------
! Algorithms:
! Interpolation of y(n) is performed in the uniform mesh n=1,2,...,
! independently of the x(n) mesh. To this purpose, function
! locate is used to find n0 such that x(n0)=x0, where x0 is the
! interpolation point. If x(n) is a uniform or logarithmic mesh,
! n0 is obtained analytically. Otherwise, bisection is used to
! locate n0, with a 6-point Lagrange interpolation for x(n).
! Interpolation in the uniform mesh n, rather than in the nonuniform
! mesh xn, is preferred because y(n) is usually better behaved than
! y(x), as well as for consistency with the algorithms used to find
! derivatives and integrals.
!
! For uniform and logarithmic meshes, the derivatives of x(n) are 
! calculated analytically. Otherwise they are approximated by a 
! 5-point formula:
! from   x(n+m) = x(n) + x'(n)*m + x''(n)*m^2/2 +
!                 x'''(n)*m^3/6 + x''''(n)*m^4/24
! we find
!   x'(n)   = ( x(n-2) - 8*x(n-1)          + 8*x(n+1) -x(n+2) )/12
!   x''(n)  = (-x(n-2) +16*x(n-1) -30*x(n) +16*x(n+1) -x(n+2) )/12
!   x'''(n) = (-x(n-2) + 2*x(n-1) -          2*x(n+1) +x(n+2) )/2
!   x''''(n)= ( x(n-2) - 4*x(n-1) + 6*x(n) - 4*x(n+1) +x(n+2  )
!
! Two algorithms are implemented for function derivatives and 
! integrals. The first one uses Lagrange interpolation to define 
! the function between mesh points. The second uses cubic splines. 
! Which method is actually used is determined by the interface 
! settings at the beginning of the module. In the Lagrange version, 
! dy/dn is obtained by the same formula as x'(n) above and we then 
! find dy/dx = (dy/dn)/(dx/dn). 
! The integral is approximated by:
!  3-point formula for first and last intervals:
!    integral_from(1)_to(2) y*dn = ( 5*y(1) + 8*y(2) - y(3) )/12
!  4-point formula for 'interior' intervals: 
!    integral_from(n)_to(n+1) f*dn = 
!       ( -y(n-1) + 13*y(n) + 13*y(n+1) - y(n+2) )/24
! In the spline version, d2y/dn2 is first obtained by imposing the
! matching of y, dy/dn, and d2y/dn2 at each mesh point. Derivatives
! and integrals are then obtained staightforwardly from the cubic
! interpolation polynomials at each point and interval:
!   dy(n)/dn = y(n) - y(n-1) + (2*d2y(n)/dn2 + d2y(n-1)/dn2)/6
!            = y(n+1) - y(n) - (d2y(n+1)/dn2 + 2*d2y(n)/dn2)/6
!            = (y(n+1) - y(n-1))/2 - (d2y(n+1)/dn2 - d2y(n-1)/dn2)/12
!   integral_from(x(n))_to(x(n+1)) y(x)*dx = 
!   integral_from(n)_to(n+1) f(n)*dn = 
!      (f(n)+f(n+1))/2 - (d2f(n)/dn2 + d2f(n+1)/dn2)/24
! with f(n)=y(x(n))*dx(n)/dn. In both versions, higher order 
! derivatives of y(x) are obtained by repeated derivation.
!
! The origin of the 1,10,1 coefs of Numerov's formula 
!    (y(x+dx)-2*y(x)+y(x-dx))/dx2 = (f(x+dx)+10*f(x)+f(x-dx))/12
! is as follows:
!   y(x+dx) = y(x) + y'(x)*dx + y''(x)*dx2/2 + 
!             y'''(x)*dx3/6 + y''''(x)*dx4/24 + O(dx5)
! Then: 
!   (y(x+dx)-2*y(x)+y(x-dx))/dx2 = y''(x) + y''''(x)*dx2/12 + O(dx4)
!                                = f(x) + f''(x)*dx2/12 + O(dx4)
! We now approximate
!   f''(x) = (f(x+dx)-2*f(x)+f(x-dx))/dx2 + O(dx2)
! to find the Numerov formula plus O(dx4) corrections.
! For variable meshes x(n), n=1,2,...,N, we define a new function
!   z(n) = y(x(n)) / (x'(n))^(1/2), where x'=dx/dn. Then, after
! some algebra, it can be shown that
!   d2z/dn2 = g(n,z) = g0(n) + g1(n)*z, where
!   g0(n) = f0(x(n))*(x')^(3/2)
!   g1(n) = f1(x(n))*(x')^2 + [3*(x'')^2/x'-2*x''']/(2*x')^2
! If x decreases with n, we can define a new variable -x and 
!   x' --> -x', what amounts to using |x'|^(1/2) and |x'|^(3/2)
! Within the numerov routine, the derivatives yp and ypp are found
! by a special algorithm, using that y''=f and
!   y(x+dx) - y(x-dx) = 2*y'(x) + y'''(x)/3 + O(dx5) = 
!                       2*y'(x) + f'(x)/3 + O(dx5) =
!                       2*y'(x) + (f(x+dx) - f(x-dx))/(6*dx) + O()
! from which y'(x) is solved. The first and last points are 
! special: taking x = xmax-dx and ignoring O(dx5) errors:
!   y'(x+dx) = y'(x) + y''(x)*dx + y'''(x)*dx2/2 + y''''(x)*dx3/6
!            = y'(x) + f(x)*dx + f'(x)*dx2/2 + f''(x)*dx3/6
!            = y'(x) + f(x)*dx + (f(x+dx)-f(x-dx))*dx/4 +
!              (f(x+dx)-2*f(x)+f(x-dx))*dx/6
!            = y'(x) + (5*f(x+dx) + 8*f(x) - f(x-dx)) * dx/12
! similarly, taking x = xmin+dx:
!   y'(x-dx) = y'(x) - (5*f(x-dx) + 8*f(x) - f(x+dx)) * dx/12
! For a variable mesh, these formulas are used to find z'(n). Then
! y'(x) is found from z' and x', using chain's rule as before
!
! Ref: J.M.Soler notes of Jan.12,1998, Nov.28,2001, Dec.4,2001,
!      Jan.16,2002, Jul.29,2002, Jun.21,2007, and Jul.5,2007.
!-----------------------------------------------------------------
! Implementation details:
! - The routines use straightforward formulas and are optimized
!   for simplicity, clarity, and ease of use, NOT for efficiency.
!-----------------------------------------------------------------
! Written by J.M.Soler. Nov.2001, Jul.2002, Jul.2007, and Jul.2008
!************************************************************************

module mesh1D

! Used module procedures
use interpolation, only: polint  ! Polynomial (Lagrange) interpolation
use interpolation, only: spline  ! Sets spline interpolation
use interpolation, only: splint  ! Finds spline interpolation

! Used module types and variables
use precision, only: dp      ! Real double precision type

implicit none

! Public module procedures
PUBLIC :: &
  derivative,        &! Returns function derivatives the same mesh points
  get_mesh,          &! Returns a previously-set 1D mesh
  get_n,             &! Returns the number of mesh points
  integral,          &! Returns the integral of a function defined in a mesh
  interpolation_local,     &! Returns interpolated values at arbitrary points
  locate,            &! Given x0, it returns real value i0 such that x(i0)=x0
  numerov,           &! Solves d2y/dx2 = f(x,y) = f0(x) + f1(x)*y
  set_interpolation, &! Sets interpolation method (lagrange|spline)
  set_mesh            ! Sets a uniform, logarithmic, or arbitrary 1D mesh

! Public types and variables: none

PRIVATE   ! Nothing is declared public beyond this point

!  integer, parameter :: dp = kind(1.d0)

! Internal parameters
  real(dp),parameter :: amax = 1.e-6_dp  ! Max. exponent coef. in log. mesh
  real(dp),parameter :: xtol = 1.e-12_dp ! Tol. for mesh type identification
  real(dp),parameter :: itol = 1.e-12_dp ! Tolerance for interpolation of x(i)

! Saved internal variables and arrays
  character(len=11),save:: mesh_type='unknown'
  character(len=10),save:: interpolation_method='spline'
  logical,  save:: defined_mesh=.false.  ! Has a mesh been defined?
  real(dp), save:: aa, b ! Log-mesh coefs: x(i)=x(1)+b*(exp(aa*(i-1))-1)
  real(dp), save:: d     ! Uniform mesh interval: x(i)=x(1)+d*(i-1)
  real(dp), save:: yp1=huge(yp1) ! dy/dx at x(1) for spline interp.
  real(dp), save:: ypn=huge(ypn) ! dy/dx at x(n) for spline interp.
  real(dp), save, dimension(:), allocatable :: &
    ivec,  &! Auxiliary vector to call polint: ivec(i)=i
    sqrxp, &! Sqrt(dx/di) at mesh points x(i). Used by numerov.
    s0,    &! (dx/di)**(3/2) at mesh points. Used by numerov.
    s1,    &! (dx/di)**2 at mesh points. Used by numerov.
    s2,    &! (3*xp2**2 - 2*xp1*xp3) / (4*xp1**2). Used by numerov.
    xi,    &! Mesh points x(i)
    xp1,   &! dx/di at mesh points
    xp2,   &! d2x/di2 at mesh points
    xp3,   &! d3x/di3 at mesh points
    xp4     ! d4x/di4 at mesh points

CONTAINS

!----------------------------------------------------------------

function get_n( xmin, xmax, dxmin, dxmax )

  implicit none
  integer              :: get_n
  real(dp), intent(in) :: xmin
  real(dp), intent(in) :: xmax
  real(dp), intent(in) :: dxmin
  real(dp), intent(in) :: dxmax

  real(dp) :: a, b

! Check signs of arguments
  if (dxmax*dxmin<=0._dp .or. (xmax-xmin)*dxmax<=0._dp) &
    stop 'get_n: ERROR: Bad arguments'

! Find required number of points
  if (dxmax==dxmin) then  ! Linear mesh
    get_n = nint( (xmax-xmin) / dxmax )
  else
    b = (xmax-xmin) * dxmin / (dxmax-dxmin)
    a = log( 1 + dxmin/b )
    get_n = 1 + nint( log(1+(xmax-xmin)/b) / a )
  endif

end function get_n

!----------------------------------------------------------------

subroutine set_mesh( n, x, xmin, xmax, a, dxndx1 )

  implicit none
  integer,            intent(in) :: n
  real(dp), optional, intent(in) :: x(n)
  real(dp), optional, intent(in) :: xmin
  real(dp), optional, intent(in) :: xmax
  real(dp), optional, intent(in) :: a
  real(dp), optional, intent(in) :: dxndx1

  integer  :: i
  real(dp) :: e, dx, dd, x1

! Check consistency of optional arguments
  if (present(x) .and. (present(xmin) .or. present(xmax) .or. &
                        present(a) .or. present(dxndx1))) &
    stop 'set_mesh: ERROR: x not compatible with other optional arguments'
  if (present(a) .and. present(dxndx1)) &
    stop 'set_mesh: ERROR: arguments a and dxndx1 are not compatible'

! Check the number of points
  if (n<3) stop 'set_mesh: ERROR: at least three points required'

! Check that mesh is monotonic
  if (present(x)) then
    if (any((x(2:n)-x(1:n-1))*(x(n)-x(1)) <= 0._dp)) &
      stop 'set_mesh: ERROR: mesh not monotonic'
  end if

! Allocate internal tables
  if (defined_mesh) deallocate( ivec, sqrxp, s0, s1, s2, &
                                xi, xp1, xp2, xp3, xp4 )
  allocate( ivec(n), sqrxp(n), s0(n), s1(n), s2(n),  &
            xi(n), xp1(n), xp2(n), xp3(n), xp4(n) )

! Set auxiliary array ivec
  forall(i=1:n) ivec(i) = i

! Set first point
  if (present(x)) then
    x1 = x(1)
  else if (present(xmin)) then
    x1 = xmin
  else
    x1 = 0
  endif

! Set mesh type
  if (present(x)) then

    ! Find temptative values of a and b (log. mesh), or d (uniform mesh)
    aa = log( (x(n)-x(n-1)) / (x(2)-x(1))) / (n-2)
    if (aa < 1.e-6_dp) then
      mesh_type = 'uniform'
      d = (x(n)-x(1)) / (n-1)
    else
      mesh_type = 'logarithmic'
      b = (x(n)-x(1)) / (exp(aa*(n-1)) - 1)
    end if

    ! Check if mesh is really uniform or logarithmic
    if (mesh_type=='uniform') then
      do i = 2,n
        dx = x(i) - x(i-1)
        if (abs(dx/d-1) > xtol) then
          mesh_type = 'numerical'
          exit
        end if
      end do
    else if (mesh_type=='logarithmic') then
      do i = 2,n
        dx = x(i) - x(i-1)
        dd = b*( exp(aa*(i-1)) - exp(aa*(i-2)) )
        if (abs(dx/dd-1) > xtol) then
          mesh_type = 'numerical'
          exit
        end if
      end do
    end if ! (mesh_type=='uniform')

  elseif (present(xmax)) then

    if (present(a)) then
      aa = a
    elseif (present(dxndx1)) then
      aa = log(dxndx1) / (n-1)
    else
      aa = 0
    endif
    if (aa < amax) then
      mesh_type = 'uniform'
      d = (xmax-x1) / (n-1)
    else ! (aa > amax)
      mesh_type = 'logarithmic'
      b = (xmax-x1) / (exp(aa*(n-1)) - 1)
    end if

  else ! (.not.present(x) .and. .not.present(xmax))
    print*, 'set_mesh: undefined mesh'
    stop
  end if ! (present(x))

! DEBUG
!  if (mesh_type=='uniform') then
!    print*, 'set_mesh: mesh_type = ', mesh_type, '   d =', d
!  else if (mesh_type=='logarithmic') then
!    print*, 'set_mesh: mesh_type = ', mesh_type, '   a, b =', aa, b
!  else
!    print*, 'set_mesh: mesh_type = ', mesh_type
!  end if
! END DEBUG

! Set the mesh points and the derivatives of x(i)
  if (mesh_type=='numerical') then
    xi = x

!   Centered 5-point derivatives for all but first/last two points
    do i = 3,n-2
      xp1(i)=( xi(i-2)- 8*xi(i-1)         + 8*xi(i+1)-xi(i+2))/12
      xp2(i)=(-xi(i-2)+16*xi(i-1)-30*xi(i)+16*xi(i+1)-xi(i+2))/12
      xp3(i)=(-xi(i-2)+ 2*xi(i-1)         - 2*xi(i+1)+xi(i+2))/2
      xp4(i)=  xi(i-2)- 4*xi(i-1)+ 6*xi(i)- 4*xi(i+1)+xi(i+2)
    enddo

!   Use first/last 5 points for derivatives of two extreme points
    xp1(1) = xp1(3) - 2*xp2(3) + 4*xp3(3)/2 - 8*xp4(3)/6
    xp1(2) = xp1(3) - 1*xp2(3) + 1*xp3(3)/2 - 1*xp4(3)/6
    xp2(1) = xp2(3) - 2*xp3(3) + 4*xp4(3)/2
    xp2(2) = xp2(3) - 1*xp3(3) + 1*xp4(3)/2
    xp3(1) = xp3(3) - 2*xp4(3)
    xp3(2) = xp3(3) - 1*xp4(3)
    xp4(1) = xp4(3)
    xp4(2) = xp4(3)
    xp1(n)   = xp1(n-2) + 2*xp2(n-2) + 4*xp3(n-2)/2 + 8*xp4(n-2)/6
    xp1(n-1) = xp1(n-2) + 1*xp2(n-2) + 1*xp3(n-2)/2 + 1*xp4(n-2)/6
    xp2(n)   = xp2(n-2) + 2*xp3(n-2) + 4*xp4(n-2)/2
    xp2(n-1) = xp2(n-2) + 1*xp3(n-2) + 1*xp4(n-2)/2
    xp3(n)   = xp3(n-2) + 2*xp4(n-2)
    xp3(n-1) = xp3(n-2) + 1*xp4(n-2)
    xp4(n)   = xp4(n-2)
    xp4(n-1) = xp4(n-2)

  elseif (mesh_type=='logarithmic') then
    do i = 1,n
      e = b * exp(aa*(i-1))
      xi(i) = x1 + e - b
      xp1(i) = aa * e
      xp2(i) = aa**2 * e
      xp3(i) = aa**3 * e
      xp4(i) = aa**4 * e
    enddo
  elseif (mesh_type=='uniform') then
    do i = 1,n
      xi(i) = x1 + (i-1) * d
      xp1(i) = d
      xp2(i) = 0
      xp3(i) = 0
      xp4(i) = 0
    enddo
  else ! ( mesh_type /= ('numerical'|'logarithmic'|'uniform') )
    stop 'set_mesh: ERROR: unknown mesh_type'
  endif ! (mesh_type=='numerical')

! Make sure that first and last points are exactly right
  if (.not.present(x)) xi(1) = x1
  if (present(xmax)) xi(n) = xmax

! Find auxiliary functions associated to the mesh
  sqrxp = abs(xp1)**0.5_dp
  s0 = abs(xp1)**1.5_dp
  s1 = xp1**2
  s2 = (3._dp*xp2**2 - 2._dp*xp1*xp3) / 4._dp / xp1**2

  defined_mesh = .true.

end subroutine set_mesh

!----------------------------------------------------------------

subroutine set_interpolation( method, ypleft, ypright )

  implicit none
  character(len=*), intent(in):: method
  real(dp),optional,intent(in):: ypleft, ypright

  if (method=='spline' .or. method=='spline' .or. &
      method=='SPLINE') then
    interpolation_method = 'spline'
  else if (method=='lagrange' .or. method=='lagrange' .or. &
           method=='LAGRANGE') then
    interpolation_method = 'lagrange'
  else
    stop 'set_interpolation: ERROR: unknown method'
  end if

  if (present(ypleft)) then 
    yp1 = ypleft
  else
    yp1 = huge(1._dp)
  end if

  if (present(ypright)) then 
    ypn = ypright
  else
    ypn = huge(1._dp)
  end if

end subroutine set_interpolation

!----------------------------------------------------------------

subroutine get_mesh( n, nx, x, dxdi, d2xdi2, d3xdi3, d4xdi4 )

  implicit none
  integer, intent(in)  :: n   ! size of argument arrays
  integer, intent(out) :: nx  ! number of points
  real(dp), optional, dimension(:), intent(out) :: &
    x, dxdi, d2xdi2, d3xdi3, d4xdi4

  if (.not.defined_mesh) stop 'get_mesh: ERROR: mesh not defined'
  nx = size(xi)
  if (present(x) .or. present(dxdi) .or. present(d2xdi2) .or. &
      present(d3xdi3) .or. present(d4xdi4)) then
    nx = max( nx, n )
  end if

  if (present(x))      x(1:nx)      = xi(1:nx)
  if (present(dxdi))   dxdi(1:nx)   = xp1(1:nx)
  if (present(d2xdi2)) d2xdi2(1:nx) = xp2(1:nx)
  if (present(d3xdi3)) d3xdi3(1:nx) = xp3(1:nx)
  if (present(d4xdi4)) d4xdi4(1:nx) = xp4(1:nx)

end subroutine get_mesh

!----------------------------------------------------------------

function locate( x0 )

  implicit none
  real(dp)             :: locate  ! Value i0 such that x(i0)=x0
  real(dp), intent(in) :: x0      ! x point

  integer :: i, i1, i2, iter, n
  real(dp):: dx, il, incr, ir, ic, xc
  logical :: found

  if (.not.defined_mesh) &
    stop 'mesh1D locate: ERROR: mesh not defined'

! Find number of mesh points
  n = size(xi)

! Trap easy cases
  if (mesh_type=='uniform') then
    locate = 1 + (x0-xi(1)) / d
    return
  else if (mesh_type=='logarithmic') then
    if (x0 <= xi(1)-b) stop 'mesh1D locate: ERROR: x0 out of range'
    locate = 1 + log(1+(x0-xi(1))/b) / aa
    return
  else if (mesh_type/='numerical') then
    stop 'mesh1D locate: ERROR: unknown mesh type'
  end if

! Find if x(i) is increasing or decreasing
  if (xi(1) < xi(n)) then
    incr = +1
  else
    incr = -1
  end if

! Find interval x(i):x(i+1) containing x0
  found = .false.
  do i = 1,n-1
    if ((x0-xi(i))*incr>=0._dp .and. (xi(i+1)-x0)*incr>=0._dp) then
      found = .true.
      exit
    end if
  end do
  if (.not.found) stop 'mesh1D locate: ERROR: x0 out of range'

! Set range of points for Lagrange interpolation
  ! Use up to 3 neighbor points on each side
  i1 = max( i-2, 1 )
  i2 = min( i+3, n )
  ! Alternatively: use always 6 points
  if (i1==1) i2 = min( 1+5, n )
  if (i2==n) i1 = max( n-5, 1 )

! Find i0 within interval x(i):x(i+1) using bisection
! This is inefficient but it is expected to be used rarely
  il = i
  ir = i+1
  do iter = 1,1000
    ! Interpolate xi(i) at midpoint of bisection interval: xc=xi(ic)
    ic = (il+ir) / 2
    call polint( ivec(i1:i2), xi(i1:i2), i2-i1+1, ic, xc, dx )
    ! Check convergence and perform bisection
    if (abs(ir-il)<itol) then
      locate = ic
      return
    else if ((xc-x0)*incr > 0._dp) then
      ir = ic
    else
      il = ic
    end if
  end do ! iter
  stop 'mesh1D locate: ERROR: bisection not converged'

end function locate

!----------------------------------------------------------------

function interpolation_local( nnew, xnew, n, y, x, dx ) &
           result(interpolation)

  implicit none
  integer,            intent(in) :: nnew
  real(dp),           intent(in) :: xnew(nnew)
  integer,            intent(in) :: n
  real(dp),           intent(in) :: y(n)
  real(dp), optional, intent(in) :: x(n)
  real(dp), optional, intent(in) :: dx
  real(dp)                       :: interpolation(nnew)

  integer  :: i, i1, i2, inew
  real(dp) :: d2ydi2(n), dy, dydi, iofxnew, ynew

! Mesh-initialization
  if (present(x)) then
    call set_mesh( n, x=x )
  elseif (present(dx)) then
    call set_mesh( n, xmax=(n-1)*dx )
  endif

! Find interpolation
  if (interpolation_method=='spline') then

    ! Set spline interpolation of y(i) in the uniform mesh i
!    call spline( 1._dp, y, n, yp1*xp1(1), ypn*xp1(n), d2ydi2 )
    call spline( 1._dp, y, n, huge(y), huge(y), d2ydi2 )

    ! Find interpolation at xnew points
    do inew = 1,nnew
      ! Find iofxnew such that x(iofxnew)=xnew
      iofxnew = locate( xnew(inew) )
      ! Make sure that iofxnew is within range [1:n)
      if (iofxnew<1-xtol .or. iofxnew>n+xtol) then
        stop 'interpolation: ERROR: xnew out of range'
      else
        iofxnew = max( iofxnew, 1._dp )
        iofxnew = min( iofxnew, n*(1-epsilon(iofxnew)) )
      end if
      ! Interpolate y(i) at iofxnew
      ! Notice: iofxnew range is [1:n) but splint expects [0:n-1)
      call splint( 1._dp, y, d2ydi2, n, iofxnew-1, ynew, dydi )
      interpolation(inew) = ynew
    end do ! inew

  else if (interpolation_method=='lagrange') then

    ! Loop on xnew points
    do inew = 1,nnew
      ! Find iofxnew such that xi(iofxnew)=xnew
      iofxnew = locate( xnew(inew) )
      i = int(iofxnew)
      ! Use up to 3 neighbor points on each side
      i1 = max( i-2, 1 )
      i2 = min( i+3, n )
      ! Alternatively: use always 6 points
      if (i1==1) i2 = min( 1+5, n )
      if (i2==n) i1 = max( n-5, 1 )
      ! Now interpolate y(iofxnew) in the uniform mesh ivec=i
      call polint( ivec(i1:i2), y(i1:i2), i2-i1+1, iofxnew, ynew, dy )
      interpolation(inew) = ynew
    end do ! inew

  else
    stop 'interpolation: ERROR: bad interpolation_method parameter'
  end if ! (interpolation_method)

end function interpolation_local

!----------------------------------------------------------------

recursive function derivative( n, y, x, dx, order ) result(dydx)

  implicit none
  integer,            intent(in) :: n
  real(dp),           intent(in) :: y(n)
  real(dp), optional, intent(in) :: x(n)
  real(dp), optional, intent(in) :: dx
  integer,  optional, intent(in) :: order
  real(dp)                       :: dydx(n)

  integer  :: i, nx, ord
  real(dp) :: dxdi(n), dydi(n), d2ydi2(n)

! Mesh-initialization
  if (present(x)) then
    call set_mesh( n, x=x )
  elseif (present(dx)) then
    call set_mesh( n, xmax=(n-1)*dx )
  endif

! Fix order of the derivative
  if (present(order)) then
    ord = order
  else
    ord = 1
  endif

! Iterate over order for order>1
  if (ord > 1) then
    dydx = derivative( n, derivative(n,y), order=ord-1 )
  elseif (ord == 1) then

    if (interpolation_method=='spline') then

      call spline( 1._dp, y, n, yp1, ypn, d2ydi2 )
!      call spline( 1._dp, y, n, yp1*xp1(1), ypn*xp1(n), d2ydi2 )

      dydi(1) = y(2) - y(1) - (2*d2ydi2(1) + d2ydi2(2)) / 6
      dydi(n) = y(n) - y(n-1) + (d2ydi2(n-1) + 2*d2ydi2(n)) / 6
      dydi(2:n-1) = (y(3:n) - y(1:n-2)) / 2 &
                  - (d2ydi2(3:n) - d2ydi2(1:n-2)) / 12

    else if (interpolation_method=='lagrange') then

!     Find derivative of y with respect to i
!     Centered 5-point formulas for all but first/last two points
      do i = 3,n-2
        dydi(i)=(y(i-2)-8*y(i-1)+8*y(i+1)-y(i+2))/12
      enddo

!     Find derivative at first/last two points
      if (n<=1) then
        dydi = 0
      else if (n==2) then
        dydi = -y(1) + y(2)
      else if (n==3) then
        dydi(1) = (-3*y(1) + 4*y(2) -   y(3)) / 2
        dydi(3) = (   y(1) - 4*y(2) + 3*y(3)) / 2
        dydi(2) = (-  y(1)          +   y(3)) / 2
      else if (n==4) then
        ! Use up to 2 neighbor points on each side
!        dydi(1) = ( -3*y(1) + 4*y(2) -   y(3)) / 2
!        dydi(4) = (    y(2) - 4*y(3) + 3*y(4)) / 2
!        dydi(2) = (- 2*y(1) - 3*y(2) + 6*y(3) -   y(4)) / 6
!        dydi(3) = (    y(1) - 6*y(2) + 3*y(3) + 2*y(4)) / 6
        ! Alternatively: use all available points
        dydi(1) = (-11*y(1) +18*y(2) - 9*y(3) + 2*y(4)) / 6
        dydi(4) = (- 2*y(1) + 9*y(2) -18*y(3) +11*y(4)) / 6
      else
        ! Use up to 2 neighbor points on each side
!        dydi(1) = ( -3*y(1)   + 4*y(2)   -   y(3)) / 2
!        dydi(n) = (    y(n-2) - 4*y(n-1) + 3*y(n)) / 2
!        dydi(2)   = (- 2*y(1)   - 3*y(2)   + 6*y(3)   -   y(4)) / 6
!        dydi(n-1) = (    y(n-3) - 6*y(n-2) + 3*y(n-1) + 2*y(n)) / 6
        ! Alternatively: use always 5 points
        dydi(1)=(-25*y(1)+48*y(2)  -36*y(3)  +16*y(4)  -3*y(5)  )/12
        dydi(n)=(+25*y(n)-48*y(n-1)+36*y(n-2)-16*y(n-3)+3*y(n-4))/12
        dydi(2)  =(-3*y(1)-10*y(2)  +18*y(3)  -6*y(4)  +y(5)    )/12
        dydi(n-1)=(+3*y(n)+10*y(n-1)-18*y(n-2)+6*y(n-3)-y(n-4)  )/12
      end if ! (n)

    else
      stop 'derivative: ERROR: bad interpolation_method parameter'
    end if ! (interpolation_method)

!   Find derivative of x with respect to i
    call get_mesh( n, nx, dxdi=dxdi )

!   Find derivative of y with respect to x
    dydx = dydi / dxdi

  endif

end function derivative

!----------------------------------------------------------------

function integral( n, y, x, dx )

  implicit none
  real(dp)                       :: integral
  integer,            intent(in) :: n
  real(dp),           intent(in) :: y(n)
  real(dp), optional, intent(in) :: x(n)
  real(dp), optional, intent(in) :: dx

  real(dp):: d2fdi2(n), f(n)

! Mesh-initialization
  if (present(x)) then
    call set_mesh( n, x=x )
  elseif (present(dx)) then
    call set_mesh( n, xmax=(n-1)*dx )
  endif

! Find integral
  if (interpolation_method=='spline') then

    f = y * xp1(1:n)
    call spline( 1._dp, f, n, yp1, ypn, d2fdi2 )
!    call spline( 1._dp, y, n, yp1*xp1(1), ypn*xp1(n), d2ydi2 )

    integral = (f(1) + f(n))/2 - (d2fdi2(1) + d2fdi2(n)) / 24 &
             + sum(f(2:n-1)) - sum(d2fdi2(2:n-1)) / 12

  else if (interpolation_method=='lagrange') then

    if (n<=1) then
      integral = 0
    else if (n==2) then
      integral = ( y(1)*xp1(1) + y(2)*xp1(2) ) / 2
    else if (n==3) then
      integral = ( y(1)*xp1(1) + 4*y(2)*xp1(2) + y(3)*xp1(3) ) / 3
    else if (n==4) then
      integral = (   3 * (y(1)*xp1(1) + y(n)  *xp1(n)  )      &
                  +  9 * (y(2)*xp1(2) + y(n-1)*xp1(n-1)) )/8
    else if (n==5) then
      integral = (   9 * (y(1)*xp1(1) + y(n)  *xp1(n)  )      &
                  + 28 * (y(2)*xp1(2) + y(n-1)*xp1(n-1))      &
                  + 22 *  y(3)*xp1(3)                    )/24
    else
      integral = (   9 * (y(1)*xp1(1) + y(n)  *xp1(n)  )      &
                  + 28 * (y(2)*xp1(2) + y(n-1)*xp1(n-1))      &
                  + 23 * (y(3)*xp1(3) + y(n-2)*xp1(n-2)) )/24 &
               + sum( y(4:n-3)*xp1(4:n-3) )
    end if

  else
    stop 'integral: ERROR: bad interpolation_method parameter'
  end if ! (interpolation_method)

end function integral

!----------------------------------------------------------------

subroutine numerov( n, y, yp, ypp, f0, f1, x, dx, norm )

  implicit none
  integer,            intent(in)    :: n
  real(dp), optional, intent(inout) :: y(n)
  real(dp), optional, intent(out)   :: yp(n)
  real(dp), optional, intent(out)   :: ypp(n)
  real(dp), optional, intent(in)    :: f0(n)
  real(dp), optional, intent(in)    :: f1(n)
  real(dp), optional, intent(in)    :: x(n)
  real(dp), optional, intent(in)    :: dx
  real(dp), optional, intent(in)    :: norm

  integer                :: i
  real(dp)               :: ynorm
  real(dp), dimension(n) :: g, g0, g1, z, zp

! Mesh-initialization
  if (present(x)) then
    call set_mesh( n, x=x )
  elseif (present(dx)) then
    call set_mesh( n, xmax=(n-1)*dx )
  endif
  if (.not.present(y)) return

! Check that y has been initialized in left boundary
  if (.not.present(f0)      & ! Homogeneous equation
      .and.abs(y(1))==0._dp &
      .and.abs(y(2))==0._dp) then
    stop 'numerov: ERROR: first two values of y needed'
  endif

! Find g0 and g1
  g0 = 0._dp
  g1 = s2
  if (present(f0)) g0 = s0 * f0
  if (present(f1)) g1 = s1 * f1 + g1

! Integrate differential equation
  z(1) = y(1) / sqrxp(1)
  z(2) = y(2) / sqrxp(2)

  do i = 2,n-1
    z(i+1) = ( (g0(i-1)+10._dp*g0(i)+g0(i+1))/12._dp +   &
               (-1._dp+g1(i-1)/12._dp)*z(i-1)        +   &
               (2._dp+(10._dp/12._dp)*g1(i))*z(i)    ) / &
                  (1._dp-g1(i+1)/12._dp)
  enddo
  y = z * sqrxp

! Find first derivative y'=dy/dx
  if (present(yp)) then
    g = g0 + g1*z
    zp(2:n-1) = (z(3:n) - z(1:n-2))/2 + (g(3:n) - g(1:n-2))/6
    zp(n) = zp(n-1) + (5*g(n) + 8*g(n-1) - g(n-2))/12
    zp(1) = zp(2)   - (5*g(1) + 8*g(2)   - g(3)  )/12
    yp = ( zp + z*xp2/(2*xp1) )/sqrxp
  endif

! Find second derivative y''=d2y/dx2
  if (present(ypp)) then
    ypp = 0
    if (present(f0)) ypp = f0
    if (present(f1)) ypp = ypp + f1*y
  endif

! Normalize solution
  if (present(norm)) then
    if (present(f0)) then
      stop 'numerov: ERROR: cannot normalize inhomogeneous solutions'
    elseif (norm <= 0._dp) then
      print*, 'numerov: ERROR: nonpositive norm =', norm
      stop
    endif
    ynorm = integral( n, y=y*y )
    y = y * sqrt(norm/ynorm)
    if (present(yp))  yp  = yp  * sqrt(norm/ynorm)
    if (present(ypp)) ypp = ypp * sqrt(norm/ynorm)
  endif

end subroutine numerov

end module mesh1D
