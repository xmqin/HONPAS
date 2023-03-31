!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
!

! It contains several routines to compute the Gaussian weights and
! abscissas for several methods of the Gauss-family

! We have the following implemented:
!  1) Gauss-Legendre (Jacobi)
!     Used for \int_a^b f(x)w(x)dx , w(x) = 1
!  2) Gauss-Gegenbauer (Jacobi)
!     Used for \int_a^b f(x)w(x)dx , w(x) = \left((b-x)(x-a)\right)^{\lambda-1/2}
!  3) Gauss-Jacobi (Jacobi)
!     Used for \int_a^b f(x)w(x)dx , w(x) = (b-x)^\alpha (x-a)^\beta
!  4) Gauss-Chebyshev (Jacobi)
!     Used for \int_a^b f(x)w(x)dx , w(x) = \frac1{\sqrt{(b-x)(x-a)}}
!  5) Gauss-Laguerre (Laguerre)
!     Used for \int_a^\infty f(x)w(x)dx , w(x) = e^{-x}
!  6) Gauss-Generalized-Laguerre (Laguerre)
!     Used for \int_a^\infty f(x)w(x)dx , w(x) = x^\alpha e^{-x}
!  7) Gauss-Hermite
!     Used for \int_{-\infty}^\infty f(x)w(x)dx , w(x) = e^{-x^2}
!  8) Gauss-Tschebysheff
!     Used for \int_a^b f(x)w(x)dx , w(x) = \frac1{\sqrt{x(1-x)}}
!  9) Tanh-Sinh
!     Used for \int_a^b f(x)w(x)dx , w(x) = 1

! All can be solved using the Golub-Welsch algorithm.
! The routines for those are suffixed with GW
! The Gauss-Legendre can be solve with two methods of recursion:
!   1) Orignal header from the code:
!
! ********************* < HEADER START > **************************
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  gauss.f: Points and weights for Gaussian quadrature                 c
!           This is the Legendre Gaussian quadrature                   c
!           The Legendre quadrature is the w(x) = 1, and the rescaling c
!           of the function by the Fermi function                      c
!								       c
!  taken from: "Projects in Computational Physics" by Landau and Paez  c 
!	       copyrighted by John Wiley and Sons, New York            c      
!                                                                      c
!  written by: Oregon State University Nuclear Theory Group            c
!	       Guangliang He & Rubin H. Landau                         c
!  supported by: US National Science Foundation, Northwest Alliance    c
!                for Computational Science and Engineering (NACSE),    c
!                US Department of Energy 	                       c
!								       c
!  comment: error message occurs if subroutine called without a main   c
!  comment: this file has to reside in the same directory as integ.c   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     rescale rescales the gauss-legendre grid points and weights
!
!     npts     number of points
!     job = 0  rescalling uniformly between (a,b)
!           1  for integral (0,b) with 50% points inside (0, ab/(a+b))
!           2  for integral (a,inf) with 50% inside (a,b+2a)
!     x, w     output grid points and weights.
!
! ********************* < HEADER END > **************************
!      This routine will be referred by: Gauss_Legendre_Rec
!
!   2) As I am not quite sure what the above routine does, I have implemented
!      the Newton-Raphson method for computing the weights and abscissas.
!      This routine will be referred by: Gauss_Legendre_NR
! 
!
! The Gauss-Chebyshev has an explicit form of the weigts and absicassas.
! Hence a specific routine can be used for this (should be encouraged)
!      Gauss_Chebyshev_Exact


module m_gauss_quad
  
  implicit none

  integer,  parameter :: dp      = selected_real_kind(p=10)
  real(dp), parameter :: Pi      = 3.14159265358979323846264338328_dp
  real(dp), parameter :: Pihalf  = 1.57079632679489655799898173427_dp 
  real(dp), parameter :: def_EPS = 1.e-12_dp

  private 

  ! Weight and abscissa creation
  public :: Gauss_Legendre_Rec
  public :: Gauss_Legendre_NR
  public :: Gauss_Legendre_GW
  public :: Gauss_Generalized_Laguerre_GW
  public :: Gauss_Laguerre_GW
  public :: Gauss_Jacobi_GW
  public :: Gauss_Gegenbauer_GW
  public :: Gauss_Hermite_GW
  public :: Gauss_Chebyshev_Exact
  public :: Gauss_Chebyshev_GW
  public :: TanhSinh_Exact

  ! Weights
  public :: w_Gauss_Chebyshev
  public :: w_Gauss_Jacobi
  public :: w_Gauss_Gegenbauer
  public :: w_Gauss_Generalized_Laguerre
  public :: w_Gauss_Laguerre
  public :: w_Gauss_Hermite

contains

! Functions for retrieving the Gauss weights and abscissas

  subroutine Gauss_Legendre_NR(N,x,w,a,b) 
    integer, intent(in) :: N
    real(dp), intent(out) :: x(N),w(N)
    real(dp), intent(in), optional :: a, b

    real(dp) :: t,t1,pp,p1,p2,p3,aj,dx

    integer :: i,j 
    real(dp), parameter :: eps = 1.0d-14

    ! equal spacing
    dx = -1._dp
    ! We use the recursion relation and the Newton-Raphson method
    do i = 1 , N
       t  = cos(Pi*(2._dp*(i-1._dp)+1._dp)/(2._dp*N)) + &
            0.27_dp/real(N,dp)*sin(Pi*dx*real(N-1,dp)/real(N+1,dp))
       dx = dx + 2._dp / real(N + 1,dp)
       t1 = huge(1._dp)
       do while ( abs(t-t1) > eps )
          p1 = t
          p2 = 1._dp
          aj = 1._dp
          ! Compute the legendre polynomials
          ! P_0(x) = 1 ; P_1(x) = x
          ! P_n(x) = \frac{2n+1}{n+1}xP_n(x) - \frac n{n+1}P_{n-1}(x)
          do j = 2 , N
             p3 = p2
             p2 = p1
             aj = aj + 1._dp
             p1 = ((2._dp*aj-1._dp)*t*p2-(aj-1._dp)*p3)/aj
          end do
          pp = real(N+1,dp)*(p2-t*p1)/(1._dp - t**2)
          t1 = t
          t  = t1 - p1/pp
       end do
       
       x(N+1-i) = t
       w(N+1-i) = 2._dp/((1._dp-t**2)*pp**2) * &
            (real(N+1,dp)/real(N,dp))**2
    end do

    call general_scale(N,x,w,a=a,b=b)
    
  end subroutine Gauss_Legendre_NR

  subroutine TanhSinh_Exact(N,x,w,a,b,p) 
    integer, intent(in) :: N
    real(dp), intent(out) :: x(N),w(N)
    real(dp), intent(in), optional :: a, b, p

    real(dp) :: lp, h, cp
    real(dp) :: u1, u2, ch
    real(dp) :: arg, cums, k
    integer :: i, j
    logical :: is_even

    ! finite spacing so that the weight becomes
    ! This is some arbitrary weight size 
    ! I have tried for various values between 2 and 5 but 4 seems
    ! to have the best properties, I am not sure whether this should
    ! be a multiple of two
    ! I have adapted this method instead as it provides some way
    ! of controlling the precision of the integral.

    lp = 1.e-2_dp
    if ( present(p) ) lp = p
    if ( lp < 0._dp ) call die('Error in precision request &
         &negative numbers can not be precioned.')

    is_even = mod(N,2) == 0
    if ( is_even ) then
       j = N / 2 + 1
    else
       j = N / 2 + 2
    end if

    k = 32._dp
    cp = lp - 1._dp
    do while ( lp > cp )
       h = 16._dp / k

       ! Weight-sum
       if ( .not. is_even ) then
          ! Add the middle abscissa-weight
          ! u1 = cosh( 0._dp ) 
          ! cosh( 0 ) == 1
          ! u2 = Pihalf * sinh( 0._dp )
          ! sinh( 0 ) == 0
          !ch = cosh( u2 ) ** 2
          ! cosh( 0 ) == 1
          ! hence w = 1

          ! for correct normalization (added a half due to 2 times later)
          cums = .5_dp
       else
          cums = 0._dp
       end if

       cp = .5_dp * h
       do i = j , N
          ! The argument for the TANH-SINH rule
          arg = cp * real(2*i-N-1,dp)
          
          u1 = cosh( arg )
          u2 = Pihalf * sinh( arg )
          ch = cosh( u2 ) ** 2

          ! for correct normalization
          cums = cums + u1 / ch

       end do

       ! we have a symmetric distribution (hence the sum of
       ! the cumultative sum divides out)
       cums = 1._dp / cums

       cp = u1 / ch * cums

       k = k + 1._dp
    end do

    cp = .5_dp * h
    do i = 1 , N
       ! The argument for the TANH-SINH rule
       arg = cp * real(2*i-N-1,dp)
          
       u1 = cosh( arg )
       u2 = Pihalf * sinh( arg )
       ch = cosh( u2 )
          
       ! Calculate the abscissa-weight pair
       x(i) = 1._dp - 1._dp / ( exp( u2 ) * ch )
       ! We immediately normalize the weight (remember that
       ! cums was the inverse sum)
       w(i) = cums * u1 / ch ** 2

    end do

    call general_scale(N,x,w,a=a,b=b)
    
  end subroutine TanhSinh_Exact



! We employ the Gauss-Chebyshev quadrature rules
! This will enable to do the integral:
! \int_a^b
! It requires a weight function
! w(x) = \frac1{\sqrt{1-x^2}}
  subroutine Gauss_Chebyshev_Exact(N,x,w,a,b,method, pure) 
    integer, intent(in)   :: N
    real(dp), intent(out) :: x(N), w(N)
    real(dp), intent(in), optional :: a, b
    integer, intent(in), optional :: method
    logical, intent(in), optional :: pure
    real(dp) :: a_arg
    integer :: lmethod, lN, lnn
    integer :: i
    
    ! Option collecting..
    lmethod = 1 ! open-method
    if ( present(method) ) lmethod = method
    ! If method is 0 it is closed
    ! if method is 1 it is open
    if ( lmethod < 0 .or. 1 < lmethod ) call die('Wrong method for &
         &Gauss-Chebyshev quadrature')

    if ( lmethod == 1 ) then
       ! calculate the Chebyshev - Open
       w(:) = Pi / real(N,dp)
       do i = 1 , N
          a_arg = .5_dp * real(2*i-1,dp) * Pi / real(N, dp)
          x(N+1-i)  = cos(a_arg)
       end do
    else if ( lmethod == 0 ) then
       i = 0
       if ( present(pure) ) then
          if ( .not. pure ) i = 1
       end if
       if ( i == 1 ) then
          ! when having the closed boundary 
          ! and we need to correct the weight, we 
          ! loose the ends, hence we need to
          ! construct it for N+2
          ! this is because abs(cos(N*Pi/N)) == abs(cos(0)) == 1
          lN = N
          lnn = lN + 1
          w(:) = Pi / real(lnn,dp)
          do i = 1 , lnn - 1
             a_arg = real(i,dp) * Pi / real(lnn, dp)
             x(N+1-i)  = cos(a_arg)
          end do
       else
          lN = N - 2
          lnn = lN + 1
          w(:) = Pi / real(lnn,dp)
          w(1) = Pi / real(2*lnn,dp)
          w(N) = Pi / real(2*lnn,dp)
          do i = 0 , lnn
             a_arg  = real(i,dp) * Pi / real(lnn, dp)
             x(N-i) = cos(a_arg)
          end do
       end if
    end if

    ! if we have a non-pure integration function, i.e.
    !  f(x) = g(x) * \sqrt{(b-x)*(x-a)}
    !  w(x) = 1/\sqrt{(b-x)*(x-a)}
    if ( present(pure) ) then
       if ( .not. pure ) then
          w(:) = w(:) / w_Gauss_Chebyshev(x)
       end if
    end if

    call general_scale(N,x,w,a=a,b=b)
       
  end subroutine Gauss_Chebyshev_Exact


! The recursive Gauss-Legendre
  subroutine Gauss_Legendre_Rec(npts,job,a,b,x,w) 
    integer npts,job,m,i,j 
    real(dp) x(npts),w(npts),a,b,xi
    real(dp) t,t1,pp,p1,p2,p3,aj
    real(dp) eps,pi,zero,two,one,half,quarter
    parameter (pi = 3.14159265358979323846264338328d0,eps = 3.0d-14)
    parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
    parameter (half=0.5d0,quarter=0.25d0)
!
! *** FIRST EXECTUABLE *************************************************
!
    m=(npts+1)/2
    do i=1,m
       t=cos(pi*(i-quarter)/(npts+half))
       t1 = t + one
       do while ( abs(t-t1).gt.eps )
          p1=one
          p2=zero
          aj=zero
          do j=1,npts
             p3=p2
             p2=p1
             aj=aj+one
             p1=((two*aj-one)*t*p2-(aj-one)*p3)/aj
          end do
          pp=npts*(t*p1-p2)/(t*t-one)
          t1=t
          t=t1-p1/pp
       end do
       x(i)=-t
       x(npts+1-i)=t
       w(i)=two/((one-t*t)*pp*pp)
       w(npts+1-i)=w(i)
    end do
!
! rescale the grid points 
    if (job.eq.0) then
!     scale to (a,b) uniformly
       do i=1,npts
          x(i)=x(i)*(b-a)/two+(b+a)/two
          w(i)=w(i)*(b-a)/two
       end do
    elseif (job.eq.1) then
! scale to (0,b) with 50% points inside (0,ab/(a+b))
       do i=1,npts
          xi=x(i)
          x(i)=a*b*(one+xi)/(b+a-(b-a)*xi)
          w(i)=w(i)*two*a*b*b/((b+a-(b-a)*xi)*(b+a-(b-a)*xi))
       end do
    elseif (job.eq.2) then
! scale to (a,inf) with 50% points inside (a,b+2a)
       do i=1,npts
          xi=x(i)
          x(i)=(b*xi+b+a+a)/(one-xi)
          w(i)=w(i)*two*(a+b)/((one-xi)*(one-xi))
       end do
    else
       call die('gauss: Wrong value of job for Gaussian quadrature')
    endif
  end subroutine Gauss_Legendre_Rec


! We employ the Gauss-Gegenbauer quadrature rule
! This will enable to do the integral:
! \int_a^b
! It requires a weight function
! w(x) = (b-x)^\alpha * (x-a)^\beta
  subroutine Gauss_Jacobi_GW(N,x,w,alpha,beta,a,b,pure,EPS)
    integer,  intent(in)  :: N
    real(dp), intent(out) :: x(N), w(N)
    real(dp), intent(in)  :: alpha, beta
    real(dp), intent(in), optional :: a, b
    logical,  intent(in), optional :: pure
    real(dp), intent(in), optional :: EPS

    integer :: N_int_mu
    real(dp) :: mu0, mu_diff, lEPS, la, lb
    real(dp), allocatable :: J(:,:)

    if ( alpha <= -1._dp ) call die('Can not perform Jacobi integration with &
         &wrong alpha, must be larger than -1.')
    if ( beta  <= -1._dp ) call die('Can not perform Jacobi integration with &
         &wrong beta, must be larger than -1.')

    la = -1._dp
    lb =  1._dp
    if ( present(a) ) la = a
    if ( present(b) ) lb = b

    lEPS = def_EPS
    if ( present(EPS) ) lEPS = EPS
    ! We first need to calculate the Jacobi weight integral
    N_int_mu = 100
    allocate(J(N_int_mu,2))
    call Gauss_Legendre_Rec(N_int_mu,0,la,lb,J(:,1),J(:,2))
    mu0 = sum(w_Gauss_Jacobi(J(:,1),alpha,beta,a=la,b=lb)*J(:,2))
    deallocate(J) ; N_int_mu = N_int_mu + 100
    do
       allocate(J(N_int_mu,2))
       call Gauss_Legendre_Rec(N_int_mu,0,la,lb,J(:,1),J(:,2))
       mu_diff = mu0
       mu0 = sum(w_Gauss_Jacobi(J(:,1),alpha,beta,a=la,b=lb)*J(:,2))
       deallocate(J)
       call adjust_Int(mu_diff,mu0,lEPS,N_int_mu)
!       print *,mu_diff,mu0
       if ( mu_diff <= lEPS ) exit
    end do

    ! Allocate for work-space (no need for memory control, we should not use that many points)
    allocate(J(0:N-1,0:N-1))

    ! Initialize the Jacobi matrix for the Golub-Welsch algorithm
    call Gauss_Jacobi_general(alpha,beta,N,J)

    ! Find the weights and abscissas by Golub-Welsch
    call solve_Golub_Welsch(N,J,x,w,mu0)

    ! clean-up
    deallocate(J)

    if ( present(pure) ) then
       if ( .not. pure ) then
          ! if we have a non-pure integration function, i.e.
          !  f(x) = g(x) / ((b-x)^\alpha * (x-a)^\beta)
          !  w(x) = (b-x)^\alpha * (x-a)^\beta
          w(:) = w(:) / w_Gauss_Jacobi(x,alpha,beta,a=a,b=b)
       end if
    end if
       
    ! The weight scaling should be corrected by the mu0
    call general_scale(N,x,w,a=a,b=b,adjust_w=.false.)!,power_slope=alpha+beta+1._dp)

  end subroutine Gauss_Jacobi_GW


! We employ the Gauss-Gegenbauer quadrature rule
! This will enable to do the integral:
! \int_a^b
! It requires a weight function
! w(x) = (1-x^2)^{\lambda-1/2}
  subroutine Gauss_Gegenbauer_GW(N,x,w,lambda,a,b,pure, EPS)
    integer,  intent(in)  :: N
    real(dp), intent(out) :: x(N), w(N)
    real(dp), intent(in)  :: lambda
    real(dp), intent(in), optional :: a, b
    logical,  intent(in), optional :: pure
    real(dp), intent(in), optional :: EPS

    integer :: N_int_mu
    real(dp) :: mu0, mu_diff, lEPS
    real(dp), allocatable :: J(:,:)

    if ( lambda <= -.5_dp ) call die('Can not perform Gegenbauer integration with &
         &wrong lambda, must be larger than -1/2. If lamda==-1/2 use Chebyshev.')
    
    lEPS = def_EPS
    if ( present(EPS) ) lEPS = EPS
    ! We first need to calculate the Jacobi weight integral
    N_int_mu = 100
    allocate(J(N_int_mu,2))
    call Gauss_Chebyshev_Exact(N_int_mu,J(:,1),J(:,2),a=a,b=b)
    mu0 = sum(w_Gauss_Gegenbauer(J(:,1),lambda+.5_dp,a=a,b=b)*J(:,2))
    deallocate(J) ; N_int_mu = N_int_mu + 100
    do
       allocate(J(N_int_mu,2))
       call Gauss_Chebyshev_Exact(N_int_mu,J(:,1),J(:,2),a=a,b=b)
       mu_diff = mu0
       mu0 = sum(w_Gauss_Gegenbauer(J(:,1),lambda+.5_dp,a=a,b=b)*J(:,2))
       deallocate(J)
       call adjust_Int(mu_diff,mu0,lEPS,N_int_mu)
!       print *,mu_diff,mu0
       if ( mu_diff <= lEPS ) exit
    end do

    ! Allocate for work-space (no need for memory control, we should not use that many points)
    allocate(J(N,N))

    ! Initialize the Jacobi matrix for the Golub-Welsch algorithm
    call Gauss_Jacobi_general(lambda-.5_dp,lambda-.5_dp,N,J)

    ! Find the weights and abscissas by Golub-Welsch
    call solve_Golub_Welsch(N,J,x,w,mu0)

    ! clean-up
    deallocate(J)

    if ( present(pure) ) then
       if ( .not. pure ) then
          ! if we have a non-pure integration function, i.e.
          !  f(x) = g(x) / (1-x^2)^{\lambda-1/2}
          !  w(x) = (1-x^2)^{\lambda-1/2}
          w(:) = w(:) / w_Gauss_Gegenbauer(x,lambda)
       end if
    end if

    ! The mu0 should have corrected the weights
    call general_scale(N,x,w,a=a,b=b,adjust_w=.false.)

  end subroutine Gauss_Gegenbauer_GW

! We employ the Gauss-Legendre quadrature rule
! This will enable to do the integral:
! \int_a^b
! It requires a weight function
! w(x) = 1
  subroutine Gauss_Legendre_GW(N,x,w,a,b) 
    integer, intent(in) :: N
    real(dp), intent(out) :: x(N), w(N)
    real(dp), intent(in), optional :: a, b

    real(dp), allocatable :: J(:,:)
    real(dp) :: mu0

    ! integration \int_a^b 1\dx (the general-scale does the same)
    mu0 = 2._dp
    if ( present(a) ) mu0 = mu0 - 1._dp - a
    if ( present(b) ) mu0 = mu0 - 1._dp + b

    ! Allocate for work-space (no need for memory control, we should not use that many points)
    allocate(J(N,N))

    ! Initialize the Jacobi matrix for the Golub-Welsch algorithm
    call Gauss_Jacobi_general(0._dp,0._dp,N,J)

    ! Find the weights and abscissas by Golub-Welsch
    call solve_Golub_Welsch(N,J,x,w,mu0)

    ! clean-up
    deallocate(J)
    
    call general_scale(N,x,w,a=a,b=b,adjust_w=.false.)

  end subroutine Gauss_Legendre_GW


! We employ the Gauss-Chebyshev quadrature rules
! This will enable to do the integral:
! \int_a^b
! It requires a weight function
! w(x) = \frac1{\sqrt{1-x^2}}
  subroutine Gauss_Chebyshev_GW(N,x,w,a,b, pure,EPS) 
    integer, intent(in)   :: N
    real(dp), intent(out) :: x(N), w(N)
    real(dp), intent(in), optional :: a, b
    logical, intent(in), optional :: pure
    real(dp), intent(in), optional :: EPS

    integer :: N_int_mu
    real(dp), allocatable :: J(:,:)
    real(dp) :: mu0, mu_diff, lEPS

    ! We first need to calculate the Chebyshev weight integral
    lEPS = def_EPS
    if ( present(EPS) ) lEPS = EPS
    N_int_mu = 100
    allocate(J(N_int_mu,2))
    call Gauss_Legendre_GW(N_int_mu,J(:,1),J(:,2),a=a,b=b)
    mu0 = sum(w_Gauss_Chebyshev(J(:,1),a=a,b=b)*J(:,2))
    deallocate(J) ; N_int_mu = N_int_mu + 100
    do
       allocate(J(N_int_mu,2))
       call Gauss_Legendre_GW(N_int_mu,J(:,1),J(:,2),a=a,b=b)
       mu_diff = mu0
       mu0 = sum(w_Gauss_Chebyshev(J(:,1),a=a,b=b)*J(:,2))
       deallocate(J)
       call adjust_Int(mu_diff,mu0,lEPS,N_int_mu)
!       print *,mu_diff,mu0
       if ( mu_diff <= lEPS ) exit
    end do

    allocate(J(N,N))
    call Gauss_Jacobi_general(-.5_dp,-.5_dp,N,J)

    ! Find the weights and abscissas by Golub-Welsch
    call solve_Golub_Welsch(N,J,x,w,mu0)

    deallocate(J)

    ! if we have a non-pure integration function, i.e.
    !  f(x) = g(x) * \sqrt{1-x^2}
    !  w(x) = 1/\sqrt{1-x^2}
    ! then we need to add the \sqrt{1-x^2} to the weight before we shift,
    ! otherwise we enter the complex realm! ;)
    if ( present(pure) ) then
       if ( .not. pure ) then
          w(:) = w(:) / w_Gauss_Chebyshev(x,a=a,b=b)
       end if
    end if

    call general_scale(N,x,w,a=a,b=b,adjust_w=.false.)
       
  end subroutine Gauss_Chebyshev_GW


! We employ the Gauss-Laguerre quadrature rule
! This will enable to do the integral:
! \int_0^\infty (it cannot easily be generalized to \int_a^\infty
! It requires a weight function
! w(x) = x^\alpha e^{-x}
  subroutine Gauss_Generalized_Laguerre_GW(N,x,w,alpha, a,pure, EPS) 
    integer,  intent(in)  :: N
    real(dp), intent(out) :: x(N), w(N)
    real(dp), intent(in)  :: alpha
    real(dp), intent(in), optional :: a
    logical,  intent(in), optional :: pure
    real(dp), intent(in), optional :: EPS

    integer :: N_int_mu
    real(dp), allocatable :: J(:,:)
    real(dp) :: mu0, mu_diff, lEPS, la
    
    ! We first need to calculate the Generalized Laguerre weight integral
    ! Luckily the Laguerre weight function is implicit in the
    ! function
    lEPS = def_EPS
    if ( present(EPS) ) lEPS = EPS
    la = 0._dp
    if ( present(a) ) la = a
    ! We first need to calculate the Laguerre weight integral
    N_int_mu = 100
    allocate(J(N_int_mu,2))
    call Gauss_Laguerre_GW(N_int_mu,J(:,1),J(:,2),a=a)
    mu0 = sum((J(:,1)-la)**alpha*J(:,2))
    deallocate(J) ; N_int_mu = N_int_mu + 100
    do
       allocate(J(N_int_mu,2))
       call Gauss_Laguerre_GW(N_int_mu,J(:,1),J(:,2),a=a)
       mu_diff = mu0
       ! As the Laguerre weight is implicit in the function
       ! we need only integrating x^alpha
       mu0 = sum((J(:,1)-la)**alpha*J(:,2))
       deallocate(J)
       call adjust_Int(mu_diff,mu0,lEPS,N_int_mu)
!       print *,mu_diff,mu0
       if ( mu_diff <= lEPS ) exit
    end do

    ! Allocate for work-space (no need for memory control, we should not use that many points)
    allocate(J(N,N))

    ! Initialize the Jacobi matrix for the Golub-Welsch algorithm
    call Gauss_Laguerre_general(alpha,N,J)

    ! Find the weights and abscissas by Golub-Welsch
    call solve_Golub_Welsch(N,J,x,w,mu0)

    ! clean-up
    deallocate(J)

    ! if we have a non-pure integration function, i.e.
    !  f(x) = g(x) * \exp{x-a}
    !  w(x) = \exp{-x+a}
    ! then we need to factor \exp{x-a} to the weight before we shift
    if ( present(pure) ) then
       if ( .not. pure ) then
          w(:) = w(:) / w_Gauss_Generalized_Laguerre(x,alpha,a=a)
       end if
    end if

    if ( present(a) ) then
       x(:) = x(:) + a
    end if
    ! potentially we could introduce a bounded integration
    ! b, this will yield a slope on the x_i

  end subroutine Gauss_Generalized_Laguerre_GW


! We employ the Gauss-Laguerre quadrature rule
! This will enable to do the integral:
! \int_a^\infty
! It requires a weight function
! w(x) = e^{-x}
  subroutine Gauss_Laguerre_GW(N,x,w,a, pure) 
    integer,  intent(in)  :: N
    real(dp), intent(out) :: x(N), w(N)
    real(dp), intent(in), optional :: a
    logical,  intent(in), optional :: pure

    real(dp), allocatable :: J(:,:)
    real(dp) :: mu0
    
    ! For normalizing the weight-function
    mu0 = 1._dp
    ! Luckily this integral is easy! ;)
    if ( present(a) ) mu0 = exp(-a)

    ! Allocate for work-space
    ! (no need for memory control, we should not use that many points)
    allocate(J(N,N))

    ! Initialize the Jacobi matrix for the Golub-Welsch algorithm
    call Gauss_Laguerre_general(0._dp,N,J)

    ! Find the weights and abscissas by Golub-Welsch
    call solve_Golub_Welsch(N,J,x,w,mu0)

    ! clean-up
    deallocate(J)

    ! if we have a non-pure integration function, i.e.
    !  f(x) = g(x) * \exp{x-a}
    !  w(x) = \exp{-x+a}
    ! then we need to factor \exp{x} to the weight before we shift
    if ( present(pure) ) then
       if ( .not. pure ) then
          w(:) = w(:) / w_Gauss_Laguerre(x,a=a)
       end if
    end if

    ! Change integration coordinates
    ! Weights are corrected by mu0
    if ( present(a) ) then
       x(:) = x(:) + a
    end if
           
  end subroutine Gauss_Laguerre_GW



! We employ the Gauss-Hermite quadrature rule
! This will enable to do the integral:
! \int_{-\infty}^\infty
! It requires a weight function
! w(x) = e^{-x^2}
  subroutine Gauss_Hermite_GW(N,x,w, pure) 
    integer, intent(in)   :: N
    real(dp), intent(out) :: x(N), w(N)
    logical, intent(in), optional :: pure

    real(dp), allocatable :: J(:,:)
    real(dp) :: mu0
    integer :: i
    
    ! For normalizing the weight-function
    ! mu0 = \int_{-\infty}^\infty \exp{-x^2}
    mu0 = sqrt(Pi)

    ! Allocate for work-space (no need for memory control, we should not use that many points)
    allocate(J(0:N-1,0:N-1))

    ! Initialize the Jacobi matrix for the Golub-Welsch algorithm
    J(:,:) = 0._dp
    do i = 1 , N - 1
       ! off-diagonal element
       J(i,i-1) = sqrt(.5_dp*real(i,dp))
       J(i-1,i) = J(i,i-1)
    end do

    ! Find the weights and abscissas by Golub-Welsch
    call solve_Golub_Welsch(N,J,x,w,mu0)

    ! clean-up
    deallocate(J)

    ! if we have a non-pure integration function, i.e.
    !  f(x) = g(x) * \exp{x^2}
    !  w(x) = \exp{-x^2}
    ! then we need to fator \exp{x^2} to the weight
    if ( present(pure) ) then
       if ( .not. pure ) then
          w(:) = w(:) / w_Gauss_Hermite(x)
       end if
    end if

  end subroutine Gauss_Hermite_GW

  subroutine solve_Golub_Welsch(N,J,x,w,mu0)
    integer, intent(in) :: N
    real(dp), intent(inout) :: J(N,N)
    real(dp), intent(out)   :: x(N), w(N)
    real(dp), intent(in)    :: mu0
    integer :: lwork, info, i, k
    real(dp), allocatable :: work(:)

    ! Calculate optimal work array size
    lwork = -1
    call dsyev('V', 'U', N, J, N, &
         x, x, lwork, info)
    lwork = nint(x(1))
    allocate(work(lwork))

    ! Solve eigen value problem
    call dsyev('V', 'U', N, J, N, &
         x, work, lwork, info)
    if ( info /= 0 ) then
       call die('Golub-Welsch algorithm solution could not be &
            &found, are you sure to have a symmetric, positive &
            &definite matrix?')
    end if
    
    ! Calculate weights
    do i = 1 , N
       w(i) = 0._dp
       do k = 1 , N
          ! it should be normalized, however, for consistency
          w(i) = w(i) + J(k,i) ** 2
       end do
       w(i) = mu0 * J(1,i) ** 2  / w(i)
    end do
    
    deallocate(work)
    
  end subroutine solve_Golub_Welsch
  
  subroutine Gauss_Jacobi_general(alpha,beta,N,J)
    real(dp), intent(in) :: alpha, beta
    integer, intent(in) :: N
    real(dp), intent(out) :: J(0:N-1,0:N-1)
    real(dp) :: a_im1, c_i, b_i, L_n
    integer :: i
    
    J(:,:) = 0._dp
    i = 0
    b_i = (beta-alpha)/(alpha+beta+2._dp)
    J(i,i) = b_i
    a_im1 = 2._dp/(alpha+beta+2._dp)
    do i = 1 , N - 1 
       L_n = 2._dp*i + alpha + beta
       b_i = (beta**2-alpha**2)/(L_n*(L_n+2._dp))
       c_i = 2._dp*(i+alpha)*(i+beta)/(L_n*(L_n+1._dp))
       if ( i > 1 ) then
          L_n = 2._dp*(i-1._dp) + alpha + beta
          a_im1 = 2._dp*i*(i+alpha+beta)/((L_n+1._dp)*(L_n+2._dp))
       end if
       J(i,i)   = b_i
       J(i-1,i) = sqrt(a_im1*c_i)
       J(i,i-1) = sqrt(a_im1*c_i)
    end do
  end subroutine Gauss_Jacobi_general

  subroutine Gauss_Laguerre_general(alpha,N,J)
    real(dp), intent(in) :: alpha
    integer, intent(in) :: N
    real(dp), intent(out) :: J(0:N-1,0:N-1)
    real(dp) :: a_im1, c_i, b_i
    integer :: i
    
    J(:,:) = 0._dp
    i = 0
    b_i = 2._dp*i+alpha+1._dp
    J(i,i) = b_i
    do i = 1 , N - 1 
       b_i = 2._dp*i+alpha+1._dp
       c_i = -(real(i,dp)+alpha)
       a_im1 = -real(i,dp)
       J(i,i) = b_i
       J(i-1,i) = sqrt(a_im1*c_i)
       J(i,i-1) = sqrt(a_im1*c_i)
    end do
  end subroutine Gauss_Laguerre_general


  ! General scaling algorithm
  subroutine general_scale(N,x,w,a,b,adjust_x,adjust_w,power_slope)
    integer, intent(in) :: N
    real(dp), intent(inout) :: x(N), w(N)
    real(dp), intent(in), optional :: a, b
    logical, intent(in), optional :: adjust_x, adjust_w
    real(dp), intent(in), optional :: power_slope
    logical :: ladjust_x, ladjust_w
    real(dp) :: la, lb, slp, slpp
    
    integer :: i
    
    la = -1._dp
    if ( present(a) ) la = a
    lb = 1._dp
    if ( present(b) ) lb = b
    ladjust_x = .true.
    ladjust_w = .true.
    if ( present(adjust_x) ) ladjust_x = adjust_x
    if ( present(adjust_w) ) ladjust_w = adjust_w

    slp = .5_dp * ( lb - la )
    slpp = slp 
    if ( present(power_slope) ) slpp = slp ** power_slope

    ! scale to (a,b) uniformly
    do i = 1 , N
       if ( ladjust_x ) x(i) = x(i)*slp+(lb+la)*.5_dp
       if ( ladjust_w ) w(i) = w(i)*slpp
    end do
  end subroutine general_scale

  subroutine adjust_Int(mu0_old, mu0, EPS, N_int_mu)
    real(dp), intent(inout) :: mu0_old
    real(dp), intent(in) :: mu0
    real(dp), intent(in) :: EPS
    integer, intent(inout) :: N_int_mu
    real(dp) :: mu_diff
    mu_diff = abs(mu0_old - mu0)
    if ( mu_diff > EPS * 1000._dp ) then
       N_int_mu = N_int_mu + 500
    else if ( mu_diff > EPS * 100._dp ) then
       N_int_mu = N_int_mu + 300
    else
       N_int_mu = N_int_mu + 100
    end if
    mu0_old = mu_diff
  end subroutine adjust_Int

! Here we contain the weight functions if ever needed
  elemental function w_Gauss_Chebyshev(x,a,b) result(w)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b
    real(dp) :: w
    if ( present(a) ) then
       if ( present(b) ) then
          w = ((b-x)*(x-a))**(-.5_dp)
       else
          w = ((1._dp-x)*(x-a))**(-.5_dp)
       end if
    else
       if ( present(b) ) then
          w = ((b-x)*(x+1._dp))**(-.5_dp)
       else
          w = (1._dp-x**2)**(-.5_dp) ! 1._dp/sqrt(1._dp-x**2)
       end if
    end if
  end function w_Gauss_Chebyshev
  
  elemental function w_Gauss_Tschebyscheff(x,b) result(w)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: b
    real(dp) :: w
    if ( present(b) ) then
       w = (x*(b-x))**(-.5_dp) ! 1._dp/sqrt(x*(1._dp-x))
    else
       w = (x*(1._dp-x))**(-.5_dp) ! 1._dp/sqrt(x*(1._dp-x))
    end if
  end function w_Gauss_Tschebyscheff

  elemental function w_Gauss_Jacobi(x,alpha,beta,a,b) result(w)
    real(dp), intent(in) :: x, alpha, beta
    real(dp), intent(in), optional :: a, b
    real(dp) :: w
    if ( present(a) ) then
       if ( present(b) ) then
          w = (b - x) ** alpha * (x-a) ** beta
       else
          w = (1._dp - x) ** alpha * (x-a) ** beta
       end if
    else
       if ( present(b) ) then
          w = (b - x) ** alpha * (x+1._dp) ** beta
       else
          w = (1._dp - x) ** alpha * (x+1._dp) ** beta
       end if
    end if
  end function w_Gauss_Jacobi

  elemental function w_Gauss_Gegenbauer(x,lambda,a,b) result(w)
    real(dp), intent(in) :: x, lambda
    real(dp), intent(in), optional :: a, b
    real(dp) :: w
    w = w_Gauss_Jacobi(x,lambda-.5_dp,lambda-.5_dp,a=a,b=b)
  end function w_Gauss_Gegenbauer

  elemental function w_Gauss_Generalized_Laguerre(x,alpha,a) result(w)
    real(dp), intent(in) :: x, alpha
    real(dp), intent(in), optional :: a
    real(dp) :: w
    w = (x-a)**alpha * exp(-x+a)
  end function w_Gauss_Generalized_Laguerre

  elemental function w_Gauss_Laguerre(x,a) result(w)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a
    real(dp) :: w
    w = exp(-x+a)
  end function w_Gauss_Laguerre

  elemental function w_Gauss_Hermite(x) result(w)
    real(dp), intent(in) :: x
    real(dp) :: w
    w = exp(-x**2)
  end function w_Gauss_Hermite

end module m_gauss_quad
