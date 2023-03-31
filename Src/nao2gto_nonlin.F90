! *** Module: nao2gto_nonlin ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Approximate Natural Atomic Orbitals by a sum of Gaussians through
!!        non-linear optimization methods
!!
!! \author Fernando Gomez Ortiz
!! \author Carlos Beltran
!! \author Yann Pouillon
!!
!! \copyright
!!      - 2019-2020 SIESTA Developers Group
!!
!! \par History
!!      - 03.2020 Imported and rewrote to be used in SIESTA [Yann Pouillon]
! *****************************************************************************
module nao2gto_nonlin

  implicit none

  private

  public :: nao2gto_gaussfit

contains

  ! ***************************************************************************
  ! *** Public routines                                                     ***
  ! ***************************************************************************

  subroutine nao2gto_gaussfit( ngfs_max, rad_pts, orb_pts, hfx_opts, &
 &                             orb_g, ierr )

    use precision, only: dp
    use nao2gto_types,  only: hfx_options_type

    implicit none

!   Arguments
    integer,  intent(in)  :: ngfs_max          ! Maximum number of Gaussians
                                               !   in the linear combination to 
                                               !   fit the NAO
    real(dp), intent(in)  :: rad_pts(:)        ! Values of the points in the
                                               !   linear mesh where the 
                                               !   radial part of the NAO
                                               !   is tabulated
    real(dp), intent(in)  :: orb_pts(:)        ! Radial part of the NAO
                                               !   (divided by r^l) 
                                               !   in the linear mesh
    type(hfx_options_type), intent(in) :: hfx_opts
                                               ! Options for the hfx calculation
    real(dp), intent(out) :: orb_g(2,ngfs_max) ! Exponents and coefficients
                                               ! Exponents and coefficients
                                               !   of the Gaussians in the 
                                               !   linear combination
                                               !   First  argument: exponent
                                               !   Second argument: coefficient
    integer,  intent(out) :: ierr              ! Error message 

!   Local variables
    integer               :: ngfs_min          ! Minimum number of Gaussians
                                               !   in the linear combination to 
                                               !   fit the NAO
    integer               :: npts_gauss        ! The Gaussians are not strictly
                                               !   confined.
                                               !   Therefore we have to extend
                                               !   the linear mesh from npts
                                               !   to npts_gauss
    integer :: i         ! Counter for the loop in the linear radial mesh
    integer :: npts      ! Number of points in the linear radial mesh where the
                         !   radial part of the numerical atomic orbital 
                         !   is given
    integer :: ngfs      ! Number of gaussians to fit the radial part 
                         !   of the numerical atomic orbital
    real(dp) :: g_err    ! Maximum difference between the radial part of the NAO
                         !   and the linear combination of the Gaussians 
                         !   at any point in the linear mesh  

    real(dp), allocatable :: t(:)     ! Value of the coordinate of the
                                      !   point in the linear mesh
    real(dp), allocatable :: ft(:)    ! Value of the radial part of the
                                      !   NAO in the linear mesh
    real(dp), allocatable :: fval(:)  ! Value of the linear combination
                                      !   of Gaussians in the linear mesh
    real(dp), allocatable :: alpha(:) ! Exponents of the gaussians 
    real(dp), allocatable :: C(:)     ! Coefficients of the gaussians
    real(dp) :: tol_g                 ! Threshold to consider
                                      !   the linear combination of
                                      !   Gaussians converged
    real(dp) :: threshold             ! Threshold to separate the 
                                      !   different exponents
    logical  :: converged             ! Has the fit converged?

! -------------------------------------------------------------------------

! Initialize some variables
    ngfs_min   = hfx_opts%min_num_gaus           ! Minimum number of Gaussians
    tol_g      = hfx_opts%tolerance_gaus         ! Tolerance for convergence
    threshold  = hfx_opts%threshold_exp_gaus     ! Separation between exponents
    npts_gauss = 2 * hfx_opts%npts_fit           ! Number of points in the mesh
    converged  = .false.

    if ( allocated(t)    ) deallocate(t)
    if ( allocated(ft)   ) deallocate(ft)
    if ( allocated(fval) ) deallocate(fval)
    allocate( t(npts_gauss)    )
    allocate( ft(npts_gauss)   )
    allocate( fval(npts_gauss) )

!!   For debugging
!    write(6,'(a,i5)')    &
! &    'nao2gto_gaussfit: ngfs_min   = ', ngfs_min
!    write(6,'(a,i5)')    &
! &    'nao2gto_gaussfit: ngfs_max   = ', ngfs_max
!    write(6,'(a,i5)')    &
! &    'nao2gto_gaussfit: npts_gauss = ', npts_gauss
!    write(6,'(a,f12.5)') &
! &    'nao2gto_gaussfit: tol_g      = ', tol_g
!    write(6,'(a,f12.5)') &
! &    'nao2gto_gaussfit: threshold  = ', threshold
!!   End debugging

    ierr       = 0
    ngfs       = ngfs_min
    npts       = size(rad_pts, 1)
    orb_g(:,:) = 0.0_dp
    g_err      = 1.0_dp

!   Check parameter consistency
    if ( size(orb_pts, 1) /= npts ) then
      ierr = -1
      return
    end if
    if ( 10*npts > 6*npts_gauss ) then
      ierr = -2
      return
    end if
    if ( ngfs > ngfs_max ) then
      ierr = -3
      return
    end if

!   Set actual function to approximate
!   t(:) is a linear grid
!   ft(:) is the radial part of the basis orbital divided by r^l
    t(:)  = 0.0_dp
    ft(:) = 0.0_dp
    t(1:npts) = rad_pts(1:npts)
    ft(1:npts) = orb_pts(1:npts)
!   The linear mesh is expanded to properly treat the tails of the Gaussians
!   The number of points is doubled with respect the linear mesh used
!   to store the radial part of the NAO
    do i = npts+1, npts_gauss
      t(i) = t(npts) + (i-npts)*(t(npts)-t(npts-1))
    end do

!!   For debugging
!    write(6,'(a,i5)')                                  &
! &    'nao2gto_gaussfit: npts = ', npts 
!    write(6,'(a,i5)')                                  &
! &    'nao2gto_gaussfit: npts_gauss = ', npts_gauss 
!    do i = 1, npts
!      write(6,'(a,2f20.12)')                            &
! &      'nao2gto_gaussfit: t, ft = ', t(i), orb_pts(i) 
!    enddo 
!!   End debugging
    

!   Find a minimum suitable fit varying the number of Gaussians
!   We start with a linear combination of ngfs_min gaussians
!   The loop is repeated until the threshold for the difference between
!   the linear combination of Gaussians and the radial part of the NAO
!   is achieved (g_err < tol_g) 
!   or the number of Gaussians that will be required is larger than
!   the maximum value allowed.
    do while ( (g_err > tol_g) .and. (ngfs <= ngfs_max) )

!     If more than one round in the loop is required, 
!     the number of gaussians has been increased
!     by one, so we have to reallocate the coefficient and the exponent arrays
      if ( allocated(alpha) ) deallocate(alpha)
      if ( allocated(C)     ) deallocate(C)
      allocate( alpha(ngfs) )
      allocate( C(ngfs)     )

      call aproxima( t, ft, ngfs, npts_gauss, hfx_opts%npts_fit, threshold, &
 &                   alpha, C, g_err, tol_g )

!     If the deviation is below the tolearance, the fitting is converged  
      if ( g_err < tol_g ) then
        converged = .true.
        exit
      endif 

!!     For debugging
!      write(6,'(a,i5,2f20.12)')                    &
! &      'nao2gto_gaussfit: ngfs, g_err, tol_g = ', &
! &         ngfs, g_err, tol_g 
!!     End debugging

!     If the deviation is above the tolerance, we increase the number
!     of gaussians in the expansion and try again
      ngfs = ngfs + 1

    end do    ! End the do-while loop

    if( .not. converged ) ngfs = ngfs - 1

!   Store the coefficients and the exponents of the fitting in the 
!   proper arrays
    do i = 1, ngfs
      orb_g(1,i) = alpha(i)
      orb_g(2,i) = C(i)
    end do

    if ( g_err < tol_g ) then
      write(6,'(a,/,a,2f12.5)') &
 &      'nao2gto_gaussfit: Fitting of the radial part of NAO to Gaussians converged:', &
 &      'nao2gto_gaussfit: error, tolerance: ', g_err, tol_g
    else
      write(6,'(a,/,a,2f12.5)') &
 &      'nao2gto_gaussfit: Fitting of the radial part of NAO to Gaussians not converged:', &
 &      'nao2gto_gaussfit: error, tolerance: ', g_err, tol_g
!      ierr = 1
    end if

    if ( allocated(alpha) ) deallocate(alpha)
    if ( allocated(C)     ) deallocate(C)
    if ( allocated(t)     ) deallocate(t)
    if ( allocated(ft)    ) deallocate(ft)
    if ( allocated(fval)  ) deallocate(fval)

  end subroutine nao2gto_gaussfit

  ! ***************************************************************************
  ! *** Private routines                                                    ***
  ! ***************************************************************************

! ***************************************************************************
!> \brief Function that computes the error of a fit with given values of 
!!        the coefficients and the exponents of the Gaussians.
!!
!! This is the function that is minimized throughout the process of finding
!! the best set of alpha exponents.
! ***************************************************************************

  function deviation_gauss_nao( t, ft, alpha, C, W, n, m )

    use precision, only: dp

    implicit none

    integer, intent(in) :: n                     ! Number of Gaussians in the 
                                                 !   linear combination
    integer, intent(in) :: m                     ! Number of points in the 
                                                 !   linear mesh, 
                                                 !   extended to deal with
                                                 !   the tail of the Gaussians
    real(dp), dimension(n), intent(in) :: alpha  ! Exponents of the Gaussians
    real(dp), dimension(n), intent(in) :: C      ! Coefficients of the Gaussians
    real(dp), dimension(m), intent(in) :: t      ! Value of the coordinate of 
                                                 !  the point in the linear mesh
                                                 !   extended to deal with the
                                                 !   tail of the Gaussians
    real(dp), dimension(m), intent(in) :: ft     ! Value of the radial part
                                                 !   of the NAO (divided by r^l)
                                                 !   in the linear mesh
    real(dp), dimension(m), intent(in) :: W      ! Weights of the different 
                                                 !   points in the mesh to 
                                                 !   compute the deviation with 
                                                 !   respect to the NAO

!   Output of the function
    real(dp) :: deviation_gauss_nao              ! Measurement of the deviation
                                                 !   between the linear 
                                                 !   combinations of Gaussians
                                                 !   and the radial part of the
                                                 !   NAO in the linear mesh

!   Internal variable
    real(dp), dimension(m)             :: fval
                                                 ! Value of the 
                                                 !   linear combination
                                                 !   of Gaussians in the 
                                                 !   linear mesh


!   We compute the values of our aproximated orbital.
    call eval_lc_gauss( C, alpha, t, n, m, fval )

!   Compute the error with respect to the original orbital.
!   The intrinsic function norm2 calculates the Euclidean vector norm 
!   i.e.
!   On the n-dimensional Euclidean space, \mathbb {R} ^{n}, 
!   the intuitive notion of length of the vector x = (x1, x2, ..., xn) 
!   is captured by the formula
!   {\displaystyle \left\|{\boldsymbol {x}}\right\|_{2}:=
!                         {\sqrt {x_{1}^{2}+\cdots +x_{n}^{2}}}.}
    deviation_gauss_nao = norm2( W(:)*abs(fval(:)-ft(:)) )**2

  end function deviation_gauss_nao


! ***************************************************************************
!> \brief This subroutine finds, for a given number of gaussians, the best
!!        fit to the atomic orbital.
!!
!! Fit: orb(r)=C_1*exp(-alpha_1*r^2)+...+C_n*exp(-alpha_n*r^2)
!! and find the best C and alpha minimizing the least square equation.
! ***************************************************************************
  subroutine aproxima( t, ft, n, m, npts, threshold,    &
 &                     alpha, C, max_diff, tol_g )

    use precision, only: dp

    implicit none

    integer, intent(in) :: n                     ! Number of Gaussians in the 
                                                 !   linear combination
    integer, intent(in) :: m                     ! Number of points in the 
                                                 !   linear mesh, 
                                                 !   extended to deal with
                                                 !   the tail of the Gaussians
    integer, intent(in) :: npts                  ! Number of points in the 
                                                 !   linear mesh
                                                 !   where the radial part of
                                                 !   the NAO is stored
    real(dp), intent(in) :: threshold            ! Threshold to separate the 
                                                 !   different exponents
    real(dp), dimension(m), intent(in)  :: t     ! Value of the coordinate of 
                                                 !   the point in the 
                                                 !   linear mesh
                                                 !   extended to deal with
                                                 !   the tail of the Gaussians
    real(dp), dimension(m), intent(in)  :: ft    ! Value of the radial part
                                                 !   of the NAO (divided by r^l)
                                                 !   in the linear mesh
    real(dp), dimension(n), intent(out) :: alpha ! Exponents of the Gaussians
    real(dp), dimension(n), intent(out) :: C     ! Coefficients of the Gaussians
    real(dp), intent(out) :: max_diff            ! Maximum difference between 
                                                 !   the radial part of the NAO
                                                 !   and the linear combination
                                                 !   of Gaussians in any point
                                                 !   of the linear mesh
    real(dp), intent(in)  :: tol_g               ! Tolerance to consider 
                                                 !   the fit converged

!   Internal variables
    real(dp), dimension(m) :: W                  ! Weights of the different 
                                                 !   points in the mesh to 
                                                 !   compute the deviation with 
                                                 !   respect to the NAO
    real(dp), dimension(m) :: fval               ! Value of the
                                                 !   linear combination
                                                 !   of Gaussians in the
                                                 !   linear mesh
    real(dp), dimension(m) :: errormioabs        ! Difference between the radial
                                                 !   part of the NAO and the 
                                                 !   linear combination of 
                                                 !   Gaussians in every point
                                                 !   of the linear mesh
    real(dp), dimension(n) :: alpha_old          ! Values of the exponents
                                                 !   in the former step
    real(dp) :: change                           ! Maximum change in the 
                                                 !   exponents from one step
                                                 !   to the next during the 
                                                 !   optimization
    real(dp) :: errormioabsmax                   ! Maximum value of errormioabs
    integer  :: i, j
    integer  :: i_steps                          ! Counter for the number of 
                                                 !   steps in the cycle for
                                                 !   the optimization of the 
                                                 !   exponents

!   Vector of weights initialized all values in 1 
!   (zero bias, all data contributes the same)
    W = 1.0_dp

!   Set up the first guess for the coefficients to zero
    C = 0.0_dp

!   Set up the first guess for the alpha exponents
!   Empirically, this choice has shown to behave quite well in all the examples
!   we have tested
    alpha(1) = 0.15_dp
    do i = 2, n
       alpha(i) = dble(i-1)
    enddo

!!   For debugging
!    write(6,'(a,i5)')                                  &
! &    'aproxima: n = ', n
!    write(6,'(a,i5)')                                  &
! &    'aproxima: m = ', m
!    do i = 1, m
!      write(6,'(a,2f20.12)')                           &
! &      'aproxima: t, ft = ', t(i), ft(i)
!    enddo 
!    do i = 1, n
!      write(6,'(a,i5,f12.5)')          &
! &      ' aproxima: i, alpha = ',      & 
! &                  i, alpha(i) 
!    enddo 
!!   End debugging

!   Alternate optimization approach of alpha and C
!   We first perform an optimization of all the exponents and coefficients
    do j = 1, n
       call optimize_j_component( j, t, ft, npts, threshold,    &
 &                                alpha, C, n, m, W )
    enddo

!   Starts the loop to reach convergence. Criteria to stop:
!   1) There is no improvement in the exponents
!   2) Surpass the maximum number of iterations (100)
!   3) The maximum difference between the orbital and the fit is below a
!      tolerance
    change         = 1.0_dp
    errormioabsmax = 1
    i_steps        = 1
    do while ( (change  > 1.d-6)  .and.     &
 &             (i_steps < 100) .and.     &
 &             (errormioabsmax>tol_g))

!      Alternate optimization untill reach convergence
       alpha_old = alpha
       do j = 1, n
          call optimize_j_component(j,t,ft,npts,threshold,alpha,C,n,m,W)
       enddo

       change = maxval( abs(alpha_old - alpha) )
       call eval_lc_gauss( C, alpha, t, n, m, fval )

!      Evaluate the difference between the radial part of the NAO
!      and the linear combination of Gaussians in every point of the linear mesh
       errormioabs    = abs(fval-ft)
!      Find the maximum of this difference
       errormioabsmax = maxval(errormioabs)

!      We update the weights. The more relative error we acummulate at a
!      certain point the more relevant it becomes for the next fit
       W(:)=(W(:)+errormioabs(:)/errormioabsmax)/2.0_dp

       i_steps = i_steps + 1
    end do

    max_diff = maxval(errormioabs)

  end subroutine aproxima

! ***************************************************************************
!> \brief Subroutine that computes the linear combination of gaussians,
!!        whose exponents and coefficients are known at the points of 
!!        a regular linear mesh.
! ***************************************************************************
  subroutine eval_lc_gauss( C, alpha, t, n, m, fval )

    use precision, only: dp

    implicit none

    integer, intent(in) :: n                    ! Number of Gaussians in the 
                                                !   linear combination
    integer, intent(in) :: m                    ! Number of points in the 
                                                !   linear mesh, 
                                                !   extended to deal with
                                                !   the tail of the Gaussians
    real(dp), dimension(n), intent(in) :: alpha ! Exponents of the Gaussians
    real(dp), dimension(n), intent(in) :: C     ! Coefficients of the Gaussians
    real(dp), dimension(m), intent(in) :: t(:)  ! Value of the coordinate of 
                                                !   the point in the linear mesh
                                                !   extended to deal with the
                                                !   tail of the Gaussians
    real(dp), dimension(m), intent(out) :: fval(:)  
                                                ! Value of the 
                                                !   linear combination
                                                !   of Gaussians in the 
                                                !   linear mesh

!   Internal variables
    integer :: i_gauss

    fval(:) = C(1) * exp(-t(:)**2*alpha(1))
    do i_gauss = 2, n
      fval(:) = fval(:) + C(i_gauss) * exp(-t(:)**2*alpha(i_gauss))
    end do

  end subroutine eval_lc_gauss

! ***************************************************************************
!> \brief Subroutine that generetes a vector of n points equally distributed
!!        in the interval (ini,fin)
! ***************************************************************************
  subroutine linspace( ini, fin, n, v )

    use precision, only: dp

    implicit none

    real(dp), intent(in) :: ini              ! First point in the interval
    real(dp), intent(in) :: fin              ! Last  point in the interval
    integer,  intent(in) :: n                ! Number of points in the interval
    real(dp), dimension(n), intent(out) :: v ! Value of the points inside the
                                             !   interval
                                             !   The initial and the final 
                                             !   points are explicitly included

!   Internal variables
    integer  :: i                            ! Counter for the points within
                                             !   the interval
    real(dp) :: esp                          ! Spacing between the points

    esp = (fin -ini) / dble(n-1)

    do i = 1, n
      v(i) = ini + esp * (i-1)
    end do

  end subroutine linspace

! ***************************************************************************
!> \brief Optimizes a function of one variable.
!!
!! To perform the minimization it evaluates 100 times the function in a
!! certain range and selects the minima. Afterwards we restrict to a range
!! near the first minima and repeat the operation. At the end, repeat the
!! process one more time.
! ***************************************************************************
  subroutine optimize_function_one_variable(pos,valor,a,b,alpha,C,W,t,ft,n,m,j)

    use precision, only: dp

    implicit none

    real(dp), intent(in) :: a         ! Initial value of the interval for
                                      !   the minimization of alpha_j
    real(dp), intent(in) :: b         ! Final value of the interval for
                                      !   the minimization of alpha_j

    integer, intent(in)  :: j         ! Index of the exponent alpha_j
                                      !   that is going to be miminzed
    integer, intent(in) :: n          ! Number of Gaussians in the 
                                      !   linear combination
    integer, intent(in) :: m          ! Number of points in the 
                                      !   linear mesh, 
                                      !   extended to deal with
                                      !   the tail of the Gaussians
    real(dp), dimension(m), intent(in)  :: t     ! Value of the coordinate of 
                                                 !   the point in the 
                                                 !   linear mesh
                                                 !   extended to deal with
                                                 !   the tail of the Gaussians
    real(dp), dimension(m), intent(in)  :: ft    ! Value of the radial part
                                                 !   of the NAO (divided by r^l)
                                                 !   in the linear mesh
    real(dp), dimension(m), intent(in) :: W      ! Weights of the different 
                                                 !   points in the linear mesh
                                                 !   to compute the
                                                 !   deviation from the NAO
    real(dp), dimension(n), intent(inout) :: alpha ! Exponents of the Gaussians
    real(dp), dimension(n), intent(inout) :: C   ! Coefficients of the Gaussians
    real(dp), intent(out) :: pos                 ! Value of alpha that minimizes
                                                 !   the error
    real(dp), intent(out) :: valor               ! Value of the error at the 
                                                 !   minimum.
                                                 !   This error refers to the
                                                 !   minimum of the Eucliden
                                                 !   vector norm between the NAO
                                                 !   and the linear combination
                                                 !   of the atomic orbitals,
                                                 !   i.e. summed over all the 
                                                 !   radial points

!   Internal variables
    integer :: indice                 ! Position in the array x where the 
                                      !   minimum is found
    integer :: i
    integer :: k                      ! Number of subdivisions in the interval
                                      !   [a,b]
    real(dp), allocatable, dimension(:) :: x   
                                      ! Value of the points inside the interval
    real(dp), allocatable, dimension(:) :: fx
                                      ! Value of the deviation with respect to
                                      !   the NAO in the linear mesh point
    real(dp)::preci

!   First grid to evaluate the function. 
!   The error function will be evaluated for 100 exponents (alphas)
!   equally spaced inside the interval 
!   [alpha_{i-1}*threshold, alpha_{i+1}/threshold],
!   [a,b] within this subroutine
    k = 100

    if ( allocated(x)  ) deallocate(x)
    allocate(x(k))
    if ( allocated(fx) ) deallocate(fx)
    allocate(fx(k))

    call linspace( a, b, k, x )

    fx    = x
    preci = (b-a)/dble(k)

    do i = 1, k
!     We evaluate the function for all possible 100 values of the target
!     alpha, for that new value of alpha compute the best C coefficients.
      alpha(j)=x(i)
      call solvemincuad( t, ft, alpha, W, n, m, C )

!     Compute the value of the error with that alpha,C election.
      fx(i) = deviation_gauss_nao( t, ft, alpha, C, W, n, m )
    end do

!   Choose the minimum among these values.
    indice = minloc(fx,1)

    if (indice==1) then
      call linspace( x(indice), x(indice)+preci, k, x )
    elseif (indice==k) then
      ! Generate the second grid to evaluate the function.
      call linspace(x(indice)-preci,x(indice),k,x)
    else
      call linspace(x(indice)-preci,x(indice)+preci,k,x)
    endif
    preci=(b-a)/dble(k**2)
    do i=1,k
      alpha(j)=x(i)
!     Second time of the minimization process.
      call solvemincuad(t,ft,alpha,W,n,m,C)
      fx(i)=deviation_gauss_nao(t,ft,alpha,C,W,n,m)
    enddo
    indice=minloc(fx,1)
    if (indice==1) then
      call linspace(x(indice),x(indice)+preci,k,x)
    elseif (indice==k) then
!     Final grid for the evaluation of the function.
      call linspace(x(indice)-preci,x(indice),k,x)
    else
      call linspace(x(indice)-preci,x(indice)+preci,k,x)
    endif
    do i=1,k
      alpha(j)=x(i)
!     Third minimization process.
      call solvemincuad(t,ft,alpha,W,n,m,C)
      fx(i)=deviation_gauss_nao(t,ft,alpha,C,W,n,m)
    enddo
    indice = minloc(fx,1)
    pos    = x(indice)
    valor  = minval(fx)

  end subroutine optimize_function_one_variable


! ***************************************************************************
!> \brief Subroutine that performs the alternate optimization.
!!
!! That is, optimizes each alpha separately with the rest of the alphas
!! frozen.
! ***************************************************************************
  subroutine optimize_j_component( j, t, ft, npts, threshold,  &
 &                                 alpha, C, n, m, W )

    use precision, only: dp

    implicit none

    integer, intent(in) :: j                     ! Index of the exponent to be
                                                 !   optimized
    integer, intent(in) :: n                     ! Number of Gaussians in the 
                                                 !   linear combination
    integer, intent(in) :: m                     ! Number of points in the 
                                                 !   linear mesh, 
                                                 !   extended to deal with
                                                 !   the tail of the Gaussians
    integer, intent(in) :: npts                  ! Number of points in the 
                                                 !   linear mesh
                                                 !   where the radial part of
                                                 !   the NAO is stored
    real(dp), intent(in) :: threshold            ! Threshold to separate the 
                                                 !   different exponents
    real(dp), dimension(m), intent(in)  :: t     ! Value of the coordinate of 
                                                 !   the point in the 
                                                 !   linear mesh
                                                 !   extended to deal with
                                                 !   the tail of the Gaussians
    real(dp), dimension(m), intent(in)  :: ft    ! Value of the radial part
                                                 !   of the NAO (divided by r^l)
                                                 !   in the linear mesh
    real(dp), dimension(n), intent(inout) :: alpha ! Exponents of the Gaussians
    real(dp), dimension(n), intent(inout) :: C   ! Coefficients of the Gaussians
    real(dp), dimension(m), intent(inout) :: W   ! Weights of the different 
                                                 !   points in the linear mesh
                                                 !   to compute the
                                                 !   deviation from the NAO

!   Internal variables
    real(dp) :: rc                               ! Cutoff radius of the 
                                                 !   NAO basis function
    real(dp) :: alpha_min                        ! Global minimal bound for the
                                                 !   first exponent
    real(dp) :: alpha_max                        ! Global maximal bound for the
                                                 !   last exponent
    real(dp) :: x                                ! Optimal value of the j-th
                                                 !   exponent
    real(dp) :: mierror2                         ! Measure of the deviation
                                                 !   between the NAO and the
                                                 !   linear combination of 
                                                 !   Gaussians at the points
                                                 !   of the linear mesh
    integer  :: index_j                          ! Counter for loops on j

!   Setup the cutoff radius of the NAO.
    rc = t(npts)

!   Set up the threshold to separate the different exponents.
!   To minimize the one dimensional function for a given $\alpha_ i$,
!   we evaluate the error function 
!   in $100$ equally spaced  points lying inside the interval 
!   $\left[ \alpha_{i-1} \cdot threshold,\frac{\alpha_{i+1}}{threshold}\right]$ 
!    threshold = 1.4_dp  

!   Set up the global minimal and maximal bounds for the 
!   first and the last exponent
!    alpha_min = 0.5_dp * 2.5_dp**2/(rc**2)
    alpha_min = 0.15_dp
    alpha_max = 0.5_dp * 499.0_dp**2/(9.0_dp*rc**2) 

    if ( j==1 ) then
!     Optimize alpha_1 in the range [alpha_min,alpha_2/threshold]
      call optimize_function_one_variable(x,mierror2,alpha_min, &
             alpha(2)/threshold,alpha,C,W,t,ft,n,m,j)

    elseif ( j==n ) then
!     Optimize alpha_n in the range [alpha_(n-1)*threshold,
!     min(alpha_max,alpha_n +1]
      call optimize_function_one_variable(x,mierror2,         &
             alpha(n-1)*threshold,min(alpha(N)+1,alpha_max), &
             alpha,C,W,t,ft,n,m,j)
    else
!     Optimize alpha_i in the range [alpha_(i-1)*threshold,
!     alpha_(i+1)/Threshold] when i=2,...,n-1
      call optimize_function_one_variable(x,mierror2, &
             alpha(j-1)*threshold,alpha(j+1)/threshold,alpha,C,W,t,ft,n,m,j)
    endif

!   Set up the optimal value of the j-th exponent 
    alpha(j) = x

!   Sort the vector of exponents in increasing order.
    call sort( n, alpha )

!   Recalculate the C coefficients for the new alpha values.
    call solvemincuad( t, ft, alpha, W, n, m, C )

!!   For debugging
!      write(6,'(a,i5)')                                          &
! &      'optimize_j_component: exponent number = ', j
!      write(6,'(a,f12.5)')                                       &
! &      'optimize_j_component: threshold  = ', threshold
!      write(6,'(a,f12.5)')                                       &
! &      'optimize_j_component: rc         = ', rc
!      write(6,'(a,f12.5)')                                       &
! &      'optimize_j_component: alpha_min  = ', alpha_min
!      write(6,'(a,f12.5)')                                       &
! &      'optimize_j_component: alpha_max  = ', alpha_max
!      write(6,'(a,f12.5)')                                       &
! &      'optimize_j_component: x          = ', x
!      write(6,'(a,f12.5)')                                       &
! &      'optimize_j_component: mierror2   = ', mierror2
!      do index_j = 1, j
!        write(6,'(a,i5,f12.5)')                                  &
! &        'optimize_j_component: j, C       = ', index_j, C(index_j)  
!        write(6,'(a,i5,f12.5)')                                  &
! &        'optimize_j_component: j, alpha   = ', index_j, alpha(index_j)  
!      enddo 
!!   End debugging

  end subroutine optimize_j_component

  ! ***************************************************************************
  !> \brief Subroutine that sorts a vector in increasing order
  ! ***************************************************************************
  subroutine sort(n,v)

    use precision, only: dp

    implicit none

    integer,intent(in)::n
    real(dp),dimension(n),intent(inout)::v
    integer::i
    real(dp),dimension(n)::aux
    logical,dimension(n)::mk

    mk=.TRUE.
    aux=v
    do i=1,n
      ! We serch the minimum of the remaining original vector
      aux(i)=minval(v,mk)

      ! We mask the values already considered.
      mk(MINLOC(v,mk))= .FALSE.
    enddo
    v=aux

  end subroutine sort

! ***************************************************************************
!> \brief Subroutine that computes the best set of coefficients, C, 
!!        for a given set of exponents, alpha, in the sense of least squares.
! ***************************************************************************
  subroutine solvemincuad( t, ft, alpha, W, n, m, C ) 

    use precision, only: dp

    implicit none

    integer, intent(in) :: n                     ! Number of Gaussians in the 
                                                 !   linear combination
    integer, intent(in) :: m                     ! Number of points in the 
                                                 !   linear mesh, 
                                                 !   extended to deal with
                                                 !   the tail of the Gaussians
    real(dp), dimension(n), intent(in) :: alpha  ! Exponents of the Gaussians
    real(dp), dimension(n), intent(inout) :: C   ! Coefficients of the Gaussians
    real(dp), dimension(m), intent(in) :: t      ! Value of the coordinate of 
                                                 !  the point in the linear mesh
                                                 !   extended to deal with the
                                                 !   tail of the Gaussians
    real(dp), dimension(m), intent(in) :: ft     ! Value of the radial part
                                                 !   of the NAO (divided by r^l)
                                                 !   in the linear mesh
    real(dp), dimension(m), intent(in) :: W      ! Weights of the different 
                                                 !   points in the mesh to 
                                                 !   compute the deviation with 
                                                 !   respect to the NAO

!   Internal variables
    real(dp), dimension(m,n) :: B
    real(dp), dimension(m*n) :: work
    real(dp), dimension(m)   :: aux
    integer :: info, i 

!   We compute the weighted matrix of coefficients
    B = 0.0_dp
    do i = 1, n
      B(:,i) = W(:) * exp(-t(:)**2*alpha(i))
    enddo

!   We compute the weighted vector of independent terms
    do i = 1, m
       aux(i) = ft(i) * W(i)
    enddo

!   Use LAPACK DGELS subroutine to compute least squares.
    call DGELS('N',m,n,1,B,m,aux,m,work,m*n,info)
    C = aux(1:n)

  end subroutine solvemincuad

end module nao2gto_nonlin
