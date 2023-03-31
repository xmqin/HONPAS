! *** Module: nao2gto_utils ***

! *****************************************************************************
!> \brief Mathematical utilities for Hartree-Fock exchange calculations.
!!
!! This module provides various utilities with a mathematical focus. It
!! contains arrays and routines to calculate the incomplete Gamma function
!! \f$F_{n}(t)\f$ for multi-center integrals over Cartesian Gaussian
!! functions.
!!
!! \author Matthias Krack
!! \author Juerg Hutter
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 01.1999 Created [Matthias Crack]
!!      - 05.2002 Adapted for use in CP2K [Juerg Hutter]
!!      - 02.2002 Added Gamma functions [Juerg Hutter]
!!      - 05.2004 Restructured and cleaned [Matthias Krack]
!!      - 08.2004 Added Euler constant (gamma) [Juerg Hutter]
!!      - 01.2016 Adapted for use in SIESTA [Xinming Qin]
!!      - 01.2018 Refactored and remodularized for SIESTA 4 [Yann Pouillon]
!!      - 01.2022 New exp_radius funtionals from cp2k for SIESTA [Xinming Qin]
! *****************************************************************************
module nao2gto_utils

  use precision, only: dp, sp
  use nao2gto_data, only: coset
  use sys, only: die

  implicit none

  private

  public :: &
&   exp_radius, &
&   exp_radius_very_extended, &
&   maxfac, fac, ifac, dfac, &
&   deallocate_md_ftable, &
&   gamma0, gamma1, &
&   fgamma_0, &
&   fgamma_ref, &
&   init_md_ftable

  !> Tables stored as constants for efficiency:
  !!
  !!   - factorial function fac
  !!   - inverse factorial function ifac
  !!   - double factorial function dfac
  !!   - gamma functions
  !!   - gamma(n) = gamma0(n) = (n - 1)!
  !!   - gamma(n + 1/2) = gamma1(n) = sqrt(pi)/2^n (2n - 1)!!
  !!
  integer, parameter :: maxfac = 30
  real(dp), parameter, dimension(0:maxfac) :: fac = (/ &
    0.10000000000000000000e+01_dp, 0.10000000000000000000e+01_dp, 0.20000000000000000000e+01_dp, &
    0.60000000000000000000e+01_dp, 0.24000000000000000000e+02_dp, 0.12000000000000000000e+03_dp, &
    0.72000000000000000000e+03_dp, 0.50400000000000000000e+04_dp, 0.40320000000000000000e+05_dp, &
    0.36288000000000000000e+06_dp, 0.36288000000000000000e+07_dp, 0.39916800000000000000e+08_dp, &
    0.47900160000000000000e+09_dp, 0.62270208000000000000e+10_dp, 0.87178291200000000000e+11_dp, &
    0.13076743680000000000e+13_dp, 0.20922789888000000000e+14_dp, 0.35568742809600000000e+15_dp, &
    0.64023737057280000000e+16_dp, 0.12164510040883200000e+18_dp, 0.24329020081766400000e+19_dp, &
    0.51090942171709440000e+20_dp, 0.11240007277776076800e+22_dp, 0.25852016738884976640e+23_dp, &
    0.62044840173323943936e+24_dp, 0.15511210043330985984e+26_dp, 0.40329146112660563558e+27_dp, &
    0.10888869450418352161e+29_dp, 0.30488834461171386050e+30_dp, 0.88417619937397019545e+31_dp, &
    0.26525285981219105864e+33_dp/)
  real(dp), parameter, dimension(0:maxfac) :: ifac = (/ &
    0.10000000000000000000e+01_dp, 0.10000000000000000000e+01_dp, 0.50000000000000000000e+00_dp, &
    0.16666666666666666667e+00_dp, 0.41666666666666666667e-01_dp, 0.83333333333333333333e-02_dp, &
    0.13888888888888888889e-02_dp, 0.19841269841269841270e-03_dp, 0.24801587301587301587e-04_dp, &
    0.27557319223985890653e-05_dp, 0.27557319223985890653e-06_dp, 0.25052108385441718775e-07_dp, &
    0.20876756987868098979e-08_dp, 0.16059043836821614599e-09_dp, 0.11470745597729724714e-10_dp, &
    0.76471637318198164759e-12_dp, 0.47794773323873852974e-13_dp, 0.28114572543455207632e-14_dp, &
    0.15619206968586226462e-15_dp, 0.82206352466243297170e-17_dp, 0.41103176233121648585e-18_dp, &
    0.19572941063391261231e-19_dp, 0.88967913924505732867e-21_dp, 0.38681701706306840377e-22_dp, &
    0.16117375710961183490e-23_dp, 0.64469502843844733962e-25_dp, 0.24795962632247974601e-26_dp, &
    0.91836898637955461484e-28_dp, 0.32798892370698379102e-29_dp, 0.11309962886447716932e-30_dp, &
    0.37699876288159056439e-32_dp/)
  real(dp), parameter, dimension(-1:2*maxfac+1) :: dfac = (/ &
    0.10000000000000000000e+01_dp, 0.10000000000000000000e+01_dp, 0.10000000000000000000e+01_dp, &
    0.20000000000000000000e+01_dp, 0.30000000000000000000e+01_dp, 0.80000000000000000000e+01_dp, &
    0.15000000000000000000e+02_dp, 0.48000000000000000000e+02_dp, 0.10500000000000000000e+03_dp, &
    0.38400000000000000000e+03_dp, 0.94500000000000000000e+03_dp, 0.38400000000000000000e+04_dp, &
    0.10395000000000000000e+05_dp, 0.46080000000000000000e+05_dp, 0.13513500000000000000e+06_dp, &
    0.64512000000000000000e+06_dp, 0.20270250000000000000e+07_dp, 0.10321920000000000000e+08_dp, &
    0.34459425000000000000e+08_dp, 0.18579456000000000000e+09_dp, 0.65472907500000000000e+09_dp, &
    0.37158912000000000000e+10_dp, 0.13749310575000000000e+11_dp, 0.81749606400000000000e+11_dp, &
    0.31623414322500000000e+12_dp, 0.19619905536000000000e+13_dp, 0.79058535806250000000e+13_dp, &
    0.51011754393600000000e+14_dp, 0.21345804667687500000e+15_dp, 0.14283291230208000000e+16_dp, &
    0.61902833536293750000e+16_dp, 0.42849873690624000000e+17_dp, 0.19189878396251062500e+18_dp, &
    0.13711959580999680000e+19_dp, 0.63326598707628506250e+19_dp, 0.46620662575398912000e+20_dp, &
    0.22164309547669977187e+21_dp, 0.16783438527143608320e+22_dp, 0.82007945326378915594e+22_dp, &
    0.63777066403145711616e+23_dp, 0.31983098677287777082e+24_dp, 0.25510826561258284646e+25_dp, &
    0.13113070457687988603e+26_dp, 0.10714547155728479551e+27_dp, 0.56386202968058350995e+27_dp, &
    0.47144007485205310027e+28_dp, 0.25373791335626257948e+29_dp, 0.21686243443194442612e+30_dp, &
    0.11925681927744341235e+31_dp, 0.10409396852733332454e+32_dp, 0.58435841445947272053e+32_dp, &
    0.52046984263666662269e+33_dp, 0.29802279137433108747e+34_dp, 0.27064431817106664380e+35_dp, &
    0.15795207942839547636e+36_dp, 0.14614793181237598765e+37_dp, 0.86873643685617511998e+37_dp, &
    0.81842841814930553085e+38_dp, 0.49517976900801981839e+39_dp, 0.47468848252659720789e+40_dp, &
    0.29215606371473169285e+41_dp, 0.28481308951595832474e+42_dp, 0.17821519886598633264e+43_dp/)
  real(dp), parameter, dimension(0:maxfac) :: gamma0 = (/ &
    0.00000000000000000000E+00_dp, 0.10000000000000000000E+01_dp, 0.10000000000000000000E+01_dp, &
    0.20000000000000000000E+01_dp, 0.60000000000000000000E+01_dp, 0.24000000000000000000E+02_dp, &
    0.12000000000000000000E+03_dp, 0.72000000000000000000E+03_dp, 0.50400000000000000000E+04_dp, &
    0.40320000000000000000E+05_dp, 0.36288000000000000000E+06_dp, 0.36288000000000000000E+07_dp, &
    0.39916800000000000000E+08_dp, 0.47900160000000000000E+09_dp, 0.62270208000000000000E+10_dp, &
    0.87178291200000000000E+11_dp, 0.13076743680000000000E+13_dp, 0.20922789888000000000E+14_dp, &
    0.35568742809600000000E+15_dp, 0.64023737057280000000E+16_dp, 0.12164510040883200000E+18_dp, &
    0.24329020081766400000E+19_dp, 0.51090942171709440000E+20_dp, 0.11240007277776076800E+22_dp, &
    0.25852016738884976640E+23_dp, 0.62044840173323943936E+24_dp, 0.15511210043330985984E+26_dp, &
    0.40329146112660563558E+27_dp, 0.10888869450418352161E+29_dp, 0.30488834461171386050E+30_dp, &
    0.88417619937397019545E+31_dp/)
 real(dp), parameter, dimension(0:maxfac) :: gamma1= (/ &
    0.17724538509055160273E+01_dp, 0.88622692545275801365E+00_dp, 0.13293403881791370205E+01_dp, &
    0.33233509704478425512E+01_dp, 0.11631728396567448929E+02_dp, 0.52342777784553520181E+02_dp, &
    0.28788527781504436100E+03_dp, 0.18712543057977883465E+04_dp, 0.14034407293483412599E+05_dp, &
    0.11929246199460900709E+06_dp, 0.11332783889487855673E+07_dp, 0.11899423083962248457E+08_dp, &
    0.13684336546556585726E+09_dp, 0.17105420683195732157E+10_dp, 0.23092317922314238412E+11_dp, &
    0.33483860987355645697E+12_dp, 0.51899984530401250831E+13_dp, 0.85634974475162063871E+14_dp, &
    0.14986120533153361177E+16_dp, 0.27724322986333718178E+17_dp, 0.54062429823350750447E+18_dp, &
    0.11082798113786903842E+20_dp, 0.23828015944641843260E+21_dp, 0.53613035875444147334E+22_dp, &
    0.12599063430729374624E+24_dp, 0.30867705405286967828E+25_dp, 0.78712648783481767961E+26_dp, &
    0.20858851927622668510E+28_dp, 0.57361842800962338401E+29_dp, 0.16348125198274266444E+31_dp, &
    0.48226969334909086011E+32_dp/)

  ! Internal constants
  real(dp), parameter  :: teps = 1.0e-13_dp

  ! Maximum n value of the tabulated F_n(t) function values
  integer, save :: current_nmax = -1

  ! F_n(t) table
  real(dp), dimension(:, :), allocatable, save :: ftable

contains

! **************************************************************************************************
!> \brief  The radius of a primitive Gaussian function for a given threshold
!>         is calculated.
!>               g(r) = prefactor*r**l*exp(-alpha*r**2) - threshold = 0
!> \param l Angular momentum quantum number l.
!> \param alpha Exponent of the primitive Gaussian function.
!> \param threshold Threshold for function g(r).
!> \param prefactor Prefactor of the Gaussian function (e.g. a contraction
!>                     coefficient).
!> \param epsabs Absolute tolerance (radius)
!> \param epsrel Relative tolerance (radius)
!> \param rlow Optional lower bound radius, e.g., when searching maximum radius
!> \return Calculated radius of the Gaussian function
!> \par Variables
!>        - g       : function g(r)
!>        - prec_exp: single/double precision exponential function
!>        - itermax : Maximum number of iterations
!>        - contract: Golden Section Search (GSS): [0.38, 0.62]
!>                    Equidistant sampling: [0.2, 0.4, 0.6, 0.8]
!>                    Bisection (original approach): [0.5]
!>                    Array size may trigger vectorization
! **************************************************************************************************
   FUNCTION exp_radius(l, alpha, threshold, prefactor, epsabs, epsrel, rlow) RESULT(radius)
      INTEGER, INTENT(IN)                           :: l
      REAL(dp), INTENT(IN)                          :: alpha, threshold, prefactor
      REAL(dp), INTENT(IN), OPTIONAL                :: epsabs, epsrel, rlow
      REAL(dp)                                      :: radius

      INTEGER, PARAMETER                            :: itermax = 5000, prec_exp = sp
      REAL(dp), PARAMETER                           :: contract(*) = (/0.38, 0.62/), &
                                                            mineps = 1.0E-12,next = 2.0, &
                                                            step = 1.0

      INTEGER                                       :: i, j
      REAL(dp)                                      :: a, d, dr, eps, r, rd, t
      REAL(dp), DIMENSION(SIZE(contract))           :: g, s


      IF (l .LT. 0) THEN
         call die("The angular momentum quantum number is negative") 
      END IF

      IF (alpha .EQ. 0.0_dp) THEN
        call die("The Gaussian function exponent is zero")
      ELSE
         a = ABS(alpha)
      END IF

      IF (threshold .NE. 0.0_dp) THEN
         t = ABS(threshold)
      ELSE
         call die("The requested threshold is zero")
      END IF

      radius = 0.0_dp
      IF (PRESENT(rlow)) radius = rlow
      IF (prefactor .EQ. 0.0_dp) RETURN

      ! MAX: facilitate early exit
      r = MAX(SQRT(0.5_dp*REAL(l, dp)/a), radius)
!      write(6,*) "exp radius" l,  a, SQRT(0.5_dp*REAL(l, dp)/a), radius

      d = ABS(prefactor); g(1) = d
      IF (l .NE. 0) THEN
         g(1) = g(1)*EXP(REAL(-a*r*r, KIND=prec_exp))*r**l
      END IF
      ! original approach may return radius=0
      ! with g(r) != g(radius)
      !radius = r
      IF (g(1) .LT. t) RETURN ! early exit

      radius = r*next + step
      DO i = 1, itermax
         g(1) = d*EXP(REAL(-a*radius*radius, KIND=prec_exp))*radius**l
         IF (g(1) .LT. t) EXIT
         r = radius; radius = r*next + step
      END DO

      ! consider absolute and relative accuracy (interval width)
      IF (PRESENT(epsabs)) THEN
         eps = epsabs
      ELSE IF (.NOT. PRESENT(epsrel)) THEN
         eps = mineps
      ELSE
         eps = HUGE(eps)
      END IF
      IF (PRESENT(epsrel)) eps = MIN(eps, epsrel*r)

      dr = 0.0_dp
      DO i = i + 1, itermax
         rd = radius - r
         ! check if finished or no further progress
         IF ((rd .LT. eps) .OR. (rd .EQ. dr)) RETURN
         s = r + rd*contract ! interval contraction
         g = d*EXP(REAL(-a*s*s, KIND=prec_exp))*s**l
         DO j = 1, SIZE(contract)
            IF (g(j) .LT. t) THEN
               radius = s(j) ! interval [r, sj)
               EXIT
            ELSE
               r = s(j) ! interval [sj, radius)
            END IF
         END DO
         dr = rd
      END DO
      IF (i .GE. itermax) THEN
         call die("Maximum number of iterations reached")
      END IF

   END FUNCTION exp_radius

! **************************************************************************************************
!> \brief computes the radius of the Gaussian outside of which it is smaller
!>      than eps
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param pab ...
!> \param o1 ...
!> \param o2 ...
!> \param ra ...
!> \param rb ...
!> \param rp ...
!> \param zetp ...
!> \param eps ...
!> \param prefactor ...
!> \param cutoff ...
!> \param epsabs ...
!> \return ...
!> \par History
!>      03.2007 new version that assumes that the Gaussians origante from spherical
!>              Gaussians
!> \note
!>      can optionally screen by the maximum element of the pab block
! **************************************************************************************************
   FUNCTION exp_radius_very_extended(la_min, la_max, lb_min, lb_max, pab, o1, o2, ra, rb, rp, &
                                     zetp, eps, prefactor, cutoff, epsabs) RESULT(radius)

      INTEGER, INTENT(IN)                           :: la_min, la_max, lb_min, lb_max
      REAL(dp), DIMENSION(:, :), OPTIONAL, POINTER  :: pab
      INTEGER, OPTIONAL                             :: o1, o2
      REAL(dp), INTENT(IN)                          :: ra(3), rb(3), rp(3), zetp, eps, &
                                                       prefactor, cutoff
      REAL(dp), OPTIONAL                            :: epsabs
      REAL(dp)                                      :: radius

      INTEGER                                       :: i, ico, j, jco, la(3), lb(3), lxa, lxb, &
                                                       lya, lyb, lza, lzb
      REAL(dp)                                      :: bini, binj, coef(0:20), epsin_local, &
                                                       polycoef(0:60), prefactor_local, &
                                                       rad_a, rad_b, s1, s2

! get the local prefactor, we'll now use the largest density matrix element of the block to screen

      epsin_local = 1.0E-2_dp
      IF (PRESENT(epsabs)) epsin_local = epsabs

      IF (PRESENT(pab)) THEN
         prefactor_local = cutoff
         DO lxa = 0, la_max
         DO lxb = 0, lb_max
            DO lya = 0, la_max - lxa
            DO lyb = 0, lb_max - lxb
               DO lza = MAX(la_min - lxa - lya, 0), la_max - lxa - lya
               DO lzb = MAX(lb_min - lxb - lyb, 0), lb_max - lxb - lyb
                  la = (/lxa, lya, lza/)
                  lb = (/lxb, lyb, lzb/)
                  ico = coset(lxa, lya, lza)
                  jco = coset(lxb, lyb, lzb)
                  prefactor_local = MAX(ABS(pab(o1 + ico, o2 + jco)), prefactor_local)
               END DO
               END DO
            END DO
            END DO
         END DO
         END DO
         prefactor_local = prefactor*prefactor_local
      ELSE
         prefactor_local = prefactor*MAX(1.0_dp, cutoff)
      END IF

      !
      ! assumes that we can compute the radius for the case where
      ! the Gaussians a and b are both on the z - axis, but at the same
      ! distance as the original a and b
      !
      rad_a = SQRT(SUM((ra - rp)**2))
      rad_b = SQRT(SUM((rb - rp)**2))

      polycoef(0:la_max + lb_max) = 0.0_dp
      DO lxa = 0, la_max
      DO lxb = 0, lb_max
         coef(0:la_max + lb_max) = 0.0_dp
         bini = 1.0_dp
         s1 = 1.0_dp
         DO i = 0, lxa
            binj = 1.0_dp
            s2 = 1.0_dp
            DO j = 0, lxb
               coef(lxa + lxb - i - j) = coef(lxa + lxb - i - j) + bini*binj*s1*s2
               binj = (binj*(lxb - j))/(j + 1)
               s2 = s2*(rad_b)
            END DO
            bini = (bini*(lxa - i))/(i + 1)
            s1 = s1*(rad_a)
         END DO
         DO i = 0, lxa + lxb
            polycoef(i) = MAX(polycoef(i), coef(i))
         END DO
      END DO
      END DO

      polycoef(0:la_max + lb_max) = polycoef(0:la_max + lb_max)*prefactor_local
      radius = 0.0_dp
      DO i = 0, la_max + lb_max
         radius = MAX(radius, exp_radius(i, zetp, eps, polycoef(i), epsin_local, rlow=radius))
      END DO

   END FUNCTION exp_radius_very_extended

! *****************************************************************************
!> \brief   Calculation of the incomplete Gamma function F(t) for multicenter
!!          integrals over Gaussian functions. f returns a vector with all
!!          \f$F_n(t)\f$ values for 0 <= n <= nmax.
!!
!! \author  Matthias Krack
!! \version 1.0
!! \date    08.01.1999,
!! \par Literature
!!      L. E. McMurchie, E. R. Davidson, J. Comp. Phys. 26, 218 (1978)
!!
!! \param nmax: Maximum n value of \f$F_{n}(t)\f$
!! \param t: argument of the incomplete Gamma function
!! \param f: The incomplete Gamma function \f$F_{n}(t)\f$
!!
!! \par History
!!      - 06.1999 Changed from a FUNCTION to a SUBROUTINE [Matthias Krack]
! *****************************************************************************
  subroutine fgamma_0(nmax, t, f)

    use units, only: pi

    implicit none

    ! Arguments
    integer , intent(in) :: nmax
    real(dp), intent(in) :: t
    real(dp), dimension(0:nmax), intent(out) :: f

    integer  :: itab, k, n
    real(dp) :: expt, g, tdelta, tmp, ttab

    ! Calculate f(t)
    if ( t < teps ) then

      ! Special case: t = 0
      do n = 0, nmax
        f(n) = 1.0_dp/real(2*n+1, dp)
      end do

    else if (t <= 12.0_dp) then

      ! 0 < t < 12 => taylor expansion
      tdelta = 0.1_dp

      ! Pretabulation of the f_n(t) function for the taylor series expansion
      if (nmax > current_nmax) then
        call init_md_ftable(nmax)
      end if

      itab = nint(t/tdelta)
      ttab = real(itab, dp)*tdelta
      f(nmax) = ftable(nmax, itab)

      tmp = 1.0_dp
      do k = 1, 6
        tmp = tmp*(ttab-t)
        f(nmax) = f(nmax)+ftable(nmax+k, itab)*tmp*ifac(k)
      end do

      expt = exp(-t)

      ! Use the downward recursion relation to generate
      ! the remaining f_n(t) values
      do n = nmax-1, 0, -1
        f(n) = (2.0_dp*t*f(n+1)+expt)/real(2*n+1, dp)
      end do

    else

      ! t > 12

      tmp = 1.0_dp/t ! reciprocal value

      if (t <= 15.0_dp) then

        ! 12 < t <= 15 => four term polynom expansion

          g = 0.4999489092_dp - 0.2473631686_dp*tmp + &
              0.321180909_dp*tmp**2 - 0.3811559346_dp*tmp**3
          f(0) = 0.5_dp*sqrt(pi*tmp) - g*exp(-t)*tmp

      else if (t <= 18.0_dp) then

!       *** 15 < t <= 18 -> Three term polynom expansion ***

          g = 0.4998436875_dp - 0.24249438_dp*tmp + 0.24642845_dp*tmp**2
          f(0) = 0.5_dp*sqrt(pi*tmp) - g*exp(-t)*tmp

      else if (t <= 24.0_dp) then

!       *** 18 < t <= 24 -> Two term polynom expansion ***

          g = 0.499093162_dp - 0.2152832_dp*tmp
          f(0) = 0.5_dp*sqrt(pi*tmp) - g*exp(-t)*tmp

      else if (t <= 30.0_dp) then

!       *** 24 < t <= 30 -> One term polynom expansion ***

          g = 0.49_dp
          f(0) = 0.5_dp*sqrt(pi*tmp) - g*exp(-t)*tmp

      else

!       *** t > 30 -> Asymptotic formula ***

          f(0) = 0.5_dp*sqrt(pi*tmp)

      end if

      if ( t > real(2*nmax+36, dp) ) then
        expt = 0.0_dp
      else
        expt = exp(-t)
      end if

      ! Use the upward recursion relation to generate
      ! the remaining f_n(t) values
      do n = 1, nmax
        f(n) = (0.5_dp*tmp)*(real(2*n-1, dp)*f(n-1)-expt)
      end do

    end if

  end subroutine fgamma_0

! *****************************************************************************
!> \brief Calculation of the incomplete Gamma function F(t) for multicenter
!!        integrals over Gaussian functions. f returns a vector with all
!!        \f$F_{n}(t)\f$ values for 0 <= n <= nmax.
!!
!! \author  Matthias Krack
!! \version 1.0
!! \date    08.01.1999
!!
!! \par Literature
!!       L. E. McMurchie, E. R. Davidson, J. Comp. Phys. 26, 218 (1978)
!!
!! \param nmax: maximum n value of \f$F_{n}(t)\f$
!! \param t: argument of the incomplete Gamma function
!! \param f: the incomplete Gamma function \f$F_{n}(t)\f$
! *****************************************************************************
  subroutine fgamma_1(nmax, t, f)

    use units, only: pi

    implicit none

    ! Arguments
    integer, intent(in)                :: nmax
    real(dp), dimension(:), intent(in) :: t
    real(dp), dimension(size(t,1),0:nmax), &
&     intent(out) :: f

    ! Local variables
    integer  :: i, itab, k, n
    real(dp) :: expt, g, tdelta, tmp, ttab

    do i=1,size(t, 1)

      !     *** calculate f(t)
      if ( t(i) < teps ) then

        ! Special case: t = 0
        do n=0,nmax
          f(i, n) = 1.0_dp/real(2*n+1, dp)
        end do

      else if ( t(i) <= 12.0_dp ) then

        ! 0 < t < 12 => taylor expansion
        tdelta = 0.1_dp

        ! Pretabulation of the f_n(t) function for the taylor series expansion
        if ( nmax > current_nmax ) then
          call init_md_ftable(nmax)
        end if

        itab = nint(t(i)/tdelta)
        ttab = real(itab, dp)*tdelta

        f(i, nmax) = ftable(nmax, itab)

        tmp = 1.0_dp
        do k=1,6
          tmp = tmp*(ttab-t(i))
          f(i, nmax) = f(i, nmax)+ftable(nmax+k, itab)*tmp*ifac(k)
        end do

        expt = exp(-t(i))

        ! Use the downward recursion relation to generate the remaining
        ! f_n(t) values
        do n=nmax-1,0,-1
          f(i, n) = (2.0_dp*t(i)*f(i, n+1)+expt)/real(2*n+1, dp)
        end do

      else

        ! t > 12
        if ( t(i) <= 15.0_dp ) then

          ! 12 < t <= 15 => four term polynom expansion
          g = 0.4999489092_dp-0.2473631686_dp/t(i)+ &
              0.321180909_dp/t(i)**2-0.3811559346_dp/t(i)**3
          f(i, 0) = 0.5_dp*sqrt(pi/t(i))-g*exp(-t(i))/t(i)

        else if ( t(i) <= 18.0_dp ) then

          ! 15 < t <= 18 => three term polynom expansion
          g = 0.4998436875_dp-0.24249438_dp/t(i)+0.24642845_dp/t(i)**2
          f(i, 0) = 0.5_dp*sqrt(pi/t(i))-g*exp(-t(i))/t(i)

        else if ( t(i) <= 24.0_dp ) then

          ! 18 < t <= 24 => two term polynom expansion
          g = 0.499093162_dp-0.2152832_dp/t(i)
          f(i, 0) = 0.5_dp*sqrt(pi/t(i))-g*exp(-t(i))/t(i)

        else if ( t(i) <= 30.0_dp ) then

          ! 24 < t <= 30 => one term polynom expansion
          g = 0.49_dp
          f(i, 0) = 0.5_dp*sqrt(pi/t(i))-g*exp(-t(i))/t(i)

        else

          ! t > 30 => asymptotic formula
          f(i, 0) = 0.5_dp*sqrt(pi/t(i))

        end if

        if ( t(i) > real(2*nmax+36, dp) ) then
          expt = 0.0_dp
        else
          expt = exp(-t(i))
        end if

        ! Use the upward recursion relation to generate the remaining
        ! f_n(t) values
        do n=1,nmax
          f(i, n) = 0.5_dp*(real(2*n-1, dp)*f(i, n-1)-expt)/t(i)
        end do

       end if

    end do

  end subroutine fgamma_1

! *****************************************************************************
!> \brief   Calculation of the incomplete Gamma function F_n(t) through a
!!          spherical Bessel function.
!!
!! This routine Calculates the incomplete Gamma function F_n(t) using a
!! spherical Bessel function expansion. fgamma_ref returns a vector
!! with all F_n(t) values for 0 <= n <= nmax. For t values greater than 50
!! an asymptotic formula is used.
!!
!! This function is expected to return accurate F_n(t) values for any
!! combination of n and t, but the calculation is slow and therefore
!! the function may only be used for a pretabulation of \f$F_{n}(t)\f$
!! values or for reference calculations.
!!
!! \author  Matthias Krack
!! \version 1.0
!! \date    07.01.1999
!!
!! \par Literature
!!        F. E. Harris, Int. J. Quant. Chem. 23, 1469 (1983)
!!
!! \par Internal variables:
!!      - expt   : Exponential term in the downward recursion of F_n(t).
!!      - factor : Prefactor of the Bessel function expansion.
!!      - p      : Product of the Bessel function quotients.
!!      - r      : Quotients of the Bessel functions.
!!      - sumterm: One term in the sum over products of Bessel functions.
!!
!! \param[in] nmax: maximum n value of F_n(t)
!! \param[in] t: argument of the incomplete Gamma function
!! \return the value of \f$F_{n}(t)\f$
! *****************************************************************************
  function fgamma_ref(nmax, t) result(f)

    use units,     only: pi

    implicit none

    ! Arguments
    integer, intent(in)  :: nmax
    real(dp), intent(in) :: t

    ! Local constants
    integer , parameter :: kmax = 50
    real(dp), parameter :: eps = epsilon(0.0_dp)

    ! Local variables
    integer  :: j, k, n
    real(dp) :: expt, factor, p, sumterm, sumtot, term
    real(dp), dimension(0:nmax)  :: f
    real(dp), dimension(kmax+10) :: r

    ! -------------------------------------------------------------------------

    ! Initialization
    f(:) = 0.0_dp

    if ( t < teps ) then

      ! Special case: t = 0 => analytic expression
      do n = 0, nmax
        f(n) = 1.0_dp/real(2*n+1, dp)
      end do

    else if ( t <= 50.0_dp ) then

      ! Initialize ratios of bessel functions
      r(kmax+10) = 0.0_dp
      do j=kmax+9,1,-1
        r(j) = -t/(real(4*j+2, dp)-t*r(j+1))
      end do
      factor = 2.0_dp*sinh(0.5_dp*t)*exp(-0.5_dp*t)/t

      do n=0,nmax

        ! Initialize iteration
        sumtot = factor/real(2*n+1, dp)
        term = 1.0_dp

        ! Begin the summation and recursion
        do k=1,kmax
          term = term*real(2*n-2*k+1, dp)/real(2*n+2*k+1, dp)

          ! Product of bessel function quotients
          p = 1.0_dp
          do j=1,k
            p = p*r(j)
          end do
          sumterm = factor*term*p*real(2*k+1, dp)/real(2*n+1, dp)

          if ( abs(sumterm) < eps ) then

            ! Iteration converged
            exit

          else if ( k == kmax ) then

            ! No convergence with kmax iterations
            call die("Error: maximum number of iterations reached")

          else

            ! Add the current term to the sum and continue the iteration
            sumtot = sumtot+sumterm

          end if

        end do

        f(n) = sumtot

      end do

    else

      ! Use asymptotic formula for t > 50
      f(0) = 0.5_dp*sqrt(pi/t)

      ! Use the upward recursion relation to generate the remaining
      ! f_n(t) values
      expt = exp(-t)

      do n=1,nmax
        f(n) = 0.5_dp*(real(2*n-1, dp)*f(n-1)-expt)/t
      end do

    end if

  end function fgamma_ref

! *****************************************************************************
! *** Internal routines                                                     ***
! *****************************************************************************

! *****************************************************************************
!> \brief   Build a table of \f$F_{n}(t)\f$ values in the range
!!          tmin <= t <= tmax with a stepsize of tdelta up to n
!!          equal to nmax.
!!
!! \author  Matthias Krack
!! \version 1.0
!! \date    11.01.1999
!!
!! \param nmax: Maximum n value of F_n(t)
!! \param tmin: Minimum t value
!! \param tmax: Maximum t value
!! \param tdelta: Difference between two consecutive t abcisses (step size)
! *****************************************************************************
  subroutine create_md_ftable(nmax, tmin, tmax, tdelta)

    implicit none

    ! Arguments
    integer , intent(in) :: nmax
    real(dp), intent(in) :: tmin, tmax, tdelta

    ! Internal variables
    integer  :: istat, itab, itabmax, itabmin, n
    real(dp) :: t

    ! -------------------------------------------------------------------------

    ! check internal consistency
    if ( current_nmax > -1 ) then
      call die("Error: an incomplete gamma function table is already allocated")
    end if

    ! Check arguments
    if ( nmax < 0 ) then
      call die("Error: a negative n value for the initialization of &
&         the incomplete gamma function is invalid")
    end if


    if ( (tmax < 0.0_dp) .or. &
&        (tmin < 0.0_dp) .or. &
&        (tdelta <= 0.0_dp) .or. &
&        (tmin > tmax) ) then
      call die('Error: invalid arguments')
    end if

    ! Create table
    n = nmax + 6
    itabmin = floor(tmin/tdelta)
    itabmax = ceiling((tmax-tmin)/tdelta)
    if ( allocated(ftable) ) then
      deallocate(ftable)
    endif
    allocate(ftable(0:n,itabmin:itabmax), stat=istat)
    if ( istat /= 0 ) then
      call die('err: memory allocate, ftable')
    end if
    ftable(:,:) = 0.0_dp

    ! Fill table
    do itab = itabmin, itabmax
      t = real(itab, dp)*tdelta
      ftable(0:n, itab) = fgamma_ref(n, t)
    end do

    ! Save initialization status to prevent double init
    current_nmax = nmax

  end subroutine create_md_ftable

! *****************************************************************************
!> \brief   Deallocate the table of F_n(t) values.
!!
!! \author  Matthias Krack
!! \version 1.0
!! \date    24.05.2004
!!
!! \note This routine is idempotent provided that the table is created with
!!       \ref create_md_ftable.
! *****************************************************************************
  subroutine deallocate_md_ftable()

    implicit none

    if ( current_nmax > -1 ) then
      if ( allocated(ftable) ) deallocate(ftable)
      current_nmax = -1
    end if

  end subroutine deallocate_md_ftable

! *****************************************************************************
!> \brief Initializes a table of incomplete Gamma function values
!!
!! This routine initializes a table of \f$F_{n}(t)\f$ values in
!! the range 0 <= t <= 12 with a stepsize of 0.1 up to n equal to
!! nmax for the Taylor series expansion used by McMurchie-Davidson.
!!
!! \author  Matthias Krack
!! \version 1.0
!! \date    10.06.1999
!!
!! \param[in] nmax: maximum n value of \f$F_{n}(t)\f$
! *****************************************************************************
  subroutine init_md_ftable(nmax)

    implicit none

    ! Arguments
    integer, intent(in) :: nmax

    ! -------------------------------------------------------------------------

    if ( nmax < 0 ) then
      call die('Error: negative n value for incomplete gamma function')
    end if

    ! Check if the current initialization is sufficient
    ! Pretabulation of the f_n(t) function  for the taylor series expansion
    if ( nmax > current_nmax ) then
      call create_md_ftable(nmax, 0.0_dp, 12.0_dp, 0.1_dp)
    end if

  end subroutine init_md_ftable

end module nao2gto_utils
