!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief All kind of helpfull little routines
!> \par History
!>      none
!> \author CJM & JGH
! *****************************************************************************
MODULE nao2gto_util

  
  USE kinds,                           ONLY: dp, sp, &
                                             dp_size
  USE mathconstants,                   ONLY: &
       dfac, fac, fourpi, oorootpi, rootpi, sqrt105, sqrt15, sqrt2, sqrt21, &
       sqrt3, sqrt35, sqrt5, sqrt7, sqrthalf

  USE atm_types,                       ONLY: coset,&
                                             indco,&
                                             nco,&
                                             nso

  use sys, only: die

  IMPLICIT NONE

  PRIVATE

  !MK sqrt* constants moved to mathconstants

  REAL(KIND=dp), PARAMETER ::  s_root1o4pi=0.5_dp*oorootpi
  REAL(KIND=dp), PARAMETER :: root4pi=2.0_dp*rootpi
  REAL(KIND=dp), PARAMETER ::  s_root3o4pi=sqrt3*s_root1o4pi
  REAL(KIND=dp), PARAMETER :: root4pio3=root4pi/sqrt3
  REAL(KIND=dp), PARAMETER :: root4pio5=root4pi/sqrt5
  REAL(KIND=dp), PARAMETER ::  s_root15o4pi=sqrt15*s_root1o4pi
  REAL(KIND=dp), PARAMETER :: root4pio15=root4pi/sqrt15
  REAL(KIND=dp), PARAMETER ::  s_root105o4pi=sqrt105*s_root1o4pi
  REAL(KIND=dp), PARAMETER :: root4pio105=root4pi/sqrt105
  REAL(KIND=dp), PARAMETER ::  s_root1o16pi=0.25_dp*oorootpi
  REAL(KIND=dp), PARAMETER :: root16pi=4.0_dp*rootpi
  REAL(KIND=dp), PARAMETER ::  s_root5o16pi=sqrt5*s_root1o16pi
  REAL(KIND=dp), PARAMETER :: root16pio5=root16pi/sqrt5
  REAL(KIND=dp), PARAMETER ::  s_2root5o16pi=2.0_dp*s_root5o16pi
  REAL(KIND=dp), PARAMETER :: root16pio5o2=root16pio5*0.5_dp
  REAL(KIND=dp), PARAMETER ::  s_3root5o16pi=3.0_dp*s_root5o16pi
  REAL(KIND=dp), PARAMETER :: root16pio5o3=root16pio5/3.0_dp
  REAL(KIND=dp), PARAMETER ::  s_18root5o16pi=18.0_dp*s_root5o16pi
  REAL(KIND=dp), PARAMETER :: root16pio5o18=root16pio5/18.0_dp
  REAL(KIND=dp), PARAMETER ::  s_2root7o16pi=2.0_dp*sqrt7*s_root1o16pi
  REAL(KIND=dp), PARAMETER :: root16pio7o2=root16pi/sqrt7*0.5_dp
  REAL(KIND=dp), PARAMETER ::  s_3root7o16pi=3.0_dp*sqrt7*s_root1o16pi
  REAL(KIND=dp), PARAMETER :: root16pio7o3=root16pi/sqrt7/3.0_dp
  REAL(KIND=dp), PARAMETER ::  s_root15o16pi=sqrt15*s_root1o16pi
  REAL(KIND=dp), PARAMETER :: root16pio15=root16pi/sqrt15
  REAL(KIND=dp), PARAMETER ::  s_3root35o16pi=sqrt5*s_3root7o16pi
  REAL(KIND=dp), PARAMETER :: root16pio35o3=root16pio7o3/sqrt5
  REAL(KIND=dp), PARAMETER ::  s_root105o16pi=0.5_dp*s_root105o4pi
  REAL(KIND=dp), PARAMETER ::  root16pio105=root4pio105*2.0_dp
  REAL(KIND=dp), PARAMETER ::  s_root1o32pi=0.25_dp*sqrthalf*oorootpi
  REAL(KIND=dp), PARAMETER ::  root32pi=root16pi*sqrt2
  REAL(KIND=dp), PARAMETER ::  s_3root5o32pi=3.0_dp*sqrt5*s_root1o32pi
  REAL(KIND=dp), PARAMETER ::  root32pio5o3=root32pi/sqrt5/3.0_dp
  REAL(KIND=dp), PARAMETER ::  s_9root5o32pi=9.0_dp*sqrt5*s_root1o32pi
  REAL(KIND=dp), PARAMETER ::  root32pio5o9=root32pi/sqrt5/9.0_dp
  REAL(KIND=dp), PARAMETER ::  s_12root5o32pi=12.0_dp*sqrt5*s_root1o32pi
  REAL(KIND=dp), PARAMETER ::  root32pio5o12=root32pi/sqrt5/12.0_dp
  REAL(KIND=dp), PARAMETER ::  s_root21o32pi=sqrt21*s_root1o32pi
  REAL(KIND=dp), PARAMETER ::  root32pio21=root32pi/sqrt21
  REAL(KIND=dp), PARAMETER ::  s_4root21o32pi=4.0_dp*s_root21o32pi
  REAL(KIND=dp), PARAMETER ::  root32pio21o4=root32pio21/4.0_dp
  REAL(KIND=dp), PARAMETER ::  s_root35o32pi=sqrt35*s_root1o32pi
  REAL(KIND=dp), PARAMETER ::  root32pio35=root32pi/sqrt35
  REAL(KIND=dp), PARAMETER ::  s_3root35o32pi=3.0_dp*s_root35o32pi
  REAL(KIND=dp), PARAMETER ::  s_9root35o32pi=9.0_dp*s_root35o32pi
  REAL(KIND=dp), PARAMETER ::  s_18root35o32pi=18.0_dp*s_root35o32pi
  REAL(KIND=dp), PARAMETER ::  s_root1o64pi=0.125_dp*oorootpi
  REAL(KIND=dp), PARAMETER ::  s_3root5o64pi=3.0_dp*sqrt5*s_root1o64pi
  REAL(KIND=dp), PARAMETER ::  s_18root5o64pi=18.0_dp*sqrt5*s_root1o64pi
  REAL(KIND=dp), PARAMETER ::  s_root1o256pi=0.0625_dp*oorootpi
  REAL(KIND=dp), PARAMETER ::  s_3root1o256pi=3.0_dp*s_root1o256pi
  REAL(KIND=dp), PARAMETER ::  s_9root1o256pi=9.0_dp*s_root1o256pi
  REAL(KIND=dp), PARAMETER ::  s_18root1o256pi=18.0_dp*s_root1o256pi
  REAL(KIND=dp), PARAMETER ::  s_24root1o256pi=24.0_dp*s_root1o256pi
  REAL(KIND=dp), PARAMETER ::  s_72root1o256pi=72.0_dp*s_root1o256pi
  REAL(KIND=dp), PARAMETER ::  s_3root35o256pi=3.0_dp*sqrt35*s_root1o256pi
  REAL(KIND=dp), PARAMETER ::  s_18root35o256pi=18.0_dp*sqrt35*s_root1o256pi

  ! *** Public subroutines ***

  PUBLIC :: exp_radius,&
            exp_radius_very_extended,&
            gauss_exponent,&
            gaussint_sph
!            trace_r_AxB

CONTAINS

! *****************************************************************************
!> \brief  The exponent of a primitive Gaussian function for a given radius
!>         and threshold is calculated.
!> \param l ...
!> \param radius ...
!> \param threshold ...
!> \param prefactor ...
!> \retval exponent ...
!> \date   07.03.1999
!> \par Variables
!>      - exponent : Exponent of the primitive Gaussian function.
!>      - l        : Angular momentum quantum number l.
!>      - prefactor: Prefactor of the Gaussian function (e.g. a contraction
!>                   coefficient).
!>      - radius   : Calculated radius of the Gaussian function.
!>      - threshold: Threshold for radius.
!> \author MK
!> \version 1.0
! *****************************************************************************
  FUNCTION gauss_exponent(l,radius,threshold,prefactor) RESULT(exponent)
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), INTENT(IN)                :: radius, threshold, prefactor
    REAL(KIND=dp)                            :: exponent

    exponent = 0.0_dp

    IF (radius < 1.0E-6_dp) RETURN
    IF (threshold < 1.0E-12_dp) RETURN

    exponent = LOG(ABS(prefactor)*radius**l/threshold)/radius**2

  END FUNCTION gauss_exponent


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
!> \brief ...
!> \param alpha ...
!> \param l ...
!> \retval gaussint_sph ...
! *****************************************************************************
  FUNCTION gaussint_sph(alpha,l)

    !  calculates the radial integral over a spherical Gaussian
    !  of the form
    !     r**(2+l) * exp(-alpha * r**2)
    !

    REAL(dp), INTENT(IN)                     :: alpha
    INTEGER, INTENT(IN)                      :: l
    REAL(dp)                                 :: gaussint_sph

    IF ((l/2)*2==l) THEN
       !even l:
       gaussint_sph=ROOTPI * 0.5_dp**(l/2+2) * dfac(l+1    )&
            /SQRT(alpha)**(l+3)
    ELSE
       !odd l:
       gaussint_sph=0.5_dp * fac((l+1    )/2) /SQRT(alpha)**(l+3)
    ENDIF

  END FUNCTION gaussint_sph

END MODULE nao2gto_util
