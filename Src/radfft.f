! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!!@LICENSE
!
C *********************************************************************
C MODULE m_radfft
C
C   Public procedures provided:
C subroutine radfft              ! Radial fast Fourier transform
C 
C   Public parameters, variables, and arrays:
C none
C
C   Used module procedures:
C use m_bessph,  only: bessph    ! Spherical Bessel functions
C use m_fft_gpfa,only: fft_gpfa_ez ! Fourier transform
C use alloc,     only: de_alloc  ! Deallocation routines
C use alloc,     only: re_alloc  ! (Re)allocation routines
C
C   Used module parameters:
C use precision, only: dp        ! Double precision real kind
C
C *********************************************************************
C SUBROUTINE RADFFT( L, NR, RMAX, F, G )
C *********************************************************************
C Makes a fast Fourier transform of a radial function.
C If function f is of the form
C   f(r_vec) = F(r_mod) * Ylm(theta,phi)
C where Ylm is a spherical harmonic with l = argument L, and
C argument F contains on input the real function F(r_mod), in a uniform
C radial grid:
C   r_mod = ir * RMAX / NR,  ir = 0,1,...,NR,
C and if g is the 3-dimensional Fourier transform of f:
C   g(k_vec) = 1/(2*pi)**(3/2) *
C              Integral( d3_r_vec * exp(-i * k_vec * r_vec) * f(r_vec) )
C then g has the form
C   g(k_vec) = (-i)**L * G(k_mod) * Ylm(theta,phi)
C where argument G contains on output the real function G(k_mod) in
C a uniform radial grid:
C   k_mod = ik * k_max / NR, ik = 0,1,...,NR,  k_max = NR*pi/RMAX
C Ref: J.M.Soler notes of 16/08/95.
C *************** INPUT ***********************************************
C INTEGER L       : Angular momentum quantum number
C INTEGER NR      : Number of radial intervals.
C                   2*NR must be an acceptable number of points for the
C                   FFT routine used.
C REAL*8  RMAX    : Maximum radius
C REAL*8  F(0:NR) : Function to be tranformed, in a radial mesh
C *************** OUTPUT **********************************************
C REAL*8  G(0:NR) : Fourier transform of F (but see point 5 below)
C *************** UNITS ***********************************************
C Units of RMAX and F are arbitrary.
C Units of k_max and G are related with those of RMAX and F in the
C   obvious way (see above).
C *************** BEHAVIOUR *******************************************
C 1) F and G may be the same physical array, i.e. it is allowed:
C      CALL RADFFT( L, NR, RMAX, F, F )
C 2) It also works in the opposite direction, but then the factor
C    multiplying the output is (+i)**L. Thus, the following two calls
C      CALL RADFFT( L, NR, RMAX, F, G )
C      CALL RADFFT( L, NR, NR*PI/RMAX, G, H )
C    make H = F
C 3) If you will divide the output by q**l, truncation errors may be
C    quite large for small k's if L and NR are large. Therefore, these
C    components are calculated by direct integration rather than FFT.
C    Parameter ERRFFT is the typical truncation error in the FFT, and
C    controls which k's are integrated directly. A good value is 1e-8.
C    If you will not divide by k**l, make ERRFFT=1.e-30.
C 4) The function F is assumed to be zero at and beyond RMAX. The last
C    point F(NR) is not used to find G, except G(0) for L=0 (see 5)
C 5) Because of the 'uncertainty principle', if f(r) is strictly zero
C    for r>RMAX, then g(k) cannot be strictly zero for k>kmax.
C    Therefore G(NR), which should be exactly zero, is used (for L=0)
C    as a 'reminder' term for the integral of G beyond kmax, to ensure 
C    that F(0)=Sum[4*pi*r**2*dr*G(IR)] (this allows to recover F(0)
C    when called again in the inverse direction). Thus, the last value
C    G(NR) should be replaced by zero for any other use. NOTICE: this
C    is commented out in this version!
C *********************************************************************
C Written by J.M.Soler. August 1996.
! Work arrays handling by Rogeli Grima, ca 2009
C *********************************************************************

      MODULE m_radfft

      USE precision, only: dp        ! Double precision real kind
      USE m_bessph,  only: bessph    ! Spherical Bessel functions
      use m_fft_gpfa,only: fft_gpfa_ez     ! 1D fast Fourier transform
      USE alloc,     only: re_alloc, de_alloc
!      USE m_timer,   only: timer_start  ! Start counting CPU time
!      USE m_timer,   only: timer_stop   ! Stop counting CPU time

      implicit none

      PUBLIC :: radfft               ! Radial fast Fourier transform
      PUBLIC :: reset_radfft         ! Deallocates work arrays

      PRIVATE

      ! Work arrays held in module to minimize reallocations
      ! Note that we avoid "automatic" arrays, which may cause stack problems
      real(dp), pointer :: GG(:)
      real(dp), pointer :: FN(:,:)
      real(dp), pointer :: P(:,:,:)
      integer           :: MAXL = -1
      integer           :: MAXNR = -1

      CONTAINS

      SUBROUTINE RADFFT( L, NR, RMAX, F, G )

      IMPLICIT NONE

C Declare argument types and dimensions -----------------------------
      INTEGER, intent(in) :: L       ! Angular momentum of function
      INTEGER, intent(in) :: NR      ! Number of radial points
      real(dp),intent(in) :: RMAX    ! Radius of last point
      real(dp),intent(in) :: F(0:NR) ! Function to Fourier-transform
      real(dp),intent(out):: G(0:NR) ! Fourier transform of F(r)
C -------------------------------------------------------------------

C ERRFFT is the typical truncation error in the FFT routine ---------
      real(dp),   PARAMETER ::    ERRFFT = 1.0E-8_dp
C -------------------------------------------------------------------

C Internal variable types and dimensions ----------------------------
      INTEGER  ::  I, IQ, IR, JR, M, MQ, N, NQ
      real(dp) ::  C, DQ, DR, FR, PI, R, RN, Q, QMAX
!!      real(dp) ::  GG(0:2*NR), FN(2,0:2*NR), P(2,0:L,0:L)
C -------------------------------------------------------------------

C Start time counter ------------------------------------------------
*     CALL TIMER_START( 'RADFFT' )
C
C     Allocate local memory 
      if (MAXL.eq.-1) nullify(P)
      if (L.GT.MAXL) then
        call re_alloc( P, 1, 2, 0, L, 0, L, 'P', 'RADFFT' )
        MAXL=L
      endif
      if (MAXNR.eq.-1) nullify(FN,GG)
      if (NR.GT.MAXNR) then
        call re_alloc( FN, 1, 2, 0, 2*NR, 'FN', 'RADFFT' )
        call re_alloc( GG, 0, 2*NR, 'GG', 'RADFFT' )
        MAXNR=NR
      endif

C Find some constants -----------------------------------------------
      PI = 4.D0 * ATAN( 1.D0 )
      NQ = NR
      DR = RMAX / NR
      DQ = PI / RMAX
      QMAX = NQ * DQ
      C = DR / SQRT( 2.D0*PI )
C -------------------------------------------------------------------

C Set up a complex polynomial such that the spherical Bessel function:
C   j_l(x) = Real( Sum_n( P(n,l) * x**n ) * exp(i*x) ) / x**(l+1)
      P(1,0,0) =  0.D0
      P(2,0,0) = -1.D0
      if (l.gt.0) then
        P(1,0,1) =  0.D0
        P(2,0,1) = -1.D0
        P(1,1,1) = -1.D0
        P(2,1,1) =  0.D0
        if (l.gt.1) then
          DO M = 2,L
          DO N = 0,M
          DO I = 1,2
            P(I,N,M) = 0.D0
            IF (N .LT. M) P(I,N,M) = P(I,N,M) + (2*M-1) * P(I,N,M-1)
            IF (N .GE. 2) P(I,N,M) = P(I,N,M) - P(I,N-2,M-2)
          ENDDO
          ENDDO
          ENDDO
        endif
      endif
C -------------------------------------------------------------------

C Initialize accumulation array -------------------------------------
      DO IQ = 0,NQ
        GG(IQ) = 0.D0
      ENDDO
C -------------------------------------------------------------------

C Iterate on terms of the j_l(q*r) polynomial -----------------------
      DO N = 0,L

C       Set up function to be fast fourier transformed
        FN(1,0) = 0.D0
        FN(2,0) = 0.D0
        DO JR = 1, 2*NR-1

          IF (JR .LT. NR) THEN
            IR = JR
            R = IR * DR
            FR = F(IR)
          ELSEIF (JR .EQ. NR) THEN
            IR = JR
            R = IR * DR
            FR = 0.D0
          ELSE
            IR = 2*NR - JR
            R = - (IR * DR)
            FR = F(IR) * (-1.D0)**L
          ENDIF

C         Find  r**2 * r**n / r**(l+1)
          RN = R**(N-L+1)

          FN(1,JR) = C * FR * RN * P(1,N,L)
          FN(2,JR) = C * FR * RN * P(2,N,L)
        ENDDO

C       Perform one-dimensional complex FFT
!
!       Only the elements from 0 to 2*NR-1 of FN are used.
!       (a total of 2*NR). The fft routine will receive a one-dimensional
!       array of size 2*NR.
!
        CALL fft_gpfa_ez( FN, 2*NR, +1 )

C       Accumulate contribution
        DO IQ = 1,NQ
          Q = IQ * DQ
          GG(IQ) = ( GG(IQ) + FN(1,IQ) ) / Q
        ENDDO

      ENDDO
C -------------------------------------------------------------------

C Special case for Q=0 ---------------------------------------------
      GG(0) = 0.D0
      IF ( L .EQ. 0 ) THEN
        DO IR = 1,NR
          R = IR * DR
          GG(0) = GG(0) + R*R * F(IR)
        ENDDO
        GG(0) = GG(0) * 2.D0 * C
      ENDIF
C -------------------------------------------------------------------

C Direct integration for the smallest Q's ---------------------------
      IF (L.EQ.0) THEN
        MQ = 0
      ELSE
        MQ = NQ * ERRFFT**(1.D0/L)
      ENDIF
      DO IQ = 1,MQ
        Q = IQ * DQ
        GG(IQ) = 0.D0
        DO IR = 1,NR
          R = IR * DR
          GG(IQ) = GG(IQ) + R*R * F(IR) * BESSPH(L,Q*R)
        ENDDO
        GG(IQ) = GG(IQ) * 2.D0 * C
      ENDDO
C -------------------------------------------------------------------

C Special case for Q=QMAX -------------------------------------------
*     IF (L.EQ.0) THEN
*       GSUM = 0.D0
*       DO IQ = 1,NQ-1
*         Q = IQ * DQ
*         GSUM = GSUM + Q*Q * GG(IQ)
*       ENDDO
*       GSUM = GSUM * 4.D0 * PI * DQ
*       GG(NQ) = (2.D0*PI)**1.5D0 * F(0) - GSUM
*       GG(NQ) = GG(NQ) / (4.D0 * PI * DQ * QMAX**2)
*     ENDIF
C -------------------------------------------------------------------

C Copy from local to output array -----------------------------------
      DO IQ = 0,NQ
        G(IQ) = GG(IQ)
      ENDDO
C -------------------------------------------------------------------

C Stop time counter ------------------------------------------------
*     CALL TIMER_STOP( 'RADFFT' )
C -------------------------------------------------------------------

      END SUBROUTINE radfft

      SUBROUTINE RESET_RADFFT( )
      implicit none
      call de_alloc( P, 'P', 'RADFFT' )
      call de_alloc( FN, 'FN', 'RADFFT' )
      call de_alloc( GG, 'GG', 'RADFFT' )
      MAXL  = -1
      MAXNR = -1
      END SUBROUTINE RESET_RADFFT

      END MODULE m_radfft

