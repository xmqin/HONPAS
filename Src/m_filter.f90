! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996-2008.
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
!******************************************************************************
! MODULE m_filter
! Handles k-filtering of strictly confined radial functions
! Written by J.M.Soler. April 2008
!******************************************************************************
!
!   PUBLIC procedures available from this module:
! function   kcPhi  : Returns planewave cutoff of a radial function
! subroutine filter : Filters out high-k components of a radial function
!
!   PUBLIC parameters, types, and variables available from this module:
! none
!
!******************************************************************************
!! subroutine filter( l, nr, r, f, kc, norm_opt, n_eigen )
!
! Filters out high-k components of a radial function
!
! Input:
!  integer  l     : Angular momentum of the function to filter
!  integer  nr    : Number of radial points
!  real(dp) r(nr) : Radial mesh, in ascending order (in bohr)
!  real(dp) kc    : Cutoff in reciprocal space (in bohr^-1)
!  integer  norm_opt : Option to renormalize (to initial norm) after filtering:
!                        0 => do not renormalize
!                        1 => renormalize f
!                        2 => renormalize f**2
! Input/output:
!  real(DP) f(nr) : Radial function to be filtered
! Output, Optional:
!  integer  n_eigen: Number of eigenvectors of the filtering kernel
!
!******************************************************************************
! function kcPhi( l, nr, r, phi, etol, kexp )
!
! Returns planewave cutoff of a radial function
!
! Input:
!  integer  l        : Angular momentum of the radial function
!  integer  nr       : Number of radial points
!  real(dp) r(nr)    : Radial mesh, in ascending order (bohr)
!  real(dp) phi(nr)  : Radial function, assumed zero beyond r(nr)
!  real(dp) etol     : Kinetic energy allowed above kc (Ry)
!
! Output:
!  real(dp) kcPhi : Planewave cutoff kc (in bohr^-1) such that
!                     etol = Integral_kc^infinity(dk*k^4*phi(k)^2) 
!                          / Integral_0^infinity(dk*k^2*phi(k)^2)
!
!******************************************************************************

MODULE m_filter

implicit none

! All public procedures (there are no public types, parameters, or variables):
PUBLIC:: &
  kcPhi, &! Returns planewave cutoff of a radial function
  filter  ! Filters out high-k components of a radial function

PRIVATE ! Nothing is declared public beyond this point

  integer, parameter :: dp = selected_real_kind(14,100)

! Internal parameters for filter subroutine
  integer, parameter:: nmesh = 128        ! Number of radial integr. points
  integer, parameter:: minj  = 10         ! Min. num. of Bessel functions
  real(dp),parameter:: njkr  = 0.65_dp    ! Num. Bess. funcs. / (kc*rc)
  real(dp),parameter:: emax  = 1.e-2_dp   ! Min. accepted eigval
  real(dp),parameter:: krtol = 1.e-12_dp  ! Tol. for roots of Bess. funcs.
!  real(dp),parameter:: k2mix = 1.e-6_dp  ! Mix weight of k2 in 'Hamiltonian'
  real(dp),parameter:: k2mix = 0.10_dp  ! Mix weight of k2 in 'Hamiltonian'

  real(dp), save :: kc_needed = 0.0_dp

CONTAINS

!------------------------------------------------------------------------------

function kcPhi( l, nr, r, phi, etol )

! Returns planewave cutoff of a radial function

  implicit none

! Arguments
  integer,  intent(in) :: l        ! Angular momentum of function to filter
  integer,  intent(in) :: nr       ! Number of radial points
  real(dp), intent(in) :: r(nr)    ! Radial mesh, in ascending order (bohr)
  real(dp), intent(in) :: phi(nr)  ! Radial function to be filtered
  real(dp), intent(in) :: etol     ! Kinetic energy allowed above kc (Ry)
  real(dp)             :: kcPhi    ! Plane wavevector cutoff (bohr^-1)

! Internal parameters
  integer, parameter:: nmesh    = 1024    ! Number of radial integr. points
  integer, parameter:: kexp     = 2       ! k exponent in integral
  real(dp),parameter:: rmaxByrc = 3._dp   ! rmax/rc

! Internal variables and arrays
  integer :: ik, ik1, ik2, ir, nrc
  real(dp):: de, dk, dr, e, f(0:nmesh), kmax, kmesh(0:nmesh), k1, k2, &
             norm, pi, rc, rmax, rmesh(0:nmesh)

! Trap bad arguments  
  if (etol<=0._dp .or. etol>=1._dp) STOP 'kcPhi: ERROR: bad tol'
  if (kexp/=0 .and. kexp/=2)        STOP 'kcPhi: ERROR: bad kexp'

! Find cutoffs for integration mesh
  pi = acos(-1._dp)
  rc = r(nr)
  rmax = rc * rmaxByrc
  dr = rmax / nmesh
  nrc = nint( rc / dr )
  dr = rc / nrc   ! Make rc one of the mesh points
  rmax = nmesh * dr
  kmax = pi / dr
  dk = kmax / nmesh

! Find integration meshes
  forall(ir=0:nmesh) rmesh(ir) = ir*dr
  forall(ik=0:nmesh) kmesh(ik) = ik*dk

! Find function at integration points
  f(0:nrc) = interpolation( phi, r, rmesh(0:nrc) )
  f(nrc+1:nmesh) = 0

! Find Fourier transform of function
  call radfft( l, nmesh, rmax, f(0:nmesh), f(0:nmesh) )

! Normalize function
  norm = sum( kmesh(:)**2 * f(:)**2 ) * dk
  if (norm==0._dp) STOP 'kcPhi: ERROR: phi=0'
  f(:) = f(:) / sqrt(norm)

! Find cutoff kc such that tol = Integral_kc^infinity dk * k^(2*kexp) * f(k)^2
  e = 0
  do ik2 = nmesh,1,-1
    ik1 = ik2 - 1
    k1 = kmesh(ik1)
    k2 = kmesh(ik2)
    dk = k2 - k1
    de = (k1**(2+kexp) * f(ik1)**2 + k2**(2+kexp) * f(ik2)**2) * dk/2
    if (e+de > etol) then
      kcPhi = k2 - dk * (etol-e)/de
      exit
    end if
    e = e + de
  end do

end function kcPhi

!------------------------------------------------------------------------------

subroutine filter( l, nr, r, f, kc, norm_opt, n_eigen)

! Filters out high-k components of a radial function
! Ref: J.M.Soler notes of 17/10/2006.
! Written by J.M.Soler. Apr.2007

  implicit none

! Arguments
  integer,  intent(in)    :: l        ! Angular momentum of function to filter
  integer,  intent(in)    :: nr       ! Number of radial points
  real(dp), intent(in)    :: kc       ! Cutoff in reciprocal space
  real(dp), intent(in)    :: r(nr)    ! Radial mesh, in ascending order
  real(dp), intent(inout) :: f(nr)    ! Radial function to be filtered
  integer,  intent(in)    :: norm_opt ! Option to renormalize after filtering
  integer,  intent(out), optional   :: n_eigen        ! Number of eigenvectors of the filtering kernel

! Internal variables and arrays
  integer :: i, ierror, ij, ij1, ij2, ik, ir, n, nj, m
  real(dp):: cij, dk, dkr, dr, f0norm, fg, fnorm, jl0, jl1, jl2, &
             k, k4j1(0:nmesh), kmesh(0:nmesh), kr, kr0, kr1, kr2, &
             pi, rc, rmesh(0:nmesh)
!             k, k4j1(0:nmesh), kmesh(0:2*nmesh), kr, kr0, kr1, kr2, &
!             pi, rc, rmesh(0:2*nmesh)
  integer, allocatable:: indx(:)
  real(dp),allocatable:: aux(:), c(:,:), e(:), fm(:), g(:,:), gr(:,:), &
                         h(:,:), jl(:,:), jlk(:,:), jlnorm(:), jlr(:), &
                         kj(:), s(:,:)

! Fix number of basis Bessel functions by empirical rule
  pi = acos(-1._dp)
  rc = r(nr)
  n = minj + nint(njkr*kc*rc)

! Find integration meshes
  dr = rc / nmesh
  dk = kc / nmesh
  do i = 0,nmesh
!  do i = 0,2*nmesh
    rmesh(i) = i*dr
    kmesh(i) = i*dk
  end do

! Allocate arrays
  allocate( aux(n), c(n,n), e(n), fm(0:nmesh), &
            g(0:nmesh,n), gr(nr,n), h(n,n), indx(n), &
            jl(0:nmesh,n), jlk(0:nmesh,n), jlnorm(n), jlr(n), kj(n), s(n,n) )
!            jl(0:nmesh,n), jlk(0:2*nmesh,n), jlnorm(n), jlr(n), kj(n), s(n,n) )

! Find roots of spherical Bessel functions 
  dkr = pi/2         ! Safe enough not to miss any root
  kr1 = dkr/2
  jl1 = bessph(l,kr1)
  nj = 0
  sampling: do
    kr2 = kr1 + dkr
    jl2 = bessph(l,kr2)
    if (jl1*jl2 < 0._dp) then  ! Root bracketted
      bisection: do
        kr0 = (kr1+kr2)/2
        if (kr2-kr1 < krtol) exit bisection
        jl0 = bessph(l,kr0)
        if (jl0*jl1 > 0._dp) then
          kr1 = kr0
          jl1 = jl0
        else
          kr2 = kr0
          jl2 = jl0
        end if
      end do bisection
      nj = nj + 1
      kj(nj) = kr0 / rc
      if (nj==n) exit sampling
    end if
    kr1 = kr2
    jl1 = jl2
  end do sampling

! Find normalized spher. Bessel funcs. at integration mesh points
  jl(:,:) = 0
  do ij = 1,n
    jlnorm(ij) = 0
    do ir = 0,nmesh
      kr = kj(ij)*rmesh(ir)
      jl(ir,ij) = bessph( l, kr )
      jlnorm(ij) = jlnorm(ij) + rmesh(ir)**2 * jl(ir,ij)**2 * dr
    end do ! ir
    jlnorm(ij) = sqrt(jlnorm(ij))
    jl(0:nmesh,ij) = jl(0:nmesh,ij) / jlnorm(ij)
  end do ! ij

! Find Fourier transform of cutted Bessel functions
! Ref: J.M.Soler notes of 16/04/2008
  do ij = 1,n
    cij = (-1._dp)**ij * 2*sqrt(rc/pi) * kj(ij)
    do ik = 0,nmesh
!    do ik = 0,2*nmesh
      k = kmesh(ik)
      if (abs(k-kj(ij)) > 100*krtol) then
        jlk(ik,ij) = cij * bessph(l,k*rc) / (k**2 - kj(ij)**2)
      else ! denominator too small
        ! Use j'(l,x)=-j(l+1,x)+l*jl(l,x)/x  and  j(l,kj(ij)*rc)=0
        k = (k + kj(ij)) / 2
        jlk(ik,ij) = cij * rc/(2*k) * &
                   ( -bessph(l+1,k*rc) + l*bessph(l,k*rc)/(k*rc) )
      end if
    end do
  end do ! ij

! Find 'Hamiltonian' (integral between kc and infinity of k**2*jl1*jl2)
  h(:,:) = 0
  do ij1 = 1,n
    h(ij1,ij1) = kj(ij1)**2 ! Integral from zero to infinity of diagonal terms
    k4j1(0:nmesh) = kmesh(0:nmesh)**4 * jlk(0:nmesh,ij1)
    do ij2 = 1,ij1
      ! Subtract integral up to kc using trapezoidal rule
!      h(ij1,ij2) = h(ij1,ij2) - sum( k4j1(:)*jlk(:,ij2) )*dk + &
!                   k4j1(nmesh)*jlk(nmesh,ij2)*dk/2
      ! Subtract integral up to kc using Simpson's rule
      h(ij1,ij2) = h(ij1,ij2) - &
        ( 4*sum(k4j1(1:nmesh-1:2)*jlk(1:nmesh-1:2,ij2)) + &
          2*sum(k4j1(2:nmesh-2:2)*jlk(2:nmesh-2:2,ij2)) + &
                k4j1(nmesh)      *jlk(nmesh,ij2)          ) * dk/3
      ! Subtract integral up to kc using higher-order routine
!      h(ij1,ij2) = h(ij1,ij2) - &
!        integral( nmesh+1, k4j1(:)*jlk(:,ij2), x=kmesh(:) )
      ! Symmetrize
      h(ij2,ij1) = h(ij1,ij2)
      write(33,*) ij1,ij2,h(ij1,ij2)
    end do ! ij2
  end do ! ij1

! Add some mixing weight of total kinetic energy to 'Hamiltonian'
  do ij = 1,n
    h(ij,ij) = (1-k2mix) * h(ij,ij) + k2mix * kj(ij)**2
    
  end do

! Initialize unity overlap matrix to call rdiag
  s = 0
  forall(i=1:n) s(i,i) = 1

! Diagonalize 'Hamiltonian'
  call filter_rdiag( h, s, n, n, n, e, c, n, 0, ierror )

! Subtract mixing weight of total kinetic energy from 'Hamiltonian'
!  do ij = 1,n
!    h(ij,ij) = ( h(ij,ij) - k2mix * kj(ij)**2 ) / (1-k2mix)
!  end do

! Print eigenvalues for debugging
!  print'(a,i4,e15.6)', ('filter: i,eigval/k2 =',i,e(i)/kj(i)**2-k2mix,i=1,n)

! Write leaked kinetic energy
!  open( unit=1, file='eigval.out' )
!  do i = 1,n
!    write(1,'(4f12.6)') kj(i)/kc, h(i,i)/kj(i)**2, &
!      sum(c(:,i)*matmul(h(:,:),c(:,i))) / kj(i)**2, e(i)/kj(i)**2
!  end do
!  close( unit=1 )

! Find how many eigenvalues are within filter tolerance
  m = 0
  do i = 1,n
    if (e(i)/kj(i)**2-k2mix > emax) exit
    m = i
  end do
  write(6,'(a)')       'Filter:    Number of eigenfunctions  of the'
  write(6,'(a,i3,i3)') '           filtering kernel (total, used)=', n, m

  if (present(n_eigen)) then
     n_eigen = m
  endif

! Find eigenvectors at integration mesh points
  g(0:nmesh,1:n) = matmul( jl(0:nmesh,1:n), c(1:n,1:n) )

! Write eigenvectors and Bessel functions
!  open( unit=1, file='bessel_k.out' )
!  do ik = 0,2*nmesh
!    write(1,'(20f12.6)') kmesh(ik)/kc, (kmesh(ik)**2*jlk(ik,i)**2,i=1,15)
!  end do
!  close( unit=1 )
!  open( unit=1, file='eigvec_k.out' )
!  do ik = 0,2*nmesh
!    write(1,'(20f12.6)') kmesh(ik)/kc, &
!                         (kmesh(ik)**2*sum(jlk(ik,:)*c(:,i))**2,i=1,15)
!  end do
!  close( unit=1 )

! Find function to be filtered at integration mesh points
  fm(0:nmesh) = interpolation( f, r, rmesh(0:nmesh) )

! Find eigenvectors at input mesh points
  do ir = 1,nr
    do ij = 1,n
      kr = kj(ij) * r(ir)
      jlr(ij) = bessph( l, kr ) / jlnorm(ij)
    end do
    gr(ir,:) = matmul( jlr(:), c(:,:) )
  end do

! Filter radial function (times r)
  f = 0
  do i = 1,n
    fg = sum( rmesh(0:nmesh)**2 * fm(:) * g(:,i) ) * dr
    if (i<=m) f(1:nr) = f(1:nr) + fg * gr(1:nr,i)
  end do

! Renormalize filtered function
  if (norm_opt==1) then
    f0norm = sum( matmul(rmesh(0:nmesh)**2*fm,g(:,1:n))*dr )
    fnorm  = sum( matmul(rmesh(0:nmesh)**2*fm,g(:,1:m))*dr )
    f = f * f0norm/fnorm
  else if (norm_opt==2) then
    f0norm = sum( (matmul(rmesh(0:nmesh)**2*fm,g(:,1:n))*dr)**2 )
    fnorm  = sum( (matmul(rmesh(0:nmesh)**2*fm,g(:,1:m))*dr)**2 )
    f = f * sqrt(f0norm/fnorm)
  else if (norm_opt/=0) then
     STOP 'filter: ERROR: invalid value of norm_opt'
  end if

! Deallocate arrays
  deallocate( aux, c, e, fm, g, gr, h, indx, jl, jlk, jlnorm, jlr, kj, s )

end subroutine filter

!------------------------------------------------------------------------------

  function interpolation( yold, xold, xnew ) result(ynew)

    ! Interpolates function yold(xold) at new points xnew
    ! Assumes that xold and xnew are in increasing order

    implicit none
    real(dp), intent(in) :: yold(:), xold(:), xnew(:)
    real(dp) :: ynew(size(xnew))

    integer :: i1, i2, inew, iold, nnew, nold

    nold = size(xold)
    nnew = size(xnew)

    iold = 1
    do inew = 1,nnew
      ! Find iold such that xold(iold) < xnew(inew) < xold(iold+1)
      do
        if (iold+1==nold .or. xold(iold+1) > xnew(inew)) exit
        iold = iold + 1
      end do
      ! Four-point Lagrange interpolation (3-point at extreme intervals)
      i1 = max(iold-1,1)
      i2 = min(iold+2,nold)
      ynew(inew) = lagrange( yold(i1:i2), xold(i1:i2), xnew(inew) )
    end do

  end function interpolation

!------------------------------------------------------------------------------

  function lagrange( yold, xold, xnew ) result(ynew)

    ! Lagrange interpolation of function yold(xold) at point xnew
    ! Based on routine polint of Numerical Recipes

    implicit none
    real(dp), intent(in) :: yold(:), xold(:), xnew
    real(dp) :: ynew

    integer :: i, i0, m, n
    real(dp):: c(size(xold)), d(size(xold)), den, dy, ho, hp, w

    n = size(xold)
    c = yold
    d = yold
    i0 = n/2
    ynew = yold(i0)
    i0 = i0 - 1
    do m = 1,n-1
      do i=1,n-m
        ho = xold(i) - xnew
        hp = xold(i+m) - xnew
        w = c(i+1) - d(i)
        den = ho - hp
        if (den==0._dp) STOP 'filter:lagrange: ERROR: den=0'
        d(i) = hp * w / den
        c(i) = ho * w / den
      end do
      if (2*i0 < n-m) then
        dy = c(i0+1)
      else
        dy = d(i0)
        i0 = i0 - 1
      endif
      ynew = ynew + dy
    end do

  end function lagrange


  ! ***************************************************************************
! subroutine rdiag(H,S,n,nm,nml,w,Z,neigvec,iscf,ierror)
!
! Simple replacement to subroutine rdiag of siesta
! J.M.Soler, April 2008
! ************************** INPUT ******************************************
! real*8 H(nml,nm)                 : Symmetric H matrix
! real*8 S(nml,nm)                 : Symmetric S matrix, ignored in this version
! integer n                        : Order of the generalized  system
! integer nm                       : Right hand dimension of H and S matrices
! integer nml                      : Left hand dimension of H and S matrices
!                                    which is greater than or equal to nm
! integer neigvec                  : No. of eigenvectors to calculate
! integer iscf                     : SCF cycle, ignored in this version
! ************************** OUTPUT *****************************************
! real*8 w(nml)                    : Eigenvalues
! real*8 Z(nml,nm)                 : Eigenvectors
! integer ierror                   : Flag indicating success code for routine
!                                  :  0 = success
!                                  : -1 = repeat call as memory is increased
!                                  :  1 = fatal error
! ***************************************************************************

subroutine filter_rdiag(H,S,n,nm,nml,w,Z,neigvec,iscf,ierror)

  implicit none

! Passed variables
  integer  :: ierror
  integer  :: iscf
  integer  :: n
  integer  :: neigvec
  integer  :: nm
  integer  :: nml
  real(dp) :: H(nml,nm)
  real(dp) :: S(nml,nm)
  real(dp) :: w(nml)
  real(dp) :: Z(nml,nm)

! Internal variables and arrays
  integer :: indx(n)
  real(dp):: aux(n), c(n,n), e(n)

! Diagonalize Hamiltonian
  c(1:n,1:n) = H(1:n,1:n)
  call tred2( c, n, n, e, aux )
  call tqli( e, aux, n, n, c )

! Order eigenvalues and eigenvectors by increasing eigval
  call ordix( e, 1, n, indx )
  call order( e, 1, n, indx )
  call order( c, n, n, indx )

! Copy eigenvectors and eigenvalues to output arrays
  w = 0
  Z = 0
  w(1:neigvec) = e(1:neigvec)
  Z(1:n,1:neigvec) = c(1:n,1:neigvec)
  ierror = 0

end subroutine filter_rdiag

      FUNCTION BESSPH (L,X) 
!
!  RETURNS THE SPHERICAL BESSEL FUNCTION JL(X).
!  REF: ABRAMOWITZ AND STEGUN, FORMULAS 10.1.2 AND 10.1.19
!  WRITTEN BY J.SOLER (JSOLER AT EMDUAM11). NOV/89.
!

      integer, intent(in)  :: L
      real(dp), intent(in) :: X
      real(dp)             :: BESSPH

      integer, parameter :: nterms = 100
      real(dp), parameter ::   zero = 0.0_dp, one = 1.0_dp, tiny=1.0e-15_dp

      integer :: i, n
      real(dp) :: switch, term, x2, sum, y, fnm1, fn, fnp1

      SWITCH=MAX(1,2*L-1)
      IF (ABS(X).LT.SWITCH) THEN
!       USE POWER SERIES
         TERM=ONE
         DO 10 I=1,L
            TERM=TERM*X/(2*I+1)
   10    CONTINUE
         X2=X*X
         SUM=ZERO
         DO 20 I=1,NTERMS
            SUM=SUM+TERM
            TERM=(-TERM)*X2/(2*I*(2*I+2*L+1))
            IF (ABS(TERM).LT.TINY) GO TO 30
   20    CONTINUE
         STOP 'BESSPH: SERIES HAS NOT CONVERGED'
   30    BESSPH=SUM
      ELSE
!       USE EXPLICIT EXPRESSIONS OR RECURRENCE RELATION
         IF (L.EQ.0) THEN
            BESSPH=SIN(X)/X
         ELSEIF (L.EQ.1) THEN
            BESSPH=(SIN(X)/X-COS(X))/X
         ELSE
            Y=ONE/X
            FNM1=SIN(X)*Y
            FN=(FNM1-COS(X))*Y
            DO 40 N=1,L-1
               FNP1=(2*N+1)*Y*FN-FNM1
               FNM1=FN
               FN=FNP1
   40       CONTINUE
            BESSPH=FN
         ENDIF
      ENDIF
      END function bessph

!
      SUBROUTINE RADFFT( L, NR, RMAX, F, G )
! *********************************************************************
! Makes a fast Fourier transform of a radial function.
! If function f is of the form
!   f(r_vec) = F(r_mod) * Ylm(theta,phi)
! where Ylm is a spherical harmonic with l = argument L, and
! argument F contains on input the real function F(r_mod), in a uniform
! radial grid:
!   r_mod = ir * RMAX / NR,  ir = 0,1,...,NR,
! and if g is the 3-dimensional Fourier transform of f:
!   g(k_vec) = 1/(2*pi)**(3/2) *
!              Integral( d3_r_vec * exp(-i * k_vec * r_vec) * f(r_vec) )
! then g has the form
!   g(k_vec) = (-i)**L * G(k_mod) * Ylm(theta,phi)
! where argument G contains on output the real function G(k_mod) in
! a uniform radial grid:
!   k_mod = ik * k_max / NR, ik = 0,1,...,NR,  k_max = NR*pi/RMAX
! Ref: J.M.Soler notes of 16/08/95.
! *************** INPUT ***********************************************
! INTEGER L       : Angular momentum quantum number
! INTEGER NR      : Number of radial intervals.
!                   2*NR must be an acceptable number of points for the
!                   FFT routine used.
! REAL*8  RMAX    : Maximum radius
! REAL*8  F(0:NR) : Function to be tranformed, in a radial mesh
! *************** OUTPUT **********************************************
! REAL*8  G(0:NR) : Fourier transform of F (but see point 5 below)
! *************** UNITS ***********************************************
! Units of RMAX and F are arbitrary.
! Units of k_max and G are related with those of RMAX and F in the
!   obvious way (see above).
! *************** BEHAVIOUR *******************************************
! 1) F and G may be the same physical array, i.e. it is allowed:
!      CALL RADFFT( L, NR, RMAX, F, F )
! 2) It also works in the opposite direction, but then the factor
!    multiplying the output is (+i)**L. Thus, the following two calls
!      CALL RADFFT( L, NR, RMAX, F, G )
!      CALL RADFFT( L, NR, NR*PI/RMAX, G, H )
!    make H = F
! 3) If you will divide the output by q**l, truncation errors may be
!    quite large for small k's if L and NR are large. Therefore, these
!    components are calculated by direct integration rather than FFT.
!    Parameter ERRFFT is the typical truncation error in the FFT, and
!    controls which k's are integrated directly. A good value is 1e-8.
!    If you will not divide by k**l, make ERRFFT=1.e-30.
! 4) The function F is assumed to be zero at and beyond RMAX. The last
!    point F(NR) is not used to find G, except G(0) for L=0 (see 5)
! 5) Because of the 'uncertainty principle', if f(r) is strictly zero
!    for r>RMAX, then g(k) cannot be strictly zero for k>kmax.
!    Therefore G(NR), which should be exactly zero, is used (for L=0)
!    as a 'reminder' term for the integral of G beyond kmax, to ensure 
!    that F(0)=Sum[4*pi*r**2*dr*G(IR)] (this allows to recover F(0)
!    when called again in the inverse direction). Thus, the last value
!    G(NR) should be replaced by zero for any other use. NOTICE: this
!    is commented out in this version!
! *********************************************************************
! Written by J.M.Soler. August 1996.
! *********************************************************************

      IMPLICIT NONE
      integer, parameter :: dp = selected_real_kind(14,100)

! Declare argument types and dimensions -----------------------------
      INTEGER           L, NR
      real(dp)          F(0:), G(0:), RMAX
! -------------------------------------------------------------------

! ERRFFT is the typical truncation error in the FFT routine ---------
      real(dp),   PARAMETER ::    ERRFFT = 1.0E-8_dp
! -------------------------------------------------------------------

! Internal variable types and dimensions ----------------------------

      INTEGER  ::  I, IQ, IR, JR, M, MQ, N, NQ
      real(dp) ::  C, DQ, DR, FR, PI, R, RN, Q, QMAX

      real(dp), allocatable      ::  GG(:), FN(:,:), P(:,:,:)

! -------------------------------------------------------------------

! Allocate local memory ---------------------------------------------

      allocate(P(1:2,0:L,0:L))
      allocate(FN(1:2,0:2*NR))
      allocate(GG(0:2*NR))

! Find some constants -----------------------------------------------
      PI = 4.D0 * ATAN( 1.D0 )
      NQ = NR
      DR = RMAX / NR
      DQ = PI / RMAX
      QMAX = NQ * DQ
      C = DR / SQRT( 2.D0*PI )
! -------------------------------------------------------------------


! Set up a complex polynomial such that the spherical Bessel function:
!   j_l(x) = Real( Sum_n( P(n,l) * x**n ) * exp(i*x) ) / x**(l+1)
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
! -------------------------------------------------------------------

! Initialize accumulation array -------------------------------------
      DO IQ = 0,NQ
        GG(IQ) = 0.D0
      ENDDO
! -------------------------------------------------------------------

! Iterate on terms of the j_l(q*r) polynomial -----------------------
      DO N = 0,L

!       Set up function to be fast fourier transformed
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

!         Find  r**2 * r**n / r**(l+1)
          RN = R**(N-L+1)

          FN(1,JR) = C * FR * RN * P(1,N,L)
          FN(2,JR) = C * FR * RN * P(2,N,L)
        ENDDO

!       Perform one-dimensional complex FFT
!
!       Only the elements from 0 to 2*NR-1 of FN are used.
!       (a total of 2*NR). Four1 will receive a one-dimensional
!       array of size 2*NR.
!
        CALL FOUR1( FN, 2*NR, +1 )

!       Accumulate contribution
        DO IQ = 1,NQ
          Q = IQ * DQ
          GG(IQ) = ( GG(IQ) + FN(1,IQ) ) / Q
        ENDDO

      ENDDO
! -------------------------------------------------------------------

! Special case for Q=0 ---------------------------------------------
      GG(0) = 0.D0
      IF ( L .EQ. 0 ) THEN
        DO IR = 1,NR
          R = IR * DR
          GG(0) = GG(0) + R*R * F(IR)
        ENDDO
        GG(0) = GG(0) * 2.D0 * C
      ENDIF
! -------------------------------------------------------------------

! Direct integration for the smallest Q's ---------------------------
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
! -------------------------------------------------------------------

! Special case for Q=QMAX -------------------------------------------
!     IF (L.EQ.0) THEN
!       GSUM = 0.D0
!       DO IQ = 1,NQ-1
!         Q = IQ * DQ
!         GSUM = GSUM + Q*Q * GG(IQ)
!       ENDDO
!       GSUM = GSUM * 4.D0 * PI * DQ
!       GG(NQ) = (2.D0*PI)**1.5D0 * F(0) - GSUM
!       GG(NQ) = GG(NQ) / (4.D0 * PI * DQ * QMAX**2)
!     ENDIF
! -------------------------------------------------------------------

! Copy from local to output array -----------------------------------
      DO IQ = 0,NQ
        G(IQ) = GG(IQ)
      ENDDO
! -------------------------------------------------------------------

      deallocate(P,GG,FN)

      end subroutine radfft

      SUBROUTINE FOUR1(DATA,NN,ISIGN)
!**********************************************************************
! Discrete Fourier transform. Modified and converted to double 
! precision from same routine in Numerical Recipes.
!**********************************************************************
! Input:
!   complex*16 DATA(NN) : Function to be Fourier transformed
!   integer    NN       : Number of points. Must be a power of 2
!   integer    ISIGN    : ISIG=+1/-1 => Direct/inverse transform
! Output:
!   complex*16 DATA(NN) : The direct Fourier transform (ISIG=+1), or
!                         NN times the inverse Fourier transf (ISIG=-1)
!**********************************************************************
      IMPLICIT NONE
      INTEGER          :: NN, ISIGN
      REAL(DP) :: DATA(2*NN)

      INTEGER          :: I, ISTEP, J, M, MMAX, N
      REAL(DP) :: TEMPI, TEMPR, THETA, WI, WPI, WPR, WR, WTEMP
      REAL(DP), PARAMETER :: TWOPI=6.28318530717959_dp,   &
                             HALF=0.5_dp, ONE=1._dp, TWO=2._dp, ZERO=0._dp

      N=2*NN
      J=1
      DO I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
        DO ! until following condition is met
          IF ((M.LT.2).OR.(J.LE.M)) EXIT
          J=J-M
          M=M/2
        END DO
        J=J+M
      END DO ! I
      MMAX=2
      DO ! until following condition is met
        IF (N.LE.MMAX) EXIT
        ISTEP=2*MMAX
        THETA=TWOPI/(ISIGN*MMAX)
        WPR=(-TWO)*SIN(HALF*THETA)**2
        WPI=SIN(THETA)
        WR=ONE
        WI=ZERO
        DO M=1,MMAX,2
          DO I=M,N,ISTEP
            J=I+MMAX
            TEMPR=WR*DATA(J)-WI*DATA(J+1)
            TEMPI=WR*DATA(J+1)+WI*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
          END DO ! I
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
        END DO ! M
        MMAX=ISTEP
      END DO ! until (N.LE.MMAX)

      END SUBROUTINE FOUR1

!
!  Helper routines for filtering package:
!
!  ordix.f order.f tqli.f tred2.f
!
      SUBROUTINE ORDIX( X, M, N, INDX )
!*******************************************************************
!Makes an index table of array X, with size N and stride M.
!Ref: W.H.Press et al. Numerical Recipes, Cambridge Univ. Press.
!Adapted by J.M.Soler from routine INDEXX of Num. Rec. May'96.
!*************** INPUT *********************************************
!REAL*8  X(M,N)   : Array with the values to be ordered
!INTEGER M, N     : Dimensions of array X
!*************** OUTPUT ********************************************
!INTEGER INDEX(N) : Array which gives the increasing order of X(1,I):
!                   X(1,INDEX(I)) .LE. X(1,INDEX(I+1)) )
!*************** USAGE *********************************************
!Example to order atomic positions X(I,IA), I=1,3, IA=1,NA by
!increasing z coordinate:
!   CALL ORDIX( X(3,1), 3, NA, INDEX )
!   CALL ORDER( X(1,1), 3, NA, INDEX )
!*******************************************************************
      IMPLICIT          NONE
      INTEGER           I, N, INDX(N), INDXT, IR, J, L, M
      DOUBLE PRECISION  X(M,N), Q

      DO 1 J=1,N
         INDX(J)=J
   1  CONTINUE
      IF (N.LE.1) RETURN
      L=N/2+1
      IR=N
   2  CONTINUE
         IF (L.GT.1) THEN
            L=L-1
            INDXT=INDX(L)
            Q=X(1,INDXT)
         ELSE
            INDXT=INDX(IR)
            Q=X(1,INDXT)
            INDX(IR)=INDX(1)
            IR=IR-1
            IF (IR.EQ.1) THEN
               INDX(1)=INDXT
               RETURN
            ENDIF
         ENDIF
         I=L
         J=L+L
   3     IF (J.LE.IR) THEN
            IF (J.LT.IR) THEN
               IF (X(1,INDX(J)).LT.X(1,INDX(J+1))) J=J+1
            ENDIF
            IF (Q.LT.X(1,INDX(J))) THEN
               INDX(I)=INDX(J)
               I=J
               J=J+J
            ELSE
               J=IR+1
            ENDIF
         GO TO 3
         ENDIF
         INDX(I)=INDXT
      GO TO 2

     END subroutine ordix


      SUBROUTINE ORDER( X, M, N, INDEX )
!*******************************************************************
!Orders array X(M,N) according to array INDEX(N), which may be
!  generated by routine ORDIX.
!Written by J.M.Soler. May'96.
!*************** INPUT *********************************************
!INTEGER M, N     : Dimensions of array X
!INTEGER INDEX(N) : Array which gives the desired order
!*************** INPUT AND OUTPUT **********************************
!REAL*8  X(M,N) : Array(s) to be ordered: Xout(I,J) = Xin(I,INDEX(J))
!*******************************************************************
      IMPLICIT          NONE
      INTEGER           I, N, INDEX(N), IORDER, ISTORE, J, M
      DOUBLE PRECISION  X(M,N), XI

      DO 40 J = 1,M
        DO 20 I = 1,N
          XI = X(J,I)
          IORDER = I
   10     CONTINUE
          ISTORE = INDEX(IORDER)
          IF (ISTORE .GT. 0) THEN
            IF (ISTORE .EQ. I) THEN
              X(J,IORDER) = XI
            ELSE
              X(J,IORDER) = X(J,ISTORE)
            ENDIF
            INDEX(IORDER) = -INDEX(IORDER)
            IORDER = ISTORE
            GOTO 10
          ENDIF
   20   CONTINUE
        DO 30 I = 1,N
          INDEX(I) = -INDEX(I)
   30   CONTINUE
   40 CONTINUE
      END subroutine order

      SUBROUTINE TQLI(D,E,N,NP,Z)

! IN COMBINATION WITH TRED2 FINDS EIGENVALUES AND EIGENVECTORS OF
! A REAL SYMMETRIC MATRIX. REF: W.H.PRESS ET AL. NUMERICAL RECIPES.

      implicit none

      integer, intent(in) :: np, n
      real(dp) D(NP),E(NP),Z(NP,NP)
      
      real(dp) :: zero, one, two
      PARAMETER (ZERO=0.0_dp,ONE=1.0_dp,TWO=2.0_dp)

      integer  :: iter, i, k, l, m
      real(dp) :: dd, g, r, s, c, p, f, b

      IF (N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
11      CONTINUE
        E(N)=ZERO
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30) STOP 'tqli: too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(TWO*E(L))
            R=SQRT(G**2+ONE)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=ONE
            C=ONE
            P=ZERO
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                C=G/F
                R=SQRT(C**2+ONE)
                E(I+1)=F*R
                S=ONE/R
                C=C*S
              ELSE
                S=F/G
                R=SQRT(S**2+ONE)
                E(I+1)=G*R
                C=ONE/R
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+TWO*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=ZERO
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END subroutine tqli


      SUBROUTINE TRED2(A,N,NP,D,E)

!HOUSEHOLDER REDUCTION OF A REAL SYMMETRIC MATRIX INTO TRIDIAGONAL FORM
!REF: W.H.PRESS ET AL. NUMERICAL RECIPES. CAMBRIDGE U.P.

      implicit none

      integer  ::  n, np
      real(dp) ::  A(NP,NP),D(NP),E(NP)

      real(dp) ::  zero, one
      PARAMETER (ZERO=0.0_dp,ONE=1.0_dp)

      integer  :: l, i, k, j
      real(dp) :: f, g, h, hh, scale

      IF(N.GT.1)THEN
        DO 18 I=N,2,-1
          L=I-1
          H=ZERO
          SCALE=ZERO
          IF(L.GT.1)THEN
            DO 11 K=1,L
              SCALE=SCALE+ABS(A(I,K))
11          CONTINUE
            IF(SCALE.EQ.ZERO)THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
12            CONTINUE
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=ZERO
              DO 15 J=1,L
                A(J,I)=A(I,J)/H
                G=ZERO
                DO 13 K=1,J
                  G=G+A(J,K)*A(I,K)
13              CONTINUE
                IF(L.GT.J)THEN
                  DO 14 K=J+1,L
                    G=G+A(K,J)*A(I,K)
14                CONTINUE
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
15            CONTINUE
              HH=F/(H+H)
              DO 17 J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 16 K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16              CONTINUE
17            CONTINUE
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
18      CONTINUE
      ENDIF
      D(1)=ZERO
      E(1)=ZERO
      DO 23 I=1,N
        L=I-1
        IF(D(I).NE.ZERO)THEN
          DO 21 J=1,L
            G=ZERO
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19          CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20          CONTINUE
21        CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=ONE
        IF(L.GE.1)THEN
          DO 22 J=1,L
            A(I,J)=ZERO
            A(J,I)=ZERO
22        CONTINUE
        ENDIF
23    CONTINUE
      RETURN
      END subroutine tred2

END MODULE m_filter
