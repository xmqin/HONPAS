!!@LICENSE
!
!******************************************************************************
! MODULE m_vv_vdwxc
! Implements the nonlocal correlation energy part of the van der Waals density
! functional of O.A.Vydrov & T.vanVoorhis, JCP 133, 244103 (2010) (VV2010):
!   Enlc = (1/2) Int Int dr1 dr2 n(r1) phi(n1,gn1,n2,gn2,r12) n(r2)
! where r12=|r2-r1|, n(r) is the electron density at point r, gn(r) its 
! gradient, and phi is a universal function defined in Eqs.(2-6) and (9).
! To be used with module m_vdwxc
! Refs: O.A.Vydrov & T.vanVoorhis, JCP 133, 244103 (2010)
!       G.Roman-Perez and J.M.Soler, PRL 103, 096102 (2009)
! Written by J.M.Soler. July 2012
!------------------------------------------------------------------------------
! Used module procedures:
!  use sys,         only: die               ! termination routine
!  use mesh1D,      only: get_mesh          ! Returns the mesh points
!  use m_radfft,    only: radfft            ! Radial fast Fourier transform
!  use mesh1D,      only: set_mesh          ! Sets a 1D mesh
!  use interpolation, only: spline            ! Sets spline interpolation
!  use interpolation, only: splint            ! Calculates spline interpolation
!------------------------------------------------------------------------------
! Used module parameters:
!   use precision,  only: dp                ! Real double precision type
!------------------------------------------------------------------------------
! Public procedures available from this module:
!   vv_vdw_beta      : Returns constant beta of eq.(13) of VV2010 JCP paper
!   vv_vdw_get_kmesh : Returns size and values of (kf,kg) mesh
!   vv_vdw_phi       : Finds phi(k) at (kf,kg) mesh points
!   vv_vdw_set_kcut  : Sets the planewave cutoff kc of the integration grid
!   vv_vdw_theta     : Finds function theta(n,gn) (eq.(8) of Roman-Soler PRL)
!------------------------------------------------------------------------------
! Public types, parameters, variables, and arrays:
!   None
!------------------------------------------------------------------------------
! Units: 
!   All lengths and energies in atomic units (Bohr, Hartree)
!******************************************************************************
! real(dp) function vv_vdw_beta()
!   Returns parameter beta=0.00497 of eq.(13) of VV2010 JCP paper
! Arguments:
!   none
! Sample usage:
!   real(dp):: beta
!   beta = vv_vdw_beta()
!------------------------------------------------------------------------------
! subroutine vdw_get_kmesh( nkf, nkg, kf, kg )
!   Returns size and values of the (kf,kg) interpolation mesh of the
!   Roman-Soler method. kf and kg are local Fermi- and gradient-wavevectors:
!   kf=(3*pi**2*n)**(1/3), kg=|grad(n)|/n
! Arguments:
!   integer,          intent(out) :: nkf    ! Number of kf mesh points
!   integer,          intent(out) :: nkg    ! Number of kg mesh points
!   real(dp),optional,intent(out) :: kf(:)  ! Mesh values of Fermi wavevect.
!   real(dp),optional,intent(out) :: kg(:)  ! Mesh values of grad(n)/n
! Sample usage:
!   integer :: nkf, nkg
!   real(dp),allocatable:: kf(:), kg(:)
!   call vv_vdw_get_kmesh( nkf, nkg )
!   allocate( kf(nkf), kg(nkg) )
!   call vv_vdw_get_kmesh( nkf, nkg, kf, kg )
! Notes:
! - If the size of arrays kf or kg is smaller than nkf and nkg, 
!   they are filled with the first size(kf) or size(kg) values
! - The size and values of the logarithmic k meshes are set by internal
!   parameters that can be changed only by editing them in this module:
!     nkf, nkg          : Number of kf and kg mesh points
!     kfcut=kfmesh(nkf) : Max. value of kf mesh
!     kgcut=kgmesh(nkg) : Max. value of kg mesh
!     dkfmaxdkfmin      : (kfmesh(nk) - kfmesh(nk-1)) / (kfmesh(2) - kfmesh(1))
!     dkgmaxdkgmin      : (kgmesh(nk) - kgmesh(nk-1)) / (kgmesh(2) - kgmesh(1))
!   Although the presently-set values have been found to yield good accuracy
!   in preliminary calculations, more tests may be required to guarantee
!   convergence in other systems. Larger nk's increase accuracy but CPU time 
!   increases between (nkf*nkg) and (nkf*nkg)**2.
!------------------------------------------------------------------------------
! subroutine vv_vdw_phi( k, phi, dphidk )
!   Finds phi(k1,k2,k) (Fourier transform of phi(k1,k2,r)) for all 
!   values of k1=(kf1,kg1) and k2=(kf2,kg2) of the (kf,kg) mesh. The mesh
!   indexes of kf and kg are combined into a single index ikfg=1,...,nkf*nkg
!   so that the size of phi and dphidk is (nkfg,nkfg), with nkfg=nkf*nkg.
!   In practice, it returns n1*n2*phi(k1,k2,k), where n1=kf1**3/(3*pi**2) and
!   n1=kf1**3/(3*pi**2) at the mesh values of kf1 and kf2. This is because
!   this function is smoother to interpolate than phi itself.
! Arguments:
!   real(dp),intent(in) :: k            ! Modulus of actual k vector
!   real(dp),intent(out):: phi(:,:)     ! phi(k1,k2,k) at given k, for all mesh
!                                       ! values of k1=(kf1,kg1) & k2=(kf2,kg2)
!   real(dp),intent(out):: dphidk(:,:)  ! dphi(k1,k2,k)/dk at given k
! Sample usage:
!   integer :: nkf, nkg
!   real(dp):: k, kcut
!   real(dp),allocatable:: phi(:,:), dphidk(:,:)
!   kcut = 10._dp ! 10 Bohr^-1 => 100 Ry (this should be adapted to your mesh)
!   call vv_vdw_set_kcut( kcut )
!   call vv_vdw_get_qmesh( nkf, nkg )
!   allocate( phi(nkf*nkg,nkf*nkg), dphidk(nkf*nkg,nkf*nkg) )
!   do k points
!     call vv_vdw_phi( k, phi, dphidk )
! Notes:
! - Requires a previous call to vv_vdw_set_kcut. Otherwise stops with an error.
! - Stops with an error message if size of array phi is smaller than
!   (nkf*nkg)**2.
!-----------------------------------------------------------------------------
! subroutine vv_vdw_set_kcut( kc )
!   Sets the reciprocal planewave cutoff kc of the integration grid, and finds 
!   the interpolation table to be used by vv_vdw_phi to obtain the vdW kernel 
!   phi at the reciprocal mesh wavevectors.
! Arguments:
!   real(dp),intent(in):: kc  ! Planewave cutoff: k>kcut => phi=0
! Sample usage:
!   real(dp):: kcut
!   kcut = 10._dp ! 10 Bohr^-1 => 100 Ry (this should be adapted to your mesh)
!   call vv_vdw_set_kcut( kcut )
!------------------------------------------------------------------------------
! subroutine vv_vdw_theta( nspin, rhos, grhos, theta, dtdrho, dtdgrho )
!   Finds the value and derivatives of function theta_ik(rho,grad_rho) 
!   of eq.(8) of the Roman-Soler PRL. ik is a combined index for the (kf,kg) 
!   interpolation mesh. 
!   In practice, it returns p_ik, rather than theta_ikgf=rho*p_ik, while
!   vv_vdw_phi returns rho1*rho2*phi, which is smoother to interpolate.
! Arguments:
!   integer, intent(in) :: nspin               ! Number of spin components
!   real(dp),intent(in) :: rhos(nspin)         ! Electron spin density
!   real(dp),intent(in) :: grhos(3,nspin)      ! Spin density gradient
!   real(dp),intent(out):: theta(nkf*nkg)      ! Function theta at mesh points
!   real(dp),intent(out):: dtdrho(nkf*nkg,nspin)     ! dtheta/drhos
!   real(dp),intent(out):: dtdgrho(3,nnkf*nkg,nspin) ! dtheta/dgrhos
! Sample usage:
!   integer,parameter:: nspin=1
!   integer,parameter:: nr=1000  ! number of mesh points in your system
!   integer :: nkf, nkg
!   real(dp):: n(nr), gn(3,nr)
!   real(dp),allocatable:: dtdn(:,:), dtdgn(:,:,:), theta(:,:)
!   call vv_vdw_get_qmesh( nkf, nkg )
!   nkfg = nkf*nkg
!   allocate( dtdn(nkfg,nr), dtdgn(3,nkfg,nr), theta(nkfg,nr) )
!   do ir ! mesh points
!     call vv_vdw_theta( nspin, n(ir), gn(:,ir), theta(:,ir), &
!                        dtdn(:,ir), dtdgn(:,:,ir) )
! Notes:
! - The values of kf(rho,grad_rho) and kg(rho,grad_rho) are saturated at
!   kfcut and kgcut (two internal parameters set below) to avoid incorrect 
!   'interpolation' of phi(k1,k2,r12) from (kf,kg) mesh points.
!******************************************************************************
! Algorithms:
!   Once we set the cutoff kcut (expressed as a wavevector), of the integration 
! mesh, to be used in evaluating the vdW nonlocal energy Enlc, we calculate and
! store (in memory) an interpolation table of phi(k1,k2,k) for all values 
! of k1=(kf1,kg1) and k2=(kf2,kg2) in the (kf,kg) mesh, and of k in a fine 
! radial grid of nk=nr points. This is done by first calculating phi(k1,k2,r) 
! in a radial grid in real space and then Fourier-transforming it to k space.
! This table is then used to interpolate phi(k1,k2,k) for any desired k.
! In practice, we interpolate (and return) rho1*rho2*phi(k1,k2,k) (which is 
! smoother than phi(k1,k2,k) itself), where rho=kf**3/(3*pi**2).
!   In order to ensure that values (kf,kg) are within the interpolation range,
! they are 'saturated' smoothly to cutoffs kfcut and kgcut. This implies an 
! approximation when either kf is very large (i.e. near the nucleus) or when 
! kg is very large, what tipically occurs in the tails of the electron density. 
! In the first case, Enlc is neglegible compared with Ex and the local part of 
! Ec. In the second case, it is neglegible because of the factor rho in the 
! integrand. Thus, the preliminary tests suggest that the presently-set values 
! of kfcut and kgcut give sufficiently accurate values.
!******************************************************************************

MODULE m_vv_vdwxc

! Used module procedures
  use sys,          only: die              ! Termination routine
  use mesh1D,       only: get_mesh         ! Returns the mesh points
  use m_radfft,     only: radfft           ! Radial fast Fourier transform
  use alloc,        only: re_alloc         ! Re-allocation routine
  use mesh1D,       only: set_mesh         ! Sets a 1D mesh
  use interpolation,only: spline           ! Sets spline interpolation
  use interpolation,only: splint           ! Calculates spline interpolation

! Used module parameters
  use precision,    only: dp               ! Real double precision type

#ifdef DEBUG_XC
  use debugXC,      only: udebug           ! File unit for debug output
!  use plot_module, only: plot
#endif /* DEBUG_XC */

  implicit none

! Called by m_vdwxc
PUBLIC ::            &
  vv_vdw_beta,       &! Returns parameter beta of the VV2010 functional
  vv_vdw_theta,      &! Finds function theta_ik(rho,grad_rho)
  vv_vdw_get_kmesh,  &! Returns size and values of (kf,kg) mesh
  vv_vdw_phi,        &! Finds and interpolates rho1*rho2*phi(k1,k2,k)
  vv_vdw_set_kcut     ! Sets the planewave cutoff kc of the integration grid

#ifdef DEBUG_XC
! Called by debugging test programs
PUBLIC :: &
  vv_vdw_phi_val,    &! Kernel phi(kf1,kf2,kg1,kg2,r12) with no interpolation
  vv_vdw_phiofr       ! Kernel phi(k1,k2,r) at tabulated k1,k2-mesh values
#endif /* DEBUG_XC */

PRIVATE  ! Nothing is declared public beyond this point

  ! Set parameters of the Vydrov-vanVoorhis functional,
  ! given in paragraph after eq.(27) of VV JCP 2010 paper
  real(dp),parameter:: vv_C = 0.0093
  real(dp),parameter:: vv_b = 5.9
  real(dp),parameter:: vv_beta = 0.00497

  ! Mesh parameters for table of phi(k1,k2,r) and its Fourier transform
  character(len=*),parameter:: kernelPrefactor='rho' ! Prefactor of the
                                             ! nonlocal kernel for interpolation
                                             ! ('kf'|'kf2'|'sqr_rho'|'rho')
  integer, parameter:: nr = 2048             ! Radial mesh points (power of 2)
  real(dp),parameter:: rcut = 100._dp        ! Radial cutoff: r>rcut => phi=0
  real(dp),parameter:: rmin = 1.e-6_dp       ! Min. radius as denominator
  integer, parameter:: nkf =  7              ! Number of Fermi wavevectors
  integer, parameter:: nkg =  5              ! Num. of grad(n)/n wavevectors
  integer, parameter:: nkfg = nkf*nkg        ! Num. of (kf,kg) mesh points
  real(dp),parameter:: kfcut = 5.0_dp        ! Max. Fermi wavevec.
  real(dp),parameter:: kgcut = 5.0_dp        ! Max. grad(n)/n
  real(dp),parameter:: dkfmaxdkfmin = 20.0_dp ! Last/first kf mesh interval
  real(dp),parameter:: dkgmaxdkgmin = 2.0_dp  ! Last/first kg mesh interval
  real(dp),parameter:: ktol = 1.e-12_dp      ! 'out of range' tolerance
  real(dp),parameter:: amin = 1.e-12_dp      ! tiny denominator to avoid /0
  real(dp),parameter:: rhoMin = 1.e-10_dp    ! Min. density for nonzero Vxc

  ! Parameters for cutoff function, used in radial Fourier transforms of phi
  integer, parameter:: ncut1 =  8      ! cutoff(x)=(1-x**ncut1)**ncut2
  integer, parameter:: ncut2 =  4

  ! Parameters for saturate function, used to enforce that q<qcut
  integer, parameter:: nsat = 15 ! h(x,xc)=1/(1/x**nsat+1/xc**nsat)**(1/nsat)
                                 ! nsat=15 approximately corresponds to 
                                 ! nsat=12 for eq.(5) of Roman-Soler PRL

  ! Private module variables and arrays
  logical, save:: kcut_set=.false.         ! Has kcut been set?
  logical, save:: kmesh_set=.false.        ! Has (kf,kg) mesh been set?
  logical, save:: phi_table_set=.false.    ! Has phi_table been set?
  integer, save:: nk                       ! Num. of radial mesh k-points
  real(dp),save:: dr                       ! r-mesh interval
  real(dp),save:: dk                       ! k-mesh interval
  real(dp),save:: kcut                     ! Planewave cutoff: k>kcut => phi=0
  real(dp),save:: kmax                     ! Max. k vector in Fourier transforms
  real(dp),save:: kfmesh(nkf)              ! Mesh points for Fermi wavevector
  real(dp),save:: kgmesh(nkg)              ! Mesh points for grad(n)/n
  real(dp),pointer,save:: &
                  phir(:,:,:)=>null(),    &! Table of phi(r)
                  phik(:,:,:)=>null(),    &! Table of phi(k)
                  d2phidr2(:,:,:)=>null(),&! Table of d^2(phi)/dr^2
                  d2phidk2(:,:,:)=>null()  ! Table of d^2(phi)/dk^2

CONTAINS

! -----------------------------------------------------------------------------

real(dp) function cutoff( x )

! Defines a smooth cutoff function that falls from y(0)=1 to y(1)=0

  implicit none
  real(dp),intent(in):: x

  if (x<=0._dp) then
    cutoff = 1
  else if (x>=1._dp) then
    cutoff = 0
  else
    cutoff = (1-x**ncut1)**ncut2
  end if

end function cutoff

!-----------------------------------------------------------------------------

subroutine iofk( kf, kg, ikf, ikg )

! Finds indexes ikf and ikg such that kfmesh(ikf) <= kf < kfmesh(ikf+1)
! and kgmesh(ikg) <= kg < kgmesh(ikg+1) in logarithmic meshes of the form
! k(ik) = b*(exp((ik-1)*a)-1)

  implicit none
  real(dp), intent(in) :: kf, kg
  integer,  intent(out):: ikf, ikg

  real(dp),save:: akf, akg, bkf, bkg
  logical, save:: first_call = .true.
  integer:: ik

  ! Mesh data initializations
  if (first_call) then
    ! Find logarithmic-mesh a & b parameters: x(j)=x(1)+b*(exp(a*(j-1))-1)
    akf = log( (kfmesh(nkf)-kfmesh(nkf-1)) / (kfmesh(2)-kfmesh(1)) ) / (nkf-2)
    akg = log( (kgmesh(nkg)-kgmesh(nkg-1)) / (kgmesh(2)-kgmesh(1)) ) / (nkg-2)
    akf = max( akf, amin )
    akg = max( akg, amin )
    bkf = (kfmesh(nkf) - kfmesh(1)) / (exp(akf*(nkf-1)) - 1)
    bkg = (kgmesh(nkg) - kgmesh(1)) / (exp(akg*(nkg-1)) - 1)
    ! Check that meshes are indeed logarithmic
    do ik = 1,nkf
      if (abs(kfmesh(ik)-kfmesh(1)-bkf*(exp(akf*(ik-1))-1))>ktol) &
         call die('vv_vdw_iofk ERROR: kfmesh not logarithmic')
    end do
    do ik = 1,nkg
      if (abs(kgmesh(ik)-kgmesh(1)-bkg*(exp(akg*(ik-1))-1))>ktol) &
         call die('vv_vdw_iofk ERROR: kgmesh not logarithmic')
    end do
    first_call = .false.
  end if

  ! Check that kf and kg are within interpolation range
  if (kf<kfmesh(1)-ktol .or. kf>kfmesh(nkf)+ktol .or. &
      kg<kgmesh(1)-ktol .or. kg>kgmesh(nkg)+ktol) then
     call die('vv_vdw_iofk ERROR: (kf,kg) out of range')
  endif

  ! Find interpolation mesh intervals
  ikf = 1 + int( log( 1 + (kf-kfmesh(1))/bkf ) / akf )
  ikg = 1 + int( log( 1 + (kg-kgmesh(1))/bkg ) / akg )
  ikf = max( 1, ikf )
  ikg = max( 1, ikg )
  ikf = min( nkf-1, ikf )
  ikg = min( nkg-1, ikg )

end subroutine iofk

! -----------------------------------------------------------------------------

subroutine kofn( n, gn, kf, kg, dkfdn, dkgdn, dkfdgn, dkgdgn )

! Finds Fermi and gradient wavevectors from density and gradient

  implicit none
  real(dp), intent(in) :: n          ! Electron density
  real(dp), intent(in) :: gn(3)      ! Density gradient
  real(dp), intent(out):: kf         ! Local Fermi wavevector
  real(dp), intent(out):: kg         ! |grad(n)|/n
  real(dp), intent(out):: dkfdn      ! dkf/dn
  real(dp), intent(out):: dkgdn      ! dkg/dn
  real(dp), intent(out):: dkfdgn(3)  ! dkf/dgn
  real(dp), intent(out):: dkgdgn(3)  ! dkg/dgn

  real(dp):: gn2, pi

! Trap exception for zero density
  if (n <= 1.e-30_dp) then
    kf = 0
    kg = 0
    dkfdn = 0
    dkfdgn = 0
    dkgdn = 0
    dkgdgn = 0
    return
  end if

! Find kf and kg
  pi = acos(-1._dp)
  kf = (3*pi**2*n)**(1._dp/3)
  gn2 = sum(gn**2)
  kg = sqrt(gn2)/n

! Find derivatives
  dkfdn = kf/n/3
  dkfdgn = 0
  dkgdn = -kg/n
  dkgdgn = kg*gn/gn2

end subroutine kofn

!-----------------------------------------------------------------------------

subroutine pofk( kf, kg, p, dpdkf, dpdkg )

! Finds the values and derivatives, at (kf,kg), of the bicubic polynomials 
! p_i(kf,kg) such that
!    y(kf,kg) = Sum_i p_i(kf,kg) * y_i
! is the bicubic spline interpolation at (kf,kg) of (any) function y with
! values y_i at mesh points (kfmesh,kgmesh)_i

  implicit none
  real(dp),intent(in) :: kf, kg ! point at which the polynomials are required
  real(dp),intent(out):: p(nkfg)      ! polynomial values at (kf,kg)
  real(dp),intent(out):: dpdkf(nkfg)  ! dp/dkf at (kf,kg)
  real(dp),intent(out):: dpdkg(nkfg)  ! dp/dkg at (kf,kg)

  integer :: ikf, ikg, ikfg, ikf0, ikg0
  real(dp):: a, b, dk, dkf0dkf, dkg0dkg, kf0, kg0
  real(dp):: pkf0(nkf), dpkfdkf0(nkf), pkg0(nkg), dpkgdkg0(nkg)
  logical, save :: first_call=.true.
  real(dp),save :: pkf(nkf,nkf), d2pkfdkf2(nkf,nkf)
  real(dp),save :: pkg(nkg,nkg), d2pkgdkg2(nkg,nkg)

! Set up spline polynomial basis
  if (first_call) then
    do ikf = 1,nkf
      pkf(:,ikf) = 0    ! ikf'th polynomial pkf(:,ikf) is one at kfmesh(ikf) 
      pkf(ikf,ikf) = 1  ! and zero at all other points
      call spline( kfmesh, pkf(:,ikf), nkf, 1.e30_dp, 1.e30_dp, &
                   d2pkfdkf2(:,ikf) )
!      call spline( kfmesh, pkf(:,ikf), nkf, 0._dp, 0._dp, d2pkfdkf2(:,ikf) )
    end do
    do ikg = 1,nkg
      pkg(:,ikg) = 0
      pkg(ikg,ikg) = 1
      call spline( kgmesh, pkg(:,ikg), nkg, 1.e30_dp, 1.e30_dp, &
                   d2pkgdkg2(:,ikg) )
!      call spline( kgmesh, pkg(:,ikg), nkg, 0._dp, 0._dp, d2pkgdkg2(:,ikg) )
    end do

!   DEBUG
    open(22,file='pkf.dat')
    do ikf = 1,nkf
      write(22,'(20e15.6)') kfmesh(ikf), pkf(:,ikf), d2pkfdkf2(:,ikf)
    end do
    close(22)
    open(22,file='pkg.dat')
    do ikg = 1,nkg
      write(22,'(20e15.6)') kgmesh(ikg), pkg(:,ikg), d2pkgdkg2(:,ikg)
    end do
    close(22)
!   END DEBUG

    first_call = .false.
  end if

! 'Saturate' (kf,kg) values to bring them to interpolation range
  call saturate( kf, kfcut, kf0, dkf0dkf )
  call saturate( kg, kgcut, kg0, dkg0dkg )

! Find interval of k mesh in which (kf0,kg0) point is included
  call iofk( kf0, kg0, ikf0, ikg0 )

! Evaluate pkf polynomials of spline basis at kf0
! The splint code is in-lined here because it is a hot point
  dk = kfmesh(ikf0+1) - kfmesh(ikf0)
  a = (kfmesh(ikf0+1) - kf0) / dk   ! dadkf0 = -1/dk
  b = (kf0 - kfmesh(ikf0)) / dk     ! dbdkf0 = +1/dk
  do ikf = 1,nkf
    pkf0(ikf) = a*pkf(ikf0,ikf) + b*pkf(ikf0+1,ikf) &
              + ((a**3-a)*d2pkfdkf2(ikf0,ikf) +     &
                 (b**3-b)*d2pkfdkf2(ikf0+1,ikf)) * dk**2/6
    dpkfdkf0(ikf) = - (pkf(ikf0,ikf) - pkf(ikf0+1,ikf)) / dk &
                    - ((3*a**2-1)*d2pkfdkf2(ikf0,ikf) -     &
                       (3*b**2-1)*d2pkfdkf2(ikf0+1,ikf)) * dk/6
  end do

! Evaluate pkg polynomials of spline basis at kg0
  dk = kgmesh(ikg0+1) - kgmesh(ikg0)
  a = (kgmesh(ikg0+1) - kg0) / dk   ! dadkg0 = -1/dk
  b = (kg0 - kgmesh(ikg0)) / dk     ! dbdkg0 = +1/dk
  do ikg = 1,nkg
    pkg0(ikg) = a*pkg(ikg0,ikg) + b*pkg(ikg0+1,ikg) &
              + ((a**3-a)*d2pkgdkg2(ikg0,ikg) +     &
                 (b**3-b)*d2pkgdkg2(ikg0+1,ikg)) * dk**2/6
    dpkgdkg0(ikg) = - (pkg(ikg0,ikg) - pkg(ikg0+1,ikg)) / dk &
                    - ((3*a**2-1)*d2pkgdkg2(ikg0,ikg) -     &
                       (3*b**2-1)*d2pkgdkg2(ikg0+1,ikg)) * dk/6
  end do

! Evaluate pkf*pkg polynomials at (kf0,kg0)
  ikfg = 0
  do ikg = 1,nkg
    do ikf = 1,nkf
      ikfg = ikfg+1
      p(ikfg) = pkf0(ikf) * pkg0(ikg)
      dpdkf(ikfg) = dpkfdkf0(ikf)*dkf0dkf * pkg0(ikg)
      dpdkg(ikfg) = pkf0(ikf) * dpkgdkg0(ikg)*dkg0dkg
    end do
  end do

end subroutine pofk

!-----------------------------------------------------------------------------

subroutine saturate( x, xc, y, dydx )

  ! Defines a function y(x,xc) = 1/(1/x^nsat+1/xc^nsat)^(1/nsat), where nsat
  ! is an integer set in the module header. This function is approximately
  ! equal to x for x<xc and it saturates to xc when x->infinity

  implicit none
  real(dp),intent(in) :: x     ! Independent variable
  real(dp),intent(in) :: xc    ! Saturation value
  real(dp),intent(out):: y     ! Function value
  real(dp),intent(out):: dydx  ! Derivative dy/dx

  real(dp):: z

  if (xc<=0._dp) then
    call die('vv_vdwxc_saturate ERROR: xc<=0')
  else if (x<0._dp) then
    call die('vv_vdwxc_saturate ERROR: x<0')
  else if (x==0._dp) then
    y = 0
    dydx = 1
  else
    z = 1/x**nsat + 1/xc**nsat
    y = 1/z**(1._dp/nsat)
    dydx = y/z/x**(nsat+1)
  end if

end subroutine saturate

!-----------------------------------------------------------------------------

subroutine saturate_inverse( y, xc, x, dxdy )

! Finds the inverse of the function defined in saturate subroutine:
!   y=1/(1/x^n+1/xc^n)^(1/n)  =>  x=1/(1/y^n-1/xc^n)^(1/n)

  implicit none
  real(dp),intent(in) :: y     ! Independent variable
  real(dp),intent(in) :: xc    ! Saturation value
  real(dp),intent(out):: x     ! Inverse function value
  real(dp),intent(out):: dxdy  ! Derivative dx/dy

  real(dp):: z

  if (xc<=0._dp) then
    call die('vv_vdwxc_saturate_inverse ERROR: xc<=0')
  else if (y<0._dp .or. y>=xc) then
    call die('vv_vdwxc_saturate_inverse ERROR: y out of range')
  else if (y==0._dp) then
    x = 0
    dxdy = 1
  else
    z = 1/y**nsat - 1/xc**nsat
    x = 1/z**(1._dp/nsat)
    dxdy = x/z/y**(nsat+1)
  end if

end subroutine saturate_inverse

!-----------------------------------------------------------------------------

subroutine set_kmesh()

! Sets mesh of q values

  implicit none
  integer :: mkf, mkg

  if (.not.kmesh_set) then
    call set_mesh( nkf, xmax=kfcut, dxndx1=dkfmaxdkfmin )
    call get_mesh( nkf, mkf, kfmesh )
    call set_mesh( nkg, xmax=kgcut, dxndx1=dkgmaxdkgmin )
    call get_mesh( nkg, mkg, kgmesh )
    kmesh_set = .true.
#ifdef DEBUG_XC
    write(udebug,'(/,a,/,(10f8.4))') 'vv_vdw_set_kmesh: kfmesh =', kfmesh
    write(udebug,'(/,a,/,(10f8.4))') 'vv_vdw_set_kmesh: kgmesh =', kgmesh
#endif /* DEBUG_XC */
  end if

end subroutine set_kmesh

! -----------------------------------------------------------------------------

subroutine set_phi_table()

! Finds and stores in memory the interpolation table (mesh points and 
! function values) for the kernel phi(k1,k2,k).

  implicit none
  character(len=*),parameter:: myName = 'vv_vdwxc/set_phi_table '
  integer :: ik, ikf1, ikf2, ikg1, ikg2, ik1, ik2, ir
  real(dp):: dkdk0, dphidk0, dphidkmax, dphidr0, dphidrmax, &
             k, kf1, kf2, kg1, kg2, pi, r(0:nr)

! Check that table was not set yet
  if (phi_table_set) return

! Check that kf, kg, r, and k meshes have been set
  if (.not.kmesh_set) call set_kmesh()
  if (.not.kcut_set) call die('vv_vdw_set_phi_table ERROR: kcut not set')
  forall(ir=0:nr) r(ir) = ir*dr
  pi = acos(-1.0_dp)

! Allocate arrays
  call re_alloc( phir,     0,nr, 1,nkfg, 1,nkfg, myName//'phir' )
  call re_alloc( phik,     0,nr, 1,nkfg, 1,nkfg, myName//'phik' )
  call re_alloc( d2phidr2, 0,nr, 1,nkfg, 1,nkfg, myName//'d2phidr2' )
  call re_alloc( d2phidk2, 0,nr, 1,nkfg, 1,nkfg, myName//'d2phidk2' )

! Loop on (k1,k2) mesh points
  do ikg2 = 1,nkg                    ! loop on kg2
    do ikf2 = 1,nkf                  ! loop on kf2
      ik2 = ikf2 + nkf*(ikg2-1)      ! combined (ikf2,ikg2) index
      do ikg1 = 1,nkg                ! loop on kg1
        do ikf1 = 1,nkf              ! loop on kf1
          ik1 = ikf1 + nkf*(ikg1-1)  ! combined (ikf1,ikg1) index
          if (ik1>ik2) cycle         ! since we will symmetrize at the end

          ! Saturated (kf,kg) values
          kf1 = kfmesh(ikf1)
          kf2 = kfmesh(ikf2)
          kg1 = kgmesh(ikg1)
          kg2 = kgmesh(ikg2)

          ! Find original (unsaturated) kf anf kg values
          ! call saturate_inverse( kfmesh(ikf1), kfcut, kf1, dkdk0 )
          ! call saturate_inverse( kfmesh(ikf2), kfcut, kf2, dkdk0 )
          ! call saturate_inverse( kgmesh(ikg1), kgcut, kg1, dkdk0 )
          ! call saturate_inverse( kgmesh(ikg2), kgcut, kg2, dkdk0 )

          ! Find kernel as a function of r12
          do ir = 0,nr
            phir(ir,ik1,ik2) = vv_vdw_phi_val( kf1, kf2, kg1, kg2, r(ir) )
          end do

          ! Kill kernel smoothly at r=rcut
          do ir = 0,nr
            phir(ir,ik1,ik2) = phir(ir,ik1,ik2) * cutoff( r(ir)/rcut )
          end do

          ! Find kernel in reciprocal space
          call radfft( 0, nr, rcut, phir(:,ik1,ik2), phik(:,ik1,ik2) )
          phik(:,ik1,ik2) = phik(:,ik1,ik2) * (2*pi)**1.5_dp

          ! Filter out above kcut
          phik(nk:nr,ik1,ik2) = 0

          ! Soft filter below kcut
          do ik = 1,nk
            k = ik * dk
            phik(ik,ik1,ik2) = phik(ik,ik1,ik2) * cutoff(k/kcut)
          end do

          ! Find filtered kernel in real space
          call radfft( 0, nr, kmax, phik(:,ik1,ik2), phir(:,ik1,ik2) )
          phir(:,ik1,ik2) = phir(:,ik1,ik2) / (2*pi)**1.5_dp

          ! Set up spline interpolation tables
          dphidr0 = 0      ! since phi(k1,k2,r) is even in r
          dphidk0 = 0      ! and therefore phi(k1,k2,k) is also even in k
          dphidrmax = 0    ! since phi->0 for r->infty
          dphidkmax = 0    ! and also when k->infty
          call spline( dr, phir(:,ik1,ik2), nr+1, dphidr0, dphidrmax, &
                       d2phidr2(:,ik1,ik2) )
          call spline( dk, phik(:,ik1,ik2), nk+1, dphidk0, dphidkmax, &
                       d2phidk2(:,ik1,ik2) )

          ! Fill symmetric elements
          phir(:,ik2,ik1) = phir(:,ik1,ik2)
          phik(:,ik2,ik1) = phik(:,ik1,ik2)
          d2phidr2(:,ik2,ik1) = d2phidr2(:,ik1,ik2)
          d2phidk2(:,ik2,ik1) = d2phidk2(:,ik1,ik2)

!#ifdef DEBUG_XC
!          if (.false. .and. ik1==ik2) then
!            print*, 'vv_vdw_set_kcut: ik1,ik2=', ik1, ik2
!            call window( 0._dp, 5._dp, -1._dp, 4._dp, 0 )
!            call axes( 0._dp, 1._dp, 0._dp, 1._dp )
!            call plot( nr+1, r, phi, phir(:,ik1,ik2) )
!            call window( 0._dp, 10._dp, -0.05_dp, 0.15_dp, 0 )
!            call axes( 0._dp, 1._dp, 0._dp, 0.05_dp )
!            call plot( nr+1, r, r**2*phi, r**2*phir(:,ik1,ik2))
!            call show()
!          end if
!#endif /* DEBUG_XC */

        end do ! ikf1
      end do ! ikg1
    end do ! ikf2
  end do ! ikg2

!#ifdef DEBUG_XC
!  open(17,file='vv_phi.table')
!  write(17,'(3a6,2a12,/,(3i6,2f15.9))') &
!    'ik1', 'ik2', 'ir', 'phi', 'd2phi/dk2', &
!!    (((ik1,ik2,ir,phir(ir,ik1,ik2),d2phidr2(ir,ik1,ik2), &
!    (((ik1,ik2,ir,phik(ir,ik1,ik2),d2phidk2(ir,ik1,ik2), &
!       ir=0,100),ik2=1,nkfg),ik1=1,nkfg)
!  close(17)
!#endif /* DEBUG_XC */

! Mark table as set
  phi_table_set = .true.

end subroutine set_phi_table

! -----------------------------------------------------------------------------

real(dp) function vv_vdw_phi_val( kf1, kf2, kg1, kg2, r12 )

! vdW energy kernel of Vydrov-vanVoorhis, eq.(2) of JCP 133, 244103 (2010)
! This subroutine calculates the 'raw' kernel directly, without interpolarions
! In practice, it returns n1*n2*phi(k1,k2,r12), which is smooth for kf->0
! Input:
!   real(dp):: kf1, kf2 ! Fermi wavevectors at points 1 and 2
!   real(dp):: kg1, kg2 ! |grad(n)|/n at points 1 and 2
!   real(dp):: r12      ! distance between points 1 and 2

! Arguments
  implicit none
  real(dp),intent(in) :: kf1, kf2, kg1, kg2, r12

! Internal variables
  real(dp):: g1, g2, kappa1, kappa2, kf1m, kf2m, n1, n2, &
             phi, pi, w01, w02, wg1, wg2, wp1, wp2

! Avoid dividing by zero when kf=0 
  kf1m = max(kf1,ktol)
  kf2m = max(kf2,ktol)

! Find kernel
  pi = acos(-1.0_dp)
  n1 = kf1m**3/(3*pi**2)          ! electron density at point 1
  n2 = kf2m**3/(3*pi**2)          ! electron density at point 2
  wp1 = sqrt(4*pi*n1)             ! local plasma frequency at point 1
  wp2 = sqrt(4*pi*n2)             ! local plasma frequency at point 2
  wg1 = sqrt(vv_C*kg1**4)         ! local band gap at point 1
  wg2 = sqrt(vv_C*kg2**4)         ! local band gap at point 2
  kappa1 = vv_b*kf1m**2/wp1       ! local VV kappa variable (eq.(9)) at point 1
  kappa2 = vv_b*kf2m**2/wp2       ! kappa variable at point 2
  w01 = sqrt(wg1**2+wp1**2/3)     ! local w0 frequency (eq.(5)) at point 1
  w02 = sqrt(wg2**2+wp2**2/3)     ! local w0 frequency at point 2
  g1 = w01*r12**2 + kappa1        ! local g variable (eq.(3)) at point 1
  g2 = w02*r12**2 + kappa2        ! local g variable at point 2
  phi = -1.5_dp/g1/g2/(g1+g2)     ! VV kernel phi (eq.(2))

  if (kernelPrefactor=='rho') then
    ! Return whole integrand of eq.(1)
    vv_vdw_phi_val = n1*n2*phi
  else if (kernelPrefactor=='kf') then
    ! Return kf1*kf2*phi
    vv_vdw_phi_val = kf1*kf2*phi
  else if (kernelPrefactor=='kf2') then
    ! Return (kf1*kf2)**2*phi
    vv_vdw_phi_val = (kf1*kf2)**2*phi
  else if (kernelPrefactor=='sqr_rho') then
    ! Find and return sqrt(n1*n2)*phi
    vv_vdw_phi_val = sqrt(n1*n2)*phi
  else
    call die('vv_vdw_phi_val ERROR: unknown kernelPrefactor')
  end if

end function vv_vdw_phi_val

! -----------------------------------------------------------------------------

real(dp) function vv_vdw_beta()

! Returns parameter beta of VV functional

  implicit none
  vv_vdw_beta = vv_beta

end function vv_vdw_beta

!-----------------------------------------------------------------------------

subroutine vv_vdw_get_kmesh( mkf, mkg, kf, kg )

! Returns size and values of (kf,kg) mesh

  implicit none
  integer,          intent(out) :: mkf    ! Number of kf mesh points
  integer,          intent(out) :: mkg    ! Number of kg mesh points
  real(dp),optional,intent(out) :: kf(:)  ! Values of kf mesh points
  real(dp),optional,intent(out) :: kg(:)  ! Values of kg mesh points
  integer:: nmax
  if (.not.kmesh_set) call set_kmesh()
  mkf = nkf
  mkg = nkg
  if (present(kf)) then
    nmax = max( nkf, size(kf) )
    kf(1:nmax) = kfmesh(1:nmax)
  end if
  if (present(kg)) then
    nmax = max( nkg, size(kg) )
    kg(1:nmax) = kgmesh(1:nmax)
  end if
end subroutine vv_vdw_get_kmesh

! -----------------------------------------------------------------------------

subroutine vv_vdw_phi( k, phi, dphidk )

! Finds by interpolation phi(k1,k2,k) (Fourier transform of phi(k1,k2,r)),
! with k1=(kf1,kg1), k2=(kf2,kg2) for all mesh values of kf (Fermi 
! wavevector) and kg (grad(n)/n). If the interpolation table does not exist,
! it is calculated in the first call to vv_vdw_phi. It requires a previous 
! call to vv_vdw_set_kc to set k mesh.

  implicit none
  real(dp),intent(in) :: k            ! Modulus of actual k vector
  real(dp),intent(out):: phi(:,:)     ! phi(k1,k2,k) at given k
                                      ! for all k1,k2 in (kf,kg) mesh
  real(dp),intent(out):: dphidk(:,:)  ! dphi(k1,k2,k)/dk at given k

  integer :: ik, ik1, ik2
  real(dp):: a, a2, a3, b, b2, b3

! Set interpolation table
  if (.not.phi_table_set) call set_phi_table()

! Check argument sizes
  if (size(phi,1)<nkfg .or. size(phi,2)<nkfg) &
    call die('vv_vdw_phi: ERROR: size(phi) too small')

! Find phi values at point k
  if (k >= kcut) then
    phi(:,:) = 0
    dphidk(:,:) = 0
  else
    ! Expand interpolation inline since this is the hottest point in VdW
    ik = int(k/dk)
    a = ((ik+1)*dk-k)/dk
    b = 1 - a
    a2 = (3*a**2 -1) * dk / 6
    b2 = (3*b**2 -1) * dk / 6
    a3 = (a**3 - a) * dk**2 / 6
    b3 = (b**3 - b) * dk**2 / 6
    do ik2 = 1,nkfg
      do ik1 = 1,ik2
!        call splint( dk, phik(:,ik1,ik2), d2phidk2(:,ik1,ik2), &
!                     nk+1, k, phi(ik1,ik2), dphidk(ik1,ik2), pr )
        phi(ik1,ik2) = a*phik(ik,ik1,ik2) + b*phik(ik+1,ik1,ik2) &
                + a3*d2phidk2(ik,ik1,ik2) + b3*d2phidk2(ik+1,ik1,ik2)
        dphidk(ik1,ik2) = (-phik(ik,ik1,ik2) &
                               +phik(ik+1,ik1,ik2) )/dk &
                - a2*d2phidk2(ik,ik1,ik2) + b2*d2phidk2(ik+1,ik1,ik2)
        phi(ik2,ik1) = phi(ik1,ik2)
        dphidk(ik2,ik1) = dphidk(ik1,ik2)
      end do
    end do
  end if

end subroutine vv_vdw_phi

!-----------------------------------------------------------------------------

subroutine vv_vdw_phiofr( r, phi )

! Finds phi(k1,k2,r) with k1=(kf1,kg1), k2=(kf2,kg2) for mesh values of
! kf's (Fermi wavevectors) and kg's (grad(n)/n)

  implicit none
  real(dp),intent(in) :: r
  real(dp),intent(out):: phi(:,:)

  integer :: ikf1, ikf2, ikg1, ikg2, ik1, ik2
  real(dp):: dphidr

  if (size(phi,1)<nkfg .or. size(phi,2)<nkfg) &
    stop 'vv_phiofr: ERROR: size(phi) too small'
  if (.not.phi_table_set) call set_phi_table()

  if (r >= rcut) then
    phi(:,:) = 0
  else
    do ik2 = 1,nkfg
      do ik1 = 1,ik2
        call splint( dr, phir(:,ik1,ik2), d2phidr2(:,ik1,ik2), &
                     nr+1, r, phi(ik1,ik2), dphidr )
        phi(ik2,ik1) = phi(ik1,ik2)
      end do ! ik1
    end do ! ik2
  end if ! (r>=rcut)

end subroutine vv_vdw_phiofr

!-----------------------------------------------------------------------------

subroutine vv_vdw_set_kcut( kc )

! Sets the reciprocal planewave cutoff kc of the integration grid, and finds 
! the interpolation table to be used by vdw_phi to obtain the vdW kernel phi
! at the reciprocal mesh wavevectors.

  implicit none
  real(dp),intent(in):: kc  ! Planewave cutoff: k>kcut => phi=0

  real(dp):: pi

  ! Set kcut and radial mesh parameters, if not already set
  if (.not. kc==kcut) then
    kcut = kc
    pi = acos(-1._dp)
    dr = rcut / nr
    dk = pi / rcut
    kmax = pi / dr
    nk = int(kcut/dk) + 1
    if (nk>nr) stop 'vv_vdw_set_kcut: ERROR: nk>nr'
    kcut_set = .true.
#ifdef DEBUG_XC
    write(udebug,'(a,5f8.3)') 'vv_vdw_set_kcut: kfcut,kgcut,rcut,kcut,kmax=', &
      kfcut, kgcut, rcut, kc, kmax
#endif /* DEBUG_XC */
  end if

  ! Set (kf,kg) mesh and phi table
  call set_kmesh()
  call set_phi_table()

end subroutine vv_vdw_set_kcut

!-----------------------------------------------------------------------------

subroutine vv_vdw_theta( nspin, rhos, grhos, theta, dtdrho, dtdgrho )

! Finds the value and derivatives of theta_i(rho,grad_rho) of eq.(8)
! of Roman-Soler PRL 2009. In practice, theta_i=p_i rather than rho*p_i,
! beacuse p_i is used here to expand rho1*rho2*phi rather than only phi.

  implicit none
  integer, intent(in) :: nspin                 ! Number of spin components
  real(dp),intent(in) :: rhos(nspin)           ! Electron spin density
  real(dp),intent(in) :: grhos(3,nspin)        ! Spin density gradient
  real(dp),intent(out):: theta(nkfg)           ! Exp. polynomials at (kf,kg)
  real(dp),intent(out):: dtdrho(nkfg,nspin)    ! dtheta(ik)/drhos
  real(dp),intent(out):: dtdgrho(3,nkfg,nspin) ! dtheta(ik)/dgrhos(ix)

  integer :: is, ix, ns
  real(dp):: rho, grho(3), dpdkf(nkfg), dpdkg(nkfg), dkfdrho, dkfdgrho(3), &
             dkgdrho, dkgdgrho(3), kf, kg, p(nkfg)

  ! Sum spin components of electron density
  ns = min(nspin,2)     ! num. of diagonal spin components (if noncollinear)
  rho = sum(rhos(1:ns))            ! local electron density
  grho = sum(grhos(:,1:ns),dim=2)  ! local density gradient

  ! If density is between threshold, return zero
  if (rho<rhoMin) then
    theta = 0
    dtdrho = 0
    dtdgrho = 0
    return
  end if

  ! Find local Fermi and gradient wavevectors, and their derivatives
  call kofn( rho, grho, kf, kg, dkfdrho, dkgdrho, dkfdgrho, dkgdgrho )

  ! Find expansion polynomials of integrand kernel of nonlocal vdW energy
  call pofk( kf, kg, p, dpdkf, dpdkg )

  ! Find theta functions and their derivatives with respect to rho and grho
  if (kernelPrefactor=='rho') then
    ! This is the right code if vv_vdw_phi_val returns rho1*rho2*phi
    theta(1:nkfg) = p(1:nkfg)
    do is = 1,ns
      dtdrho(1:nkfg,is) = dpdkf(1:nkfg)*dkfdrho + dpdkg(1:nkfg)*dkgdrho
      do ix = 1,3
        dtdgrho(ix,1:nkfg,is) = dpdkf(1:nkfg)*dkfdgrho(ix) &
                              + dpdkg(1:nkfg)*dkgdgrho(ix)
      end do
    end do
  else if (kernelPrefactor=='kf') then
    ! This is the right code if vv_vdw_phi_val returns kf1*kf2*phi
    kf = max(kf,ktol)
    theta(1:nkfg) = p(1:nkfg)*rho/kf
    do is = 1,ns
      dtdrho(1:nkfg,is) = (dpdkf(1:nkfg)*dkfdrho &
                          +dpdkg(1:nkfg)*dkgdrho)*rho/kf &
                        + p(1:nkfg)*(1/kf-rho/kf**2*dkfdrho)
      do ix = 1,3
        dtdgrho(ix,1:nkfg,is) = (dpdkf(1:nkfg)*dkfdgrho(ix) &
                                +dpdkg(1:nkfg)*dkgdgrho(ix))*rho/kf
      end do
    end do
  else if (kernelPrefactor=='kf2') then
    ! This is the right code if vv_vdw_phi_val returns (kf1*kf2)**2*phi
    kf = max(kf,ktol)
    theta(1:nkfg) = p(1:nkfg)*rho/kf**2
    do is = 1,ns
      dtdrho(1:nkfg,is) = (dpdkf(1:nkfg)*dkfdrho &
                          +dpdkg(1:nkfg)*dkgdrho)*rho/kf**2 &
                        + p(1:nkfg)*(1/kf**2-2*rho/kf**3*dkfdrho)
      do ix = 1,3
        dtdgrho(ix,1:nkfg,is) = (dpdkf(1:nkfg)*dkfdgrho(ix) &
                                +dpdkg(1:nkfg)*dkgdgrho(ix))*rho/kf**2
      end do
    end do
  else if (kernelPrefactor=='sqr_rho') then
    ! This is the right code if vv_vdw_phi_val returns sqrt(rho1*rho2)*phi
    rho = max(rho,1.e-12_dp)  ! to avoid division by zero
    theta(1:nkfg) = p(1:nkfg) * sqrt(rho)
    do is = 1,ns
      dtdrho(1:nkfg,is) = ( dpdkf(1:nkfg)*dkfdrho &
                          + dpdkg(1:nkfg)*dkgdrho ) * sqrt(rho) &
                        + p(1:nkfg) * 0.5/sqrt(rho)
      do ix = 1,3
        dtdgrho(ix,1:nkfg,is) = dpdkf(1:nkfg)*dkfdgrho(ix) * sqrt(rho) &
                              + dpdkg(1:nkfg)*dkgdgrho(ix) * sqrt(rho)
      end do
    end do
  else
    call die('vv_vdw_theta ERROR: unknown kernelPrefactor')
  end if ! (kernelPrefactor=='rho')

end subroutine vv_vdw_theta

END MODULE m_vv_vdwxc
