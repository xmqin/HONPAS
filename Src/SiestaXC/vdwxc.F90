!!@LICENSE
!
!******************************************************************************
! MODULE m_vdwxc
! Implements the nonlocal correlation energy part of the van der Waals density
! functional of M.Dion et al, PRL 92, 246401 (2004):
!   Enlc = (1/2) Int Int dr1 dr2 rho(r1) phi(q1*r12,q2*r12) rho(r2)
! with r12=|r2-r1|, q1=q0(rho(r1),grad_rho(r1)), q2=q0(rho(r2),grad_rho(r2)),
! where rho(r) is the electron density at point r, grad_rho(r) its gradient,
! and q0(rho,grad_rho) and phi(d1,d2) are universal functions defined in 
! Eqs.(11-12) and (14-16) of Dion et al.
! Using separate module m_vv_vdwxc, it implements also the similar functional
! of Vydrov and Van Voorhis. See that module for more information.
! Refs: M.Dion et al, PRL 92, 246401 (2004)
!       J.Klimes et al, JPCM 22, 022201 (2009)
!       V.R.Cooper, PRB 81, 161104(R) (2010)
!       K.Lee et al, PRB 82, 081101 (2010)
!       O.A.Vydrov and T.VanVoorhis, JCP 133, 244103 (2010)
!       G.Roman-Perez and J.M.Soler, PRL 103, 096102 (2009)
! Written by J.M.Soler. July 2007 - April 2010 and July 2012.
!------------------------------------------------------------------------------
! Used module procedures:
!  use sys,         only: die               ! termination routine
!  use mesh1D,      only: derivative        ! Derivative of a function in a mesh
!  use m_ggaxc      only: ggaxc             ! General GGA XC routine
!  use m_ldaxc,     only: ldaxc             ! General LDA XC routine
!  use mesh1D,      only: integral          ! Integral of a function in a mesh
!  use mesh1D,      only: get_mesh          ! Returns the mesh points
!  use mesh1D,      only: get_n             ! Returns the number of mesh points
!  use m_radfft,    only: radfft            ! Radial fast Fourier transform
!  use mesh1D,      only: set_interpolation ! Sets interpolation method
!  use mesh1D,      only: set_mesh          ! Sets a 1D mesh
!  use interpolation, only: spline            ! Sets spline interpolation
!  use interpolation, only: splint            ! Calculates spline interpolation
!  use m_vv_vdwxc,  only: vv_vdw_beta       ! Parameter beta of VV2010 functionl
!  use m_vv_vdwxc,  only: vv_vdw_theta      ! Func. theta of VV2010 functional
!  use m_vv_vdwxc,  only: vv_vdw_get_kmesh  ! Size and values of (kf,kg) mesh
!  use m_vv_vdwxc,  only: vv_vdw_phi        ! Interpolate rho1*rho2*phi(k1,k2,k)
!  use m_vv_vdwxc,  only: vv_vdw_set_kcut   ! Sets the planewave cutoff kc
!------------------------------------------------------------------------------
! Used module parameters:
!   use precision,  only: dp                ! Real double precision type
!------------------------------------------------------------------------------
! Public procedures available from this module:
!   vdw_decusp    : Energy due to the softening of the VdW kernel cusp
!   vdw_get_qmesh : Returns size and values of q-mesh
!   vdw_localxc   : LDA/GGA xc part apropriate for the used vdW flavour
!   vdw_phi       : Finds and interpolates phi(q1,q2,k)
!   vdw_set_author: Sets the vdW functional flavour (author initials)
!   vdw_set_kcut  : Sets the planewave cutoff kc of the integration grid
!   vdw_theta     : Finds function theta_q(rho,grad_rho)
!------------------------------------------------------------------------------
! Public types, parameters, variables, and arrays:
!   None
!------------------------------------------------------------------------------
! Units: 
!   All lengths and energies in atomic units (Bohr, Hartree)
!******************************************************************************
! subroutine vdw_decusp( nspin, rhos, grhos, eps, dedrho, dedgrho )
!   Finds the local energy correction due to the softening of the VdW kernel
!   cusp, defined as eps(rho,grad_rho) =
!   (1/2) * rho * Int 4*pi*r**2*dr * (phi(q1*r,q2*r) - phi_soft(q1*r,q2*r))
!   where q1=q2=q0(rho,grad_rho). Notice that grad_rho is included in the value 
!   of q0 at the origin but not in the change of rho(r) in the integrand.
!   phi_soft(d1,d2) is a softened version of the nonlocal VdW kernel phi(d1,d2)
!   (Eq.(14) of Dion et al) in which the logarithmic divergence at d1=d2=0 is
!   substituted by a smooth analytic function of the form defined in vdw_phi
! Arguments:
!   integer, intent(in) :: nspin            ! Number of spin components
!   real(dp),intent(in) :: rhos(nspin)      ! Electron spin density
!   real(dp),intent(in) :: grhos(3,nspin)   ! Spin density gradient
!   real(dp),intent(out):: eps              ! Energy correction per electron
!   real(dp),intent(out):: dedrho(nspin)    ! d(rho*eps)/d(rho)
!   real(dp),intent(out):: dedgrho(3,nspin) ! d(rho*eps)/d(grad_rho)
! Notes:
! - Requires a previous call to vdw_set_kcut. Otherwise stops with an error msg.
!------------------------------------------------------------------------------
! subroutine vdw_exchng( iRel, nSpin, D, GD, epsX, dEXdD, dEXdGD )
!   Finds the exchange energy density and its derivatives, using the GGA
!   functional apropriate for the previously-set vdW functional flavour
! Arguments:
!   integer, intent(in) :: iRel            ! Relativistic exchange? 0=no, 1=yes
!   integer, intent(in) :: nSpin           ! Number of spin components
!   real(dp),intent(in) :: D(nSpin)        ! Local electron (spin) density
!   real(dp),intent(in) :: GD(3,nSpin)     ! Gradient of electron density
!   real(dp),intent(out):: epsX            ! Exchange energy per electron
!   real(dp),intent(out):: dEXdD(nSpin)    ! dEx/dDens, Ex=Int(dens*epsX)
!   real(dp),intent(out):: dEXdGD(3,nSpin) ! dEx/dGrad(Dens)
! Sample usage:
!   integer,parameter:: iRel=0, nSpin=2
!   real(dp):: D(nSpin), dEXdD(nSpin), dEXdGD(3,nSpin), epsX, GD(3,nSpin)
!   call vdw_set_author('DRSLL')
!   do r points
!     Find D and GD at r
!     call vdw_exchng( iRel, nSpin, D, GD, epsX, dEXdD, dEXdGD )
!     Ex = Ex + dVolume*sum(D)*epsX
!   end do
!------------------------------------------------------------------------------
! subroutine vdw_get_qmesh( n, q )
!   Returns size and values of q-mesh
! Arguments:
!   integer,          intent(out) :: n      ! Number of q mesh points
!   real(dp),optional,intent(out) :: q(:)   ! Values of q mesh points
! Sample usage:
!   integer :: nq
!   real(dp):: kcut
!   real(dp),allocatable:: qmesh(:)
!   kcut = 10._dp ! 10 Bohr^-1 => 100 Ry (this should be adapted to your mesh)
!   call vdw_set_kcut( kcut )
!   call vdw_get_qmesh( nq )
!   allocate( qmesh(nq) )
!   call vdw_get_qmesh( nq, qmesh )
! Notes:
! - Requires a previous call to vdw_set_kcut. Otherwise stops with an error msg.
! - If the size of array q is smaller than that of the stored qmesh, it is 
!   filled with the first size(q) values of qmesh
! - The size and values of the logarithmic q mesh are set by internal
!   parameters that can be changed only by editing them in this module:
!     nq             : Number of q mesh points
!     qcut=qmesh(nq) : Max. value of q mesh
!     dqmaxdqmin     : (qmesh(nq) - qmesh(nq-1)) / (qmesh(2) - qmesh(1))
!   Although the presently-set values have been found to yield good accuracy
!   in preliminary calculations, more tests may be required to guarantee
!   convergence in other systems. The value of nq is particularly important:
!   larger nq increases accuracy but CPU time increases between nq and nq**2.
!------------------------------------------------------------------------------
! subroutine vdw_phi( k, phi, dphidk )
!   Finds phi_soft(q1,q2,k) (Fourier transform of phi_soft(q1,q2,r)) for all 
!   values of q1 and q2 in qmesh. The q's are local wavevectors defined in 
!   eqs.(11-12), and phi_soft is a smoothed version of Eq.(14) of Dion et al,
!   defined as phi_soft(d1,d2) = phi_soft(d,a) = phi0 + phi2*d**2 + phi4*d**4, 
!   where phi0 is a parameter, d=sqrt(d1**2+d2**2), a=atan(d2/d1), and phi2, 
!   phi4 are chosen so that phi_soft(d,a) matches phi(d,a) in value and slope
!   at d=dsoft (another parameter).
! Arguments:
!   real(dp),intent(in) :: k            ! Modulus of actual k vector
!   real(dp),intent(out):: phi(:,:)     ! phi(q1,q2,k) at given k
!                                       ! for all q1,q2 in qmesh
!   real(dp),intent(out):: dphidk(:,:)  ! dphi(q1,q2,k)/dk at given k
! Sample usage:
!   integer :: nq
!   real(dp):: k, kcut
!   real(dp),allocatable:: phi(:,:), dphidk(:,:)
!   kcut = 10._dp ! 10 Bohr^-1 => 100 Ry (this should be adapted to your mesh)
!   call vdw_set_kcut( kcut )
!   call vdw_get_qmesh( nq )
!   allocate( phi(nq,nq), dphidk(nq,nq) )
!   do k points
!     call vdw_phi( k, phi, dphidk )
! Notes:
! - Requires a previous call to vdw_set_kcut. Otherwise stops with an error msg.
! - Stops with an error message if size of array phi is smaller than nq*nq.
!-----------------------------------------------------------------------------
! subroutine vdw_set_author( author )
!   Sets the functional flavour (author initials) and subsequent parameters
! Arguments:
!   character(len=*),intent(in):: author ! Functnl flavour 
!                                     ('DRSLL'|'LMKLL'|'KBM'|'C09'|'BH'|'VV')
! Notes:
! - If vdw_set_author is not called, author='DRSLL' is set by default
! - Stops with an error message if author has not an allowed value
!-----------------------------------------------------------------------------
! subroutine vdw_set_kcut( kc )
!   Sets the reciprocal planewave cutoff kc of the integration grid, and finds 
!   the interpolation table to be used by vdw_phi to obtain the vdW kernel phi
!   at the reciprocal mesh wavevectors.
! Arguments:
!   real(dp),intent(in):: kc  ! Planewave cutoff: k>kcut => phi=0
! Notes:
! - An interpolation table to calculate the VdW kernel phi is read from file 
!   'vdw_kernel.table'. If the file does not exist, the table is calculated
!   and the file written when vdw_set_kcut is called.
!------------------------------------------------------------------------------
! subroutine vdw_theta( nspin, rhos, grhos, theta, dtdrho, dtdgrho )
!   Finds the value and derivatives of theta_i(rho,grad_rho) = rho*p_i(q0), 
!   where q0(rho,grad_rho) is the local wavevector defined in eqs.(11-12) 
!   of Dion et al, PRL 92, 246401 (2004). 
!   p_i(q0) are the cubic polynomials such that
!     y(q0) = Sum_i p_i(q0) * y_i
!   is the cubic spline interpolation at q0 of (any) function y(q) with
!   values y_i at mesh points qmesh_i
! Arguments:
!   integer, intent(in) :: nspin               ! Number of spin components
!   real(dp),intent(in) :: rhos(nspin)         ! Electron spin density
!   real(dp),intent(in) :: grhos(3,nspin)      ! Spin density gradient
!   real(dp),intent(out):: theta(nq)           ! Expansion of rho*q in qmesh
!   real(dp),intent(out):: dtdrho(nq,nspin)    ! dtheta(iq)/drhos
!   real(dp),intent(out):: dtdgrho(3,nq,nspin) ! dtheta(iq)/dgrhos(ix)
! Notes:
! - Requires a previous call to vdw_set_kcut
! - The value of q0(rho,grad_rho) is saturated at qmesh(nq)=qcut (a parameter)
!   to avoid incorrect 'interpolation' of phi(q1,q2,r12) from qmesh points.
!******************************************************************************
! Algorithms:
!   Although Eqs.(14-16) of Dion et al provide a straightforward definition of
! the nonlocal VdW kernel phi(d1,d2), its direct evaluation with the required 
! accuracy turns out to be too time consuming. In addition, the kernel has a 
! logarithmic divergence at d1=d2=0, which prevents its direct numerical 
! Fourier transform. Therefore, an elaborate layered procedure is followed:
!   A smoothed version of phi(d1,d2) is defined as 
!     phi_soft(d1,d2) = phi_soft(d,a) = phi0 + phi2*d**2 + phi4*d**4, 
! where phi0 is a parameter, d=sqrt(d1**2+d2**2), a=atan(d2/d1), and phi2, 
! phi4 are chosen so that phi_soft(d,a) matches phi(d,a) in value and slope
! at d=dsoft (another parameter).
!   The difference in energy between phi_soft and the right phi (called DEcusp) 
! is approximated in a local-density approximation as eps(rho,grad_rho) =
!   (1/2) * rho * Int 4*pi*r**2*dr * (phi(q1*r,q2*r) - phi_soft(q1*r,q2*r))
! where q1=q2=q0(rho,grad_rho). Notice that grad_rho is included in the value 
! of q0 at the origin but not in the change of rho(r) in the integrand. This
! could be improved in the future, although preliminary tests suggest that
! the corrections would be quite small.
!   A one-dimensional table of phi(d,d) is first calculated and stored as a
! function of d. Then, phi(d1,d2) is calculated as phi(dmax,dmax)+dphi(d1,d2),
! with dmax=max(d1,d2). The reason is that dphi converges better for large d,
! and therefore it requires smaller integration limits in Eq.(14) of Dion et al,
! while phi(dmax,dmax) is interpolated from the stored table.
!   Using the above methods, a two-dimensional table of phi_soft(d1,d2) is
! calculated and stored for all values of d1 and d2 in a logarithmic dmesh
! of size nd with a cutoff dcut (two internal parameters). This table is 
! written on disk file 'vdw_kernel.table' for reuse in subsequent runs.  
! If the file exists already, the table is simply read and stored in memory.
!   Once we set the cutoff kcut (expressed as a wavevector), of the integration 
! mesh to be used in evaluating the VdW nonlocal energy Enlc, we calculate and
! store (in memory) another interpolation table of phi_soft(q1,q2,k) for all
! values of q1 and q2 in qmesh, and of k in a fine radial grid of nk=nr points.
! This is done by first calculating phi_soft(q1,q2,r)=phi_soft(q1*r,q2*r) in a
! radial grid in real space and then Fourier-transforming it to k space.
! This table is then used to interpolate phi_soft(q1,q2,k) for any desired k.
!   In order to ensure that values of q are within the interpolation range,
! they are 'saturated' smoothly to a cutoff qcut=qmesh(nq). This implies an 
! approximation when either rho is very large (i.e. near the nucleus) or when 
! (grad_rho/rho)**2/kF -> infinity, what tipically occurs in the tails of
! the electron density. In the first case, Enlc is neglegible compared with
! Ex and the local part of Ec. In the second case, it is neglegible because of
! the factor rho in the integrand. Thus, the preliminary tests suggest that the
! presently-set value of qcut gives sufficiently accurate values.
!******************************************************************************

MODULE m_vdwxc

! Used module procedures
  use sys,         only: die               ! termination routine
  use mesh1D,      only: derivative        ! Derivative of a function in a mesh
  use mesh1D,      only: integral          ! Integral of a function in a mesh
  use mesh1D,      only: get_mesh          ! Returns the mesh points
  use mesh1D,      only: get_n             ! Returns the number of mesh points
  use m_ggaxc,     only: ggaxc             ! General GGA XC routine
  use m_ldaxc,     only: ldaxc             ! General LDA XC routine
  use m_radfft,    only: radfft            ! Radial fast Fourier transform
  use alloc,       only: re_alloc          ! Re-allocation routine
  use mesh1D,      only: set_interpolation ! Sets interpolation method
  use mesh1D,      only: set_mesh          ! Sets a 1D mesh
  use interpolation,only: spline           ! Sets spline interpolation
  use interpolation,only: splint           ! Calculates spline interpolation
  use m_vv_vdwxc,  only: vv_vdw_beta       ! Parameter beta of VV2010 functional
  use m_vv_vdwxc,  only: vv_vdw_theta      ! Func. theta of VV2010 functional
  use m_vv_vdwxc,  only: vv_vdw_get_kmesh  ! Size and values of (kf,kg) mesh
  use m_vv_vdwxc,  only: vv_vdw_phi        ! Interpolates rho1*rho2*phi(k1,k2,k)
  use m_vv_vdwxc,  only: vv_vdw_set_kcut   ! Sets the planewave cutoff kc

! Used module parameters
  use precision,   only: dp                ! Real double precision type

#ifdef DEBUG_XC
  use m_vv_vdwxc,  only: vv_vdw_phiofr     ! Interpolates rho1*rho2*phi(k1,k2,r)
  use debugXC,     only: udebug     ! File unit for debug output
!  use plot_module, only: plot
#endif /* DEBUG_XC */

  implicit none

! Called by xc routines
PUBLIC ::         &
  vdw_decusp,     &! Energy due to the softening of the VdW kernel cusp
  vdw_get_qmesh,  &! Returns size and values of q-mesh
  vdw_localxc,    &! LDA/GGA exchange-corr apropriate for the used vdW flavour
  vdw_phi,        &! Finds and interpolates phi(q1,q2,k)
  vdw_set_author, &! Sets the vdW functional flavour (author initials)
  vdw_set_kcut,   &! Sets the planewave cutoff kc of the integration grid
  vdw_theta        ! Finds function theta_q(rho,grad_rho)

#ifdef DEBUG_XC
! Called by debugging test programs
PUBLIC ::     &
  phiofr,     &! Finds the kernel phi(q1,q2,r) at tabulated q-mesh values
  phi_interp, &! Finds soft phi(d1,d2) kernel by interpolation of phi_table
  phi_soft,   &! Finds phi(d1,d2) softened near d1=d2=0
  phi_val,    &! Finds kernel phi(d1,d2) by direct integration
  pofq,       &! Finds polynomials p(q) from cubic spline interpolation
  qofrho       ! Finds the local wavevector parameter q(rho,grad_rho)
#endif /* DEBUG_XC */

PRIVATE  ! Nothing is declared public beyond this point

!  integer, parameter:: dp = kind(1.d0)

  ! Precision parameters for the integral defining phi in routine phi_val
  real(dp),parameter:: acutmin = 10.0_dp   ! Upper integration limit
  real(dp),parameter:: acutbyd  = 30.0_dp  ! Upper integration limit / d
  real(dp),parameter:: damax = 0.5_dp      ! Max. integration interval
  real(dp),parameter:: damin = 1.e-2_dp    ! Min. integration interval
  real(dp),parameter:: dabyd = 0.1_dp      ! Min. integration interval / d

  ! Precision parameter for the integral in routine dphi
  real(dp),parameter:: ashortbyd = 2.5_dp  ! Shorter integration limit / d

  ! Parameters for phi_soft. Some recommended pairs of values are
  ! (dsoft,phi0_soft)=(0.5,0.8)|(0.7|0.6)|(1.0|0.4)|(1.5,0.22)|(2.0,0.12)
  real(dp),parameter:: dsoft = 1.0_dp       ! Softening matching radius
  real(dp),parameter:: phi0_soft = 0.40_dp  ! phi_soft(0,0) (depends on dsoft)
  real(dp),parameter:: dd = 0.01_dp         ! Delta(d) for derivatives

  ! Mesh parameters for table phi(d1,d2)
  integer, parameter:: nd = 20               ! Number of d mesh points
  real(dp),parameter:: dcut = 30.0_dp        ! Max. value of d mesh
  real(dp),parameter:: ddmaxddmin = 20.0_dp  ! Last d mesh interval / first one

  ! Use routine dphi for better efficiency in setting phi table?
  logical,parameter:: use_dphi =.true. !(.true.=>efficiency|.false.=>accuracy)

  ! Set derivation methods to use for interpolation table
  character(len=*),parameter:: deriv_method = 'numeric'  !('numeric'|'interp')
  character(len=*),parameter:: interp_method= 'Spline' !('Lagrange'|'Spline')

  ! Parameters to find numerical derivatives of phi, to set interpolation
  real(dp),parameter:: ddbyd = 0.01_dp       ! Delta to find phi derivs / d
  real(dp),parameter:: ddmin = 0.001_dp      ! Min. delta to find phi derivs

  ! Mesh parameters for table of phi(q1,q2,r) and its Fourier transform
  integer, parameter:: nr = 2048             ! Radial mesh points (power of 2)
  integer, parameter:: mq = 30               ! Total number of q mesh points
!  integer, parameter:: nq = mq-1             ! Effective number of q mesh points
  real(dp),parameter:: qcut = 5.0_dp         ! Max. value of q mesh
  real(dp),parameter:: dqmaxdqmin = 20.0_dp  ! Last q mesh interval / first one
  real(dp),parameter:: rcut = 100._dp        ! Radial cutoff: r>rcut => phi=0
  real(dp),parameter:: rmin = 1.e-6_dp       ! Min. radius as denominator
  real(dp),parameter:: rsoft = 0.0_dp        ! Soften kernel in r<rsoft

  ! Parameters for cutoff function, used in radial Fourier transforms of phi
  integer, parameter:: ncut1 =  8      ! cutoff(x)=(1-x**ncut1)**ncut2
  integer, parameter:: ncut2 =  4

  ! Parameters for saturate function, used to enforce that q<qcut
  integer, parameter:: nsat  = 12      ! xsat(x)=1-exp(-sum_n=1:nsat x**n/n)

  ! Parameters for saturate_inverse function
  real(dp),parameter:: xmaxbyxc = 1.5_dp   ! qmax/qcut
  real(dp),parameter:: ytol = 1.e-15_dp     ! Tol. for saturated q

  ! Private module variables and arrays
  character(len=5),save:: vdw_author='DRSLL' ! Functional 'flavour' name
  real(dp),save:: dmesh(nd)                ! Mesh points for phi(d1,d2) table
  real(dp),save:: qmesh(mq)                ! Mesh points for phi(q1,q2,r)
  real(dp),save:: phi_table(0:3,0:3,nd,nd) ! Coefs. for bicubic interpolation
  logical, save:: phi_table_set=.false.    ! Has phi_table been set?
  logical, save:: qmesh_set=.false.        ! Has qmesh been set?
  logical, save:: kcut_set=.false.         ! Has kcut been set?
  real(dp),save:: dr                       ! r-mesh interval
  real(dp),save:: dk                       ! k-mesh interval
  real(dp),save:: kcut                     ! Planewave cutoff: k>kcut => phi=0
  integer, save:: nk                       ! # k points within kcut
  real(dp),save:: zab=-0.8491_dp           ! Parameter of the vdW functional
  real(dp),pointer,save:: &
                  phir(:,:,:)=>null(),    &! Table of phi(r)
                  phik(:,:,:)=>null(),    &! Table of phi(k)
                  d2phidr2(:,:,:)=>null(),&! Table of d^2(phi)/dr^2
                  d2phidk2(:,:,:)=>null()  ! Table of d^2(phi)/dk^2

!  real(dp),save:: dqmaxdqmin, qcut

CONTAINS

! -----------------------------------------------------------------------------

SUBROUTINE bcucof( n1, n2, x1, x2, y, dydx1, dydx2, d2ydx1dx2, c )
! Finds coefficients for bicubic interpolation
! Adapted from Numerical recipes
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n1, n2
  REAL(dp),INTENT(IN) :: x1(n1)
  REAL(dp),INTENT(IN) :: x2(n2)
  REAL(dp),INTENT(IN) :: y(n1,n2)
  REAL(dp),INTENT(IN) :: dydx1(n1,n2)
  REAL(dp),INTENT(IN) :: dydx2(n1,n2)
  REAL(dp),INTENT(IN) :: d2ydx1dx2(n1,n2)
  REAL(dp),INTENT(OUT):: c(0:3,0:3,n1,n2)

  INTEGER  :: i1, i11, i12, i13, i14, i2, i21, i22, i23, i24
  REAL(dp) :: dx1, dx2, wt(16,16), z(16)
  DATA wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,&
    8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,&
    2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,&
    2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,0,1,-2,1,5*0,&
    -3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0,&
    -1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,&
    -1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/

! Set coefs. for i1<n1 and i2<n2
  do i2 = 1,n2-1
    do i1 = 1,n1-1
      dx1 = x1(i1+1) - x1(i1)
      dx2 = x2(i2+1) - x2(i2)
      i11 = i1
      i12 = i1+1
      i13 = i1+1
      i14 = i1
      i21 = i2
      i22 = i2
      i23 = i2+1
      i24 = i2+1
      z( 1) = y(i11,i21)
      z( 2) = y(i12,i22)
      z( 3) = y(i13,i23)
      z( 4) = y(i14,i24)
      z( 5) = dydx1(i11,i21) * dx1
      z( 6) = dydx1(i12,i22) * dx1
      z( 7) = dydx1(i13,i23) * dx1
      z( 8) = dydx1(i14,i24) * dx1
      z( 9) = dydx2(i11,i21) * dx2
      z(10) = dydx2(i12,i22) * dx2
      z(11) = dydx2(i13,i23) * dx2
      z(12) = dydx2(i14,i24) * dx2
      z(13) = d2ydx1dx2(i11,i21) * dx1 * dx2
      z(14) = d2ydx1dx2(i12,i22) * dx1 * dx2
      z(15) = d2ydx1dx2(i13,i23) * dx1 * dx2
      z(16) = d2ydx1dx2(i14,i24) * dx1 * dx2
      z = matmul(wt,z)
      c(0:3,0:3,i1,i2) = reshape(z,(/4,4/),order=(/2,1/))
    end do ! i1
  end do ! i2

! Set c for i1=n1 and i2=n2 (valid only at the very border) using
! sum_i,j c(i,j,n1,i2) * 0**i * x2**j = sum_i,j c(i,j,n1-1,i2) * 1**i * x2**j
! sum_i,j c(i,j,i1,n2) * x1**i * 0**j = sum_i,j c(i,j,i1,n2-1) * x1**i * 1**j
! sum_i,j c(i,j,n1,n2) * 0**i * 0**j = sum_i,j c(i,j,n1-1,n2-1) * 1**i * 1**j
  c(:,:,n1,:) = 0
  c(:,:,:,n2) = 0
  c(0,:,n1,1:n2-1) = sum(c(:,:,n1-1,1:n2-1),dim=1)
  c(:,0,1:n1-1,n2) = sum(c(:,:,1:n1-1,n2-1),dim=2)
  c(0,0,n1,n2) = sum(c(:,:,n1-1,n2-1))

END SUBROUTINE bcucof

! -----------------------------------------------------------------------------

function cutoff( x )

  implicit none
  real(dp),intent(in):: x
  real(dp)           :: cutoff

  if (x<=0._dp) then
    cutoff = 1
  else if (x>=1._dp) then
    cutoff = 0
  else
    cutoff = (1-x**ncut1)**ncut2
  end if

end function cutoff

!-----------------------------------------------------------------------------

real(dp) function dphi( d1, d2 )

! Finds phi(d1,d2) - phi(dmax,dmax), with dmax=max(d1,d2)

  real(dp),intent(in) :: d1, d2

  integer  :: ia, ib, n, nmesh, nshort
  real(dp) :: a, acut, b, da1, dan, deff, dmax, dmin, gamma, &
              pi, t, t0, w
  real(dp),allocatable :: amesh(:), c(:), dphida(:), dphidb(:), &
                          s(:), nu0(:), nu1(:), nu2(:)

! Find integration mesh
  dmax = max( abs(d1), abs(d2) )
  dmin = min( abs(d1), abs(d2) )
  deff = dmax
!  deff = sqrt(d1**2+d2**2)
  acut = max( acutmin, acutbyd * deff )
  da1 = min( damax, dabyd * deff )
  da1 = max( da1, damin )
  dan = damax
  n = get_n( 0._dp, acut, da1, dan )
  allocate( amesh(n), c(n), dphida(n), dphidb(n), s(n), &
            nu0(n), nu1(n), nu2(n) )
  call set_mesh( n, xmax=acut, dxndx1=dan/da1 )
  call get_mesh( n, nmesh, amesh )

! Find limit of shorter mesh
  nshort = n
  do ia = n-1,1,-1
    if (amesh(ia) > ashortbyd*deff) nshort = ia
  end do

! Find cos(a), sin(a), nu1(a), and nu2(a)
  pi = acos(-1._dp)
  gamma = 4*pi/9
  do ia = 1,n
    a = amesh(ia)
    c(ia) = cos(a)
    s(ia) = sin(a)
    if (a==0._dp) then
      nu0(ia) = dmax**2 / 2 / gamma
      nu1(ia) = d1**2 / 2 / gamma
      nu2(ia) = d2**2 / 2 / gamma
    else ! (a/=0)
      if (d1==0._dp) then
        nu1(ia) = a*a / 2
      else
        nu1(ia) = a*a / 2 / (1-exp(-gamma*(a/d1)**2))
      end if
      if (d2==0._dp) then
        nu2(ia) = a*a / 2
      else
        nu2(ia) = a*a / 2 / (1-exp(-gamma*(a/d2)**2))
      end if
      if (dmax<=0._dp) then
        nu0(ia) = a*a / 2
      else
        nu0(ia) = a*a / 2 / (1-exp(-gamma*(a/dmax)**2))
      end if
    end if ! (a==0)
  end do

! Make integral on variable a
  dphida(1) = 0
  do ia = 2,nshort
    a = amesh(ia)

    ! Make integral on variable b
    dphidb(1) = 0
    do ib = 2,n
      b = amesh(ib)

      w = 2*( (3-a*a)*b*c(ib)*s(ia) + (3-b*b)*a*c(ia)*s(ib) &
            + (a*a+b*b-3)*s(ia)*s(ib) - 3*a*b*c(ia)*c(ib) )/(a*b)**3

      t = 0.5_dp * ( 1/(nu1(ia)+nu1(ib)) + 1/(nu2(ia)+nu2(ib)) ) &
                 * ( 1/(nu1(ia)+nu2(ia))/(nu1(ib)+nu2(ib)) &
                   + 1/(nu1(ia)+nu2(ib))/(nu2(ia)+nu1(ib)) )

      t0 = 0.5_dp * ( 1/(nu0(ia)+nu0(ib)) + 1/(nu0(ia)+nu0(ib)) ) &
                  * ( 1/(nu0(ia)+nu0(ia))/(nu0(ib)+nu0(ib)) &
                    + 1/(nu0(ia)+nu0(ib))/(nu0(ia)+nu0(ib)) )

      dphidb(ib) = a*a * b*b * w * (t-t0)

    end do ! ib
    dphida(ia) = integral( n, dphidb ) - integral( ia, dphidb )

  end do ! ia
  dphi = 2/pi**2 * 2*integral( nshort, dphida )

  deallocate( amesh, c, dphida, dphidb, s, nu0, nu1, nu2 )

#ifdef DEBUG_XC
!  print'(a,2f12.6,i8,f12.6)', 'dphi: d1,d2, na, phi =', d1, d2, n, dphi
#endif /* DEBUG_XC */

end function dphi

! -----------------------------------------------------------------------------

function dphi_fast( d1, d2 )

! Finds phi(d1,d2)-phi(dmax,dmax), with dmax=max(d1,d2), by 
!  - Direct integration if d < dsoft, where d=sqrt(d1**2+d2**2)
!  - Interpolation of phi_table if d > dsoft

  implicit none
  real(dp),intent(in) :: d1, d2
  real(dp)            :: dphi_fast

  real(dp):: d, dmax

  if (.not.phi_table_set) call set_phi_table()

  d = sqrt( d1**2 + d2**2 )
  dmax = max( d1, d2 )

  if (d < dsoft) then
    dphi_fast = dphi( d1, d2 )
  else
    dphi_fast = phi_interp( d1, d2 ) - phi_interp( dmax, dmax )
  end if

end function dphi_fast

! -----------------------------------------------------------------------------

function dphi_soft( d1, d2 )

! Finds phi_soft(d1,d2)-phi_soft(dmax,dmax), with dmax=max(d1,d2)

  implicit none
  real(dp),intent(in) :: d1, d2
  real(dp)            :: dphi_soft

  real(dp):: d, dmax

  d = sqrt( d1**2 + d2**2 )
  dmax = max( d1, d2 )

  if (d < dsoft) then
    dphi_soft = phi_soft( d1, d2 ) - phi_soft( dmax, dmax )
  else
    dphi_soft = dphi( d1, d2 )
  end if

end function dphi_soft

!-----------------------------------------------------------------------------

integer function iofd( d )

! Finds index i such that dmesh(i) <= d < dmesh(i+1)

  implicit none
  real(dp), intent(in) :: d

  real(dp),parameter :: amin = 1.e-12_dp
  real(dp),save:: a, b
  logical, save:: first_call = .true.

  if (first_call) then
    a = log( (dmesh(nd)-dmesh(nd-1)) / (dmesh(2)-dmesh(1)) ) / (nd-2)
    a = max( a, amin )
    b = (dmesh(2) - dmesh(1)) / (exp(a) - 1)
    first_call = .false.
  end if

  iofd = 1 + log( 1 + (d-dmesh(1))/b ) / a
  iofd = max( 1, iofd )
  iofd = min( nd-1, iofd )

end function iofd

!-----------------------------------------------------------------------------

integer function iofq( q )

! Finds index i such that qmesh(i) <= q < qmesh(i+1)

  implicit none
  real(dp), intent(in) :: q

  real(dp),parameter :: amin = 1.e-12_dp
  real(dp),save:: a, b
  logical, save:: first_call = .true.

  if (first_call) then
    a = log( (qmesh(mq)-qmesh(mq-1)) / (qmesh(2)-qmesh(1)) ) / (mq-2)
    a = max( a, amin )
    b = (qmesh(2) - qmesh(1)) / (exp(a) - 1)
    first_call = .false.
  end if

  iofq = 1 + log( 1 + (q-qmesh(1))/b ) / a
  iofq = max( 1, iofq )
  iofq = min( mq-1, iofq )

end function iofq

! -----------------------------------------------------------------------------

function phi_fast( d1, d2 )

! Finds hard phi(d1,d2) kernel by 
!  - Direct integration if d < dsoft, where d=sqrt(d1**2+d2**2)
!  - Interpolation of phi_table if d > dsoft

  implicit none
  real(dp),intent(in) :: d1, d2
  real(dp)            :: phi_fast

  real(dp):: d

  if (.not.phi_table_set) call set_phi_table()

  d = sqrt( d1**2 + d2**2 )

  if (d < dsoft) then
    phi_fast = phi_val( d1, d2 )
  else
    phi_fast = phi_interp( d1, d2 )
  end if

end function phi_fast

! -----------------------------------------------------------------------------

function phi_interp( d1, d2 )

! Finds soft phi(d1,d2) kernel by interpolation of phi_table

  implicit none
  real(dp),intent(in) :: d1, d2
  real(dp)            :: phi_interp

  integer :: i1, i2, id1, id2
  real(dp):: dd1(0:3), dd2(0:3)

  if (.not.phi_table_set) call set_phi_table()

  if (d1>=dcut .or. d2>=dcut) then
    phi_interp = 0
    return
  end if

  id1 = iofd( d1 )
  id2 = iofd( d2 )

  dd1(0) = 1
  dd2(0) = 1
  dd1(1) = (d1 - dmesh(id1)) / (dmesh(id1+1) - dmesh(id1))
  dd2(1) = (d2 - dmesh(id2)) / (dmesh(id2+1) - dmesh(id2))
  dd1(2) = dd1(1)**2
  dd2(2) = dd2(1)**2
  dd1(3) = dd1(1)**3
  dd2(3) = dd2(1)**3

  phi_interp = 0
  do i2 = 0,3
    do i1 = 0,3
      phi_interp = phi_interp + phi_table(i1,i2,id1,id2) * dd1(i1) * dd2(i2)
    end do
  end do

end function phi_interp

! -----------------------------------------------------------------------------

subroutine phiofr( r, phi )

! Finds phi(q1,q2,r) = phi(q1*r,q2*r) with q1=qmesh(iq1), q2=qmesh(iq2), 
! by interpolation from phi_table
! Notice that q is a density parameter, related to the Fermi wavevector

  implicit none
  real(dp),intent(in) :: r
  real(dp),intent(out):: phi(:,:)

  integer :: iq1, iq2
  real(dp):: dphidr

  ! Trap VV-version exception
#ifdef DEBUG_XC
  if (vdw_author=='VV') then
    call vv_vdw_phiofr( r, phi )
    return
  end if
#endif /* DEBUG_XC */

  if (size(phi,1)<mq .or. size(phi,2)<mq) &
    stop 'phiofr: ERROR: size(phi) too small'
  if (.not.qmesh_set) call set_qmesh()

  if (r >= rcut) then
    phi(:,:) = 0
  else
    do iq2 = 1,mq
      do iq1 = 1,iq2
        ! Use unfiltered kernel
!        phi(iq1,iq2) = phi_interp( qmesh(iq1)*r, qmesh(iq2)*r )
        ! Use filtered kernel
        call splint( dr, phir(:,iq1,iq2), d2phidr2(:,iq1,iq2), nr+1, r, &
                     phi(iq1,iq2), dphidr )
        phi(iq2,iq1) = phi(iq1,iq2)
      end do ! iq1
    end do ! iq2
  end if ! (r>=rcut)

end subroutine phiofr

!-----------------------------------------------------------------------------

function phi_soft( d1, d2 )

! Finds phi(d1,d2) softened near d1=d2=0

  implicit none
  real(dp),intent(in) :: d1, d2
  real(dp)            :: phi_soft

  real(dp):: d, d1m, d1p, d2m, d2p, dphidd, &
             phi, phi0, phi2, phi4, phim, phip

  d = sqrt( d1**2 + d2**2 )

  if (d<=0._dp) then
    phi_soft = phi0_soft
  else if (d > dsoft) then
    phi_soft = phi_val( d1, d2 )
  else ! (0<d<dsoft)

    d1p = d1/d * (dsoft+dd)
    d2p = d2/d * (dsoft+dd)
    d1m = d1/d * (dsoft-dd)
    d2m = d2/d * (dsoft-dd)
    phip = phi_val( d1p, d2p )
    phim = phi_val( d1m, d2m )
    phi = (phip + phim) / 2
    dphidd = (phip - phim) / (2*dd)
!    phi0 = phi - dphidd*dsoft/2
    phi0 = phi0_soft
    phi2 = (4*(phi-phi0) - dphidd*dsoft) / (2*dsoft**2)
    phi4 = (2*(phi0-phi) + dphidd*dsoft) / (2*dsoft**4)
    phi_soft = phi0 + phi2*d**2 + phi4*d**4
#ifdef DEBUG_XC
!    print'(a,5f8.3)', 'phi_soft: d,delta,phi0,phi2,phi4=', &
!      d, (d1-d2)/(d1+d2), phi0, phi2, phi4
#endif /* DEBUG_XC */

  end if ! (d<=0)

end function phi_soft

! -----------------------------------------------------------------------------

real(dp) function phi_val( d1, d2 )

! Finds kernel phi by direct integration of Eq.(14) of 
! Dion et al, PRL 92, 246401 (2004)

  real(dp),intent(in) :: d1, d2

  integer  :: ia, ib, n, nmesh
  real(dp) :: a, acut, b, da1, dan, deff, dmax, dmin, gamma, &
              pi, t, w
  real(dp),allocatable :: amesh(:), c(:), dphida(:), dphidb(:), &
                          s(:), nu1(:), nu2(:)

! Find integration mesh
  dmax = max( abs(d1), abs(d2) )
  dmin = min( abs(d1), abs(d2) )
!  deff = dmax
  deff = sqrt(d1**2+d2**2)
  acut = max( acutmin, acutbyd * deff )
  da1 = min( damax, dabyd * deff )
  da1 = max( da1, damin )
  dan = damax
  n = get_n( 0._dp, acut, da1, dan )
  allocate( amesh(n), c(n), dphida(n), dphidb(n), s(n), nu1(n), nu2(n) )
  call set_mesh( n, xmax=acut, dxndx1=dan/da1 )
  call get_mesh( n, nmesh, amesh )

#ifdef DEBUG_XC
!  print'(a,i6,/,(10f8.3))', 'phi_val: size, amesh =', n, amesh
#endif /* DEBUG_XC */

! Find cos(a), sin(a), nu1(a), and nu2(a)
  pi = acos(-1._dp)
  gamma = 4*pi/9
  do ia = 1,n
    a = amesh(ia)
    c(ia) = cos(a)
    s(ia) = sin(a)
    if (a==0._dp) then
      nu1(ia) = d1**2 / 2 / gamma
      nu2(ia) = d2**2 / 2 / gamma
    else ! (a/=0)
      if (d1==0._dp) then
        nu1(ia) = a*a / 2
      else
        nu1(ia) = a*a / 2 / (1-exp(-gamma*(a/d1)**2))
      end if
      if (d2==0._dp) then
        nu2(ia) = a*a / 2
      else
        nu2(ia) = a*a / 2 / (1-exp(-gamma*(a/d2)**2))
      end if
    end if ! (a==0)
  end do

! Make integral on variable a
  dphida(1) = 0
  do ia = 2,n
    a = amesh(ia)

    ! Make integral on variable b
    dphidb(1) = 0
    do ib = 2,n
      b = amesh(ib)

      w = 2*( (3-a*a)*b*c(ib)*s(ia) + (3-b*b)*a*c(ia)*s(ib) &
            + (a*a+b*b-3)*s(ia)*s(ib) - 3*a*b*c(ia)*c(ib) )/(a*b)**3

      t = 0.5_dp * ( 1/(nu1(ia)+nu1(ib)) + 1/(nu2(ia)+nu2(ib)) ) &
                 * ( 1/(nu1(ia)+nu2(ia))/(nu1(ib)+nu2(ib)) &
                   + 1/(nu1(ia)+nu2(ib))/(nu2(ia)+nu1(ib)) )

      dphidb(ib) = a*a * b*b * w * t

    end do ! ib
    dphida(ia) = integral( n, dphidb )

  end do ! ia
  phi_val = 2/pi**2 * integral( n, dphida )

  deallocate( amesh, c, dphida, dphidb, s, nu1, nu2 )

#ifdef DEBUG_XC
!  print'(a,2f12.6,i8,f12.6)', 'phi_val: d1,d2, na, phi =', d1, d2, n, phi_val
#endif /* DEBUG_XC */

end function phi_val

! -----------------------------------------------------------------------------

subroutine pofq( q0, p0, dp0dq0 )

! Finds the values and derivatives, at q0, of the cubic polynomials 
! p_i(q0) such that
!    y(q0) = Sum_i p_i(q0) * y_i
! is the cubic spline interpolation at q0 of (any) function y(q) with
! values y_i at mesh points qmesh_i

  implicit none
  real(dp),intent(in) :: q0
  real(dp),intent(out):: p0(mq)
  real(dp),intent(out):: dp0dq0(mq)

  integer :: iq, iq0
  real(dp):: a, b, dq
  logical, save :: first_call=.true.
  real(dp),save :: p(mq,mq), d2pdq2(mq,mq)

! Set up spline polynomial basis
  if (first_call) then
    p = 0
    do iq = 1,mq
      p(iq,iq) = 1
      call spline( qmesh, p(:,iq), mq, 1.e30_dp, 1.e30_dp, d2pdq2(:,iq) )
!      call spline( qmesh, p(:,iq), mq, 0._dp, 0._dp, d2pdq2(:,iq) )
    end do
    first_call = .false.
  end if

! Find interval of qmesh in which q0 is included
  if (q0>qmesh(mq)) then   ! q0 out of range
    p0 = 0
    dp0dq0 = 0
    return
  end if
  iq0 = iofq( q0 )

! Evaluate polynomials of spline basis
  dq = qmesh(iq0+1) - qmesh(iq0)
  a = (qmesh(iq0+1) - q0) / dq   ! dadq0 = -1/dq
  b = (q0 - qmesh(iq0)) / dq     ! dbdq0 = +1/dq
  do iq = 1,mq
    p0(iq) = a*p(iq0,iq) + b*p(iq0+1,iq) &
      + ((a**3-a)*d2pdq2(iq0,iq) + (b**3-b)*d2pdq2(iq0+1,iq)) * dq**2/6
    dp0dq0(iq) = - (p(iq0,iq) - p(iq0+1,iq)) / dq &
      - ((3*a**2-1)*d2pdq2(iq0,iq) - (3*b**2-1)*d2pdq2(iq0+1,iq)) * dq/6
  end do

end subroutine pofq

!-----------------------------------------------------------------------------

subroutine qofrho( rho, grho, q, dqdrho, dqdgrho )

! Finds the local wavevector parameter q0 defined in eqs.(11-12) of
! M.Dion et al, PRL 92, 246401 (2004)

  implicit none
  real(dp), intent(in) :: rho        ! Electron density
  real(dp), intent(in) :: grho(3)    ! Density gradient
  real(dp), intent(out):: q          ! Local wave vector parameter q0
  real(dp), intent(out):: dqdrho     ! d_q/d_rho
  real(dp), intent(out):: dqdgrho(3) ! d_q/d_grho

  character(len=*),parameter :: author = 'PW92'  ! Perdew-Wang'92 for LDA
  integer,         parameter :: irel   = 0       ! Non-relativistic exchange
  integer,         parameter :: nspin  = 1       ! Unpolarized electron gas
  real(dp):: decdrho, dexdrho, dkfdrho, dq0dgrho(3), dq0dgrho2, dq0drho, &
             dqdq0, dvxdrho(nspin), dvcdrho(nspin), ex, ec, grho2, &
             kf, pi, q0, rhos(nspin), vx(nspin), vc(nspin)

! Trap exception for zero density
  if (rho <= 1.e-15_dp) then
    q = qcut
    dqdrho = 0
    dqdgrho = 0
    return
  end if

  pi = acos(-1._dp)
  kf = (3*pi**2 * rho)**(1._dp/3)

! Find exchange and correlation energy densities
  rhos(1) = rho
  call ldaxc( author, irel, nspin, rhos, ex, ec, vx, vc, dvxdrho, dvcdrho )

! Find q
  grho2 = sum(grho**2)
  q0 = ( 1 + ec/ex - zab/9 * grho2 / (2*kf*rho)**2 ) * kf

! Find derivatives
  dkfdrho = kf / (3*rho)
  dexdrho = (vx(1) - ex) / rho  ! Since vx = d(rho*ex)/d_rho = ex + rho*dexdrho
  decdrho = (vc(1) - ec) / rho
  dq0drho = ( decdrho/ex - ec/ex**2*dexdrho + 2 * zab/9 * grho2 / &
              (2*kf*rho)**3 * (2*dkfdrho*rho + 2*kf) ) * kf &
            + q0/kf * dkfdrho
  dq0dgrho2 = -(zab/9) / (2*kf*rho)**2 * kf
  dq0dgrho(:) = dq0dgrho2 * 2*grho(:)  ! Since d(vector**2)/d_vector = 2*vector

! Saturate q to qcut smoothly
  call saturate( q0, qcut, q, dqdq0 )
  dqdrho = dqdq0 * dq0drho
  dqdgrho = dqdq0 * dq0dgrho

end subroutine qofrho

!-----------------------------------------------------------------------------

subroutine saturate( x, xc, y, dydx )

  ! Defines a function y(x) = xc * (1 - exp(-Sum_n=1:nsat (x/xc)**n/n))
  ! It is approx. equal to x for x<xc and it saturates to xc when x->infinity

  implicit none
  real(dp),intent(in) :: x     ! Independent variable
  real(dp),intent(in) :: xc    ! Saturation value
  real(dp),intent(out):: y     ! Function value
  real(dp),intent(out):: dydx  ! Derivative dy/dx

  integer :: n
  real(dp):: dpdx, p

  if (nsat >= 100) then
    if (x < xc) then
      y = x
      dydx = 1
    else ! (x >= xc)
      y = xc
      dydx = 0
    end if ! (x < xc)
  else ! (nsat < 100)
!    This is the straightforward polynomial evaluation
!    p = x/xc
!    dpdx = 1/xc
!    do n = 2,nsat
!      p = p + (x/xc)**n / n
!      dpdx = dpdx + (x/xc)**(n-1) / xc
!    end do
!   And this is a more accurate way to evaluate it
    p = (x/xc)/nsat
    dpdx = 1._dp/xc
    do n = nsat-1,1,-1
      p = (p + 1._dp/n) * x/xc
      dpdx = (dpdx*x + 1) / xc
    end do
    y = xc * (1 - exp(-p))
    dydx = xc * dpdx * exp(-p)
  end if ! (nsat >= 100)

end subroutine saturate

!-----------------------------------------------------------------------------

subroutine saturate_inverse( y, xc, x, dydx )

! Finds the inverse of the function defined in saturate subroutine

  implicit none
  real(dp),intent(in) :: y     ! Independent variable
  real(dp),intent(in) :: xc    ! Saturation value
  real(dp),intent(out):: x     ! Inverse function value
  real(dp),intent(out):: dydx  ! Derivative dy/dx

  real(dp):: x1, x2, yx

  if (y<0._dp .or. y>xc) stop 'vdw:saturate_inverse: y out of range'
  x1 = 0
  x2 = xmaxbyxc * xc
  do
    x = (x1+x2)/2
    call saturate( x, xc, yx, dydx )
    if (abs(y-yx)<ytol) then
      return
    else if (yx < y) then
      x1 = x
    else
      x2 = x
    end if
  end do

end subroutine saturate_inverse

!-----------------------------------------------------------------------------

subroutine set_phi_table()

! Finds and writes in disk the interpolation table (mesh points and function 
! values) for the kernel phi(d1,d2). If the table file already exists, it is
! simply read and stored in memory.

  implicit none

  logical :: file_found
  integer :: id, id1, id2, nmesh
  real(dp):: d, d1, d1m, d1p, d2, d2m, d2p, dd, &
             dphidd1(nd,nd), dphidd2(nd,nd), d2phidd1dd2(nd,nd), &
             phi(nd,nd), phi1, phim, phip, phimm, phimp, phipm, phipp

! Read file with table, if present
  inquire( file='vdw_kernel.table', exist=file_found )
  if (file_found) then
    open( unit=1, file='vdw_kernel.table', form='unformatted' )
    read(1,end=1) nmesh
    if (nmesh==nd) then
      read(1,end=1) dmesh
      read(1,end=1) phi_table
    end if
    close( unit=1 )
    phi_table_set = .true.
    return
  end if
1 continue ! Come here if end-of-file found

! Set d-mesh points
  call set_mesh( nd, xmax=dcut, dxndx1=ddmaxddmin )
  call get_mesh( nd, nmesh, dmesh )
#ifdef DEBUG_XC
  write(udebug,'(a,/,(10f8.4))') 'set_phi_table: dmesh =', dmesh
#endif /* DEBUG_XC */

! Find function at mesh points
  do id1 = 1,nd
    d1 = dmesh(id1)
    phi1 = phi_soft( d1, d1 )
    phi(id1,id1) = phi1
    do id2 = 1,id1-1
      d2 = dmesh(id2)
      d = sqrt( d1**2 + d2**2 )
      if (d < dsoft) then
        phi(id1,id2) = phi_soft( d1, d2 )
      else
        if (use_dphi) then ! Use dphi for better efficiency
          phi(id1,id2) = phi1 + dphi( d1, d2 )
        else ! Use only phi_val, to eliminate uncertainties
          phi(id1,id2) = phi_val( d1, d2 )
        end if
      end if
      phi(id2,id1) = phi(id1,id2)
    end do ! id2
  end do ! id1

#ifdef DEBUG_XC
!  open( unit=44, file='phi.out' )
!  do id2 = 1,nd
!    write(44,'(/,(f12.6))') phi(:,id2)
!  end do
!  close( unit=44 )
#endif /* DEBUG_XC */

  if (deriv_method == 'numeric') then
!    print*, 'set_phi_table: Using numerical derivatives'

    ! Find derivatives at mesh points
     do id1 = 1,nd
      d1 = dmesh(id1)
      dd = ddbyd * d1
      dd = max( dd, ddmin )
!      d1 = max( d1, dd )
      d1m = d1 - dd
      d1p = d1 + dd
      phim = phi_soft( d1m, d1m )
      phip = phi_soft( d1p, d1p )
      do id2 = 1,id1
        d2  = dmesh(id2)
!        d2 = max( d2, dd )
        d = sqrt( d1**2 + d2**2 )
        d1m = d1 - dd
        d1p = d1 + dd
        d2m = d2 - dd
        d2p = d2 + dd

        if (d < dsoft) then
          phimm = phi_soft( d1m, d2m )
          phipm = phi_soft( d1p, d2m )
          phipp = phi_soft( d1p, d2p )
          phimp = phi_soft( d1m, d2p )
        else ! (d>dsoft)
          if (use_dphi) then
            phimm = phim + dphi( d1m, d2m )
            phipm = phip + dphi( d1p, d2m )
            phipp = phip + dphi( d1p, d2p )
            if (id1==id2) then
              phimp = phip + dphi( d1m, d2p )
            else
              phimp = phim + dphi( d1m, d2p )
            end if
          else ! (.not.use_dphi)
            phimm = phi_val( d1m, d2m )
            phipm = phi_val( d1p, d2m )
            phipp = phi_val( d1p, d2p )
            phimp = phi_val( d1m, d2p )
          end if ! (use_dphi)
        end if ! (d<dsoft)

        dphidd1(id1,id2)     = (phipp+phipm-phimp-phimm) / (4*dd)
        dphidd2(id1,id2)     = (phipp-phipm+phimp-phimm) / (4*dd)
        d2phidd1dd2(id1,id2) = (phipp-phipm-phimp+phimm) / (2*dd)**2
        dphidd1(id2,id1)     = dphidd2(id1,id2)
        dphidd2(id2,id1)     = dphidd1(id1,id2)
        d2phidd1dd2(id2,id1) = d2phidd1dd2(id1,id2)

      end do ! id2
    end do ! id1

  else if (deriv_method == 'interp') then
!    print*, 'set_phi_table: Using interpolation for derivatives'

    ! Reset mesh, which has been changed by phi_val
    call set_mesh( nd, xmax=dcut, dxndx1=ddmaxddmin )

    ! Set interpolation method
    if (interp_method=='Lagrange') then
      call set_interpolation( 'Lagrange' )
    else if (interp_method=='Spline') then
!      call set_interpolation( 'Spline', huge(phi), huge(phi) )
      call set_interpolation( 'Spline', 0._dp, 0._dp )
    else
      stop 'set_phi_val: ERROR: Unknown interp_method'
    end if

    ! Find first partial derivative d_phi/d_d1
    do id = 1,nd
      dphidd1(:,id) = derivative( nd, phi(:,id) )
      dphidd2(id,:) = dphidd1(:,id)
    end do

    ! Find second cross partial derivative d_phi/d_d1/d_d2
    do id = 1,nd
      d2phidd1dd2(id,:) = derivative( nd, dphidd1(id,:) )
    end do

    ! Symmetrize d_phi/d_d1/d_d2
    do id2 = 2,nd
      do id1 = 1,id2-1
        d2phidd1dd2(id1,id2) = (d2phidd1dd2(id1,id2) + &
                                d2phidd1dd2(id2,id1)) / 2
        d2phidd1dd2(id2,id1) = d2phidd1dd2(id1,id2) 
      end do
    end do

  else
    stop 'set_phi_table: ERROR: Unknown deriv_method'
  end if ! (deriv_method)

! Make values and derivatives strictly zero when d1=dmax or d2=dmax
  phi(:,nd) = 0
  phi(nd,:) = 0
  dphidd1(:,nd) = 0
  dphidd1(nd,:) = 0
  dphidd2(:,nd) = 0
  dphidd2(nd,:) = 0
  d2phidd1dd2(:,nd) = 0
  d2phidd1dd2(nd,:) = 0

! Make dphi(d1,d2)/dd1=0 for d1=0 and dphi(d1,d2)/dd2=0 for d2=0
  dphidd1(1,:) = 0
  dphidd2(:,1) = 0
  d2phidd1dd2(:,1) = 0
  d2phidd1dd2(1,:) = 0

#ifdef DEBUG_XC
! Print values and derivatives for debugging
!  print'(a,/,2a10,4a15)', &
!   'set_phi_table:', 'd1', 'd2', 'phi', 'dphi/dd1', 'dphi/dd2', 'd2phi/dd1dd2'
!  do id1 = 1,nd
!    do id2 = id1,nd
!      d1 = dmesh(id1)
!      d2 = dmesh(id2)
!      print'(2f10.6,4e15.6)', d1, d2, phi(id1,id2), &
!        dphidd1(id1,id2), dphidd2(id1,id2), d2phidd1dd2(id1,id2)
!    end do
!  end do
#endif /* DEBUG_XC */

! Set up bicubic interpolation coefficients
  call bcucof( nd, nd, dmesh, dmesh, phi, dphidd1, dphidd2, d2phidd1dd2, &
               phi_table )

! Save phi_table in file
  open( unit=1, file='vdw_kernel.table', form='unformatted' )
  write(1) nd
  write(1) dmesh
  write(1) phi_table
  close( unit=1 )

#ifdef DEBUG_XC
! Save phi_table also in formatted file
  open( unit=1, file='vdw_kernel.table.formatted', form='formatted' )
  write(1,*) nd
  write(1,*) dmesh
  do id2 = 1,nd
    do id1 = 1,nd
      write(1,*) phi_table(:,:,id1,id2)
    end do
  end do
  close( unit=1 )
#endif /* DEBUG_XC */

! Mark table as set
  phi_table_set = .true.

end subroutine set_phi_table

! -----------------------------------------------------------------------------

subroutine set_qmesh()

! Sets mesh of q values

  implicit none
  integer :: nmesh

  call set_mesh( mq, xmax=qcut, dxndx1=dqmaxdqmin )
  call get_mesh( mq, nmesh, qmesh )
  qmesh_set = .true.

#ifdef DEBUG_XC
  write(udebug,'(/,a,/,(10f8.4))') 'vdw:set_qmesh: qmesh =', qmesh
#endif /* DEBUG_XC */

end subroutine set_qmesh

! -----------------------------------------------------------------------------

subroutine vdw_decusp( nspin, rhos, grhos, eps, dedrho, dedgrho )

! Finds the local energy correction due to the softening of the VdW kernel
! cusp, defined as DEcusp(rho,grad_rho) = 
!   (1/2) * rho**2 * Int 4*pi*r**2*dr * (phi(q1*r,q2*r) - phi_soft(q1*r,q2*r))
! where q1=q2=q0(rho,grad_rho). Notice that grad_rho is included in the value 
! of q0 at the origin but not in the change of rho(r) in the integrand.

  implicit none
  integer, intent(in) :: nspin            ! Number of spin components
  real(dp),intent(in) :: rhos(nspin)      ! Electron spin density
  real(dp),intent(in) :: grhos(3,nspin)   ! Spin density gradient
  real(dp),intent(out):: eps              ! Energy correction per electron
  real(dp),intent(out):: dedrho(nspin)    ! d(rho*eps)/d(rho)
  real(dp),intent(out):: dedgrho(3,nspin) ! d(rho*eps)/d(grad_rho)

  logical, save:: initialized=.false.
  real(dp),save:: table(mq,mq)
  integer :: iq1, iq2, ir, is, ix, ns
  real(dp):: dptpdq, dpdq(mq), dqdrho, dqdgrho(3), grho(3), p(mq), &
             ptp, phi22, phi(mq,mq), pi, pt(mq), q, r, rho

  ! Trap VV-version exception
  if (vdw_author=='VV') then
    eps = vv_vdw_beta()
    dedrho = eps
    dedgrho = 0
    return
  end if

  if (.not.initialized) then

    if (.not.phi_table_set) call set_phi_table()
    if (.not.kcut_set) stop 'vdw_decusp: ERROR: kcut has not been set'

    pi = acos(-1._dp)
    table = 0
    do ir = 1,nr
      r = ir * dr
      call phiofr( r, phi )
      do iq2 = 1,mq
        do iq1 = 1,iq2
          table(iq1,iq2) = table(iq1,iq2) - 2*pi*dr * r**2 * phi(iq1,iq2)
        end do ! iq1
      end do ! iq2
    end do

    do iq2 = 2,mq
      do iq1 = 1,iq2-1
        table(iq2,iq1) = table(iq1,iq2)
      end do
    end do

    initialized = .true.

  end if ! (.not.initialized)

  ns = min(nspin,2)
  rho = sum(rhos(1:ns))
  do ix = 1,3
    grho(ix) = sum(grhos(ix,1:ns))
  end do

  call qofrho( rho, grho, q, dqdrho, dqdgrho )
  call pofq( q, p, dpdq )

  pt = matmul(p,table)
  ptp = sum( pt * p )
  dptpdq = 2 * sum( pt * dpdq )
  eps = rho * ptp
  dedrho(:) = 2 * rho * ptp + rho**2 * dptpdq * dqdrho
  do ix = 1,3
    dedgrho(ix,:) = rho**2 * dptpdq * dqdgrho(ix)
  end do

end subroutine vdw_decusp

!-----------------------------------------------------------------------------

subroutine vdw_get_qmesh( n, q )

! Returns size and values of q-mesh

  implicit none
  integer,          intent(out) :: n     ! Number of qmesh points
  real(dp),optional,intent(out) :: q(:)  ! Values of qmesh points
  integer:: nmax, nkf, nkg

  ! Trap VV-version exception
  if (vdw_author=='VV') then
    call vv_vdw_get_kmesh( nkf, nkg )
    n = nkf*nkg
    if (present(q)) &
      call die('vdw_get_qmesh: ERROR q-mesh not available for author=VV')
    return
  end if

  if (.not.qmesh_set) call set_qmesh()
  n = mq
  if (present(q)) then
    nmax = max( mq, size(q) )
    q(1:nmax) = qmesh(1:nmax)
  end if
end subroutine vdw_get_qmesh

! -----------------------------------------------------------------------------

subroutine vdw_localxc( iRel, nSpin, D, GD, epsX, epsC, &
                        dEXdD, dECdD, dEXdGD, dECdGD )

! Finds the (semi)local part of the correlation energy density and its 
! derivatives, using the GGA functional apropriate for the previously-set
! vdW functional flavour

  implicit none
  integer, intent(in) :: iRel            ! Relativistic exchange? 0=no, 1=yes
  integer, intent(in) :: nSpin           ! Number of spin components
  real(dp),intent(in) :: D(nSpin)        ! Local electron (spin) density
  real(dp),intent(in) :: GD(3,nSpin)     ! Gradient of electron density
  real(dp),intent(out):: epsX            ! Local exchange energy per electron
  real(dp),intent(out):: epsC            ! Local correlation energy per electron
  real(dp),intent(out):: dEXdD(nSpin)    ! dEx/dDens, Ex=Int(dens*epsX)
  real(dp),intent(out):: dECdD(nSpin)    ! dEc/dDens, Ec=Int(dens*epsC)
  real(dp),intent(out):: dEXdGD(3,nSpin) ! dEx/dGrad(Dens)
  real(dp),intent(out):: dECdGD(3,nSpin) ! dEc/dGrad(Dens)

! Internal variables and arrays
  real(dp):: epsAux, dEdDaux(nSpin), dEdGDaux(3,nSpin)
  real(dp):: dVXdD(nSpin,nSpin), dVCdD(nSpin,nSpin)

! Initialize output
  epsX = 0
  epsC = 0
  dEXdD = 0
  dECdD = 0
  dEXdGD = 0
  dECdGD = 0

! Call the appropriate GGA functional.
! Use aux arrays to avoid overwritting the wrong ones
  if (vdw_author=='DRSLL' .or. vdw_author=='drsll' .or. &
      vdw_author=='DF1' .or. vdw_author=='df1') then
    ! Dion et al, PRL 92, 246401 (2004)
    ! GGA exchange and LDA correlation (we choose PW92 for the later)
    call GGAxc( 'revPBE', iRel, nSpin, D, GD, epsX, epsAux, &
                 dEXdD, dEdDaux, dEXdGD, dEdGDaux )
    call LDAxc( 'PW92', iRel, nSpin, D, epsAux, epsC,  &
                dEdDaux, dECdD, dVXdD, dVCdD )
  else if (vdw_author=='LMKLL' .or. vdw_author=='lmkll' .or. &
           vdw_author=='DF2' .or. vdw_author=='df2') then
    ! Lee et al, PRB 82, 081101 (2010)
    call GGAxc( 'PW86R', iRel, nSpin, D, GD, epsX, epsAux, &
                 dEXdD, dEdDaux, dEXdGD, dEdGDaux )
    call LDAxc( 'PW92', iRel, nSpin, D, epsAux, epsC,  &
                dEdDaux, dECdD, dVXdD, dVCdD )
  else if (vdw_author=='KBM' .or. vdw_author=='kbm') then
    ! optB88-vdW of Klimes et al, JPCM 22, 022201 (2009)
    call GGAxc( 'B88KBM', iRel, nSpin, D, GD, epsX, epsAux, &
                 dEXdD, dEdDaux, dEXdGD, dEdGDaux )
    call LDAxc( 'PW92', iRel, nSpin, D, epsAux, epsC,  &
                dEdDaux, dECdD, dVXdD, dVCdD )
  else if (vdw_author=='C09' .or. vdw_author=='c09') then
    ! C09x-vdWc of Cooper, PRB 81, 161104(R) (2010)
    call GGAxc( 'C09', iRel, nSpin, D, GD, epsX, epsAux, &
                 dEXdD, dEdDaux, dEXdGD, dEdGDaux )
    call LDAxc( 'PW92', iRel, nSpin, D, epsAux, epsC,  &
                dEdDaux, dECdD, dVXdD, dVCdD )
  else if (vdw_author=='BH' .or. vdw_author=='bh') then
    ! Berland and Hyldgaard, PRB 89, 035412 (2014)
    call GGAxc( 'BH', iRel, nSpin, D, GD, epsX, epsAux, &
                dEXdD, dEdDaux, dEXdGD, dEdGDaux )
    call LDAxc( 'PW92', iRel, nSpin, D, epsAux, epsC,  &
                dEdDaux, dECdD, dVXdD, dVCdD )
  else if (vdw_author=='VV' .or. vdw_author=='vv') then
    ! Vydrov and VanVoorhis, JCP 133, 244103 (2010)
    ! GGA for both exchange and correlation, but with different flavours
    call GGAxc( 'PW86R', iRel, nSpin, D, GD, epsX, epsAux, &
                dEXdD, dEdDaux, dEXdGD, dEdGDaux )
    call GGAxc( 'PBE', iRel, nSpin, D, GD, epsAux, epsC, &
                dEdDaux, dECdD, dEdGDaux, dECdGD )
  else
    stop 'vdw_exchng ERROR: unknown author'
  end if

end subroutine vdw_localxc

!-----------------------------------------------------------------------------

subroutine vdw_phi( k, phi, dphidk )

! Finds by interpolation phi(q1,q2,k) (Fourier transform of phi(q1,q2,r)) 
! for all values of q1 and q2 in qmesh. If the interpolation table does not
! exist, it is calculated in the first call to vdw_phi. It requires a 
! previous call to vdw_set_kc to set qmesh.

  implicit none
  real(dp),intent(in) :: k            ! Modulus of actual k vector
  real(dp),intent(out):: phi(:,:)     ! phi(q1,q2,k) at given k
                                      ! for all q1,q2 in qmesh
  real(dp),intent(out):: dphidk(:,:)  ! dphi(q1,q2,k)/dk at given k

  integer :: ik, iq1, iq2
  real(dp):: a, a2, a3, b, b2, b3

  ! Trap VV-version exception
  if (vdw_author=='VV') then
    call vv_vdw_phi( k, phi, dphidk )
    return
  end if

  if (.not.kcut_set) stop 'vdw_phi: ERROR: kcut must be previously set'

! Check argument sizes
  if (size(phi,1)<mq .or. size(phi,2)<mq) &
    stop 'vdw_phi: ERROR: size(phi) too small'

! Find phi values at point k
  if (k >= kcut) then
    phi(:,:) = 0
  else
    ! Expand interpolation inline since this is the hottest point in VdW
    ik = k/dk
    a = ((ik+1)*dk-k)/dk
    b = 1 - a
    a2 = (3*a**2 -1) * dk / 6
    b2 = (3*b**2 -1) * dk / 6
    a3 = (a**3 - a) * dk**2 / 6
    b3 = (b**3 - b) * dk**2 / 6
    do iq2 = 1,mq
      do iq1 = 1,iq2
!        call splint( dk, phik(:,iq1,iq2), d2phidk2(:,iq1,iq2), nk+1, k, &
!                     phi(iq1,iq2), dphidk(iq1,iq2), pr )
        phi(iq1,iq2) = a*phik(ik,iq1,iq2) + b*phik(ik+1,iq1,iq2) &
                + a3*d2phidk2(ik,iq1,iq2) + b3*d2phidk2(ik+1,iq1,iq2)
        dphidk(iq1,iq2) = (-phik(ik,iq1,iq2) + phik(ik+1,iq1,iq2)) / dk &
                   - a2*d2phidk2(ik,iq1,iq2) + b2*d2phidk2(ik+1,iq1,iq2)
        phi(iq2,iq1) = phi(iq1,iq2)
        dphidk(iq2,iq1) = dphidk(iq1,iq2)
      end do
    end do
  end if

end subroutine vdw_phi

!-----------------------------------------------------------------------------

subroutine vdw_set_author( author )

! Sets the functional flavour (author initials) and subsequent parameters

  implicit none
  character(len=*),intent(in):: author ! Functional flavour 
                                   ! ('DRSLL'|'LMKLL'|'KBM'|'C09'|'BH'|'VV')

  if (author=='DRSLL') then
    ! Dion et al, PRL 92, 246401 (2004)
    zab = -0.8491_dp
  else if (author=='LMKLL') then
    ! Lee et al, PRB 82, 081101 (2010)
    zab = -1.887_dp
  else if (author=='KBM') then
    ! optB88-vdW of Klimes et al, JPCM 22, 022201 (2009)
    zab = -0.8491_dp
  else if (author=='C09') then
    ! Cooper, PRB 81, 161104(R) (2010)
    zab = -0.8491_dp
  else if (author=='BH') then
    ! Berland and Hyldgaard, PRB 89, 035412 (2014)
    zab = -0.8491_dp
  else if (author=='VV') then
    ! Vydrov and Van Voorhis, JCP 133, 244103 (2010)
  else
    stop 'vdw_set_author: ERROR: author not known'
  end if
  vdw_author = author

end subroutine vdw_set_author

!-----------------------------------------------------------------------------

subroutine vdw_set_kcut( kc )

! Sets the reciprocal planewave cutoff kc of the integration grid, and finds 
! the interpolation table to be used by vdw_phi to obtain the vdW kernel phi
! at the reciprocal mesh wavevectors.

  implicit none
  real(dp),intent(in):: kc  ! Planewave cutoff: k>kcut => phi=0

  character(len=*),parameter:: myName = 'vdw_set_kcut '
  integer :: ik, iq1, iq2, ir, nrs
  real(dp):: dphids, dphidk0, dphidkmax, dphidr0, dphidrmax, dqdq0, &
             k, kmax, phi(0:nr), phi0, phi2, phis, pi, q1, q2, r(0:nr), rs

  ! Trap VV-version exception
  if (vdw_author=='VV') then
    call vv_vdw_set_kcut( kc )
    return
  end if

  if (kc == kcut) return   ! Work alredy done
  if (.not.qmesh_set) call set_qmesh()

  ! Allocate arrays
  call re_alloc( phir,     0,nr, 1,mq, 1,mq, myName//'phir' )
  call re_alloc( phik,     0,nr, 1,mq, 1,mq, myName//'phik' )
  call re_alloc( d2phidr2, 0,nr, 1,mq, 1,mq, myName//'d2phidr2' )
  call re_alloc( d2phidk2, 0,nr, 1,mq, 1,mq, myName//'d2phidk2' )

  pi = acos(-1._dp)
  dr = rcut / nr
  dk = pi / rcut
  kmax = pi / dr
  nk = int(kc/dk) + 1
  nrs = nint(rsoft/dr)
  rs = nrs * dr
  if (nk>nr) stop 'vdw_set_kcut: ERROR: nk>nr'

#ifdef DEBUG_XC
  write(udebug,'(a,5f8.3)') 'vdw_set_kcut: dcut,qcut,rcut,kcut,kmax=', &
    dcut, qcut, rcut, kc, kmax
#endif /* DEBUG_XC */

  ! For each pair of values q1 and q2
  do iq2 = 1,mq
    do iq1 = 1,iq2
!      print*, 'vdw_set_kcut: iq1,iq2=', iq1, iq2

!     Saturated q values
!      q1 = qmesh(iq1)
!      q2 = qmesh(iq2)

!     Find original (unsaturated) q values
      call saturate_inverse( qmesh(iq1), qcut, q1, dqdq0 )
      call saturate_inverse( qmesh(iq2), qcut, q2, dqdq0 )

      ! Find kernel in real space
      do ir = 0,nr
        r(ir) = ir * dr
        phir(ir,iq1,iq2) = phi_interp( q1*r(ir), q2*r(ir) )
      end do
      phi(:) = phir(:,iq1,iq2)

      ! Change kernel near origin to a parabola matching at rs
      if (nrs>0) then
        phis = phi(nrs)
        dphids = (phi(nrs+1) - phi(nrs-1)) / (2*dr)
        phi0 = phis - dphids * rs/2
        phi2 = dphids / (2*rs)
        do ir = 0,nrs
          phir(ir,iq1,iq2) = phi0 + phi2 * r(ir)**2
        end do
      end if ! (nrs>0)

      ! Kill kernel smoothly at r=rcut
      do ir = 0,nr
        phir(ir,iq1,iq2) = phir(ir,iq1,iq2) * cutoff( r(ir)/rcut )
      end do

      ! Optimized filter (warning: inaccurate for very large kcut*rcut)
!      call filter( 0, nr+1, r(:), phir(:,iq1,iq2), kc, 0 )

      ! Find kernel in reciprocal space
      call radfft( 0, nr, rcut, phir(:,iq1,iq2), phik(:,iq1,iq2) )
      phik(:,iq1,iq2) = phik(:,iq1,iq2) * (2*pi)**1.5_dp

      ! Filter out above kcut
      phik(nk:nr,iq1,iq2) = 0

      ! Soft filter below kcut
      do ik = 1,nk
        k = ik * dk
        phik(ik,iq1,iq2) = phik(ik,iq1,iq2) * cutoff(k/kc)
      end do

      ! Find filtered kernel in real space
      call radfft( 0, nr, kmax, phik(:,iq1,iq2), phir(:,iq1,iq2) )
      phir(:,iq1,iq2) = phir(:,iq1,iq2) / (2*pi)**1.5_dp

      ! Set up spline interpolation tables
      dphidr0 = 0
      dphidrmax = 0
      dphidk0 = 0
      dphidkmax = 0
      call spline( dr, phir(:,iq1,iq2), nr+1, dphidr0, dphidrmax, &
                   d2phidr2(:,iq1,iq2) )
      call spline( dk, phik(:,iq1,iq2), nk+1, dphidk0, dphidkmax, &
                   d2phidk2(:,iq1,iq2) )

      ! Fill symmetric elements
      phir(:,iq2,iq1) = phir(:,iq1,iq2)
      phik(:,iq2,iq1) = phik(:,iq1,iq2)
      d2phidr2(:,iq2,iq1) = d2phidr2(:,iq1,iq2)
      d2phidk2(:,iq2,iq1) = d2phidk2(:,iq1,iq2)

!      if (.false. .and. iq1==iq2) then
!        print*, 'vdw_set_kcut: iq1,iq2=', iq1, iq2
!        call window( 0._dp, 5._dp, -1._dp, 4._dp, 0 )
!        call axes( 0._dp, 1._dp, 0._dp, 1._dp )
!        call plot( nr+1, r, phi, phir(:,iq1,iq2) )
!        call window( 0._dp, 10._dp, -0.05_dp, 0.15_dp, 0 )
!        call axes( 0._dp, 1._dp, 0._dp, 0.05_dp )
!        call plot( nr+1, r, q1*q2*r**2*phi, q1*q2*r**2*phir(:,iq1,iq2) )
!        call show()
!      end if

    end do ! iq1
  end do ! iq2

!  print'(a,/,(2i6,f12.6))', 'vdw_set_kcut: iq1, iq2, phir(0,iq1,iq2) =', &
!    ((iq1,iq2,phir(0,iq1,iq2),iq1=2,iq2),iq2=2,mq)

  kcut = kc
  kcut_set = .true.

end subroutine vdw_set_kcut

!-----------------------------------------------------------------------------

subroutine vdw_theta( nspin, rhos, grhos, theta, dtdrho, dtdgrho )

! Finds the value and derivatives of theta_i(rho,grad_rho) = rho*p_i(q0), 
! where q0(rho,grad_rho) is the local wavevector defined in eqs.(11-12) 
! of Dion et al. p_i(q0) are the cubic polynomials such that
!     y(q0) = Sum_i p_i(q0) * y_i
! is the cubic spline interpolation at q0 of (any) function y(q) with
! values y_i at mesh points qmesh_i

  implicit none
  integer, intent(in) :: nspin               ! Number of spin components
  real(dp),intent(in) :: rhos(nspin)         ! Electron spin density
  real(dp),intent(in) :: grhos(3,nspin)      ! Spin density gradient
  real(dp),intent(out):: theta(:)            ! Expansion of rho*q in qmesh
  real(dp),intent(out):: dtdrho(:,:)         ! dtheta(iq)/drhos
  real(dp),intent(out):: dtdgrho(:,:,:)      ! dtheta(iq)/dgrhos(ix)
!  real(dp),intent(out):: theta(mq)           ! Expansion of rho*q in qmesh
!  real(dp),intent(out):: dtdrho(mq,nspin)    ! dtheta(iq)/drhos
!  real(dp),intent(out):: dtdgrho(3,mq,nspin) ! dtheta(iq)/dgrhos(ix)

  integer :: is, ix, ns
  real(dp):: rho, grho(3), dpdq(mq), dqdrho, dqdgrho(3), p(mq), q

  ! Trap VV-version exception
  if (vdw_author=='VV') then
    call vv_vdw_theta( nspin, rhos, grhos, theta, dtdrho, dtdgrho )
    return
  end if

  ns = min(nspin,2)
  rho = sum(rhos(1:ns))
  do ix = 1,3
    grho(ix) = sum(grhos(ix,1:ns))
  end do

  call qofrho( rho, grho, q, dqdrho, dqdgrho )
  call pofq( q, p, dpdq )

  theta(1:mq) = rho * p(1:mq)
  dtdrho(:,:) = 0
  dtdgrho(:,:,:) = 0
  do is = 1,ns
    dtdrho(1:mq,is) = p(1:mq) + rho * dpdq(1:mq) * dqdrho
    do ix = 1,3
      dtdgrho(ix,1:mq,is) = rho * dpdq(1:mq) * dqdgrho(ix)
    end do
  end do

end subroutine vdw_theta

END MODULE m_vdwxc
