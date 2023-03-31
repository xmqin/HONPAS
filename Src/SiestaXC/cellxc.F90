!!@LICENSE
!
! *******************************************************************
! subroutine cellXC( irel, cell, nMesh, lb1, ub1, lb2, ub2, lb3, ub3, 
!    .               nSpin, dens, Ex, Ec, Dx, Dc, stress, Vxc, dVxcdD )
! *******************************************************************
! Finds total exchange-correlation energy and potential in a
!   periodic cell.
! This version implements the Local (spin) Density Approximation and
!   the Generalized-Gradient-Aproximation with the 'explicit mesh 
!   functional' approach of White & Bird, PRB 50, 4954 (1994).
! Gradients are 'defined' by numerical derivatives, using 2*nn+1 mesh
!   points, where nn is a parameter defined below
! Ref: L.C.Balbas et al, PRB 64, 165110 (2001)
! Wrtten by J.M.Soler using algorithms developed by 
!   L.C.Balbas, J.L.Martins and J.M.Soler, Dec.1996 - Aug.1997
! Parallel version written by J.Gale. June 1999.
! Argument dVxcdD added by J.Junquera. November 2000.
! Adapted for multiple functionals in the same run by J.Gale 2005
! Van der Waals functional added by J.M.Soler, Jan.2008, as explained in
!   G.Roman-Perez and J.M.Soler, PRL 103, 096102 (2009)
! ************************* INPUT ***********************************
! integer  irel        : Relativistic exchange? (0=>no, 1=>yes)
! real(dp) cell(3,3)   : Unit cell vectors cell(ixyz,ivector)
! integer  nMesh(3)    : Total mesh divisions of each cell vector
! integer  lb1,lb2,lb3 : Lower bounds of arrays dens, Vxc, dVxcdD
! integer  ub1,ub2,ub3 : Upper bounds of arrays dens, Vxc, dVxcdD
! integer  nSpin       : nSpin=1 => unpolarized; nSpin=2 => polarized;
!                        nSpin>2 => non-collinear polarization
! real(grid_p) dens(lb1:ub1,lb2:ub2,lb3:ub3,nSpin) : Total (nSpin=1) or 
!                        spin (nSpin=2) electron density at mesh points
!                        For non-collinear polarization, the density
!                        matrix is given by: dens(1)=D11, dens(2)=D22,
!                        dens(3)=Real(D12), dens(4)=Im(D12)
! ************************* OUTPUT **********************************
! real(dp) Ex             : Total exchange energy per unit cell
! real(dp) Ec             : Total correlation energy per unit cell
! real(dp) Dx             : IntegralOf( rho * (eps_x - v_x) ) in unit cell
! real(dp) Dc             : IntegralOf( rho * (eps_c - v_c) ) in unit cell
! real(dp) stress(3,3)    : xc contribution to the stress tensor, in unit
!                           cell, assuming constant density (not charge),
!                           i.e. r->r' => rho'(r') = rho(r)
!                           For plane-wave and grid (finite diff) basis
!                           sets, density rescaling gives an extra term
!                           (not included) (Dx+Dc-Ex-Ec)/cell_volume for
!                           the diagonal elements of stress. For other
!                           basis sets, the extra term is, in general:
!                           IntegralOf(v_xc * d_rho/d_strain) / cell_vol
! real(grid_p) Vxc(lb1:ub1,lb2:ub2,lb3:ub3,nSpin) : (Spin) xc potential
! ************************* OPTIONAL OUTPUT *************************
! real(grid_p) dVxcdD(lb1:ub1,lb2:ub2,lb3:ub3,nSpin*nSpin) : Derivatives
!                           of xc potential respect to charge density
!                           Available only for LDA
! ************************ UNITS ************************************
! Distances in atomic units (Bohr).
! Densities in atomic units (electrons/Bohr**3)
! Energy unit depending of parameter EUnit below
! Stress in EUnit/Bohr**3
! ************************ USAGE ************************************
! With the prototype module xcmod below, you must make a previous call
!     CALL setXC( nFunc, XCfunc, XCauth, XCweightX, XCweightC )
! before calling cellXC for the first time
!
! A typical serial program call is:
!
!   use precision, only: dp, grid_p
!   use xcmod,     only: setXC
!   integer  :: nMesh(3), nSpin
!   real(dp) :: cell(3,3), Dc, Dx, Ec, Ex, stress(3,3), 
!   real(grid_p),allocatable :: dens(:,:,:,:), Vxc(:,:,:,:)
!     Find nSpin, cell(:,:), and nMesh(:)
!   allocate( dens(nMesh(1),nMesh(2),nMesh(3),nSpin), &
!              Vxc(nMesh(1),nMesh(2),nMesh(3),nSpin)) )
!     Find dens(:,:,:,:) at all mesh points
!   call setXC( 1, (/'GGA'/), (/'PBE'/), (/1._dp/), (/1._dp/) )
!   call cellXC( 0, cell, nMesh, 1,nMesh(1), 1,nMesh(2), 1,nMesh(3), &
!                nSpin, dens, Ex, Ex, Dx, Dc, stress, Vxc )
!
! A typical parallel program call is:
!
!   use precision, only: dp, grid_p
!   use xcmod,     only: setXC
!   integer  :: iSpin, myBox(2,3), nMesh(3), nSpin
!   real(dp) :: cell(3,3), Dc, Dx, Ec, Ex, stress(3,3), 
!   real(grid_p),allocatable :: dens(:,:,:,:), Vxc(:,:,:,:)
!     Find nSpin, cell(:,:), nMesh(:), and myBox(:,:)
!   allocate( dens(myBox(1,1):myBox(2,1),        &
!                  myBox(1,2):myBox(2,2),        &
!                  myBox(1,3):myBox(2,3),nSpin), &
!              Vxc(myBox(1,1):myBox(2,1),        &
!                  myBox(1,2):myBox(2,2),        &
!                  myBox(1,3):myBox(2,3),nSpin) )
!   do i3 = myBox(1,3),myBox(2,3)
!   do i2 = myBox(1,2),myBox(2,2)
!   do i1 = myBox(1,1),myBox(2,1)
!     do iSpin = 1,nSpin
!       dens(i1,i2,i3,iSpin) = (spin)density at point (i1,i2,i3)
!     end do
!   end do
!   end do
!   end do
!   call setXC( 1, (/'GGA'/), (/'PBE'/), (/1._dp/), (/1._dp/) )
!   call cellXC( 0, cell, nMesh, myBox(1,1), myBox(2,1), &
!                                myBox(1,2), myBox(2,2), &
!                                myBox(1,3), myBox(2,3), &
!                nSpin, dens, Ex, Ex, Dx, Dc, stress, Vxc )
!
! IMPORTANT: arrays dens, Vxc, and dVxcdD may be alternatively 
! allocated and initialized with indexes beginning in 0 or 1, 
! or even use a single-index array, e.g.:
!   real(grid_p),allocatable :: dens(:,:), Vxc(:,:)
!   myMesh(1:3) = myBox(2,:) - myBox(1,:) + 1
!   myPoints = myMesh(1)*myMesh(2)*myMesh(3)
!   allocate( dens(myPoints,nSpin), Vxc(myPoints,nSpin) )
!   iPoint = 0
!   do i3 = myBox(1,3),myBox(2,3)
!   do i2 = myBox(1,2),myBox(2,2)
!   do i1 = myBox(1,1),myBox(2,1)
!     iPoint = iPoint + 1
!     do iSpin = 1,nSpin
!       dens(iPoint,iSpin) = (spin)density at point (i1,i2,i3)
!     end do
!   end do
!   end do
!   end do
! but the call to cellXC must still be as given above, i.e. the
! arguments lb1,ub1,... must be the lower and upper bounds of the
! mesh box stored by each processor (not the allocated array bounds).
! However, the arrays size MUST be (ub1-lb1+1)*(ub2-lb2+1)*(ub3-lb3+1)
! The processor mesh boxes must not overlap, and they must cover the 
! entire unit cell mesh. This is NOT checked inside cellXC.
!
! ********* BEHAVIOUR ***********************************************
! - Stops and prints a warning if functl is not one of LDA, GGA, or VDW
! - The output values of Ex, Ec, Dx, Dc, and stress, are integrals over
!   the whole unit cell, not over the mesh box of the local processor
! - Since the exchange and correlation part is usually a small fraction
!   of a typical electronic structure calculation, this routine has
!   been coded with emphasis on simplicity and functionality, not in
!   efficiency.
! ********* DEPENDENCIES ********************************************
! Routines called: 
!   GGAXC, LDAXC, meshKcut, RECLAT, TIMER, VOLCEL
! Modules used in serial:
!   precision : defines parameters 'dp' and 'grid_p' for real kinds
!   xcmod     : sets and gets the xc functional(s) used
!   vdW       : povides routines for the Van der Waals functional
!   mesh3D    : provides routines to handle mesh arrays distributed
!               among processors
!   sys       : provides the stopping subroutine 'die'
! Additional modules used in parallel:
!   mpi_siesta
! ********* PRECISION MODULE PROTOTYPE ******************************
! The following module will set the required real types
!
! module precision
!   integer,parameter:: dp     = kind(1.d0)
!   integer,parameter:: grid_p = kind(1.d0)
! end module precision
!
! ********* XCMOD MODULE PROTOTYPE **********************************
! The following module will set the xc functional(s) directly, when
! calling routine setxc (rather than reading them from the datafile,
! as done in siesta):
!
! module xcmod
!   use precision, only: dp
!   implicit none
!   integer, parameter :: maxFunc = 10
!   integer,           save :: nXCfunc
!   character(len=10), save :: XCauth(MaxFunc), XCfunc(MaxFunc)
!   real(dp),          save :: XCweightX(MaxFunc), XCweightC(MaxFunc)
! contains
!   subroutine setXC( n, func, auth, wx, wc )
!     implicit none
!     integer,         intent(in):: n       ! Number of functionals
!     character(len=*),intent(in):: func(n) ! Functional name labels
!     character(len=*),intent(in):: auth(n) ! Functional author labels
!     real(dp),        intent(in):: wx(n)   ! Functl weights for exchng
!     real(dp),        intent(in):: wc(n)   ! Functl weights for correl
!     nXCfunc = n
!     XCfunc(1:n) = func(1:n)
!     XCauth(1:n) = auth(1:n)
!     XCweightX(1:n) = wx(1:n)
!     XCweightC(1:n) = wc(1:n)
!   end subroutine setXC
! end module xcmod
!
! Sample usage:
!   use precision, only: dp
!   use xcmod, only: setxc
!   call setXC( 1, (/'GGA'/), (/'PBE'/), (/1._dp/), (/1._dp/) )
!   call cellXC( ... )
!
! Allowed functional/author values:
! XCfunc: 
!   'LDA' or 'LSD' => Local density approximation
!            'GGA' => Generalized gradients approx.
!            'VDW' => Van der Waals functional
! XCauth:
!     'CA' or 'PZ' => LSD Perdew & Zunger, PRB 23, 5075 (1981)
!           'PW91' => GGA Perdew & Wang, JCP, 100, 1290 (1994) 
!           'PW92' => LSD Perdew & Wang, PRB, 45, 13244 (1992). This is
!                     the local density limit of the next:
!            'PBE' => GGA Perdew, Burke & Ernzerhof, PRL 77, 3865 (1996)
!           'RPBE' => GGA Hammer, Hansen & Norskov, PRB 59, 7413 (1999)
!         'revPBE' => GGA Zhang & Yang, PRL 80,890(1998)
!            'LYP' => GGA Becke-Lee-Yang-Parr (see subroutine blypxc)
!             'WC' => GGA Wu-Cohen (see subroutine wcxc)
!         'PBESOL' => GGA Perdew et al, PRL, 100, 136406 (2008)
!           'AM05' => GGA Mattsson & Armiento, PRB, 79, 155101 (2009)
!    'PBE(JsJrLO)' => GGA Reparametrizations of the PBE functional by
!   'PBE(JsJrHEG)' => GGA   L.S.Pedroza et al, PRB 79, 201106 (2009) and
!    'PBE(GcGxLO)' => GGA   M.M.Odashima et al, JCTC 5, 798 (2009)
!   'PBE(GcGxHEG)' => GGA using 4 different combinations of criteria
!          'DRSLL' => VDW Dion et al, PRL 92, 246401 (2004)
!          'LMKLL' => VDW K.Lee et al, PRB 82, 081101 (2010)
!            'KBM' => VDW optB88-vdW of J.Klimes et al, JPCM 22, 022201 (2010)
!            'C09' => VDW V.R. Cooper, PRB 81, 161104 (2010)
!             'BH' => VDW K. Berland and Per Hyldgaard, PRB 89, 035412 (2014)
!             'VV' => VDW Vydrov-VanVoorhis, JCP 133, 244103 (2010)
! *******************************************************************

MODULE m_cellXC

  implicit none

  PUBLIC:: cellXC  ! Exchange and correlation in a periodic unit cell

  PRIVATE ! Nothing is declared public beyond this point

CONTAINS ! nothing else but public routine cellXC

SUBROUTINE cellXC( irel, cell, nMesh, lb1, ub1, lb2, ub2, lb3, ub3, &
                   nSpin, dens, Ex, Ec, Dx, Dc, stress, Vxc, dVxcdD, &
                   keep_input_distribution )

  ! Module routines
  use mesh3D,  only: addMeshData   ! Accumulates a mesh array
  use alloc,   only: alloc_default ! Sets (re)allocation defaults
  use mesh3D,  only: associateMeshTask ! Associates mesh tasks & distr.
  use mesh3D,  only: copyMeshData  ! Copies a box of a mesh array
  use alloc,   only: de_alloc      ! Deallocates arrays
  use sys,     only: die           ! Termination routine
  use fftr,    only: fftk2r        ! Inverse real Fourier transform
  use fftr,    only: fftr2k        ! Direct real Fourier transform
  use mesh3D,  only: fftMeshDistr  ! Sets/gets distribution for FFTs
  use mesh3D,  only: freeMeshDistr ! Frees a mesh distribution ID
  use moreParallelSubs, only: miscAllReduce ! Adds variables from all processes
  use xcmod,   only: getXC         ! Returns the XC functional to be used
  use m_ggaxc, only: ggaxc         ! General GGA XC routine
  use m_ldaxc, only: ldaxc         ! General LDA XC routine
  use m_chkgmx,only: meshKcut      ! Returns the planewave cutoff of a mesh
  use mesh3D,  only: myMeshBox     ! Returns the mesh box of my processor
  use parallel,only: parallel_init ! Initializes nodes variables
#ifdef DEBUG_XC
  use moreParallelSubs, only: nodeString ! Returns a string with my node number
  use m_vdwxc, only: qofrho        ! For debugging only
#endif /* DEBUG_XC */
  use alloc,   only: re_alloc      ! Reallocates arrays
  use cellsubs,only: reclat        ! Finds reciprocal unit cell vectors
  use mesh3D,  only: sameMeshDistr ! Finds if two mesh distr. are equal
  use mesh3D,  only: setMeshDistr  ! Defines a new mesh distribution
#ifdef DEBUG_XC
  use debugXC, only: setDebugOutputUnit ! Sets udebug variable
#endif /* DEBUG_XC */
  use m_timer, only: timer_get     ! Returns counted times
  use m_timer, only: timer_start   ! Starts counting time
  use m_timer, only: timer_stop    ! Stops counting time
  use m_vdwxc, only: vdw_localxc   ! Local LDA/GGA xc apropriate for vdW flavour
  use m_vdwxc, only: vdw_decusp    ! Cusp correction to VDW energy
  use m_vdwxc, only: vdw_get_qmesh ! Returns q-mesh for VDW integrals
  use m_vdwxc, only: vdw_phi       ! Returns VDW functional kernel
  use m_vdwxc, only: vdw_set_kcut  ! Fixes k-cutoff in VDW integrals
  use m_vdwxc, only: vdw_theta     ! Returns VDW theta function
  use cellsubs,only: volcel        ! Finds volume of unit cell

  ! Module types and variables
  use alloc,     only: allocDefaults ! Derived type for allocation defaults
  use precision, only: dp            ! Double precision real type
  use precision, only: gp=>grid_p    ! Real precision type of mesh arrays
  use parallel,  only: nodes         ! Number of processor nodes
#ifdef DEBUG_XC
  use debugXC,   only: udebug        ! Output file unit for debug info
#endif /* DEBUG_XC */

  implicit none

  ! Argument types and dimensions
  integer, intent(in) :: irel        ! Relativistic exchange? (0=>no, 1=>yes)
  real(dp),intent(in) :: cell(3,3)   ! Unit cell vectors cell(ixyz,ivector)
  integer, intent(in) :: nMesh(3)    ! Total mesh divisions of each cell vector
  integer, intent(in) :: lb1,lb2,lb3 ! Lower bounds of dens, Vxc, and dVxcdD
  integer, intent(in) :: ub1,ub2,ub3 ! Upper bounds of dens, Vxc, and dVxcdD
  integer, intent(in) :: nSpin       ! Number of spin components
  real(gp),intent(in) :: dens(0:ub1-lb1,0:ub2-lb2,0:ub3-lb3,1:nSpin)
                                     ! (spin) Density at mesh points
  real(dp),intent(out):: Ex          ! Total exchange energy in unit cell
  real(dp),intent(out):: Ec          ! Total correlation energy in unit cell
  real(dp),intent(out):: Dx          ! IntegralOf(rho*(eps_x-v_x)) in unit cell
  real(dp),intent(out):: Dc          ! IntegralOf(rho*(eps_c-v_c)) in unit cell
  real(dp),intent(out):: stress(3,3) ! xc contribution to stress in unit cell
  real(gp),intent(out):: Vxc(0:ub1-lb1,0:ub2-lb2,0:ub3-lb3,1:nSpin) 
                                     ! (spin) xc potential
  real(gp),intent(out),optional:: &  ! dVxc/dDens (LDA only)
                         dVxcdD(0:ub1-lb1,0:ub2-lb2,0:ub3-lb3,1:nSpin**2) 
  logical, intent(in),optional::  keep_input_distribution

  ! Fix the order of the numerical derivatives
  ! nn is the number of points used in each coordinate and direction,
  ! i.e. a total of 6*nn neighbour points is used to find the gradients
  integer,parameter  :: nn = 3

  ! Fix energy unit:  Eunit=1.0 => Hartrees,
  !                   Eunit=0.5 => Rydbergs,
  !                   Eunit=0.03674903 => eV
  real(dp),parameter :: Eunit = 0.5_dp

  ! Fix the maximum allowed dispersion in CPU-time among different processors
  real(dp),parameter :: maxUnbalance = 0.10_dp

  ! Fix density threshold below which it will be taken as zero
  real(dp),parameter :: Dmin = 1.0e-15_dp

  ! Fix density threshold below which we make Vxc=0
  real(dp),parameter :: Dcut = 1.0e-9_dp

  ! Fix a minimum value of k vectors to avoid division by zero
  real(dp),parameter :: kmin = 1.0e-15_dp

  ! Fix the maximum number of functionals to be combined
  integer, parameter :: maxFunc = 10

  ! Subroutine name
  character(len=*),parameter :: myName = 'cellXC '
  character(len=*),parameter :: errHead = myName//'ERROR: '

  ! Static internal variables
  integer,save::   &
     io2my=-1,     &! ID of communications task between ioDistr and myDistr
     ioDistr=-1,   &! ID of input mesh distribution
     kDistr=-1,    &! ID of mesh distribution used for FFTs
     left2my(3)=-1,&! ID of commun. task from left-border boxes and myBox
     my2io=-1,     &! ID of commun. task between myDistr and ioDistr
     my2left(3)=-1,&! ID of commun. task from myBox to left-border boxes
     my2rght(3)=-1,&! ID of commun. task from myBox to right-border boxes
     myDistr=-1,   &! ID of mesh distrib. used internally in this routine
     oldMesh(3)=-1,&! Total number of mesh points in previous call
     rght2my(3)=-1  ! ID of commun. task from right-border boxes and myBox
  real(dp),save::  &
     myTime=1,     &! CPU time in this routine and processor in last iteration
     timeAvge=1,   &! Average over processors of CPU time
     timeDisp=huge(1.0_dp)  ! Dispersion over processors of CPU time

  ! Internal pointers for dynamical allocation
  real(dp), pointer:: &
     dphidk(:,:)=>null(), dtdgd(:,:,:)=>null(), dtdd(:,:)=>null(), &
     dudk(:)=>null(), phi(:,:)=>null(), tk(:)=>null(), tr(:)=>null(), &
     ur(:)=>null(), uk(:)=>null()
  real(gp), pointer:: &
     tvac(:)=>null(), tvdw(:,:)=>null(), uvdw(:,:)=>null(), &
     Vj(:)=>null(), workload(:,:,:)=>null()
  real(gp), dimension(:,:,:,:), pointer:: &
     Dleft=>null(), Dleft1=>null(), Dleft2=>null(), Dleft3=>null(), &
     Drght=>null(), Drght1=>null(), Drght2=>null(), Drght3=>null(), &
     myDens=>null(), mydVxcdD=>null(), myVxc=>null(), &
     tq=>null(), uq=>null(), &
     Vleft=>null(), Vleft1=>null(), Vleft2=>null(), Vleft3=>null(), &
     Vrght=>null(), Vrght1=>null(), Vrght2=>null(), Vrght3=>null()
  logical, pointer:: &
     nonempty(:,:,:)=>null()

  ! Other internal variables and arrays
  integer :: &
     boxLeft(2,3), boxRght(2,3), &
     i1, i2, i3, ic, ii1, ii2, ii3, ik, in, &
     ioBox(2,3), ip, iq, is, ix, iuvdw,  &
     jj(3), jn, jp, js, jx, kBox(2,3), kMesh(3), kPoints, ks,  &
     l11, l12, l13, l21, l22, l23, &
     m11, m12, m13, m21, m22, m23, maxPoints, mesh(3), &
     myBox(2,3), myMesh(3), myOldDistr, myPoints, &
     ndSpin, nf, nonemptyPoints, nPoints, nq, ns, nXCfunc, &
     r11, r12, r13, r21, r22, r23
  real(dp):: &
     comTime, D(nSpin), dedk, dEcdD(nSpin), dEcdGD(3,nSpin), &
     dEcidDj, dEcuspdD(nSpin), dEcuspdGD(3,nSpin),  &
     dExdD(nSpin), dExdGD(3,nSpin), dExidDj, &
     dGdM(-nn:nn), dGidFj(3,3,-nn:nn), Dj(nSpin), &
     dMdX(3,3), DV, dVol, Dtot, dXdM(3,3), &
     dVcdD(nSpin*nSpin), dVxdD(nSpin*nSpin), &
     EcuspVDW, Enl, epsC, epsCusp, epsNL, epsX, f1, f2, &  
     GD(3,nSpin), k, kcell(3,3), kcut, kvec(3),  &
     stressVDW(3,3), sumTime, sumTime2, totTime, VDWweightC, volume, &
     XCweightC(maxFunc), XCweightVDW, XCweightX(maxFunc)
#ifdef DEBUG_XC
  integer :: iip, jjp, jq
  real(dp):: rmod, rvec(3)
  integer,save:: myIter=0
#endif /* DEBUG_XC */
  logical :: &
     GGA, GGAfunctl, VDW, VDWfunctl
  character(len=20):: &
     XCauth(maxFunc), XCfunc(maxFunc)
  character(len=80):: &
     errMsg
  type(allocDefaults):: &
     prevAllocDefaults

  logical :: keep_input_distr

#ifdef DEBUG_XC
  ! Variables for debugging
  real(dp):: GDtot(3), q, dqdD, dqdGD(3)
#endif /* DEBUG_XC */

  ! Make sure that variables in parallel module are set
  call parallel_init()

  ! Start time counter
  call timer_start( myName )

  keep_input_distr = .false.
  if (present(keep_input_distribution)) then
     keep_input_distr =   keep_input_distribution
  endif
     
#ifdef DEBUG_XC
  ! Initialize udebug variable
  call setDebugOutputUnit()
#endif /* DEBUG_XC */

  ! Get the functional(s) to be used
  call getXC( nXCfunc, XCfunc, XCauth, XCweightX, XCweightC )

  ! Set routine name for allocations
  call alloc_default( old=prevAllocDefaults, copy=.true., shrink=.true. )

  ! Find number of diagonal spin components
  ndSpin = min( nSpin, 2 )

  ! Find total number of mesh points in unit cell
  nPoints = nMesh(1) * nMesh(2) * nMesh(3)

  ! Find the cell volume and the volume per mesh point
  volume = volcel( cell )
  dVol = volume / nPoints

  ! Set GGA and VDW switches
  GGA = .false.
  VDW = .false.
  XCweightVDW = 0
  do nf = 1,nXCfunc
    if ( XCfunc(nf).eq.'LDA' .or. XCfunc(nf).eq.'lda' .or. &
         XCfunc(nf).eq.'LSD' .or. XCfunc(nf).eq.'lsd' ) then
      cycle ! nf loop
    else if ( XCfunc(nf).eq.'GGA' .or. XCfunc(nf).eq.'gga') then
      GGA = .true.
    else if ( XCfunc(nf).eq.'VDW' .or. XCfunc(nf).eq.'vdw') then
      GGA = .true.
      VDW = .true.
      XCweightVDW = XCweightVDW + XCweightC(nf)
    else
      write(errMsg,*) errHead//'Unknown functional ', XCfunc(nf)
      call die( trim(errMsg) )
    endif
  enddo

  ! Check argument dVxcdD
  if (present(dVxcdD) .and. GGA) &
    call die(errHead//'dVxcdD available only for LDA')

  ! Find my mesh box in I/O distribution of mesh points
  ioBox(1,1)=lb1;  ioBox(1,2)=lb2;  ioBox(1,3)=lb3
  ioBox(2,1)=ub1;  ioBox(2,2)=ub2;  ioBox(2,3)=ub3
  myMesh(:) = ioBox(2,:) - ioBox(1,:) + 1
  myPoints = product(myMesh)

  ! Get ID of the I/O distribution of mesh points
  call setMeshDistr( ioDistr, nMesh, ioBox )

  if (keep_input_distr) then

     myDistr = ioDistr

  else

  ! If nMesh has changed, use input distribution also initially as myDistr
  if (any(nMesh/=oldMesh)) call setMeshDistr( myDistr, nMesh, ioBox )
  oldMesh = nMesh

  ! Find new mesh distribution, if previous iteration was too unbalanced
  if (nodes>1 .and. timeDisp/timeAvge>maxUnbalance) then
    ! Find my node's mesh box using myDistr
    call myMeshBox( nMesh, myDistr, myBox )
    ! Allocate arrays for the density and workload per point in my box
    call re_alloc( myDens, myBox(1,1),myBox(2,1), &
                           myBox(1,2),myBox(2,2), &
                           myBox(1,3),myBox(2,3), 1,ndSpin, myName//'myDens' )
    call re_alloc(workload,myBox(1,1),myBox(2,1), &
                           myBox(1,2),myBox(2,2), &
                           myBox(1,3),myBox(2,3), myName//'workload' )
    ! Copy density to my box 
    call associateMeshTask( io2my, ioDistr, myDistr )
    call copyMeshData( nMesh, ioDistr, dens(:,:,:,1:ndSpin), &
                              myBox, myDens(:,:,:,1:ndSpin), io2my )
!    call copyMeshData( nMesh, ioDistr, dens(:,:,:,1:ndSpin), &
!                              myBox, myDens(:,:,:,1:ndSpin) )
    ! Find initial expected workload (zero if dens=0, one otherwise)
    do i3 = myBox(1,3),myBox(2,3)
    do i2 = myBox(1,2),myBox(2,2)
    do i1 = myBox(1,1),myBox(2,1)
      Dtot = sum( myDens(i1,i2,i3,1:ndSpin) )
      if (Dtot>Dmin) then
         workload(i1,i2,i3) = 1
      else
         workload(i1,i2,i3) = 0
      end if
    end do
    end do
    end do
    ! Modify workload according to previous CPU time
    workload = workload * myTime/timeAvge
!    workload = workload * (myTime/myPoints) / (nodes*timeAvge/nPoints)
    ! Find new distribution
    call setMeshDistr( myDistr, nMesh, wlDistr=myDistr, workload=workload )
    ! Deallocate temporary arrays
    call de_alloc( workload, myName//'workload' )
    call de_alloc( myDens,   myName//'myDens' )
#ifdef DEBUG_XC
    myIter = myIter+1
    call myMeshBox( nMesh, myDistr, myBox )
    write(udebug,'(a,i6,1x,3i5,3(1x,2i5))') &
      myName//'Iter,nMesh,myBox=', myIter, nMesh, myBox
#endif /* DEBUG_XC */
  end if ! (nodes>1 .and. timeDisp/timeAvge>maxUnbalance)

  endif
#ifdef DEBUG_XC
  ! Keep input distribution for the time being
!  myDistr = ioDistr
#endif /* DEBUG_XC */

  ! Find the box of mesh points that belong to this node
  call myMeshBox( nMesh, myDistr, myBox )

  ! Find local number of mesh points along each axis
  myMesh(:) = myBox(2,:) - myBox(1,:) + 1
  myPoints = product(myMesh)

  ! Allocate myDens and myVxc arrays if either:
  ! - The parallel mesh distribution is changed internally (even with LDA)
  ! - We need finite differences (GGA) and have a distributed mesh
  if (.not.sameMeshDistr(myDistr,ioDistr) .or. (GGA .and. myDistr/=0)) then

    m11=myBox(1,1); m12=myBox(1,2); m13=myBox(1,3)
    m21=myBox(2,1); m22=myBox(2,2); m23=myBox(2,3)
    ns = nSpin  ! Just a shorter name
    call re_alloc( myDens, m11,m21, m12,m22, m13,m23, 1,ns, myName//'myDens' )
    call re_alloc( myVxc,  m11,m21, m12,m22, m13,m23, 1,ns, myName//'myVxc'  )
    if (present(dVxcdD)) &
      call re_alloc( mydVxcdD, m11,m21, m12,m22, m13,m23, 1,ns**2, &
                     myName//'mydVxcdD'  )
    ! Allocate arrays for density and potential in neighbor regions
    if (GGA) then
      l11=m11-nn;     l12=m12-nn;     l13=m13-nn
      l21=m11-1;      l22=m12-1;      l23=m13-1
      r11=m21+1;      r12=m22+1;      r13=m23+1
      r21=m21+nn;     r22=m22+nn;     r23=m23+nn
      call re_alloc( Dleft1, l11,l21, m12,m22, m13,m23, 1,ns, myName//'Dleft1' )
      call re_alloc( Drght1, r11,r21, m12,m22, m13,m23, 1,ns, myName//'Drght1' )
      call re_alloc( Dleft2, m11,m21, l12,l22, m13,m23, 1,ns, myName//'Dleft2' )
      call re_alloc( Drght2, m11,m21, r12,r22, m13,m23, 1,ns, myName//'Drght2' )
      call re_alloc( Dleft3, m11,m21, m12,m22, l13,l23, 1,ns, myName//'Dleft3' )
      call re_alloc( Drght3, m11,m21, m12,m22, r13,r23, 1,ns, myName//'Drght3' )
      call re_alloc( Vleft1, l11,l21, m12,m22, m13,m23, 1,ns, myName//'Vleft1' )
      call re_alloc( Vrght1, r11,r21, m12,m22, m13,m23, 1,ns, myName//'Vrght1' )
      call re_alloc( Vleft2, m11,m21, l12,l22, m13,m23, 1,ns, myName//'Vleft2' )
      call re_alloc( Vrght2, m11,m21, r12,r22, m13,m23, 1,ns, myName//'Vrght2' )
      call re_alloc( Vleft3, m11,m21, m12,m22, l13,l23, 1,ns, myName//'Vleft3' )
      call re_alloc( Vrght3, m11,m21, m12,m22, r13,r23, 1,ns, myName//'Vrght3' )
    end if ! (GGA)
  end if ! (myDistr/=ioDistr .or. GGA)

  ! Redistribute dens data
  if (associated(myDens)) then
    call associateMeshTask( io2my, ioDistr, myDistr )
    call copyMeshData( nMesh, ioDistr, dens, myBox, myDens, io2my )
!    call copyMeshData( nMesh, ioDistr, dens, myBox, myDens )
  end if

  ! If mesh arrays are distributed, find density of neighbor regions
  if (GGA .and. myDistr/=0) then ! Distributed Dens data
    ! Find density of neighbor regions
    do ic = 1,3  ! Loop on cell axes
      if (ic==1) then
        Dleft => Dleft1
        Drght => Drght1
      else if (ic==2) then 
        Dleft => Dleft2
        Drght => Drght2
      else ! (ic==3)
        Dleft => Dleft3
        Drght => Drght3
      end if ! (ic==1)
      boxLeft = myBox
      boxLeft(1,ic) = myBox(1,ic)-nn
      boxLeft(2,ic) = myBox(1,ic)-1
      boxRght = myBox
      boxRght(1,ic) = myBox(2,ic)+1
      boxRght(2,ic) = myBox(2,ic)+nn
      call associateMeshTask( my2left(ic), myDistr )
      call associateMeshTask( my2rght(ic), myDistr )
      call copyMeshData( nMesh, myDistr, myDens, boxLeft, Dleft, my2left(ic) )
      call copyMeshData( nMesh, myDistr, myDens, boxRght, Drght, my2rght(ic) )
!      call copyMeshData( nMesh, myDistr, myDens, boxLeft, Dleft )
!      call copyMeshData( nMesh, myDistr, myDens, boxRght, Drght )
    end do ! ic
  end if ! (GGA .and. myDistr/=0)

  ! Find Jacobian matrix dx/dmesh and gradient coeficients dGrad_i/dDens_j
  if (GGA) then

    ! Find mesh unit vectors dx/dmesh
    do ic = 1,3
      dXdM(:,ic) = cell(:,ic) / nMesh(ic)
    enddo

    ! Find reciprocal mesh vectors dmesh/dx
    call reclat( dXdM, dMdX, 0 )

    ! Find weights of numerical derivation from Lagrange interp. formula
    do in = -nn,nn
      f1 = 1.0_dp
      f2 = 1.0_dp
      do jn = -nn,nn
        if (jn/=in .and. jn/=0) f1 = f1 * (0  - jn)
        if (jn/=in)             f2 = f2 * (in - jn)
      enddo
      dGdM(in) = f1 / f2
    enddo
    dGdM(0) = 0.0_dp

    ! Find the weights for the derivative d(gradF(i))/d(F(j)) of
    ! the gradient at point i with respect to the value at point j
    do in = -nn,nn
      do ic = 1,3
        dGidFj(:,ic,in) = dMdX(:,ic) * dGdM(in)
      enddo
    enddo

  endif ! (GGA)

  ! Initialize output
  Ex = 0.0_dp
  Ec = 0.0_dp
  Dx = 0.0_dp
  Dc = 0.0_dp
  stress(:,:) = 0.0_dp
  Vxc(:,:,:,:) = 0.0_gp
  if (present(dVxcdD)) dVxcdD(:,:,:,:) = 0.0_gp

  ! VdW initializations -------------------------------------------------------
  if (VDW) then

    ! Find mask of nonempty points
    call re_alloc( nonempty, 0,myMesh(1)-1, 0,myMesh(2)-1, 0,myMesh(3)-1, &
                   myName//'nonempty' )
    if (associated(myDens)) then
      do i3 = myBox(1,3),myBox(2,3)
      do i2 = myBox(1,2),myBox(2,2)
      do i1 = myBox(1,1),myBox(2,1)
        Dtot = sum( myDens(i1,i2,i3,1:ndSpin) )
        if ( Dtot >= Dmin ) then
          nonempty(i1-myBox(1,1),i2-myBox(1,2),i3-myBox(1,3)) = .true.
        else
          nonempty(i1-myBox(1,1),i2-myBox(1,2),i3-myBox(1,3)) = .false.
        end if
      end do
      end do
      end do
    else
      do i3 = 0, myMesh(3)-1
      do i2 = 0, myMesh(2)-1
      do i1 = 0, myMesh(1)-1
        Dtot = sum( dens(i1,i2,i3,1:ndSpin) )
        if ( Dtot >= Dmin ) then
          nonempty(i1,i2,i3) = .true.
        else
          nonempty(i1,i2,i3) = .false.
        end if
      end do
      end do
      end do
    end if
    nonemptyPoints = count( nonempty )

    ! Find distribution of k-mesh points
    call fftMeshDistr( nMesh, kDistr )

    ! Find my node's box of k-mesh points
    call myMeshBox( nMesh, kDistr, kBox )
    kMesh(:) = kBox(2,:) - kBox(1,:) + 1
    kPoints = product( kMesh )

    ! Find a size large enough for both real and reciprocal points
    maxPoints = max( nonemptyPoints, kPoints )
    mesh = max( myMesh, kMesh )

    ! Allocate VdW arrays
    call vdw_get_qmesh( nq )
    call re_alloc( dudk,  1,nq,                   myName//'dudk' )
    call re_alloc( dtdd,  1,nq, 1,nSpin,          myName//'dtdd' )
    call re_alloc( dtdgd,  1,3,    1,nq, 1,nSpin, myName//'dtdgd')
    call re_alloc( dphidk,1,nq,    1,nq,          myName//'dphidk')
    call re_alloc( phi,   1,nq,    1,nq,          myName//'phi'  )
    call re_alloc( tk,    1,nq,                   myName//'tk'   )
    call re_alloc( tr,    1,nq,                   myName//'tr'   )
    call re_alloc( tvac,  1,nq,                   myName//'tvac' )
    call re_alloc( ur,    1,nq,                   myName//'ur'   )
    call re_alloc( uk,    1,nq,                   myName//'uk'   )
    call re_alloc( tvdw, 1,maxPoints, 1,nq,       myName//'tvdw' )
    call re_alloc( tq, 0,mesh(1)-1, 0,mesh(2)-1, 0,mesh(3)-1, 1,1, &
                   myName//'tq' )

    ! Assign another name to these arrays (but be careful!!!)
    uvdw => tvdw
    uq   => tq

    ! Set mesh cutoff to filter VdW kernel
    kcut = meshKcut( cell, nMesh )
    call vdw_set_kcut( kcut )

    ! Find vacuum value of theta
    D = 0._dp
    GD = 0._dp
    call vdw_theta( nSpin, D, GD, tr, dtdd, dtdgd )
    tvac = tr

#ifdef DEBUG_XC
!    call timer_start( 'cellXC1' )
#endif /* DEBUG_XC */

    ! Loop on mesh points to find theta_q(r)
    ip = 0
    do i3 = 0,myMesh(3)-1
    do i2 = 0,myMesh(2)-1
    do i1 = 0,myMesh(1)-1

      ! Skip point if density=0
      if (.not.nonempty(i1,i2,i3)) cycle
      ip = ip + 1

      ! Find mesh indexes relative to cell origin
      ii1 = i1 + myBox(1,1)  
      ii2 = i2 + myBox(1,2)
      ii3 = i3 + myBox(1,3)

      ! Find density at this point. Notice that mesh indexes of dens and myDens 
      ! are relative to box and cell origins, respectively
      if (associated(myDens)) then
        D(:) = myDens(ii1,ii2,ii3,:)
      else
        D(:) = dens(i1,i2,i3,:)
      end if

      ! Avoid negative densities
      D(1:ndSpin) = max( D(1:ndSpin), 0._dp )

      !  Find gradient of density at this point
      call getGradDens( ii1, ii2, ii3, GD )   ! This subr. is contained below

      ! Find expansion of theta(q(r)) for VdW
      call vdw_theta( nSpin, D, GD, tr, dtdd, dtdgd )
      tvdw(ip,1:nq) = tr(1:nq)

#ifdef DEBUG_XC
!      ! Write q(r) for debugging
!      if (i3==myBox(1,3)) then
!        if (i1==myBox(1,1) .and. i2==myBox(1,2)) then
!          open(unit=47,file='qvdw.out')
!        else if (i1==myBox(2,1) .and. i2==myBox(2,2)) then
!          close(unit=47)
!        end if
!        Dtot = sum(D(1:ndSpin))
!        do ix = 1,3
!          GDtot(ix) = sum(GD(ix,1:ndSpin))
!        end do
!        call qofrho( Dtot, GDtot, q, dqdD, dqdGD )
!        if (i1==myBox(1,1)) write(47,*) ' '
!        write(47,*) q
!      end if
#endif /* DEBUG_XC */

    enddo ! i1
    enddo ! i2
    enddo ! i3  (End of loop on mesh points to find theta_q(r))

#ifdef DEBUG_XC
!    call timer_stop( 'cellXC1' )
!    call timer_start( 'cellXC2' )
!    call timer_start( 'cellXC2.1' )
#endif /* DEBUG_XC */

    ! Fourier-tranform theta_iq(r)
    do iq = 1,nq
      ! Unpack nonempty points into a full-mesh array
      ! Last index (=1) is just because fftr2k expects a rank 4 array
!      tq(0:myMesh(1)-1,0:myMesh(2)-1,0:myMesh(3)-1,1) = &  ! Slower!!!
!        unpack( tvdw(1:nonemptyPoints,iq), mask=nonempty, field=tvac(iq) )
      ip = 0
      do i3 = 0,myMesh(3)-1
      do i2 = 0,myMesh(2)-1
      do i1 = 0,myMesh(1)-1
        if (nonempty(i1,i2,i3)) then
          ip = ip+1
          tq(i1,i2,i3,1) = tvdw(ip,iq)
        else
          tq(i1,i2,i3,1) = tvac(iq)
        end if
      end do
      end do
      end do
      ! Peform Fourier transform
      call fftr2k( nMesh, myDistr, tq(:,:,:,1:1) )
      ! Store all reciprocal space points
!      tvdw(1:kPoints,iq) = &   ! Slower!!!
!        reshape( tq(0:kMesh(1)-1,0:kMesh(2)-1,0:kMesh(3)-1,1), (/kPoints/) )
      ik = 0
      do i3 = 0,kMesh(3)-1
      do i2 = 0,kMesh(2)-1
      do i1 = 0,kMesh(1)-1
        ik = ik+1
        tvdw(ik,iq) = tq(i1,i2,i3,1)
      end do
      end do
      end do
    end do ! iq

#ifdef DEBUG_XC
!    call timer_stop( 'cellXC2.1' )
!    call timer_start( 'cellXC2.2' )
#endif /* DEBUG_XC */

    ! Find reciprocal unit vectors
    call reclat( cell, kcell, 1 )

    ! Initialize nonlocal parts of VdW correlation energy and stress
    Enl = 0.0_dp
    EcuspVDW = 0.0_dp
    stressVDW = 0.0_dp

    ! Loop on k-mesh points
    ik = 0
    do i3 = 0,kMesh(3)-1   ! Mesh indexes relative to my k-box origin
    do i2 = 0,kMesh(2)-1
    do i1 = 0,kMesh(1)-1
      ik = ik + 1

      ! Find k vector
      ii1 = i1 + kBox(1,1)  ! Mesh indexes relative to reciprocal cell origin
      ii2 = i2 + kBox(1,2)
      ii3 = i3 + kBox(1,3)
      if (ii1 > nMesh(1)/2) ii1 = ii1 - nMesh(1)
      if (ii2 > nMesh(2)/2) ii2 = ii2 - nMesh(2)
      if (ii3 > nMesh(3)/2) ii3 = ii3 - nMesh(3)
      kvec(:) = kcell(:,1)*ii1 + kcell(:,2)*ii2 + kcell(:,3)*ii3
      k = sqrt( sum(kvec**2) )

      if (k<kcut) then

        ! Find Fourier transform of VdW kernel phi(r,r')
        call vdw_phi( k, phi, dphidk )

        ! Find Fourier transform of Int_dr'*phi(r,r')*rho(r')
        ! Warning: tvdw and uvdw are the same array
        tk(1:nq) = tvdw(ik,1:nq)
        uk(1:nq) = matmul( tk(1:nq), phi(1:nq,1:nq) )
        uvdw(ik,1:nq) = uk(1:nq)

#ifdef DEBUG_XC
!        ! Find contribution to 0.5*Int_dr*Int_dr'*rho(r)*phi(r,r')*rho(r')
!        ! Factor 0.5 in the integral cancels with a factor 2 required 
!        ! because tk and uk contain only the real or imaginary parts
!        ! of the Fourier components (see fftr2k) 
!        Enl = Enl + volume * sum(uk*tk)
#endif /* DEBUG_XC */

        ! Find contribution to stress from change of k vectors with strain
        if (k > kmin) then  ! Avoid k=0 (whose contribution is zero)
          dudk(1:nq) = matmul( tk(1:nq), dphidk(1:nq,1:nq) )
          ! See note above on cancelation of factors 0.5 and 2 in Enl
          dedk = sum(dudk*tk) * (volume / k)
          do jx = 1,3
            do ix = 1,3
              stressVDW(ix,jx) = stressVDW(ix,jx) - dedk * kvec(ix) * kvec(jx)
            end do
          end do
        end if ! (k>kmin)

      else
        uvdw(ik,1:nq) = 0.0_gp
      end if ! (k<kcut)

    end do ! i1
    end do ! i2
    end do ! i3  End of loop on k-mesh points

#ifdef DEBUG_XC
!    call timer_stop( 'cellXC2.2' )
!    call timer_start( 'cellXC2.3' )
#endif /* DEBUG_XC */

#ifdef DEBUG_XC
!    print'(a,3f12.6)','cellXC: Ex,Ec,Enl (eV) =', &
!      Ex/0.03674903_dp, Ec/0.03674903_dp, Enl/0.03674903_dp
#endif /* DEBUG_XC */

    ! Fourier-tranform u_q(k) back to real space
    do iq = 1,nq
!      uq(0:kMesh(1)-1,0:kMesh(2)-1,0:kMesh(3)-1,1) = &    ! Slower!!!
!        reshape( uvdw(1:kPoints,iq), kMesh )
      ik = 0
      do i3 = 0,kMesh(3)-1
      do i2 = 0,kMesh(2)-1
      do i1 = 0,kMesh(1)-1
        ik = ik+1
        uq(i1,i2,i3,1) = uvdw(ik,iq)
      end do
      end do
      end do
      call fftk2r( nMesh, myDistr, uq(:,:,:,1:1) )
!      uvdw(1:nonemptyPoints,iq) = &    ! Slower!!!
!        pack( uq(0:myMesh(1)-1,0:myMesh(2)-1,0:myMesh(3)-1,1), mask=nonempty )
      ip = 0
      do i3 = 0,myMesh(3)-1
      do i2 = 0,myMesh(2)-1
      do i1 = 0,myMesh(1)-1
        if (nonempty(i1,i2,i3)) then
          ip = ip+1
          uvdw(ip,iq) = uq(i1,i2,i3,1)
        end if
      end do
      end do
      end do
    end do

#ifdef DEBUG_XC
!    call timer_stop( 'cellXC2.3' )
!    call timer_stop( 'cellXC2' )
#endif /* DEBUG_XC */

#ifdef DEBUG_XC
!    ! Re-initialize Enl if it has been calculated in reciprocal space
!    Enl = 0.0_dp
#endif /* DEBUG_XC */

  end if ! (VDW) End of VdW initializations------------------------------------

#ifdef DEBUG_XC
!  call timer_start( 'cellXC3' )
#endif /* DEBUG_XC */

  ! Loop on mesh points -------------------------------------------------------
  ip = 0
  do i3 = 0,myMesh(3)-1   ! Mesh indexes relative to my box origin
  do i2 = 0,myMesh(2)-1
  do i1 = 0,myMesh(1)-1
    ii1 = i1 + myBox(1,1) ! Mesh indexes relative to cell origin
    ii2 = i2 + myBox(1,2)
    ii3 = i3 + myBox(1,3)

    ! Find density at this point. Notice that mesh indexes of dens and myDens 
    ! are relative to box and cell origins, respectively
    if (associated(myDens)) then
      D(:) = myDens(ii1,ii2,ii3,:)
    else
      D(:) = dens(i1,i2,i3,:)
    end if

    ! Skip point if density=0
    Dtot = sum(D(1:ndSpin))
    if (Dtot < Dmin) cycle ! i1 loop on mesh points
    ip = ip + 1

    ! Avoid negative densities
    D(1:ndSpin) = max( D(1:ndSpin), 0._dp )

    ! Find gradient of density at this point
    if (GGA) call getGradDens( ii1, ii2, ii3, GD )

    ! Loop over all functionals
    do nf = 1,nXCfunc

      ! Is this a VDW or GGA?
      if (XCfunc(nf).eq.'VDW' .or. XCfunc(nf).eq.'vdw') then
        VDWfunctl = .true.
        GGAfunctl = .true.
      else if (XCfunc(nf).eq.'GGA' .or. XCfunc(nf).eq.'gga') then
        VDWfunctl = .false.
        GGAfunctl = .true.
      else
        VDWfunctl = .false.
        GGAfunctl = .false.
      endif

      ! Find exchange and correlation energy densities and their 
      ! derivatives with respect to density and density gradient
      if (VDWfunctl) then

        ! Local exchange-corr. part from the apropriate LDA/GGA functional
        call vdw_localxc( irel, nSpin, D, GD, epsX, epsC, &
                          dExdD, dEcdD, dExdGD, dEcdGD )

#ifdef DEBUG_XC
!        ! Select only non local correlation energy and potential
!        epsX = 0
!        epsC = 0
!        dExdD = 0
!        dEcdD = 0
!        dExdGD = 0
!        dEcdGD = 0
#endif /* DEBUG_XC */

        ! Local cusp correction to nonlocal VdW energy integral
        call vdw_decusp( nSpin, D, GD, epsCusp, dEcuspdD, dEcuspdGD )

#ifdef DEBUG_XC
!        ! Select only non local correlation energy and potential
!        epsCusp = 0
!        dEcuspdD = 0
!        dEcuspdGD = 0
#endif /* DEBUG_XC */

        ! Find expansion of theta(q(r)) for VdW
        call vdw_theta( nSpin, D, GD, tr, dtdd, dtdgd )

        ! Add nonlocal VdW energy contribution and its derivatives
        Dtot = sum(D(1:ndSpin))
        ur(1:nq) = uvdw(ip,1:nq)
        epsNL = epsCusp + 0.5_dp*sum(ur*tr)/(Dtot+tiny(Dtot))
        epsC = epsC + epsNL
        do is = 1,nSpin
          dEcdD(is) = dEcdD(is) + dEcuspdD(is) + sum(ur(1:nq)*dtdd(1:nq,is))
          do ix = 1,3
            dEcdGD(ix,is) = dEcdGD(ix,is) + dEcuspdGD(ix,is) + &
                            sum(ur(1:nq)*dtdgd(ix,1:nq,is))
          end do
        end do

        ! Sum nonlocal VdW contributions for debugging
        EcuspVDW = EcuspVDW + dVol * Dtot * epsCusp
        Enl = Enl + dVol * Dtot * epsNL

#ifdef DEBUG_XC
!        if (i1==0 .and. i2==0 .and. i3==0) then
!          open( unit=33, file='epsNL'//trim(nodeString()), form='formatted' )
!          write(33,'(3f12.6,i6)') (cell(:,ix),nMesh(ix),ix=1,3)
!        end if
!        write(33,'(3i6,3e15.6)') ii1, ii2, ii3, Dtot, epsNL, epsX+epsC
!        ! Alternatively, write x,y,z instead of i1,i2,i3
!        write(33,'(3f12.6,2e15.6)') &
!          ii1*cell(:,1)/nMesh(1), &
!          ii2*cell(:,2)/nMesh(2), &
!          ii3*cell(:,3)/nMesh(3), Dtot, epsNL
#endif /* DEBUG_XC */

      else if (GGAfunctl) then
        call ggaxc( XCauth(nf), irel, nSpin, D, GD, &
                    epsX, epsC, dExdD, dEcdD, dExdGD, dEcdGD )
      else ! (.not.VDWfunctl .and. .not.GGAfunctl)
        call ldaxc( XCauth(nf), irel, nSpin, D, &
                    epsX, epsC, dExdD, dEcdD, dVxdD, dVcdD )
      endif ! (VDWfunctl)

      ! Scale return values by weight for this functional
      epsX = XCweightX(nf)*epsX
      epsC = XCweightC(nf)*epsC
      dExdD(:) = XCweightX(nf)*dExdD(:)
      dEcdD(:) = XCweightC(nf)*dEcdD(:)
      if (GGAfunctl) then
        dExdGD(:,:) = XCweightX(nf)*dExdGD(:,:)
        dEcdGD(:,:) = XCweightC(nf)*dEcdGD(:,:)
      endif
      if (present(dVxcdD)) then
        dVxdD(:) = XCweightX(nf)*dVxdD(:)
        dVcdD(:) = XCweightC(nf)*dVcdD(:)
      endif

      ! Add contributions to exchange-correlation energy and its
      ! derivatives with respect to density at this point
      Ex = Ex + dVol * Dtot * epsX
      Ec = Ec + dVol * Dtot * epsC
      Dx = Dx + dVol * Dtot * epsX
      Dc = Dc + dVol * Dtot * epsC
      Dx = Dx - dVol * sum(D(:)*dExdD(:))
      Dc = Dc - dVol * sum(D(:)*dEcdD(:))
      if (associated(myVxc)) then
        myVxc(ii1,ii2,ii3,:) = myVxc(ii1,ii2,ii3,:) + dExdD(:) + dEcdD(:)
      else  ! Add directly to output array Vxc
        Vxc(i1,i2,i3,:) = Vxc(i1,i2,i3,:) + dExdD(:) + dEcdD(:)
      end if ! (associated(myVxc))
      if (present(dVxcdD)) then
        if (associated(mydVxcdD)) then
          mydVxcdD(ii1,ii2,ii3,:) = mydVxcdD(ii1,ii2,ii3,:) + dVxdD(:)+dVcdD(:)
        else
          dVxcdD(i1,i2,i3,:) = dVxcdD(i1,i2,i3,:) + dVxdD(:) + dVcdD(:)
        end if
      end if ! (present(dVxcdD))

      ! Add contributions to exchange-correlation potential
      ! with respect to density at neighbor points
      if (GGAfunctl) then
        if (myDistr==0) then   ! dens data not distributed
          do ic = 1,3          ! Loop on cell axes
            do in = -nn,nn     ! Loop on finite difference index
              ! Find mesh indexes of neighbor point
              jj(1) = ii1
              jj(2) = ii2
              jj(3) = ii3
              jj(ic) = modulo( jj(ic)+in, nMesh(ic) )
              ! Add contributions from dE/dGradDi * dGradDi/dDj
              ! Notice: for myDistr==0, box and cell origins are equal
              do is = 1,nSpin  ! Loop on spin component
                dExidDj = sum( dExdGD(:,is) * dGidFj(:,ic,in) )
                dEcidDj = sum( dEcdGD(:,is) * dGidFj(:,ic,in) )
                Dx = Dx - dVol * dens(jj(1),jj(2),jj(3),is) * dExidDj
                Dc = Dc - dVol * dens(jj(1),jj(2),jj(3),is) * dEcidDj
                Vxc(jj(1),jj(2),jj(3),is) = Vxc(jj(1),jj(2),jj(3),is) + &
                                            dExidDj + dEcidDj
              end do ! is
            end do ! in
          end do ! ic
        else ! (myDistr/=0)      Distributed dens data
          do ic = 1,3          ! Loop on cell axes
            if (ic==1) then
              Dleft => Dleft1
              Drght => Drght1
              Vleft => Vleft1
              Vrght => Vrght1
            else if (ic==2) then 
              Dleft => Dleft2
              Drght => Drght2
              Vleft => Vleft2
              Vrght => Vrght2
            else ! (ic==3)
              Dleft => Dleft3
              Drght => Drght3
              Vleft => Vleft3
              Vrght => Vrght3
            end if ! (ic==1)
            do in = -nn,nn     ! Loop on finite difference index
              ! Find indexes jj(:) of neighbor point
              jj(1) = ii1  ! Warning: jj(:)=(/ii1,ii2,ii3/) is VERY slow!!!
              jj(2) = ii2
              jj(3) = ii3
              jj(ic) = jj(ic) + in
              ! Point Dj and Vj to the apropriate array
              if (jj(ic)<myBox(1,ic)) then ! Left neighbor region
                Dj =  Dleft(jj(1),jj(2),jj(3),1:nSpin)
                Vj => Vleft(jj(1),jj(2),jj(3),1:nSpin)
              else if (jj(ic)>myBox(2,ic)) then ! Right region
                Dj =  Drght(jj(1),jj(2),jj(3),1:nSpin)
                Vj => Vrght(jj(1),jj(2),jj(3),1:nSpin)
              else ! j within myBox
                Dj = myDens(jj(1),jj(2),jj(3),1:nSpin)
                Vj => myVxc(jj(1),jj(2),jj(3),1:nSpin)
              end if
              ! Add contributions from dE/dGradDi * dGradDi/dDj
              do is = 1,nSpin  ! Loop on spin component
                dExidDj = sum( dExdGD(:,is) * dGidFj(:,ic,in) )
                dEcidDj = sum( dEcdGD(:,is) * dGidFj(:,ic,in) )
                Dx = Dx - dVol * Dj(is) * dExidDj
                Dc = Dc - dVol * Dj(is) * dEcidDj
                Vj(is) = Vj(is) + dExidDj + dEcidDj
              end do ! is
            end do ! in
          end do ! ic
        end if ! (myDistr==0)
      end if ! (GGAfunctl)

      ! Add contribution to stress due to change in gradient of density
      ! originated by the deformation of the mesh with strain
      if (GGAfunctl) then
        do jx = 1,3
          do ix = 1,3
            do is = 1,nSpin
              stress(ix,jx) = stress(ix,jx) - dVol * GD(ix,is) * &
                               ( dExdGD(jx,is) + dEcdGD(jx,is) )
            enddo
          enddo
        enddo
      endif ! (GGAfunctl)

    enddo ! nf (End of loop over functionals)

  enddo ! i1
  enddo ! i2
  enddo ! i3  (End of loop over mesh points)-----------------------------------

#ifdef DEBUG_XC
  close( unit=33 )
#endif /* DEBUG_XC */

  ! If mesh arrays are distributed, add Vxc contribution from neighbor regions
  if (GGA .and. myDistr/=0) then ! Distributed Vxc data
    ! Add neighbor regions contribution to myVxc array
    do ic = 1,3  ! Loop on cell axes
      if (ic==1) then
        Vleft => Vleft1
        Vrght => Vrght1
      else if (ic==2) then 
        Vleft => Vleft2
        Vrght => Vrght2
      else ! (ic==3)
        Vleft => Vleft3
        Vrght => Vrght3
      end if ! (ic==1)
      boxLeft = myBox
      boxLeft(1,ic) = myBox(1,ic)-nn
      boxLeft(2,ic) = myBox(1,ic)-1
      boxRght = myBox
      boxRght(1,ic) = myBox(2,ic)+1
      boxRght(2,ic) = myBox(2,ic)+nn
      call associateMeshTask( left2my(ic), myDistr )
      call associateMeshTask( rght2my(ic), myDistr )
      call addMeshData( nMesh, boxLeft, Vleft, myDistr, myVxc, left2my(ic) )
      call addMeshData( nMesh, boxRght, Vrght, myDistr, myVxc, rght2my(ic) )
!      call addMeshData( nMesh, boxLeft, Vleft, myDistr, myVxc )
!      call addMeshData( nMesh, boxRght, Vrght, myDistr, myVxc )
    end do ! ic
  end if ! (GGA .and. myDistr/=0)

  ! Make Vxc=0 if VDWfunctl and Dens<Dcut, to avoid singularities
  if (VDWfunctl) then
    do i3 = 0,myMesh(3)-1   ! Mesh indexes relative to my box origin
    do i2 = 0,myMesh(2)-1
    do i1 = 0,myMesh(1)-1
      ii1 = i1 + myBox(1,1) ! Mesh indexes relative to cell origin
      ii2 = i2 + myBox(1,2)
      ii3 = i3 + myBox(1,3)
      if (associated(myDens)) then
        Dtot = sum( myDens(ii1,ii2,ii3,1:ndSpin) )
      else
        Dtot = sum( dens(i1,i2,i3,1:ndSpin) )
      end if
      if (Dtot<Dcut) then
        if (associated(myVxc)) then
          myVxc(ii1,ii2,ii3,:) = 0
        else
          Vxc(i1,i2,i3,:) = 0
        end if
        if (present(dVxcdD)) then
          if (associated(mydVxcdD)) then
            mydVxcdD(ii1,ii2,ii3,:) = 0
          else
            dVxcdD(i1,i2,i3,:) = 0
          end if
        end if ! (present(dVxcdD))
      end if ! (Dtot<Dcut)
    end do
    end do
    end do
  end if ! (VDWfunctl)

  ! Copy Vxc data to output arrays
  if (associated(myVxc)) then  ! Distributed Vxc array
    if (sameMeshDistr(ioDistr,myDistr)) then ! Just copy myVxc to output array
      Vxc = myVxc
      if (present(dVxcdD)) dVxcdD = mydVxcdD
    else ! Redistribution required
      ! Copy myVxc and dVxcdD to output box
      call associateMeshTask( my2io, myDistr )
      call copyMeshData( nMesh, myDistr, myVxc, ioBox, Vxc, my2io )
!      call copyMeshData( nMesh, myDistr, myVxc, ioBox, Vxc )
      if (present(dVxcdD)) &
        call copyMeshData( nMesh, myDistr, mydVxcdD, ioBox, dVxcdD, my2io )
!        call copyMeshData( nMesh, myDistr, mydVxcdD, ioBox, dVxcdD )
    end if ! (sameMeshDistr(ioDistr,myDistr))
  end if ! (associated(myVxc))

#ifdef DEBUG_XC
!  call timer_stop( 'cellXC3' )
#endif /* DEBUG_XC */

#ifdef DEBUG_XC
!  ! Some printout for debugging
!  if (VDW) then
!    print'(a,f12.6)', &
!         'cellXC: EcuspVDW (eV) =', EcuspVDW/0.03674903_dp
!    print'(a,3f12.6)','cellXC: Ex,Ecl,Ecnl (eV) =', &
!         Ex/0.03674903_dp, (Ec-Enl)/0.03674903_dp, Enl/0.03674903_dp
!  end if
#endif /* DEBUG_XC */

  ! Add contribution to stress from the change of volume with strain
  forall(ix=1:3) stress(ix,ix) = stress(ix,ix) + Ex + Ec

  ! Add contribution to stress from change of k vectors with strain
  if (VDW) stress = stress + stressVDW * XCweightVDW

  ! Divide by volume to get correct stress definition (dE/dStrain)/Vol
  stress = stress / volume

  ! Divide by energy unit
  Ex = Ex / Eunit
  Ec = Ec / Eunit
  Dx = Dx / Eunit
  Dc = Dc / Eunit
  Vxc = Vxc / Eunit
  stress = stress / Eunit
  if (present(dVxcdD)) dVxcdD = dVxcdD / Eunit

  ! Deallocate VDW-related arrays
  if (VDW) then
    call de_alloc( tq,       myName//'tq' )
    call de_alloc( tvdw,     myName//'tvdw' )
    call de_alloc( uk,       myName//'uk' )
    call de_alloc( ur,       myName//'ur' )
    call de_alloc( tvac,     myName//'tvac' )
    call de_alloc( tr,       myName//'tr' )
    call de_alloc( tk,       myName//'tk' )
    call de_alloc( phi,      myName//'phi' )
    call de_alloc( dphidk,   myName//'dphidk' )
    call de_alloc( dtdgd,    myName//'dtdgd' )
    call de_alloc( dtdd,     myName//'dtdd' )
    call de_alloc( dudk,     myName//'dudk' )
    call de_alloc( nonempty, myName//'nonempty' )
  end if

  ! Deallocate GGA-related arrays
  if (associated(myDens)) then
    if (GGA) then
      call de_alloc( Vrght3, myName//'Vrght3' )
      call de_alloc( Vleft3, myName//'Vleft3' )
      call de_alloc( Vrght2, myName//'Vrght2' )
      call de_alloc( Vleft2, myName//'Vleft2' )
      call de_alloc( Vrght1, myName//'Vrght1' )
      call de_alloc( Vleft1, myName//'Vleft1' )
      call de_alloc( Drght3, myName//'Drght3' )
      call de_alloc( Dleft3, myName//'Dleft3' )
      call de_alloc( Drght2, myName//'Drght2' )
      call de_alloc( Dleft2, myName//'Dleft2' )
      call de_alloc( Drght1, myName//'Drght1' )
      call de_alloc( Dleft1, myName//'Dleft1' )
    end if
    if (present(dVxcdD)) call de_alloc( mydVxcdD, myName//'mydVxcdD'  )
    call de_alloc( myVxc,  myName//'myVxc'  )
    call de_alloc( myDens, myName//'myDens' )
  end if ! (associated(myDens))

  ! Restore previous allocation defaults
  call alloc_default( restore=prevAllocDefaults )

  ! Stop time counter
  call timer_stop( myName )

  ! Get local calculation time (excluding communications)
  call timer_get( myName, lastTime=totTime, lastCommTime=comTime )
  myTime = totTime - comTime
  sumTime  = myTime
  sumTime2 = myTime**2

  ! Add integrated magnitudes from all processors
  call miscAllReduce( 'sum', Ex, Ec, Dx, Dc, sumTime, sumTime2, a2=stress )

  ! Find average and dispersion of CPU time
  timeAvge = sumTime / nodes
  timeDisp = sqrt( max( sumTime2/nodes - timeAvge**2, 0._dp ) )

#ifdef DEBUG_XC
  write(udebug,'(a,3f12.6,/)') &
    myName//'My CPU time, avge, rel.disp =', &
    myTime, timeAvge, timeDisp/timeAvge
#endif /* DEBUG_XC */

CONTAINS !---------------------------------------------------------------------

  subroutine getGradDens( ii1, ii2, ii3, GD )

  ! Finds the density gradient at one mesh point

  ! Arguments
  integer, intent(in) :: ii1, ii2, ii3  ! Global mesh point indexes
  real(dp),intent(out):: GD(3,nSpin) ! Density gradient

  ! Variables and arrays accessed from parent subroutine:
  !   dens, Dleft1, Dleft2, Dleft3, 
  !   Drght1, Drght2, Drght3, DGiDFj,
  !   myBox, myDistr, nn, nSpin

  ! Local variables and arrays
  integer :: ic, in, is, jj(3)
  real(dp):: Dj(nSpin)
  real(gp),pointer:: Dleft(:,:,:,:), Drght(:,:,:,:)

  GD(:,:) = 0
  if (myDistr==0) then   ! dens data not distributed
    do ic = 1,3          ! Loop on cell axes
      do in = -nn,nn     ! Loop on finite difference index
        ! Find index jp of neighbor point
        jj(1) = ii1  ! Warning: jj(:)=(/ii1,ii2,ii3/) is VERY slow!!!
        jj(2) = ii2
        jj(3) = ii3
        jj(ic) = modulo( jj(ic)+in, nMesh(ic) )
        ! Find contribution of density at j to gradient at i
        do is = 1,nSpin
          do ix = 1,3  ! Warning: GD(:,is)=GD(:,is)+... is slower!!!
            GD(ix,is) = GD(ix,is) + DGiDFj(ix,ic,in) * &
                                    dens(jj(1),jj(2),jj(3),is)
          end do ! ix
        end do ! is
      end do ! in
    end do ! ic
  else                   ! Distributed dens data
    do ic = 1,3          ! Loop on cell axes
      if (ic==1) then
        Dleft => Dleft1
        Drght => Drght1
      else if (ic==2) then 
        Dleft => Dleft2
        Drght => Drght2
      else ! (ic==3)
        Dleft => Dleft3
        Drght => Drght3
      end if ! (ic==1)
      do in = -nn,nn     ! Loop on finite difference index
        ! Find indexes jj(:) of neighbor point
        jj(1) = ii1
        jj(2) = ii2
        jj(3) = ii3
        jj(ic) = jj(ic) + in
        ! Find density Dj at neighbor point j
        if (jj(ic)<myBox(1,ic)) then ! Left neighbor region
          Dj(1:nSpin) = Dleft(jj(1),jj(2),jj(3),1:nSpin)
        else if (jj(ic)>myBox(2,ic)) then ! Right region
          Dj(1:nSpin) = Drght(jj(1),jj(2),jj(3),1:nSpin)
        else ! j within myBox
          Dj(1:nSpin) = myDens(jj(1),jj(2),jj(3),1:nSpin)
        end if
        ! Find contribution of density at j to gradient at i
        do is = 1,nSpin  ! Loop on spin component
          do ix = 1,3    ! Warning: GD(:,is)=GD(:,is)+... is slower!!!
            GD(ix,is) = GD(ix,is) + DGiDFj(ix,ic,in) * Dj(is)
          end do ! ix
        end do ! is
      end do ! in
    end do ! ic
  end if ! (myDistr==0)

  end subroutine getGradDens

END SUBROUTINE cellXC

END MODULE m_cellXC
