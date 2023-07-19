module main_vars
  use precision, only: dp, sp
  use subs, only: ival
  use m_getopts
  use units, only: Ang, pi

  implicit none

  public 

  integer :: ierr, klb, it, is, k, nw
  integer :: nao, ia, iz, ko, nkp, nsp, nen, mxwf, io, ie
  integer :: nnao, ik, is0, iw, iw0, i1, i2, i3, i4
  integer :: no_s, nspin, nh, im, ii, io2
  integer :: ncb, nln, il, nwd, n, l, i, noc0, ic, io1, isr, iov
  integer :: j, m, idos, lorb, naoatx

  integer :: nspecies, na_u, no_u
  integer, allocatable :: no(:), nquant(:,:), lquant(:,:), zeta(:,:)
  integer, allocatable :: iaorb(:), iphorb(:)
  character(len=20), allocatable :: label(:)

  real(dp) :: qtot, temp_in_file, dm, alfa, vvv
  real(dp) :: qcos, qsin, w0, want
  real(dp) :: min_energy, max_energy, e_step, energy, weight, efermi
  real(dp) :: low_e, high_e, eigval
  integer  :: intdos_u, number_of_wfns

  real(dp) :: minimum_spec_energy = -huge(1.0_dp)
  real(dp) :: maximum_spec_energy = huge(1.0_dp)

  real(dp)            ::  smear = 0.5        ! Units of energy are eV 
  integer             ::  npts_energy = 200

  integer :: nshmx
  integer, parameter :: ncbmx=20
  integer, parameter :: nlwmx=30

  character :: sflnm*50, taux*100, wrd(nlwmx)*20, cx*20
  integer :: mpr_u=50, wk_u=51
  integer :: out_u=70, wfs_u=72, hs_u=73
  integer :: stt_u=60, tab_u=61
  logical   :: wk_x, wfs_x, hs_x, tab_x

  ! OUT file
  integer, allocatable :: isa(:)
  integer, allocatable :: za(:), zc(:), zn(:), zl(:), zx(:), zz(:)
  real(dp), allocatable :: zval(:)

  ! WFS file
  integer, allocatable :: nwf(:,:)
  real(dp),    allocatable :: pk(:,:)
  real(dp), allocatable ::   ados(:,:), ww(:)
  real(dp), allocatable ::   intdos(:), intebs(:)

  ! HS file
  integer, allocatable :: numh(:), listhptr(:), listh(:)  
  integer, allocatable :: indxuo(:)
  real(dp),    allocatable :: hamilt(:,:), Sover(:), xij(:,:), dij(:)

  real(dp),    allocatable :: wk(:)
  real(SP),  allocatable :: wf(:,:)       ! Note single precision

  ! MPR file
  character :: what*4, tit(ncbmx)*30
  logical   :: dos, coop
  real(dp)  :: dtc(ncbmx,2)

  integer              :: noc(ncbmx,2)
  integer, allocatable :: koc(:,:,:)
  logical, allocatable :: orb_mask(:,:,:)

  ! RESULTS
  real(dp),  allocatable :: coop_vals(:,:,:)
  real(dp),  allocatable :: cohp_vals(:,:,:)
  real(dp),  allocatable :: pdos_vals(:,:,:)
  logical,   allocatable :: ref_mask(:)

  !
  logical :: gamma
  real(dp) :: r_dummy(3), dummy_weight, ztot
  integer  :: idummy
  !
  character(len=200) :: opt_arg, mflnm, ref_line
  character(len=10)  :: opt_name 
  integer :: nargs, iostat, n_opts, nlabels, iorb, ikb
  integer :: nkb, nkp_wfs

  logical :: debug    = .false.
  logical :: simple_dos = .true.
  logical :: ref_line_given = .false.
  logical :: energies_only = .false.

end module main_vars
!------------------------------------------------------------------
