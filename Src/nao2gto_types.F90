! *** Module: nao2gto_types ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Management of NAO2GTO data
!!
!! This module defines data structures to handle NAO2GTO options and
!! prescreening tolerances of ERIs in HONPAS. The specifications are read
!! from the FDF input file.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 03.2016 Created [Xinming Qin]
!!      - 01.2018 Refactored for SIESTA 4 [Yann Pouillon]
!!      - 01.2022 Optimized on SIESTA v4.1.5 [Xinming Qin]
! *****************************************************************************
module nao2gto_types

  use precision, only: dp
  use atm_types, only: maxn_orbnl, maxnorbs

  implicit none

  private

  ! HFX potential types
  !
  !> Use the full Coulomb potential to compute Hartree-Fock XC contributions
  integer, parameter, public :: do_hfx_potential_coulomb = 1

  !> Use the short-range potential to compute Hartree-Fock XC contributions
  integer, parameter, public :: do_hfx_potential_short   = 2

  !> Truncate the Coulomb potential to compute Hartree-Fock XC contributions
  integer, parameter, public :: do_hfx_potential_truncated = 3


  !
  ! Boundaries
  !

  !> Limit to approximate log(0) when evaluating Hartree-Fock gradients
  real(dp), parameter, public :: log_zero = -1000.0_dp

  !> Minimum log value for the Powell optimization method
  real(dp), parameter, public :: powell_min_log = -20.0_dp

  !
  ! Gaussian-Type Orbitals
  !

  ! NAO2GTO: NAO = sum RGTOs = sum CGTOs
  ! RGTOs  : Real Spherical Harmonic GTOs or slater-type GTOs
  ! CGTOs  : Cartesian GTOs

  ! Number of CGTOs larger than maxnorbs (PAOs of an atom)
  integer, parameter, public :: maxnorbs_cphi = 120

  ! Max numbers of CGTOs and RGTOs
  integer, parameter, public :: maxn_contract = 10  !   Maximum number of
                                                    !   Gaussians that enter in
                                                    !   the linear combination
                                                    !   to expand a given NAO
  integer, parameter, public :: ncon_max = 60

  ! Support Max 6 GTOs for d orbitals, 10 GTOs for f orbitals
  integer, parameter, public :: l_max = 3

  !> \brief Data type to store Hartree-Fock exchange options that can be
  !!        read from SIESTA input files
  !!
  !! \par Default values from FDF:
  !!      - DM_trunc       = .true. (use sparse DM to screen ERIs)
  !!      - dump_fit_data  = .true. (dump fit data to nao2gto_fit.yml)
  !!      - farfield       = .true. (far-near field screening)
  !!      - Dynamic_parallel = .ture.  1-level Master-Worker dynamic parallel for ERIs and HFX
  !!      - frag_size  = 10000    Number of batched ERIs of dynamic paralle per request-send task
  !!      - npts_fit       = NTBMAX (number of data points to fit orbitals)
  !!      - min_num_gaus   = 3    (minimum number of gaussians in the expansion)
  !!      - max_num_gaus   = 6    (maximum number of gaussians in the expansion)
  !!      - potential_type = 1 (1/r/erfc(wr)/r ...)
  !!      - omega          = 0.11
  !!      - cutoff_radius  = ???
  !!      - eps_far        = 1.0e-6  (far-field screening tolerance)
  !!      - eps_pairlist   = 1.0e-6  (build shell pair-list tolerance)
  !!      - eps_schwarz    = 1.0e-6  (Schwarz tolerance)
  !!      - eps_stored     = 1.0e-6  (stored ERIs tolerance)
  !!      - threshold_exp_gaus = 1.4e0 (threshold to separate the consecutive exp)
  !!      - tolerance_gaus = 1.0e-3  (threshold to consider convergence in fit)
  type, public :: hfx_options_type
    logical   ::  DM_trunc = .false.
    logical   ::  is_fitted_nao =.true.
    logical   ::  dump_fit_data = .false.
    logical   ::  farfield = .false.
    logical   ::  Dynamic_parallel = .false.
    integer   ::  frag_size  = 10000
    integer   ::  npts_fit = -1
    integer   ::  potential_type = -1
    integer   ::  max_num_gaus_s   =  6
    integer   ::  max_num_gaus_p   =  5
    integer   ::  max_num_gaus_d   =  4
    real(dp)  ::  omega = 0.11_dp
    real(dp)  ::  cutoff_radius = 1.0d-6
    real(dp)  ::  eps_farfield = 1.0d-6
    real(dp)  ::  eps_pairlist = 1.0d-6
    real(dp)  ::  eps_schwarz = 1.0d-6
    real(dp)  ::  eps_stored = 1.0d-6
    real(dp)  ::  threshold_exp_gaus = 1.4_dp
    real(dp)  ::  tolerance_gaus = 1.d-3
    real(dp)  ::  gto_eps = 1.d-5
  end type hfx_options_type

  !> \brief Data type to store information about orbital pairs
  type, public :: pair_list_element_type

    integer  :: pair(2)          ! Indices of the two neighbour orbitals
                                 !   in the supercell
    integer  :: nl_index         ! Index that labels the neighbour pair of
                                 !   shells in such a way that will allow 
                                 !   to determine whether a quartet will be
                                 !   computed or not
!     Note that currently LIBINT has a very important restriction on the
!     angular momentum ordering of the functions in shell quartets that it
!     can handle.
!     LIBINT can evaluate a shell quartet (ab|cd) if
!     \lambda(a) \ge \lambda(b),
!     \lambda(c) \ge \lambda(d),
!     and \lambda(c) + \lambda(d) \ge \lambda(a) + \lambda(b).
!     If one needs to compute a quartet that doesn’t conform the rule, e.g.
!     of type (pf|sd), permutational symmetry of integrals can be utilized
!     to compute such quartet
!     (pq|rs)=(pq|sr)=(qp|rs)=(qp|sr)=(rs|pq)=(rs|qp)=(sr|pq)=(sr|qp)
!     In the case of (pf|sd) shell quartet, one computes quartet (ds|fp)
!     instead, and then permutes function indices back to obtain the
!     desired (pf|sd).
    real(dp) :: r1(3), r2(3)     ! Centers of the two neighbour orbitals 
                                 !   in the supercell
    real(dp) :: dist2            ! Square of the distance between the two
                                 !   neighbour orbitals

  end type

  !> \brief Data type to store a list of orbital-pair information
  type, public :: pair_list_type
    type(pair_list_element_type), dimension(:), allocatable :: element
    integer :: nelement = 0
  end type pair_list_type

!> Data type to store information about screening coefficients
  type, public :: hfx_screen_coeff_type
    real(dp) :: x(2)
  end type hfx_screen_coeff_type

  !> \brief Data type to point NAO2GTO routines to relevant SIESTA data
  type, public :: hfx_system_type

    real(dp) :: cell(3,3)
    real(dp) :: cell_r(3,3)

    integer , pointer :: maxnh
    integer , pointer :: na
    integer , pointer :: norb
    integer , pointer :: nspin
    integer , pointer :: nua
    integer , pointer :: nuo
    integer , pointer :: nuotot
    integer , pointer :: iaorb(:)
    integer , pointer :: indxua(:)
    integer , pointer :: iphorb(:)
    integer , pointer :: isa(:)
    integer , pointer :: listh(:)
    integer , pointer :: listhptr(:)
    integer , pointer :: nsc(:)
    integer , pointer :: numh(:)
    real(dp), pointer :: xa(:,:)

  end type hfx_system_type

                    ! ------------------------------------ !

  ! sphi: Slater-type orbital (spherical)
  ! cphi: Cartesian orbital
  type, public :: gto_info_type

    ! Added by Honghui Shang
    integer, dimension(maxn_orbnl)  ::  orbnl_contract      
                                        ! Number of Gaussians to expand the
                                        !  radial part of a NAO of a given shell
    integer, dimension(maxn_orbnl)  ::  orbnl_index_cphi
                                        ! Index where the first cartesian 
                                        !  gaussian function of a given shell 
                                        !  appears
    integer, dimension(maxn_orbnl)  ::  orbnl_index_sphi
                                        ! Index where the first spherical 
                                        !  gaussian function of a given shell 
                                        !  appears
    real(dp)                        ::  kind_radius             
                                        ! Largest radii (see pgf_radius below)
                                        !   between all the 
                                        !   Gaussians that expand any NAO 
                                        !   orbital of any shell for a given
                                        !   species
    real(dp), dimension(maxn_contract,maxn_orbnl) :: orbnl_zeta
                                        ! Exponent of the Gaussians in the 
                                        !   expansion
    real(dp), dimension(maxn_contract,maxn_orbnl) :: orbnl_coefficient
                                        ! Coefficient of the Gaussian in the
                                        !   expansion
    real(dp), dimension(maxn_contract,maxn_orbnl) :: pgf_radius
                                        ! Point where
                                        !   g(r)=coeff * r**l * exp(-zeta*r**2)
                                        !   is smaller than a given threshold, 
                                        !   1.d-5 in nao2gto_transfer 
                                        !   subroutine
    real(dp), dimension(maxn_orbnl) ::   shell_radius
                                        ! Largest value of pgf_radius  
                                        !   between all the Gaussians that
                                        !   expand a NAO of a given shell

    ! Added by Xinming Qin
    ! zeta and coeff of the most diffuse GTO. Normalized adjoint GTO
    ! NAO = sum Cm*GTOm

    ! minimum value of the exponent of the gaussian (zeta), the most diffuse gto         
    real(dp), dimension(maxn_orbnl) :: orbnl_adjoined_zeta

    ! the corresponding coefficient
    !real(dp), dimension(maxn_orbnl) :: orbnl_diff_coeff

    ! Added by Xinming Qin
    ! Max contraction coefficient: C=cartesian, R=spherical (radial)
    ! shell : same n and l,  
    ! CGTO (nco*k) ---> RGTO (nso*k) ---> PAO (nso)
    ! c2s(nco, nso)   Fit: Sum_k D_k*RGTO_k,  k = 0, M 
    ! RGTO (nso) = Transpose[c2s(nco, nso)] * CGTO (nco), matrix * vector 
    ! RGTO (io)   = sum c2s(1:nco, io) * CGTO(1:nco)
    !
    ! ERI (nsoa*nsob*nsoc*nsod),  nso PAO once !
    ! Contribution of each CGTO shell (nco) to PAO shell (nso) ERIs : 
    !
    !     Sum_k D_k* c2s_k(nco,nso)*GTO_k(nco) ), k=1, M
    !
    ! Actually, both D_k and normorlized coeff have been involved in c2s_k (nco, nso), 
    ! We calc the trans matrix sphi(nco*M, nso) !
    !
    ! So max contribution coeff of echo CGTO shell is:
    ! max( sum [sphi(1: nco, i)],  i = 1, nso ) !
    real(dp), dimension(maxn_orbnl) ::  orbnl_contraction_coeff

    !----------------cphi's data------------------------!
    integer :: norbs_cphi           ! Number of Cartessian Gaussian functions 
                                    !   per atom in the basis set
    integer,  dimension(maxnorbs_cphi) :: orb_cphi_lx
                                    ! Exponent of the x-coordinate in the 
                                    !   cartesian gaussian function
    integer,  dimension(maxnorbs_cphi) :: orb_cphi_ly
                                    ! Exponent of the y-coordinate in the 
                                    !   cartesian gaussian function
    integer,  dimension(maxnorbs_cphi) :: orb_cphi_lz
                                    ! Exponent of the z-coordinate in the 
                                    !   cartesian gaussian function
    real(dp), dimension(maxnorbs_cphi) :: norm_cphi
    integer,  dimension(maxnorbs_cphi) :: orb_index_cphi
                                    ! Index of the atomic shell to which a
                                    !   Cartesian Gaussian function
                                    !   is associated
    integer,  dimension(maxnorbs_cphi) :: orb_n_cphi
                                    ! Principal quantum number of the 
                                    !   Cartesian Gaussian function
    integer,  dimension(maxnorbs_cphi) :: orb_l_cphi
                                    ! Angular quantum number of the 
                                    !   Cartesian Gaussian function
    real(dp), dimension(ncon_max,maxnorbs_cphi) :: cphi
                                    ! Product of the coefficient of the 
                                    !   primitive Gaussian in the expansion
                                    !   of the radial part of the NAO times the
                                    !   normalization function of the 
                                    !   Cartessian Gaussian function
    !---------------end cphi's data---------------------!

    !---------------sphi's data-------------------------!
    ! In SIESTA, all orbitals are actually spherical
    real(dp),dimension(ncon_max,maxnorbs) :: sphi
                                    ! Transformation of cphi to Spherical 
                                    !   Gaussian functions
    !---------------end sphi's data---------------------!

  end type gto_info_type

!  integer, dimension(:,:), pointer, public, save   ::  D2Sindx => null()

  !> \brief Data type to store a chained list of 4-center integral parameters
  type, public :: eri_link_type

    real(dp) :: gto_eri
    integer  :: four_index(4)
    type(eri_link_type), pointer :: next

  end type eri_link_type

  type(eri_link_type), pointer, public :: head
  type(eri_link_type), pointer, public :: tail
  type(eri_link_type), pointer, public :: ptr

  !> \brief Data type to store relevant cell parameters
  type, public :: hfx_cell_type
    integer :: nsc(3)
    real(dp) :: cell(3)
    real(dp) :: cell_r(3)
  end type hfx_cell_type

end module nao2gto_types
