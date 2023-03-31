! *** Module: nao2gto_data ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief Management of NAO2GTO data
!!
!! This module defines global variables required by the NAO2GTO routines.
!! All of them should ideally be substituted by a design ensuring a
!! seamless data flow throughout the whole program.
!!
!! \author Yann Pouillon
!!
!! eri_prescreen: Following Sec. III A 1 of Ref. \cite Shang:2011, 
!! we are going to compute
!! \f$ \left( \mu \nu \vert \mu \nu \right)_{\rm SR} \f$, where according
!! to the heading of the Section, and Eq. (14) of Ref. \cite Shang:2011,
!! corresponds to
!!
!! \f{eqnarray*}{
!!   \left( \mu \nu \vert \mu \nu \right)_{\rm SR} \equiv
!!   \left(\phi_{\mu}^{\vec{0}} \phi_{\nu}^{\vec{N}} \vert
!!         \phi_{\mu} ^{\vec{0}} \phi_{\nu}^{\vec{N}}\right)_{\rm SR}
!!   = \int \int \frac{\phi_{\mu}(\vec{r}-\vec{R}_{\mu})
!!                     \phi_{\nu}(\vec{r}-\vec{N}-\vec{R}_{\nu}) \:\:
!!                     {\rm erfc}(\omega \vert \vec{r} - \vec{r}^{\prime}\vert)
!!                     \:\: \phi_{\mu}(\vec{r}^{\prime}-\vec{R}_{\mu})
!!                     \phi_{\nu}(\vec{r}^{\prime}-\vec{N}-\vec{R}_{\nu})
!!   }
!!     {\vert \vec{r} - \vec{r}^{\prime} \vert} d^{3}r d^{3}r^{\prime},
!! \f}
!!
!! These matrix elements are computed in calc_prescreen_eri.
!!
!! \par History
!!      - 01.2018 Created [Yann Pouillon]
! *****************************************************************************
module nao2gto_data

  use precision, only: dp
  use nao2gto_libint, only: Libderiv_t, Libint_t
  use nao2gto_types

  implicit none

  private

  integer, public :: hfx_call_counter = 0     ! Number of times that the 
                                              !   subroutine setup_hfx is called

  real(dp), public :: log10_eps_schwarz = 0.0_dp
                                              ! Decimal logarithm of the 
                                              !   Schwarz tolerance.
                                              !   Computed in init_prescreen_eri

  integer, pointer, public, save   ::  D2Sindx(:,:) => null()

  real(dp), pointer, public :: eri_prescreen(:) => null()

  type(Libint_t)        , public :: hfx_libint
  type(Libderiv_t)      , public :: hfx_libderiv
  type(hfx_options_type), public :: hfx_options
  type(hfx_system_type) , public :: hfx_system

  type(pair_list_type), public :: list_ij     ! List of pair of shells that
                                              !    are preselected according 
                                              !    to the Schwarz screening
                                              !    Orbital i is within the
                                              !    home unit cell,
                                              !    Orbital j might be anywhere
                                              !    in the auxiliary supercell
  type(pair_list_type), public :: list_kl     ! List of pair of shells that
                                              !    are preselected according
                                              !    to the Schwarz screening
                                              !    Orbitals k and l might
                                              !    be anywhere in the 
                                              !    auxiliary supercell

  integer, pointer, public :: subshell(:) => null()
                                          ! It has the dimension of the
                                          !  number of orbitals in the supercell
                                          !  and points to the index of the
                                          !  "different" angular momentum shell.
                                          !  For instance, bulk Si
                                          !  (two atom per unit cell)
                                          !  and a DZP basis, we have
                                          !  ten different shells in the 
                                          !  unit cell,  
                                          !  and (ten*number of repetitions on 
                                          !  the unit cell in the supercell)
                                          !  different shells in the supercell.
                                          !  The value of the subshell for each
                                          !  orbital would be
                                          !  Si #1, 3s, first  zeta, subshell=1
                                          !  Si #1, 3s, second zeta, subshell=2
                                          !  Si #1, 3p, first  zeta, subshell=3
                                          !  Si #1, 3p, second zeta, subshell=4
                                          !  Si #1, 3d, first  zeta, subshell=5
                                          !  Si #2, 3s, first  zeta, subshell=6
                                          !  Si #2, 3s, second zeta, subshell=7
                                          !  Si #2, 3p, first  zeta, subshell=8
                                          !  Si #2, 3p, second zeta, subshell=9
                                          !  Si #2, 3d, first  zeta, subshell=10
                                          !  Si #3, 3s, first  zeta, subshell=11
                                          !  ...



  type(hfx_screen_coeff_type), dimension(:,:,:,:,:,:), pointer, public :: &
    pair_dist_radii_pgf => null()
                                          !  Fit the of the precomputed values
                                          !     of the center of the product 
                                          !     densities as a function of the 
                                          !     inter Gaussian distance to a 
                                          !     two-parameter function of the 
                                          !     kind x(1)*rab^2 + x(2)
                                          !     Described in Appendix 5C of 
                                          !     Ref. \cite Guidon:2009

  type(hfx_screen_coeff_type), dimension(:,:,:,:,:,:), pointer, public :: &
    sfc_pgf => null()
  type(hfx_screen_coeff_type), dimension(:,:,:,:), pointer, public :: &
    sfc_shell => null()
  type(hfx_screen_coeff_type), dimension(:,:), pointer, public :: &
    sfc_kind => null()

!  logical, pointer, public :: um_cut(:,:) => null()

                    ! ------------------------------------ !

  ! CO: Cartesian GTOS, SO: Spherical GTOs
  !   1s+3p+6d     = 10
  !   1s+3p+6d+10f = 20
  !   ncosum(lmax)=sum_l (l+1)(l+2)/2 ; l=0,l_max
  !
  ! Consider magnetic quantum number m for RGTOs
  !
  ! nco(l): number of Cartesian Gaussians for a given angular momentum
  ! nso(l): number of Spherical Harmonic Gaussians for a given angular momentum
  ! ncosum(l)   : accumulated number of Cartesian Gaussian functions up to a 
  !               given angular momentum
  ! co(lx,ly,lz): matrix that transforms the three exponents on the cartesian
  !               components of a Cartesian Gaussian function (lx, ly, lz) 
  !               into a single index that identifies the Cartesian Gaussian
  !               function of a given angular momentum
  !               For a given l, co runs between 1 and (l + 1)*(l + 2)/2
  !               For instance:
  !               - l = 0 => lx = ly = lz = 0 
  !                     => only one Cartesian Gaussian 
  !                     co( 0, 0, 0 ) = 1       coset( 0, 0, 0 ) = 1
  !               - l = 1 => there are three Cartesian Gaussians, with
  !                     lx = 1, ly = 0, lz = 0 
  !                     lx = 0, ly = 1, lz = 0 
  !                     lx = 0, ly = 0, lz = 1 
  !                     co( 1, 0, 0 ) = 1       coset( 1, 0, 0 ) = 2
  !                     co( 0, 1, 0 ) = 2       coset( 0, 1, 0 ) = 3
  !                     co( 0, 0, 1 ) = 3       coset( 0, 0, 1 ) = 4 
  !               - l = 2 => there are six Cartesian Gaussians, with
  !                     lx = 2, ly = 0, lz = 0 
  !                     lx = 1, ly = 1, lz = 0 
  !                     lx = 1, ly = 0, lz = 1 
  !                     lx = 0, ly = 2, lz = 0 
  !                     lx = 0, ly = 1, lz = 1 
  !                     lx = 0, ly = 0, lz = 2 
  !                     co( 2, 0, 0 ) = 1       coset( 2, 0, 0 ) = 5
  !                     co( 1, 1, 0 ) = 2       coset( 1, 1, 0 ) = 6
  !                     co( 1, 0, 1 ) = 3       coset( 1, 0, 1 ) = 7
  !                     co( 0, 2, 0 ) = 4       coset( 0, 2, 0 ) = 8
  !                     co( 0, 1, 1 ) = 5       coset( 0, 1, 1 ) = 9
  !                     co( 0, 0, 2 ) = 6       coset( 0, 0, 2 ) = 10
  ! coset(lx,ly,lz): matrix that transforms the three exponents on the cartesian
  !                  into a single index that identifies the Cartesian Gaussian
  !                  function within the total number of Cartesian Gaussian 
  !                  functions in the basis (see example above)
  ! indco:           inverse of coset. This is defined only up to f-orbitals 
  !                     (l=3)
  !                  In this array we enter with an integer number 
  !                  from 1 to the total number of cartesian Gaussian functions
  !                  (second index of indco), and it returns the lx, ly, and lz,
  !                  depending on the value we introduce in the first index

  integer, save, target, public :: &
    nco(0:l_max), &
    nso(0:l_max), &
    ncosum(-1:l_max), &
    co(0:l_max,0:l_max,0:l_max), &
    coset(0:l_max,0:l_max,0:l_max), &
    indco(3,20)

                    ! ------------------------------------ !

  ! Information about GTOs, to use along species_info
  type(gto_info_type), target, allocatable, save, public :: hfx_gtos(:)

end module nao2gto_data
