! *** Module: nao2gto_io ***

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! *****************************************************************************
!> \brief I/O module for Gaussian-based Hartree-Fock exchange
!!
!!  This module bridges the input file of SIESTA with the NAO2GTO routines,
!!  which calculate the Hartree-Fock exchange interaction using Gaussians.
!!
!! \note
!!      This file currently works with a version of Libint configured for
!!      LIBINT_MAX_AM=5 and LIBINT_MAX_AM1=4 only.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \copyright
!!      - 2010-2018 SIESTA Developers Group
!!
!! \par History
!!      - 11.2017 Reviewed for inclusion in SIESTA [Xinming Qin]
!!      - 01.2018 Brought together from separate files [Yann Pouillon]
! *****************************************************************************
module nao2gto_io

  implicit none

  private

  public :: &
    nao2gto_dump_system, &
    nao2gto_transfer

contains

  ! ***************************************************************************
  ! *** Public routines                                                     ***
  ! ***************************************************************************

  ! ***************************************************************************
  !> \brief Displays a summary of a hfx_system_type data structure
  !!
  !! \param[in] hfx_sys: data structure to display
  ! ***************************************************************************
  subroutine nao2gto_dump_system(hfx_sys)

    use nao2gto_types, only: hfx_system_type

    implicit none

    ! Arguments
    type(hfx_system_type), intent(in) :: hfx_sys

    ! -------------------------------------------------------------------------

    if ( associated(hfx_sys%maxnh) ) then
      write(*, fmt='(2X,A,": ",I8)') "maxnh", hfx_sys%maxnh
    else
      write(*, fmt='(2X,A,": ",A8)') "maxnh", "null"
    end if
    if ( associated(hfx_sys%na) ) then
      write(*, fmt='(2X,A,": ",I8)') "na", hfx_sys%na
    else
      write(*, fmt='(2X,A,": ",A8)') "na", "null"
    end if
    if ( associated(hfx_sys%norb) ) then
      write(*, fmt='(2X,A,": ",I8)') "norb", hfx_sys%norb
    else
      write(*, fmt='(2X,A,": ",A8)') "norb", "null"
    end if
    if ( associated(hfx_sys%nspin) ) then
      write(*, fmt='(2X,A,": ",I8)') "nspin", hfx_sys%nspin
    else
      write(*, fmt='(2X,A,": ",A8)') "nspin", "null"
    end if
    if ( associated(hfx_sys%nua) ) then
      write(*, fmt='(2X,A,": ",I8)') "nua", hfx_sys%nua
    else
      write(*, fmt='(2X,A,": ",A8)') "nua", "null"
    end if
    if ( associated(hfx_sys%nuo) ) then
      write(*, fmt='(2X,A,": ",I8)') "nuo", hfx_sys%nuo
    else
      write(*, fmt='(2X,A,": ",A8)') "nuo", "null"
    end if
    if ( associated(hfx_sys%nuotot) ) then
      write(*, fmt='(2X,A,": ",I8)') "nuotot", hfx_sys%nuotot
    else
      write(*, fmt='(2X,A,": ",A8)') "nuotot", "null"
    end if
    if ( associated(hfx_sys%iaorb) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "iaorb", size(hfx_sys%iaorb, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "iaorb", "null"
    end if
    if ( associated(hfx_sys%indxua) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "indxua", size(hfx_sys%indxua, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "indxua", "null"
    end if
    if ( associated(hfx_sys%iphorb) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "iphorb", size(hfx_sys%iphorb, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "iphorb", "null"
    end if
    if ( associated(hfx_sys%isa) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "isa", size(hfx_sys%isa, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "isa", "null"
    end if
    if ( associated(hfx_sys%listh) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "listh", size(hfx_sys%listh, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "listh", "null"
    end if
    if ( associated(hfx_sys%listhptr) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "listhptr", size(hfx_sys%listhptr, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "listhptr", "null"
    end if
    if ( associated(hfx_sys%nsc) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "nsc", size(hfx_sys%nsc, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "nsc", "null"
    end if
    if ( associated(hfx_sys%numh) ) then
      write(*, fmt='(2X,A,": array(",I8,")")') "numh", size(hfx_sys%numh, 1)
    else
      write(*, fmt='(2X,A,": ", A8)') "numh", "null"
    end if
    if ( associated(hfx_sys%xa) ) then
      write(*, fmt='(2X,A,": array(",I8,",",I8")")') "xa", &
        size(hfx_sys%xa, 1), size(hfx_sys%xa, 2)
    else
      write(*, fmt='(2X,A,": ", A8)') "xa", "null"
    end if

  end subroutine nao2gto_dump_system

! *****************************************************************************
!> \brief Reads coefficients and builds the corresponding spherical Gaussians
!!
!! This routine reads the coefficients of Numerical Atomic Orbitals (NAO) onto
!! Gaussian-Type Orbitals (GTO) from the FDF input of SIESTA and builds the
!! corresponding spherical Gaussians.
!!
!! This is called from siesta_init, and only by the input/output node.
!!
!! \author Honghui Shang
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 12.2010 Imported and connected to atm_transfer [Honghui Shang]
!!      - 12.2013 Modified for input from FDF format [Xinming Qin]
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[in] l_max: Maximum angular momentum
!! \param[in] co: matrix that transforms the three exponents on the cartesian
!!                components of a Cartesian Gaussian function (lx, ly, lz)
!!                into a single index that identifies the Cartesian Gaussian
!!                function of a given angular momentum
!!                For a given l, co runs between 1 and (l + 1)*(l + 2)/2
!! \param[out] orbtramat: Coefficients of the linear transformation
!!                from cartesian Gaussian functions
!!                to real spherical harmonic Gaussian functions
! *****************************************************************************
  subroutine nao2gto_transfer(gtos, hfx_opts)
    use precision, only: dp
    use units,     only: pi
    use alloc,     only: re_alloc, de_alloc
    use atm_types, only: maxnorbs, nspecies
    use atm_types, only: maxn_orbnl, species, species_info
    use atmfuncs,  only: lofio,mofio,labelfis, rcut
    use atomlist,  only: rmaxo
    use basis_specs, only: label2species
    use chemical,  only: species_label
    use m_io,      only: io_assign
    use radial,    only: rad_get
    use sys,       only: die
    use fdf
    use parsing
    !FIXME: restore after debugging
    !use nao2gto_data, only: co, coset, indco, nco, ncosum, nso
    use nao2gto_types, only: maxn_contract, maxnorbs_cphi, l_max, ncon_max
    use nao2gto_data
    use nao2gto_nonlin
    use nao2gto_types, only: gto_info_type, hfx_options_type
    use nao2gto_utils, only: dfac, fac, exp_radius, exp_radius_very_extended
    use nao2gto_transform, only: calc_c2s_matrix, cphi2sphi, orbtramat_type
    use nao2gto_wrappers, only: nao2gto_libint_dump, nao2gto_libderiv_dump

    implicit none

    ! Arguments
    type(gto_info_type), intent(inout)  :: gtos(nspecies)
    type(hfx_options_type), intent(out) :: hfx_opts

    ! Local variables
    character(len=132) :: msg
    logical :: do_fit     ! If true, the fitting of the NAO into a linear
                          !   expansions of Gaussians will be done within SIESTA
                          ! If false, the exponents and coefficients of the 
                          !   Gaussian expansion will be read from the 
                          !   NAO2GTO block in the FDF file
    logical :: do_stop    ! Variable that determines if fitting is incomplete
                          !   If true at the end of the subroutine, the program
                          !   will die with an error message
    integer  :: errno, fit_fd, isp, io, iorbrd, inlz, l, m, n, jnlz, &
&     ipgf, irad, itry, jpgf, npts, lx, ly, lz, nco_max, nso_max, &
&     ico, ifdfblk, i_cphi, ren, rel, rezeta, recoeffs, npt
    integer  :: num_cphi(maxnorbs)
    real(dp) :: expzet, fnorm, zeta, zetb, prefac, gcca, gccb, &
&     gauss_resid, dummy, fit_rad, fit_orb, fit_amin, fit_amax, fit_sig3, &
&     new_sig3
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    type(species_info), pointer :: spp => null()

    logical, dimension(:),      pointer :: chk_coeffs => null()
    integer, dimension(:,:),    pointer :: hfx_contract => null()
!                                          Actual number of Gaussians in the
!                                          expansion of a given shell of NAO
!                                          First index: shell of NAO
!                                          Second index: atomic species
!                 
    real(dp), dimension(:,:),   pointer :: fit_trial => null()
!                                          Temporary array to store the
!                                          exponents and coefficients of the
!                                          expansions of the Gaussians as they
!                                          come from nao2gto_gaussfit
    real(dp), dimension(:,:,:), pointer :: nao2gto_zeta => null() 
!                                          Exponents of the Gaussians that 
!                                          fit a particular NAO 
!                                          First index:  number of the Gaussian 
!                                              in the expansion
!                                          Second index: shell of NAO
!                                          Third  index: atomic species
    real(dp), dimension(:,:,:), pointer :: nao2gto_coefficient => null()
!                                          Coefficients of the Gaussians that
!                                          fit a particular NAO 
!                                          First index:  number of the Gaussian 
!                                              in the expansion
!                                          Second index: shell of NAO
!                                          Third  index: atomic species
    real(dp), dimension(:),     pointer :: rad_pts => null()
    real(dp), dimension(:),     pointer :: orb_pts => null()
    real(dp), dimension(:),     pointer :: fval => null()
    real(dp), dimension(:),     pointer :: orb_dr => null()

    type(orbtramat_type), dimension(:), pointer :: orbtramat => null()
    real(dp):: orbital_rcut, rmax
   

    ! -------------------------------------------------------------------------

    write(*,'(A,/)') "************************ Begin: HYBRID XC INITIALIZATION **********************"

! ------------------------------------------------------------------------------
!> Step 1: Initialization of nco, nso, co, coset, orbtramat, indco,
!!         hfx_contract, nao2gto_zeta, and nao2gto_coefficient
!!         The different variables to run an EXX calculation are read
!!         from the FDF input file calling the subroutine nao2gto_read_options
! ------------------------------------------------------------------------------
    ncosum(-1:l_max) = 0

!   Loop on angular momentum
    do l = 0, l_max

!     Compute the number of Cartesian Gaussians for a given angular momentum
      nco(l) = (l + 1)*(l + 2)/2

!     Compute the number of Spherical Harmonic Gaussians for a given angular
!     momentum
      nso(l) = 2*l + 1

!     Compute the accumulated number of Cartesian Gaussian functions up to a 
!     given angular momentum
      ncosum(l) = ncosum(l-1) + nco(l)

!!     For debugging
!      write(6,'(a,6i5)')                                 &
! &      'nao2gto_transfer: l, nco, nso, ncosum = ',      &
! &      l, nco(l), nso(l), ncosum(l) 
!!     End debugging

    end do

!   Compute the matrices that tranform the exponents of the cartesian components
!   of a Cartesian Gaussian function into a single index
!   co(lx,ly,lz): matrix that transforms the three exponents on the cartesian
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
!   coset(lx,ly,lz): matrix that transforms the three exponents on the cartesian
!                  into a single index that identifies the Cartesian Gaussian
!                  function within the total number of Cartesian Gaussian
!                  functions in the basis (see example above)
    do lx = 0, l_max
      do ly = 0, l_max
        do lz = 0, l_max
          l = lx + ly + lz
          if ( l > l_max ) cycle
          co(lx,ly,lz) = 1 + (l - lx)*(l - lx + 1)/2 + lz
          coset(lx,ly,lz) = ncosum(l-1) + co(lx,ly,lz)
!!         For debugging
!          write(6,'(a,6i5)')                                 &
! &          'nao2gto_transfer: lx, ly, lz, l, co, coset = ', &
! &            lx, ly, lz, l, co(lx,ly,lz), coset(lx,ly,lz) 
!!         End debugging
        end do
      end do
    end do

!   Define the inverse of coset.
!   This is the array indco
!   In this array we enter with an integer number from 1 to the total number
!   of cartesian Gaussian functions (second index of indco),
!   and it returns the lx, ly, and lz, depending on the value we introduce
!   in the first index
    indco(:,:) = 0
    do l = 0, l_max
      do lx = 0, l
        do ly = 0, l-lx
          lz = l - lx - ly
          indco(1:3,coset(lx,ly,lz)) = (/lx,ly,lz/)
!!         For debugging
!          write(6,'(a,8i5)')                                    &
! &          'nao2gto_transfer: lx, ly, lz, l, coset, indco = ', &
! &            lx, ly, lz, l, coset(lx,ly,lz), indco(1:3,coset(lx,ly,lz))
!!         End debugging
        end do
      end do
    end do

!   Compute the coefficients of the linear transformation
!   from cartesian Gaussian functions
!   to real spherical harmonic Gaussian functions.
!   One cannot use re_alloc with orbtramat, because it is a vector
!   of structured types
    allocate(orbtramat(0:l_max))
    do l=0,l_max
      nco_max = nco(l)
      nso_max = nso(l)
      write(msg,'("orbtramat(",I1,")%c2s")') l
      nullify(orbtramat(l)%c2s)
      call re_alloc(orbtramat(l)%c2s, 1, nso_max, 1, nco_max, &
&       name=trim(msg), routine='nao2gto_transfer')
    enddo
    call calc_c2s_matrix(l_max, co, orbtramat)

!   Allocate the variable that will contain the actual number of Gaussians
!   that fit a given shell of NAO
    call re_alloc(hfx_contract, 1, maxn_orbnl, 1, nspecies, &
&     "nao2gto_transfer")

!   Allocate the variable that will contain the exponents of the Gaussians
!   that fit a NAO
    call re_alloc(nao2gto_zeta, 1, maxn_contract, 1, maxn_orbnl, &
&     1, nspecies, "nao2gto_transfer")

!   Allocate the variable that will contain the coefficients of the Gaussians
!   that fit a NAO
    call re_alloc(nao2gto_coefficient, 1, maxn_contract, 1, maxn_orbnl, &
&     1, nspecies, "nao2gto_transfer")

!   Initialize the variables
    hfx_contract(:,:)          = 0
    nao2gto_zeta(:,:,:)        = 0.0_dp
    nao2gto_coefficient(:,:,:) = 0.0_dp

!   Read Hartree-Fock exchange parameters from FDF input file
    call nao2gto_read_options(hfx_opts)

! ------------------------------------------------------------------------------
!> Step 2: Read Hartree-Fock exchange parameters from FDF input file
!!
!! \note See ldau_specs.f to understand how to read blocks
! ------------------------------------------------------------------------------

    do_fit = .true.

!   If present, read the exponents and the coefficients of the Gaussians 
!   from the block NAO2GTO in the FDF file
    if (fdf_block('NAO2GTO', bfdf)) then

      do_fit = .false.

      ifdfblk = 1
      do while(fdf_bline(bfdf,pline))

        !> Step 2.a(nofit): Read species
        if ( .not. fdf_bmatch(pline,'ni') ) then
          write(msg,'(A," (line ",I4,")")') &
&           'Wrong format in NAO2GTO', ifdfblk
          call die(trim(msg))
        endif
        isp = label2species(fdf_bnames(pline,1))
        if (isp .eq. 0) then
          write(*,'(a,1x,a)') &
&           'WRONG species symbol in NAO2GTO:', &
&         trim(fdf_bnames(pline,1))
          call die()
        endif
        spp => species(isp)
        ifdfblk = ifdfblk + 1

        !> Step 2.b(nofit): Prepare variables that will tell whether
        !! all coefficients for all orbitals have been read
        call re_alloc(chk_coeffs, 1, spp%n_orbnl, "nao2gto_transfer")
        chk_coeffs(:) = .false.

        !> Step 2.c(nofit): Read data for each (n, l, zeta) triplet
        !! \note Coefficients can be provided in any order.
        do iorbrd=1,spp%n_orbnl

          !> Step 2.c.1(nofit): Read information about which orbital the
          !! coefficients correspond to
          if (.not. fdf_bline(bfdf, pline)) &
&           call die('Not enough information on the Gaussian expansion')
          if (fdf_bmatch(pline,'iiii')) then
            ren = fdf_bintegers(pline,1)
            rel = fdf_bintegers(pline,2)
            rezeta = fdf_bintegers(pline,3)
            recoeffs = fdf_bintegers(pline,4)
!          write(6,*) iorbrd, ren, rel, rezeta, recoeffs
          else
            write(msg,'(A," (line ",I4,")")') &
&             'Wrong format in NAO2GTO', ifdfblk
            call die(trim(msg))
          endif

          !> Step 2.c.2(nofit): Translate n, l, zeta into an orbital index
!          write(6,*)"ttt", ren, rel, rezeta
          inlz = orb_nlz_to_index(spp, ren, rel, rezeta)
!          write(6,*)"ttt", ren, rel, rezeta, inlz
          if ( inlz == -1 ) then
            write(msg,'(A,3(1X,I2),A," (line ",I4,")")') &
&             'Could not find orbital(', ren, rel, rezeta, ') in SIESTA', &
&             ifdfblk
            call die(trim(msg))
          end if
          hfx_contract(inlz,isp) = recoeffs
          ifdfblk = ifdfblk + 1

          !> Step 2.c.3(nofit): Check for duplicates
          if ( .not. chk_coeffs(inlz) ) then
            chk_coeffs(inlz) = .true.
          else
            write(msg,'(3(A),3(1X,I2),A," (line ",I4,")")') &
&             'Duplicate coefficients for ', trim(spp%symbol), ' orbital(', &
&             ren, rel, rezeta,') in NAO2GTO block', ifdfblk
            call die(trim(msg))
          end if

          !> Step 2.c.3(nofit): Read Gaussian coefficients
          do jnlz=1,hfx_contract(inlz,isp)
            if (.not. fdf_bline(bfdf, pline)) &
&               call die('Not enough information on the Gaussian expansion')
            if (fdf_bmatch(pline,'vv')) then
              nao2gto_zeta(jnlz,inlz,isp) = fdf_bvalues(pline,1)
              nao2gto_coefficient(jnlz,inlz,isp) = fdf_bvalues(pline,2)
            else
              write(msg,'(A," (line ",I4,")")') &
&                 'Wrong format in NAO2GTO', ifdfblk
              call die(trim(msg))
            endif
            ifdfblk = ifdfblk + 1
          enddo

        enddo ! iorbrd=1,spp%n_orbnl

        !> Step 2.d(nofit): Check that all coefficients for the
        !! current species have actually been read
        do iorbrd=1,spp%n_orbnl
          if ( .not. chk_coeffs(iorbrd) ) then
            write(msg,'(A,3(1X,I2),A," (line ",I4,")")') &
&             'Missing coefficients for orbital(', ren, rel, rezeta, &
&             ') in NAO2GTO block', ifdfblk
            call die(trim(msg))
          end if
        end do
        call de_alloc(chk_coeffs, 'chk_coeffs', 'nao2gto_transfer')

      enddo ! while fdf_bline(bfdf,pline)

    else    ! If the fitting will be done within SIESTA

      !> Step 2.a(fit): Prepare data stores
      call re_alloc(rad_pts, 1, hfx_opts%npts_fit, 'rad_pts', &
        'nao2gto_transfer')
      call re_alloc(orb_pts, 1, hfx_opts%npts_fit, 'orb_pts', &
        'nao2gto_transfer')
!     Allocate the temporal variable that will contain the 
!     coefficients of the Gaussians that fit a NAO
      call re_alloc(fit_trial, 1, 2, 1, hfx_opts%max_num_gaus, 'nao2gto_transfer')

!     Loop over atomic species
      do isp = 1, nspecies

        spp => species(isp)
        spp%label = species_label(isp)

        inlz = 0

        do io = 1, spp%norbs
!       Identify the angular momentum quantum number
        l = spp%orb_l(io)
        m = spp%orb_m(io)
!       Identify the magnetic quantum number
!       Only the first time that an orbital of a given shell appears is
!       considered
        if ( m .ne. -l ) cycle   ! not a normal orbital

        inlz = inlz + 1
!       Store in the variable gto the number of Gaussians in the expansion
!       Loop over different NAO shells in a given atomic species
!        do inlz = 1, spp%n_orbnl
          !> Step 2.b(fit): Translate rad_func representation of the orbital
          do irad=1,hfx_opts%npts_fit
            rad_pts(irad) = spp%orbnl(inlz)%cutoff * real(irad-1, dp) / &
&                   (real(hfx_opts%npts_fit, dp) - 1)
            call rad_get(spp%orbnl(inlz), rad_pts(irad), orb_pts(irad), dummy)
          enddo
          orb_pts(hfx_opts%npts_fit) = 0.0_dp

          !> Step 2.c(fit): Perform Gaussian fitting and filter results
          fit_trial(:,:) = 0.0_dp
          if( l == 0)  hfx_opts%max_num_gaus = 5
          if( l > 0)  hfx_opts%max_num_gaus = 4
          if (spp%orbnl_ispol(inlz)) hfx_opts%max_num_gaus = 3
          write(6,*) inlz, l, spp%orbnl_ispol(inlz), hfx_opts%max_num_gaus

          call nao2gto_gaussfit( hfx_opts%max_num_gaus, rad_pts, orb_pts, &
 &                               hfx_opts, fit_trial, errno )
          ipgf = 0
          do jnlz = 1, hfx_opts%max_num_gaus
!!           For debugging
            write(6,'(a,3i5,3f12.5,i5,l5)')                                 &
 &            'nao2gto_transfer: species, orbital, gauss, exp, coef, epsilon, error = ', &
 &            isp, inlz, jnlz, fit_trial(1,jnlz), fit_trial(2,jnlz),     &
 &            epsilon(1.0_dp), errno, abs(fit_trial(1,jnlz)) < epsilon(1.0_dp)
!!           End debugging
            if ( abs(fit_trial(1,jnlz)) < epsilon(1.0_dp) ) exit
            ipgf = jnlz
          end do
          if ( ipgf == 0 ) errno = 10

          if ( errno == 0 ) then
              hfx_contract(inlz,isp) = ipgf
              do jnlz = 1, ipgf
                nao2gto_zeta(jnlz,inlz,isp)        = fit_trial(1,jnlz)
                nao2gto_coefficient(jnlz,inlz,isp) = fit_trial(2,jnlz)
              end do
          else
            write(*, fmt='("GAUSSIAN FITTING ERROR ",I4.4)') errno
          end if

!!         For debugging
!          write(6,'(a,3i5)')                                          &
! &          'nao2gto_transfer: isp, inlz, hfx_contract(inlz,isp) = ', &
! &           isp, inlz, hfx_contract(inlz,isp) 
!          do jnlz = 1, hfx_contract(inlz,isp)
!            write(6,'(a,3i5,2f12.5)')                                   &
! &            'nao2gto_transfer: isp, inlz, jnlz, zeta, coeff = ',      &
! &             isp, inlz, jnlz, nao2gto_zeta(jnlz,inlz,isp),            &
! &             nao2gto_coefficient(jnlz,inlz,isp) 
!          enddo 
!!         End debugging

        end do   ! End loop over different atomic orbital shells for a specie

!!       For debugging
!        write(6,'(a)')' End of the fitting of the atomic orbitals'
!        call die()
!!       End debugging

      end do   ! End loop over atomic species

      !> Step 2.d(fit): Clean-up the mess
      call de_alloc( rad_pts,   'rad_pts',   'nao2gto_transfer' )
      call de_alloc( orb_pts,   'orb_pts',   'nao2gto_transfer' )
      call de_alloc( fit_trial, 'fit_trial', 'nao2gto_transfer' )

    endif   ! NAO2GTO block

! ------------------------------------------------------------------------------
!> Step 3: Process input data for each species
! ------------------------------------------------------------------------------
    do isp=1,nspecies
      spp => species(isp)
      inlz = 0

      !> Step 3.a: Store GTOs for normal orbitals (m == -l)
      !!
      !! \bug Did not check the use of negative Z numbers
      !!      in original implementation

!     Loop on all the atomic orbitals of a given species
      do io = 1, spp%norbs
!       Identify the angular momentum quantum number
        l = spp%orb_l(io)
!       Identify the magnetic quantum number
        m = spp%orb_m(io)

!       Only the first time that an orbital of a given shell appears is considered
        if ( m .ne. -l ) cycle   ! not a normal orbital

        inlz = inlz + 1
!       Store in the variable gto the number of Gaussians in the expansion
        gtos(isp)%orbnl_contract(inlz) = hfx_contract(inlz,isp)
        do ipgf=1,gtos(isp)%orbnl_contract(inlz)
!         Store in the variable gto the exponents of the Gaussians
          gtos(isp)%orbnl_zeta(ipgf,inlz) = nao2gto_zeta(ipgf,inlz,isp)
!         Store in the variable gto the coefficients of the Gaussians
          gtos(isp)%orbnl_coefficient(ipgf,inlz) = nao2gto_coefficient(ipgf,inlz,isp)
        enddo
!       Store in the variable gto the minimum value of the exponent of the gaussian (zeta),
!       (i.e. the most diffuse gto)
        gtos(isp)%orbnl_adjoined_zeta(inlz) = &
&         minval(gtos(isp)%orbnl_zeta(1:gtos(isp)%orbnl_contract(inlz),inlz))
!!       For debugging
!        write(6,'(a,3i5,f12.5)')                                              &
! &       'nao2gto_transfer: isp, inlz, orbnl_contract, orbnl_adjoined_zeta = ',&
! &         isp, inlz, gtos(isp)%orbnl_contract(inlz),                          &
! &         gtos(isp)%orbnl_adjoined_zeta(inlz)
!!       End debugging
      enddo ! endloop on the orbitals of a given species.

      !> Step 3.b: Add cutoff radius for primitive GTOs and contracted
      !! GTOs, actually PAO. PAO itself has cutoff radius in siesta, we
      !! use the shell cutoff just for comparison
      !! (added by Xinming Qin, Oct. 2018).
      gtos(isp)%pgf_radius(1:maxn_contract,1:maxn_orbnl) = 0.0_dp
      gtos(isp)%shell_radius(1:maxn_orbnl) = 0.0_dp
      inlz = 0
      do io=1,spp%norbs
        orbital_rcut = rcut(isp,io)
        l = spp%orb_l(io)
        m = spp%orb_m(io)
        if ( m /= -l ) cycle   ! Not a normal orbital
        inlz = inlz+1
        do ipgf=1,gtos(isp)%orbnl_contract(inlz)
          gcca = nao2gto_coefficient(ipgf,inlz,isp)
          zeta = nao2gto_zeta(ipgf,inlz,isp)
!         Compute the point where 
!         Gaussian(r) = coefficient * r**l * exp(-zeta*r**2) is smaller than
!         a given threshold, here 1.d-5
!          gtos(isp)%pgf_radius(ipgf,inlz)= exp_radius(l, zeta, 1.d-10, gcca)
!          if(zeta.le.0.1) then
!             gtos(isp)%pgf_radius(ipgf,inlz)= exp_radius(l, zeta, 0.5d-3, 1.0_dp )
!          else
            gtos(isp)%pgf_radius(ipgf,inlz)= exp_radius(l, zeta, hfx_opts%gto_eps, 1.0_dp ) ! gcca) !1.0_dp)
!            gtos(isp)%pgf_radius(ipgf,inlz)= exp_radius(l, zeta, 1.0d-4, 1.0_dp ) !
!          endif
!!         For debugging
          write(6,'(a,3i5,f12.5)')                                           &
 &          'nao2gto_transfer: isp, inlz, ipgf, pgf_radius = ',              &
 &            isp, inlz, ipgf, gtos(isp)%pgf_radius(ipgf,inlz)                
!!         End debugging
        enddo
        gtos(isp)%shell_radius(inlz) = &
&         maxval(gtos(isp)%pgf_radius(1:gtos(isp)%orbnl_contract(inlz),inlz))
!!       For debugging
        write(6,'(a,2i5,2f12.5)')                                              &
 &        'nao2gto_transfer: isp, inlz, shell_radius = ',                     &
 &          isp, inlz, gtos(isp)%shell_radius(inlz), orbital_rcut                  
!       End debugging
      enddo
      gtos(isp)%kind_radius = maxval(gtos(isp)%shell_radius(1:maxn_orbnl))
      !rmax = 
!!     For debugging
      write(6,'(a,1i5,f12.5)')                                              &
 &      'nao2gto_transfer: isp, kind_radius = ',                 &
 &        isp, gtos(isp)%kind_radius
!!     End debugging

      !> Step 3.c: Compute the number of cphi coefficients for each orbital
      !! and store their total number in i_cphi
      num_cphi(:) = 0
      i_cphi = 1

!     Loop on all the atomic orbitals of a given species
      do io=1,spp%norbs

!!       For debugging
!        write(6,'(a,2i5)')'nao2gto_transfer: io, orb_m = ', &
! &        io, spp%orb_m(io)
!!       End debugging

        if ( spp%orb_m(io) .lt. 2 ) then
          num_cphi(i_cphi) = io
          i_cphi = i_cphi + 1
        else if ( spp%orb_m(io) .eq. 2 ) then
          num_cphi(i_cphi) = io
          num_cphi(i_cphi+1) = io
          i_cphi = i_cphi + 2
        else if ( spp%orb_m(io) .eq. 3 ) then
          num_cphi(i_cphi) = io
          num_cphi(i_cphi+1) = io
          num_cphi(i_cphi+2) = io
          num_cphi(i_cphi+3) = io
          i_cphi = i_cphi + 4
        endif
      enddo
      i_cphi = i_cphi - 1

!!     For debugging
!      write(6,'(a,i5)')'nao2gto_transfer: i_cphi = ', i_cphi 
!      do io = 1, i_cphi
!        write(6,'(a,2i5)')                                        &
! &        'nao2gto_transfer: i_cphi, num_cphi = ',                &
! &        io, num_cphi(io)
!      enddo 
!!     End debugging

      !> Step 3.d: Initialize indices for cartesian coefficients
      !!
      !! \note Here, we only consider the lmax <= 3 case (s, p, d, and f)
!     Set up the number of Cartesian Gaussian functions per atomic species
!     If the maximum angular momentum is smaller than l = 2 
!     (i.e. only s and p orbitals) the number of Cartesian Gaussian functions
!     is the same as the number of real spherical harmonics, 
!     (i.e. equal to the traditional number of NAO.
!     However if the maximum angular momentum is equal or larger than l = 2
!     (i.e. we are considering d or f orbitals) 
!     the number of Cartesian Gaussian functions is larger than the number
!     of real spherical harmonics:
!     Six CG instead of 5 real spherical harmonics for the d-shell
!     Ten CG instead of 7 real spherical harmonics for the f-shell
!  
      if ( spp%lmax_basis .lt. 2 ) then
        gtos(isp)%norbs_cphi = spp%norbs   !< Coefficients for s and p orbitals
      else
        gtos(isp)%norbs_cphi = i_cphi      !< Coefficients for d and f orbitals
      endif

      do i_cphi = 1,gtos(isp)%norbs_cphi
        io = num_cphi(i_cphi)
!       Identify the principal quantum number of the Cartesian Gaussian function
        gtos(isp)%orb_n_cphi(i_cphi) = spp%orb_n(io)
!       Identify the angular quantum number of the Cartesian Gaussian function
        gtos(isp)%orb_l_cphi(i_cphi) = spp%orb_l(io)
!       Identify the shell to which this Cartesian Gaussian function is 
!       associated
        gtos(isp)%orb_index_cphi(i_cphi) = spp%orb_index(io)
!!       For debugging
!        write(*,'(A,A,5(1X,A,"=",I4))') &
!          "[DEBUG][orb_index] ", trim(spp%label), &
!          "norbs_cphi", gtos(isp)%norbs_cphi, &
!          "io", io, &
!          "n", spp%orb_n(io), &
!          "l", spp%orb_l(io), &
!          "orb_index", spp%orb_index(io)
!!       End debugging
      enddo

      io = 0
!     Loop over all the different shells of NAO
!     i.e. For Si with a DZP basis, spp%n_orbnl = 5
!     1. First  zeta for the 3s shell (l = 0)
!     2. Second zeta for the 3s shell (l = 0)
!     3. First  zeta for the 3p shell (l = 1)
!     4. Second zeta for the 3p shell (l = 1)
!     5. First  zeta for the 3d shell (l = 2)
      do inlz = 1,spp%n_orbnl
        l = spp%orbnl_l(inlz)

!       Indentify the indices where the 
!       first cartesian gaussian function (gtos(isp)%orbnl_index_cphi(inlz)) 
!       or the 
!       first spherical gaussian function (gtos(isp)%orbnl_index_sphi(inlz)) 
!       appear
        if ( inlz .eq. 1 ) then
          gtos(isp)%orbnl_index_cphi(inlz) = 1
          gtos(isp)%orbnl_index_sphi(inlz) = 1
        else
          gtos(isp)%orbnl_index_cphi(inlz) = gtos(isp)%orbnl_index_cphi(inlz-1) &
&                                    + nco(spp%orbnl_l(inlz-1))
          gtos(isp)%orbnl_index_sphi(inlz) = gtos(isp)%orbnl_index_sphi(inlz-1) &
&                                    + nso(spp%orbnl_l(inlz-1))
        endif

!!       For debugging
!        write(6,'(a,4i5)')                                                  &
! &        'nao2gto_transfer: inlz, l, gtos(isp)%orbnl_index_cphi(inlz), gtos(isp)%orbnl_index_sphi(inlz) = ', &
! &        inlz, l, gtos(isp)%orbnl_index_cphi(inlz),                        &
! &        gtos(isp)%orbnl_index_sphi(inlz)
!!       End debugging

!       Loop over all the cartesian gaussian functions of a given shell
        do ico=ncosum(l-1)+1,ncosum(l)
          io = io + 1
!         Identify the exponents lx, ly, and lz of the cartesian gaussian
!         functions
          gtos(isp)%orb_cphi_lx(io) = indco(1,ico)
          gtos(isp)%orb_cphi_ly(io) = indco(2,ico)
          gtos(isp)%orb_cphi_lz(io) = indco(3,ico)
!!         For debugging
!          write(6,'(a,6i5)')                                                 &
! &          'nao2gto_transfer: inlz, l, ico, lx, ly, lz = ',                 &
! &          inlz, l, ico, gtos(isp)%orb_cphi_lx(io), gtos(isp)%orb_cphi_ly(io),&
! &          gtos(isp)%orb_cphi_lz(io)
!!         End debugging
        enddo
      enddo

      !> Step 3.e: Compute norms of cartesian GTOs
!     Loop over all the cartesian gaussian orbitals of a given atomic species
      do io = 1, gtos(isp)%norbs_cphi
!       Indentify the angular momentum
        l = gtos(isp)%orb_l_cphi(io)
        expzet = 0.5_dp*REAL(2*l + 3,dp)
        fnorm = 0.0_dp
!       Indentify the shell to which this cartesian gaussian function is 
!       associated
        inlz = gtos(isp)%orb_index_cphi(io)
!!       For debugging
!        write(6,'(a,2i5,f12.5,i5)')                                          &
! &        'nao2gto_transfer: io, l, expzet, inlz  = ',                       &
! &                           io, l, expzet, inlz  
!!       End debugging
!       Loop over all the primitive gaussians that enter in the linear 
!       combination of the radial part of the NAO for that given shell
        do ipgf = 1, gtos(isp)%orbnl_contract(inlz)
!         Identify the coefficient of the primitive gaussian
          gcca = gtos(isp)%orbnl_coefficient(ipgf,inlz)
!         Identify the exponent of the primitive gaussian
          zeta = gtos(isp)%orbnl_zeta(ipgf,inlz)
          do jpgf=1,gtos(isp)%orbnl_contract(inlz)
            gccb = gtos(isp)%orbnl_coefficient(jpgf,inlz)
            zetb = gtos(isp)%orbnl_zeta(jpgf,inlz)
            fnorm = fnorm + gcca*gccb/((zeta + zetb)**expzet)
          end do
        end do

        fnorm = (0.5_dp**l)*(pi**1.5_dp)*fnorm
        lx = gtos(isp)%orb_cphi_lx(io)
        ly = gtos(isp)%orb_cphi_ly(io)
        lz = gtos(isp)%orb_cphi_lz(io)
        prefac = dfac(2*lx - 1)*dfac(2*ly - 1)*dfac(2*lz - 1)
        gtos(isp)%norm_cphi(io) = 1.0_dp/dsqrt(prefac*fnorm)
      enddo

      !> Step 3.f: Compute cphi and sphi
      gtos(isp)%cphi(1:ncon_max,1:gtos(isp)%norbs_cphi) = 0.0_dp
      gtos(isp)%sphi(1:ncon_max,1:spp%norbs) = 0.0_dp

!!     For debugging
!      write(6,'(a,3i5)')                                       &
! &      'nao2gto_transfer: ncon_max, norbs_cphi, norbs = ',    &
! &       ncon_max, gtos(isp)%norbs_cphi, spp%norbs
!!     End debugging

!     Loop over all the Cartesian Gaussian functions of a particular species
      do io = 1, gtos(isp)%norbs_cphi

!       Identify the index of the shell to which this Cartesian Gaussian 
!       function is associated
        inlz = gtos(isp)%orb_index_cphi(io)

!       Identify the exponents of the cartesian coordinates x, y, and z
!       in the Cartessian Gaussian function
        lx = gtos(isp)%orb_cphi_lx(io)
        ly = gtos(isp)%orb_cphi_ly(io)
        lz = gtos(isp)%orb_cphi_lz(io)
!!       For debugging
!        write(6,'(a,5i5)')                                       &
! &        'nao2gto_transfer: io, inlz, lx, ly, lz = ',           &
! &        io, inlz, lx, ly, lz 
!!       End debugging

!       Loop over all the primitive Gaussians that expand this particular
!       Cartessian Gaussian function
        do jnlz = 1, gtos(isp)%orbnl_contract(inlz)
!!         For debugging
!          write(6,'(a,6i5)')                                                 &
! &          'nao2gto_transfer: io, inlz, orb_l_cphi, jnlz, co, counter = ',  &
! &          io, inlz, gtos(isp)%orb_l_cphi(io), jnlz, co(lx,ly,lz),          &
! &          (jnlz-1)*nco(gtos(isp)%orb_l_cphi(io))+co(lx,ly,lz)
!!         End debugging

!         Compute the coefficient of any primitive Gaussian that enters
!         in the expansion as the product of the coefficient times
!         the normalization function
!         The index of cphi runs between 1 and 
!         (number of Cartesian Gaussians for this shell) times 
!         (number of primitive Gaussians to expand the radial part of the
!         NAO for this shell).
!         For instance, for the l=2 shell, and if there are 
!         6 primitive Gaussians to expand the corresponding NAO, then
!         this index runs between 1 and 36
          gtos(isp)%cphi((jnlz-1)*nco(gtos(isp)%orb_l_cphi(io))+co(lx,ly,lz),io) =  &
&           gtos(isp)%orbnl_coefficient(jnlz,inlz) * gtos(isp)%norm_cphi(io)
        enddo
      enddo

!     Compute the transformation of the coefficients from Cartesian Gaussian
!     functions to Spherical Gaussian functions
      call cphi2sphi(ncon_max, gtos(isp)%norbs_cphi, spp%norbs, l_max, &
&                    spp%n_orbnl, spp%orbnl_l, nco, nso, &
&                     gtos(isp)%orbnl_index_cphi, gtos(isp)%orbnl_index_sphi, &
&                     gtos(isp)%cphi, gtos(isp)%sphi, orbtramat)

      inlz = 0
      do io = 1, spp%norbs
        l = spp%orb_l(io)
        m = spp%orb_m(io)
        if ( m .ne. -l ) cycle     ! not a normal orbital
        inlz = inlz + 1
!!       For debugging
!        do jnlz = io, io+nso(l)-1
!          write(6,'(a,7i5,36f12.5)')                                       &
! &        'nao2gto_transfer: isp, io, l, m, inlz, jnlz, index1, sphi = ',  &
! &                           isp, io, l, m, inlz, jnlz,                    &
! &                           nco(l)*gtos(isp)%orbnl_contract(inlz),        &
! &                           gtos(isp)%sphi(1:nco(l)*gtos(isp)%orbnl_contract(inlz),jnlz)
!        enddo
!!       End debugging
        gtos(isp)%orbnl_contraction_coeff(inlz) =  &
&         maxval([(sum(abs(gtos(isp)%sphi(1:nco(l)*gtos(isp)%orbnl_contract(inlz),jnlz))), &
&           jnlz=io,io+nso(l)-1)])
!!       For debugging
!        write(6,'(a,5f12.5)') &
! &        'nao2gto_transfer: sum(abs(gtos(isp)%sphi))', &
! &        (sum(abs(gtos(isp)%sphi(1:nco(l)*gtos(isp)%orbnl_contract(inlz),jnlz))), jnlz=io,io+nso(l)-1)
!        write(6,'(a,f12.5)') &
! &        'nao2gto_transfer: orbnl_contraction_coeff', &
! &        gtos(isp)%orbnl_contraction_coeff(inlz)
!!       End debugging
        if(.not.hfx_opts%is_fitted_nao) then
         expzet = 0.5_dp*REAL(2*l+3, dp)
         prefac = fac(2*l+2)*SQRT(pi)/2._dp**REAL(2*l+3, dp)/fac(l+1)
         fnorm = 0.0_dp
         do ipgf =1,gtos(isp)%orbnl_contract(inlz)
             gcca = gtos(isp)%orbnl_coefficient(ipgf,inlz)
             zeta = gtos(isp)%orbnl_zeta(ipgf,inlz)
             do jpgf=1, gtos(isp)%orbnl_contract(inlz)
                gccb = gtos(isp)%orbnl_coefficient(jpgf,inlz)
                zetb = gtos(isp)%orbnl_zeta(jpgf,inlz)
                fnorm = fnorm+gcca*gccb*prefac/(zeta+zetb)**expzet
             enddo
         enddo
         fnorm = SQRT(fnorm)
         !write(6,*) "Fiited normalization factor",isp inlz, fnorm
!         write(6,*) "fnorm1", io, l, nso(l), nso(l)-1
         gtos(isp)%sphi(:,io:io+nso(l)-1) =  gtos(isp)%sphi(:,io:io+nso(l)-1)*fnorm
        endif
      enddo

    enddo   ! is=1,nspecies

! ------------------------------------------------------------------------------
!> Step 4: Report about fitting parameters
! ------------------------------------------------------------------------------
!    call nao2gto_print_info(gtos, hfx_opts)

    if ( hfx_opts%dump_fit_data ) then

      call io_assign(fit_fd)
      open(unit=fit_fd, file="nao2gto_fit.yml", form="formatted", &
&       position="rewind", action="write", status="unknown")
      write(unit=fit_fd, fmt='(A/"---")') "%YAML 1.1"

      write(unit=fit_fd, fmt='(/A,":")') "hfx_options"
      write(unit=fit_fd, fmt='(2X,A,": ",I2)') "potential_type", &
&       hfx_opts%potential_type
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "omega", &
&       hfx_opts%omega
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "eps_farfield", &
&       hfx_opts%eps_farfield
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "eps_pairlist", &
&       hfx_opts%eps_pairlist
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "eps_schwarz", &
&       hfx_opts%eps_schwarz
      write(unit=fit_fd, fmt='(2X,A,": ",E20.8)') "eps_stored", &
&       hfx_opts%eps_stored
      if ( hfx_opts%DM_trunc ) then
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "DM_trunc", "True"
      else
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "DM_trunc", "False"
      end if
      if ( hfx_opts%dump_fit_data ) then
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "dump_fit_data", "True"
      else
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "dump_fit_data", "False"
      end if

      if ( hfx_opts%farfield ) then
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "farfield", "True"
      else
        write(unit=fit_fd, fmt='(2X,A,": ",A)') "farfield", "False"
      end if

      do isp=1,nspecies
        spp => species(isp)
        spp%label = species_label(isp)

        do inlz=1,spp%n_orbnl
          write(unit=fit_fd, fmt='(/"---"//A,":")') "orbital"
          write(unit=fit_fd, fmt='(2X,A,": ",A)') "species", &
&           trim(adjustl(spp%label))
          write(unit=fit_fd, fmt='(2X,A,": ",I4)') "qn_n", spp%orbnl_n(inlz)
          write(unit=fit_fd, fmt='(2X,A,": ",I4)') "qn_l", spp%orbnl_l(inlz)
          write(unit=fit_fd, fmt='(2X,A,": ",I4)') "zeta", spp%orbnl_z(inlz)
          if ( spp%orbnl_ispol(inlz) ) then
            write(unit=fit_fd, fmt='(2X,A,": ",A)') "polarized", "True"
          else
            write(unit=fit_fd, fmt='(2X,A,": ",A)') "polarized", "False"
          end if
          write(unit=fit_fd, fmt='(2X,A,": ",E15.5)') &
&           "population", spp%orbnl_pop(inlz)

          if ( do_fit ) then
            write(unit=fit_fd, fmt='(2X,A,": ",A)') "do_fit", "True"
          else
            write(unit=fit_fd, fmt='(2X,A,": ",A)') "do_fit", "False"
          end if

          write(unit=fit_fd, fmt='(2X,A,":")') "raw_data"
          do irad=1,hfx_opts%npts_fit
            fit_rad = spp%orbnl(inlz)%cutoff * real(irad-1, dp) / &
&               real(hfx_opts%npts_fit, dp)
            call rad_get(spp%orbnl(inlz), fit_rad, fit_orb, dummy)
            write(unit=fit_fd, fmt='(4X,"- [",E20.8,", ",E20.8,"]")') &
&               fit_rad, fit_orb
          enddo

          write(unit=fit_fd, fmt='(2X,A,":")') "trial"
          do jnlz=1,hfx_contract(inlz,isp)
            write(unit=fit_fd, fmt='(4X,"- [",E24.8,", ",E24.8,"]")') &
&             gtos(isp)%orbnl_zeta(jnlz,inlz), &
              gtos(isp)%orbnl_coefficient(jnlz,inlz)
          end do
        end do
      end do

      write(unit=fit_fd, fmt='(/A)') "..."
      close(unit=fit_fd)

    end if   ! hfx_opts%dump_fit_data

! ------------------------------------------------------------------------------
!> Step 5: Check the completeness of the fitting
! ------------------------------------------------------------------------------
    do_stop = .false.
    do isp = 1, nspecies
      spp => species(isp)
      spp%label = species_label(isp)
      do inlz = 1, spp%n_orbnl
        if ( hfx_contract(inlz,isp) == 0 ) then
          write(msg, fmt='(A," (",A,", n=",I2,", l=",I2,", zeta=",I2,")")') &
&           "Gaussian fitting of the orbitals failed for", &
&           trim(spp%label), spp%orbnl_n(inlz), spp%orbnl_l(inlz), spp%orbnl_z(inlz)
          write(*, fmt='(A)') msg
          do_stop = .true.
        end if
      end do
    end do


! ------------------------------------------------------------------------------
!> Step 6: Replace original orbitals with fitted NAOs from GTOs
! Add  by  Xinming, 20220328
!-------------------------------------------------------------------------------

    call nao2gto_print_info(gtos, hfx_opts)

    if(hfx_opts%is_fitted_nao) then
     call nao2gto_replace_nao(gtos, hfx_opts)
    endif

!> Step 7: Properly free memory
! ------------------------------------------------------------------------------
    do l=0,l_max
      write(msg,'("orbtramat(",I1,")%c2s")') l
      call de_alloc(orbtramat(l)%c2s, name=trim(msg), &
&       routine='nao2gto_transfer')
      nullify(orbtramat(l)%c2s)
    enddo
    deallocate(orbtramat)
    nullify(orbtramat)

    call de_alloc(hfx_contract, name="hfx_contract", &
&     routine="nao2gto_transfer")
    call de_alloc(nao2gto_zeta, name="nao2gto_zeta", &
&     routine="nao2gto_transfer")
    call de_alloc(nao2gto_coefficient, name="nao2gto_coefficient", &
&     routine="nao2gto_transfer")

! ------------------------------------------------------------------------------
!> Step 7: Stop program if the fitting is incomplete
! ------------------------------------------------------------------------------
    if ( do_stop ) then
      write(msg, fmt='(A)') &
&       "The fitting procedure was not successfully completed"
      call die(trim(msg))
    end if

    call nao2gto_libint_dump(hfx_libint)
    call nao2gto_libderiv_dump(hfx_libderiv)

    write(*,'(A,/)') "************************ End: HYBRID XC INITIALIZATION ************************"

  end subroutine nao2gto_transfer

  ! ***************************************************************************
  ! *** Private routines                                                    ***
  ! ***************************************************************************

  ! ***************************************************************************
  !> \brief Displays the values of the specified Hartree-Fock exchange data
  !!        structure
  !!
  !! \param[in] gtos: data structure storing Gaussian fitting information
  !! \param[in] hfx_opts: data structure storing Hartree-Fock exchange options
  ! ***************************************************************************
  subroutine nao2gto_print_info(gtos, hfx_opts)
    use precision, only: dp
    use units, only: pi
    use nao2gto_utils, only: fac
    use atm_types, only: nspecies, species, species_info
    use chemical,  only: species_label
    use parallel, only: IOnode
#ifdef BSC_CELLXC
    use bsc_xcmod, only: nXCfunc, XCfunc, XCauth
#else
    use siestaXC, only: getXC
#endif /* BSC_CELLXC */
    use nao2gto_types, only: gto_info_type, hfx_options_type
    use nao2gto_types, only: do_hfx_potential_coulomb
    use nao2gto_types, only: do_hfx_potential_short
    use nao2gto_types, only: do_hfx_potential_truncated
    use radial,    only: rad_get
    use alloc,     only: re_alloc, de_alloc

    implicit none

    ! Arguments
    type(gto_info_type), intent(in) :: gtos(nspecies)
    type(hfx_options_type), intent(in) :: hfx_opts

    ! Local variables
    integer :: inlz, isp, l, jnlz, nf, lun2, n_series, zet, ispol, npts, ip, irad
    integer :: Nwrite
    real(dp) :: dummy
    real(dp), allocatable :: fval(:)           ! Value of the linear combination
                                               !   of Gaussians in the 
                                               !   linear mesh
    integer :: ipgf, jpgf
!    real(dp) :: dummy
    real(dp) :: fnorm, gcca, gccb, zeta, zetb, expzet, prefac

    real(dp), dimension(:), pointer   :: rad_pts => null()
                                               ! Values of the points in the
                                               !   linear mesh where the
                                               !   radial part of the NAO
                                               !   is tabulated
    real(dp), dimension(:), pointer   :: orb_pts => null()
                                               ! Radial part of the NAO
                                               !   (divided by r^l)
                                               !   in the linear mesh
    character*30 fileid
    type(species_info), pointer :: spp => null()
#ifndef BSC_CELLXC
    character(len=20) :: XCfunc(10), XCauth(10)
    integer :: nXCfunc
#endif /* BSC_CELLXC */

    ! -------------------------------------------------------------------------

    if ( .not. IOnode ) return

    npts = hfx_opts%npts_fit           ! Number of points in the mesh

    write(*,'(/,A,/)') "nao2gto_print_info: NAO2GTO fitting information -------------------------------"

    ! Mark the beginning of the FDF data for automatic extraction
    write(*,'("# %%%",1X,A,/)') "HYBRID XC FDF BEGIN"

    ! Display XC parameters
#ifndef BSC_CELLXC
    call getXC(nXCfunc, XCfunc, XCauth)
#endif /* BSC_CELLXC */
    do nf=1,nXCfunc
      select case(hfx_opts%potential_type)
        case(do_hfx_potential_coulomb)
          write(*, &
            '("#",1X,"XC parameters with a ",A," EXX potential")') &
            "full Coulomb"
        case(do_hfx_potential_short)
          write(*, &
            '("#",1X,"XC parameters with a ",A," EXX potential")') &
            "short-range"
        case(do_hfx_potential_truncated)
          write(*, &
            '("#",1X,"XC parameters with a ",A," EXX potential")') &
            "truncated Coulomb"
        case default
          write(*, &
            '("#",1X,"XC parameters with a ",A," EXX potential")') &
            "*UNKNOWN*"
      end select
      write(*,'("XC.functional",3X,A)') trim(XCfunc(nf))
      write(*,'("XC.authors",6X,A)') trim(XCauth(nf))
    enddo

    ! Display the NAO2GTO block
    write(*,'(/,"#",1X,A,/,"#")') &
      "The NAO2GTO block is structured as following:"
    write(*,'("#",5X,A)') "species1 norbs_species1"
    write(*,'("#",5X,A)') "n1 l1 zeta1 ngaussians_orbital1"
    write(*,'("#",5X,A)') "exponent_gaussian1 coefficient_gaussian1"
    write(*,'("#",5X,A)') "exponent_gaussian2 coefficient_gaussian2"
    write(*,'("#",5X,A)') "exponent_gaussian3 coefficient_gaussian3"
    write(*,'("#",5X,A)') "..."
    write(*,'("#",5X,A)') "n2 l2 zeta2 ngaussians_orbital2"
    write(*,'("#",5X,A)') "exponent_gaussian1 coefficient_gaussian1"
    write(*,'("#",5X,A)') "..."
    write(*,'("#",5X,A)') "species2 norbs_species2"
    write(*,'("#",5X,A,/,"#")') "..."
    write(*,'(A)') "%block NAO2GTO"
    do isp=1,nspecies
      spp => species(isp)
      spp%label = species_label(isp)

      write(*,'(A,1X,I3)') trim(spp%label), spp%n_orbnl
      do inlz=1,spp%n_orbnl
        write(*,'(I1,3(1X,I2))') spp%orbnl_n(inlz), spp%orbnl_l(inlz), &
          spp%orbnl_z(inlz), gtos(isp)%orbnl_contract(inlz)
        do jnlz=1,gtos(isp)%orbnl_contract(inlz)
          write(*,'(E18.8,1X,E18.8)') gtos(isp)%orbnl_zeta(jnlz,inlz), &
            gtos(isp)%orbnl_coefficient(jnlz,inlz)
        enddo
      enddo
    enddo
    write(*,'(A)') "%endblock NAO2GTO"

    ! Display Hartree-Fock options
    write(*,'(/,"#",1X,A)') "Hartree-Fock exchange options"
    write(*,'(A,9X,L12)') "HFX.UseFittedNAOs", hfx_opts%is_fitted_nao
    write(*,'(A,9X,L12)') "HFX.TruncateDM", hfx_opts%DM_trunc
    write(*,'(A,8X,L12)') "HFX.DumpFitData", hfx_opts%dump_fit_data
    write(*,'(A,11X,L12)') "HFX.FarField", hfx_opts%farfield
    write(*,'(A,6X,I12)') "HFX.FitDataPoints", hfx_opts%npts_fit
    write(*,'(A,2X,E12.3)') "HFX.FarFieldTolerance", hfx_opts%eps_farfield
    write(*,'(A,2X,E12.3)') "HFX.PairListTolerance", hfx_opts%eps_pairlist
    write(*,'(A,3X,E12.3)') "HFX.SchwarzTolerance", hfx_opts%eps_schwarz
    write(*,'(A,1X,E12.3)') "HFX.StoreERIsTolerance", hfx_opts%eps_stored
    write(*,'(A,14X,E12.3)') "HFX.Omega", hfx_opts%omega
    write(*,'(A,14X,I5)') "HFX.MinimumNumberGaussians", hfx_opts%min_num_gaus
    write(*,'(A,14X,I5)') "HFX.MaximumNumberGaussians", hfx_opts%max_num_gaus
    write(*,'(A,14X,E12.3)') "HFX.SeparationExponents", hfx_opts%threshold_exp_gaus
    write(*,'(A,14X,E12.3)') "HFX.ToleranceFit", hfx_opts%tolerance_gaus
    write(*,'(A,14X,E12.3)') "HFX.GaussianEPS", hfx_opts%gto_eps
    write(*,'(A,9X,L12)') "HFX.Dynamic_parallel", hfx_opts%Dynamic_parallel
    write(*,'(A,8X,L12)') "HFX.FragSize", hfx_opts%frag_size


    ! Mark the end of the FDF data for automatic extraction
    write(*,'(/,"# %%%",1X,A)') "HYBRID XC FDF END"

    write(*, '(/,a,/)') "nao2gto_print_info: END -------------------------------------------------------"
 
    Nwrite = hfx_opts%npts_fit
!   Dump the linear combinations of the Gaussians at the linear point mesh
    call re_alloc(rad_pts, 1, Nwrite, 'rad_pts', &
                  'nao2gto_print_info')
    call re_alloc(orb_pts, 1, Nwrite, 'orb_pts', &
                  'nao2gto_print_info')
    if ( allocated(fval) ) deallocate(fval)
    allocate( fval(Nwrite) )
    orb_pts = 0.0_dp

    do isp=1,nspecies
      spp => species(isp)
      spp%label = species_label(isp)
      write(*,'(A,1X,I3)') trim(spp%label), spp%n_orbnl
      n_series = 0
      do inlz=1,spp%n_orbnl

        do irad=1,hfx_opts%npts_fit
          rad_pts(irad) = spp%orbnl(inlz)%cutoff * real(irad-1, dp) / &
&                 (real(hfx_opts%npts_fit, dp) - 1)
          call rad_get(spp%orbnl(inlz), rad_pts(irad), orb_pts(irad), dummy)
        enddo
        orb_pts(hfx_opts%npts_fit) = 0.0_dp

!        do irad=hfx_opts%npts_fit+1, Nwrite
!            rad_pts(irad) = spp%orbnl(inlz)%cutoff * real(irad-1, dp) / &
!&                 (real(hfx_opts%npts_fit, dp) - 1)
!        enddo

        if (spp%orbnl_ispol(inlz)) ispol = 1
          zet = spp%orbnl_z(inlz)
          if (zet == 1) n_series = n_series + 1
          write(fileid,'(a,i1,a,i1,a,a)') "ORB.GAUSS.S", n_series, ".", &
 &             zet, ".", trim(spp%label)
          call io_assign(lun2)
          open(unit=lun2,file=fileid,status='replace',                  &
 &             form='formatted')
!          write(lun2,'("# ",5X,A)')                                     &
! &          " species label, n, l, zeta, ngaussians_orbital"
 !         write(lun2,'("# ",5X,A,2x, I1,3(1X,I2))')                     &
 !&            trim(spp%label), spp%orbnl_n(inlz), spp%orbnl_l(inlz),    &
 !             spp%orbnl_z(inlz), gtos(isp)%orbnl_contract(inlz)
 !         write(lun2,'("#",5X,A)') "exponent_gaussian  coefficient_gaussian"
 !         do jnlz=1,gtos(isp)%orbnl_contract(inlz)
 !           write(lun2,'("#",E18.8,1X,E18.8)')gtos(isp)%orbnl_zeta(jnlz,inlz), &
 !             gtos(isp)%orbnl_coefficient(jnlz,inlz)
 !         enddo

        l = spp%orbnl_l(inlz)
         expzet = 0.5_dp*REAL(2*l+3, dp)
         prefac = fac(2*l+2)*SQRT(pi)/2._dp**REAL(2*l+3, dp)/fac(l+1)
        fnorm = 0.0_dp
        do ipgf =1,gtos(isp)%orbnl_contract(inlz)
             gcca = gtos(isp)%orbnl_coefficient(ipgf,inlz)
             zeta = gtos(isp)%orbnl_zeta(ipgf,inlz)
             do jpgf=1, gtos(isp)%orbnl_contract(inlz)
                gccb = gtos(isp)%orbnl_coefficient(jpgf,inlz)
                zetb = gtos(isp)%orbnl_zeta(jpgf,inlz)
                fnorm = fnorm+gcca*gccb*prefac/(zeta+zetb)**expzet
             enddo
        enddo
        fnorm = 1.0_dp/SQRT(fnorm)

 !         write(lun2,'("#",5X,A)')     &
 !&           "r, NAO/r^l, cGTOs/r^l, normalized-cGTO/r^l"

          fval(:) = 0.0_dp

          do jnlz = 1, gtos(isp)%orbnl_contract(inlz)
            fval(:) = fval(:) + gtos(isp)%orbnl_coefficient(jnlz,inlz) *  &
 &                    exp(-rad_pts(:)**2*gtos(isp)%orbnl_zeta(jnlz,inlz))
          end do
          do ip = 1, Nwrite
            write(lun2,'(4g26.16)')rad_pts(ip),orb_pts(ip), fval(ip), fval(ip)*fnorm
          enddo 
         call io_close(lun2)
      enddo 
    enddo 

    call de_alloc(rad_pts, 'rad_pts', 'nao2gto_print_info')
    call de_alloc(orb_pts, 'orb_pts', 'nao2gto_print_info')
    deallocate( fval )

  end subroutine  nao2gto_print_info

! *****************************************************************************
!> \brief Reads Hartree-Fock exchange parameters from the FDF input
!!
!! This routine reads the values of the Hartree-Fock exchange parameters
!! from the FDF input of SIESTA and strores them into the specified data
!! structure. It also sets the default values for parameters not present
!! in the input file.
!!
!! \author Xinming Qin
!! \author Yann Pouillon
!!
!! \par History
!!      - 12.2017 Adapted to SIESTA 4 [Yann Pouillon]
!!      - 01.2018 Refactored and documented [Yann Pouillon]
!!
!! \param[out] hfx_opts: data structure storing Hartree-Fock exchange
!!                            parameters
! *****************************************************************************
  subroutine nao2gto_read_options(hfx_opts)

    use fdf, only: fdf_get
#ifdef BSC_CELLXC
    use bsc_xcmod, only: nXCfunc, XCfunc, XCauth
#else
    use siestaXC, only: getXC
#endif /* BSC_CELLXC */
    use atmparams, only: NTBMAX
    use nao2gto_types, only: hfx_options_type
    use nao2gto_types, only: do_hfx_potential_coulomb
    use nao2gto_types, only: do_hfx_potential_short
    use nao2gto_types, only: do_hfx_potential_truncated

    implicit none

    ! Arguments
    type(hfx_options_type), intent(out) :: hfx_opts

    ! Local variables
    integer :: nf
#ifndef BSC_CELLXC
    integer :: nXCfunc
    character(len=20) :: XCfunc(10), XCauth(10)
#endif /* BSC_CELLXC */

    ! -------------------------------------------------------------------------

#ifndef BSC_CELLXC
    call getXC(nXCfunc, XCfunc, XCauth)
#endif /* BSC_CELLXC */

    hfx_opts%DM_trunc = fdf_get("HFX.TruncateDM", .true.)
    hfx_opts%dump_fit_data = fdf_get("HFX.DumpFitData", .true.)
    hfx_opts%farfield = fdf_get("HFX.FarField", .true.)
    hfx_opts%npts_fit = fdf_get("HFX.FitDataPoints", NTBMAX)
    hfx_opts%min_num_gaus = fdf_get("HFX.MinimumNumberGaussians", 3)
    hfx_opts%max_num_gaus = fdf_get("HFX.MaximumNumberGaussians", 6)
    hfx_opts%threshold_exp_gaus = fdf_get("HFX.SeparationExponents", 1.4)
    hfx_opts%tolerance_gaus = fdf_get("HFX.ToleranceFit", 1.d-3)
    hfx_opts%is_fitted_nao = fdf_get("HFX.UseFittedNAOs", .true.)
    hfx_opts%gto_eps = fdf_get( "HFX.GaussianEPS", 1.0d-5)
    hfx_opts%Dynamic_parallel = fdf_get("HFX.Dynamic_parallel", .false.)
    hfx_opts%frag_size = fdf_get("HFX.FragSize", 10000)

    do nf = 1,nXCfunc
      if ( ((XCauth(nf).eq.'pbe0') .or. (XCauth(nf).eq.'PBE0')) &
&          .and. (XCfunc(nf).eq.'GGA') ) then
        hfx_opts%potential_type = do_hfx_potential_coulomb
      else
        if ( ((XCauth(nf).eq.'hse06') .or. (XCauth(nf).eq.'HSE06')) &
&            .and. (XCfunc(nf).eq.'GGA') ) then
          hfx_opts%potential_type = do_hfx_potential_short
        else
          hfx_opts%potential_type = do_hfx_potential_truncated
        endif
      endif
    enddo

    hfx_opts%omega        = fdf_get('HFX.Omega',              0.11d0)
    hfx_opts%eps_farfield = fdf_get('HFX.FarFieldTolerance',  1.0d-6)
    hfx_opts%eps_pairlist = fdf_get('HFX.PairListTolerance',  1.0d-6)
    hfx_opts%eps_schwarz  = fdf_get('HFX.SchwarzTolerance',   1.0d-6)
    hfx_opts%eps_stored   = fdf_get('HFX.StoreERIsTolerance', 1.0d-6)

  end subroutine nao2gto_read_options

  ! ***************************************************************************
  !> \brief Finds the orbital index of a (n, l, zeta) triplet
  ! ***************************************************************************
  function orb_nlz_to_index(spp, orb_n, orb_l, orb_z) result(orb_index)

    use atm_types, only: species_info

    implicit none

    ! Arguments
    type(species_info), intent(in) :: spp
    integer, intent(in) :: orb_n
    integer, intent(in) :: orb_l
    integer, intent(in) :: orb_z

    ! Local variables
    integer :: orb_index
    integer :: iorb

    orb_index = -1

    do iorb=1,spp%n_orbnl
      if ( spp%orbnl_n(iorb) == orb_n ) then
        if ( spp%orbnl_l(iorb) == orb_l ) then
          if ( spp%orbnl_z(iorb) == orb_z ) then
            orb_index = iorb
            exit
          end if
        end if
      end if
    end do

  end function orb_nlz_to_index

  ! ***************************************************************************
  !> \brief Replace original orbitals with fitted NAOs, Xinming 20220404
  ! ***************************************************************************

  subroutine nao2gto_replace_nao(gtos, hfx_opts)
    use precision, only: dp
    use units, only: pi
    use atm_types, only: nspecies, species, species_info
!    use chemical,  only: species_label
!    use parallel, only: IOnode
    use nao2gto_types, only: gto_info_type, hfx_options_type
    use radial,     only: rad_alloc, rad_get, rad_setup_d2
    use radial,  only: reset_rad_func
    use radial,     only: rad_func
    use alloc,     only: re_alloc, de_alloc
    use nao2gto_utils, only: fac 

    implicit none

    ! Arguments
    type(gto_info_type), intent(in) :: gtos(nspecies)
    type(hfx_options_type), intent(in) :: hfx_opts

    ! Local variables
    integer :: inlz, isp, l, jnlz, nf, lun2, n_series, ispol, npts, ip, irad
    integer :: ipgf, jpgf
!    real(dp) :: dummy
!    real(dp), allocatable :: fval(:)           ! Value of the linear combination
                                               !   of Gaussians in the
                                               !   linear mesh
    real(dp) :: fnorm, gcca, gccb, zeta, zetb, expzet, prefac
    real(dp), dimension(:), pointer   :: rad_pts => null()
                                               ! Values of the points in the
                                               !   linear mesh where the
                                               !   radial part of the NAO
                                               !   is tabulated
    real(dp), dimension(:), pointer   :: orb_pts => null()
                                               ! Radial part of the NAO
                                               !   (divided by r^l)
                                               !   in the linear mesh
    type(species_info), pointer :: spp => null()
    type(rad_func), pointer      :: func => null()


    call re_alloc(rad_pts, 1, hfx_opts%npts_fit, 'rad_pts', &
                  'nao2gto_replace_nao')

    do isp=1,nspecies
      spp => species(isp)
      n_series = 0
      do inlz = 1,spp%n_orbnl

        l = spp%orbnl_l(inlz)
         expzet = 0.5_dp*REAL(2*l+3, dp)
         prefac = fac(2*l+2)*SQRT(pi)/2._dp**REAL(2*l+3, dp)/fac(l+1)
        fnorm = 0.0_dp
        do ipgf =1,gtos(isp)%orbnl_contract(inlz)
             gcca = gtos(isp)%orbnl_coefficient(ipgf,inlz)
             zeta = gtos(isp)%orbnl_zeta(ipgf,inlz)
             do jpgf=1, gtos(isp)%orbnl_contract(inlz)
                gccb = gtos(isp)%orbnl_coefficient(jpgf,inlz)
                zetb = gtos(isp)%orbnl_zeta(jpgf,inlz)
                fnorm = fnorm+gcca*gccb*prefac/(zeta+zetb)**expzet
             enddo
        enddo
        fnorm = 1.0_dp/SQRT(fnorm)

        func => spp%orbnl(inlz)
!         call reset_rad_func(func)
!         spp%orbnl(inlz)%cutoff = gtos(isp)%shell_radius(inlz)
!         spp%orbnl(inlz)%delta = spp%orbnl(inlz)%cutoff/(dble(hfx_opts%npts_fit-1)+1.0d-20)
         func%cutoff = gtos(isp)%shell_radius(inlz)
         func%n = hfx_opts%npts_fit
         func%delta =func%cutoff/(dble(func%n-1)+1.0d-20)
!        write(6,*) "cutoff 1 :: ", func%cutoff, func%n, func%delta
!        write(6,*) "cutoff 1 :: ", gtos(isp)%shell_radius(inlz),spp%orbnl(inlz)%cutoff
!        write(6,*) "n, delta 1 :: ",  spp%orbnl(inlz)%n, spp%orbnl(inlz)%delta

!        orb_pts(hfx_opts%npts_fit) = 0.0_dp

        do irad=1,func%n
           rad_pts(irad) = func%cutoff * real(irad-1, dp) / &
&                 (real(func%n, dp) - 1)
!          call rad_get(spp%orbnl(inlz), rad_pts(irad), orb_pts(irad), dummy)

        enddo

        func%f(:) = 0.0_dp

        do jnlz=1,gtos(isp)%orbnl_contract(inlz)
!           write(6,'("#",E18.8,1X,E18.8)')gtos(isp)%orbnl_zeta(jnlz,inlz), &
!                   gtos(isp)%orbnl_coefficient(jnlz,inlz)
           func%f(:) = func%f(:) + gtos(isp)%orbnl_coefficient(jnlz,inlz) *  &
 &                    exp(-rad_pts(:)**2*gtos(isp)%orbnl_zeta(jnlz,inlz))
        enddo
       
        func%f(func%n) = 0.0_dp       
        func%f  =  func%f * fnorm
        call rad_setup_d2(func,yp1=0.0_dp,ypn=huge(1.0_dp))

       enddo

     enddo

    call de_alloc(rad_pts, 'rad_pts', 'nao2gto_replace_nao')

end subroutine nao2gto_replace_nao


end module nao2gto_io
