! 
! This file is part of the HONPAS package.
!
! Coded by xmqin, 12. 03, 2016.
!
!
      module hfx_types
!
!=======================================================================
!
!     This module defines data structures to handle the specification
!     of the NAO2GTO set and prescreening tolerances of ERIs in HONPAS. 
!     This specification is read from the fdf file by routine 'nao2gto_read'.
!
      use kinds,     only : dp, int_8
!      use alloc,         only : re_alloc, de_alloc
      use m_fdf_global,  only :  fdf_global_get
      use parallel,      only :  Node, Nodes
      use xcmod,         only : nXCfunc, XCfunc, XCauth

      IMPLICIT NONE
      !          hfx_pgf_image, hfx_pgf_list, &
       !         hfx_pgf_product_list, hfx_cell_type

      PUBLIC :: read_hfx_info, print_hfx_info

!      INTEGER,DIMENSION(:,:), POINTER, PUBLIC, SAVE :: D2SIndx
      integer, pointer, save, public   ::  D2Sindx(:,:) 
      real(dp), pointer, save, public  ::  eri_prescreen(:) 

      INTEGER(int_8), PARAMETER, PRIVATE  :: one = 1_int_8
      REAL(dp), PARAMETER, PUBLIC         :: log_zero = -1000.0_dp
      REAL(dp), PARAMETER, PUBLIC   :: powell_min_log = -20.0_dp
      REAL(KIND=dp), DIMENSION(0:10), &
      PARAMETER, PUBLIC                      :: mul_fact = (/1.0_dp,&
                                                             1.1781_dp,&
                                                             1.3333_dp,&
                                                             1.4726_dp,&
                                                             1.6000_dp,&
                                                             1.7181_dp,&
                                                             1.8286_dp,&
                                                             1.9328_dp,&
                                                             2.0317_dp,&
                                                             2.1261_dp,&
                                                             2.2165_dp/)


      INTEGER, PARAMETER, PUBLIC    ::    do_hfx_potential_coulomb=1,&
                                          do_hfx_potential_short=2,&
                                          do_hfx_potential_truncated=3

      INTEGER, SAVE                 :: init_t_c_g0_lmax = -1
      integer, pointer, save, public :: subshell(:)

      type, public :: hfx_input_parameter
        integer   ::  potential_type = -1 != 1 !! 1/r/erfc(wr)/r ...
        logical   ::  far_field = .true.    ! Far-near field screening ?
        logical   ::  DM_trunc = .false.
        logical   ::  eri_on_disk = .false.
        logical   ::  on_the_fly  = .false.
        logical   ::  is_fitted_nao =.true.
        logical   ::  dump_fit_data = .false.
        logical   ::  parallel = .false.
        integer   ::  fragsize  = 10000
        integer   ::  npts_fit = -1
        integer   ::  min_num_gaus   =  3
        integer   ::  max_num_gaus   =  6
        real(dp)  ::  omega = -1.0_dp
        real(dp)  ::  cutoff_radius = -1.0_dp
        real(dp)  ::  eps_far = -1.0_dp
        real(dp)  ::  eps_pairlist = -1.0_dp
        real(dp)  ::  eps_schwarz = -1.0_dp
        real(dp)  ::  eps_stored = -1.0_dp
        real(dp)  ::  threshold_exp_gaus = 1.4_dp
        real(dp)  ::  tolerance_gaus = 1.d-3
        real(dp)  ::  gto_eps = 1.d-5
      end type hfx_input_parameter

!      type(hfx_input_parameter), save, public :: hfx_parameter

      TYPE cell_type
       INTEGER                           :: ref_count, id_nr
       LOGICAL                           :: orthorhombic
       REAL(KIND = dp)                   :: deth
       INTEGER, DIMENSION(3)             :: perd
       REAL(KIND = dp), DIMENSION(3,3)   :: hmat,h_inv
      END TYPE cell_type

! *****************************************************************************
      TYPE cell_p_type
        TYPE(cell_type),POINTER :: cell
      END TYPE cell_p_type


      TYPE hfx_cell_type
        INTEGER                                  :: nsc(3)
        REAL(dp)                                 :: cell(3)
        REAL(dp)                                 :: cell_r(3)
      END TYPE

!----------------------  pair list type  -----------------------------------
      TYPE, PUBLIC :: hfx_screen_coeff_type
        REAL(dp)                                 :: x(2)
      END TYPE
! --------------------------------------------------------------------------

!----------------------pair list type -----------------------------
      TYPE, public :: pair_list_element_type
        INTEGER, DIMENSION(2) :: pair
        INTEGER               :: nl_index
        REAL(dp)              :: r1(3),r2(3)
        REAL(KIND=dp)         :: vec(3)
        REAL(dp)              :: dist2
      END TYPE

      TYPE, public :: pair_list_type
        TYPE(pair_list_element_type), DIMENSION(:), ALLOCATABLE :: element
        INTEGER :: nelement
      END TYPE pair_list_type

      TYPE(pair_list_type), save :: list_ij, list_kl

!-----------------------------------------------------------------

      type, public :: ERI_link
        real(dp) :: gto_eri
        integer, dimension(4) :: four_index
        type(ERI_link), pointer :: next
      end type
      type(ERI_link), pointer, public :: head
      type(ERI_link), pointer, public :: tail
      type(ERI_link), pointer, public :: ptr

      
!----------------------------------------------------------------
      type(hfx_screen_coeff_type), &
      dimension(:, :, :, :, :, :), pointer, save  :: &
                                   pair_dist_radii_pgf, &
                                   sfc_pgf
      type(hfx_screen_coeff_type), &
      dimension(:, :, :, :), pointer, save  :: &
                                   sfc_shell
      type(hfx_screen_coeff_type), &
      dimension(:, :), pointer, save  :: &
                                   sfc_kind
!---------------------------------------------------------------
      PRIVATE
      PUBLIC :: list_ij, list_kl, pair_dist_radii_pgf, &
                sfc_pgf, sfc_shell, sfc_kind

! ---------------------------------------------------------------------


     CONTAINS

      subroutine  read_hfx_info(hfx_parameter)
        type(hfx_input_parameter) :: hfx_parameter
        integer nf
     
      do nf = 1,nXCfunc
       if((XCauth(nf).eq.'pbe0').or.(XCauth(nf).eq.'PBE0') &  
              .and.(XCfunc(nf).eq.'GGA')) then
           hfx_parameter%potential_type = 1 !do_hfx_potential_coulomb
       elseif((XCauth(nf).eq.'hse06').or.(XCauth(nf).eq.'HSE06') &
             .and.(XCfunc(nf).eq.'GGA')) then
           hfx_parameter%potential_type = 2 !do_hfx_potential_short
           call fdf_global_get(hfx_parameter%omega, "Omega", 0.11d0)
       endif
      enddo

     call fdf_global_get(hfx_parameter%on_the_fly,"On_the_fly", .false.)
     call fdf_global_get(hfx_parameter%eps_schwarz,'HFX.SchwarzTolerance',1.0d-6)
     call fdf_global_get(hfx_parameter%eps_stored,'HFX.StoreERIsTolerance',1.0d-6)
     call fdf_global_get(hfx_parameter%eps_pairlist,'HFX.PairListTolerance',1.0d-6)
     call fdf_global_get(hfx_parameter%eps_far,'HFX.FarFieldTolerance',1.0d-6)
     call fdf_global_get(hfx_parameter%far_field,"HFX.FarField", .true.)
     call fdf_global_get(hfx_parameter%eri_on_disk,'HFX.StoreERIs', .false.)
     call fdf_global_get(hfx_parameter%DM_trunc,'HFX.TruncateDM', .true.)
     call fdf_global_get(hfx_parameter%parallel,'HFX.Dynamic_parallel', .true.)
     call fdf_global_get(hfx_parameter%fragsize,'HFX.FragSize', 10000)
     call fdf_global_get(hfx_parameter%dump_fit_data, 'HFX.DumpFitData',.true.)
     call fdf_global_get(hfx_parameter%npts_fit, 'HFX.FitDataPoints',500)
     call fdf_global_get(hfx_parameter%min_num_gaus, &
                         'HFX.MinimumNumberGaussians', 3)
     call fdf_global_get(hfx_parameter%max_num_gaus, &
                         'HFX.MaximumNumberGaussians', 6)
     call fdf_global_get(hfx_parameter%threshold_exp_gaus, &
                         'HFX.SeparationExponents', 1.4_dp)
     call fdf_global_get(hfx_parameter%tolerance_gaus, &
                         'HFX.ToleranceFit', 1.d-3)
     call fdf_global_get(hfx_parameter%is_fitted_nao, &
                         'HFX.UseFittedNAOs', .true.)
     call fdf_global_get(hfx_parameter%gto_eps, "HFX.GaussianEPS", 1.0d-5)

      end subroutine  read_hfx_info

      subroutine  print_hfx_info(hfx_parameter)
       type(hfx_input_parameter) :: hfx_parameter
       integer  nf
      if(node.eq.0) then
       write(6,'(/, a)') &
       " Input, default information for hybrid DFT calculations"   
      do nf = 1,nXCfunc
       if((XCauth(nf).eq.'pbe0').or.(XCauth(nf).eq.'PBE0') &
              .and.(XCfunc(nf).eq.'GGA')) then
        write(6,'(/, a)') " Hybrid DFT :: PBE0"
           hfx_parameter%potential_type = do_hfx_potential_coulomb
       elseif((XCauth(nf).eq.'hse06').or.(XCauth(nf).eq.'HSE06') &
             .and.(XCfunc(nf).eq.'GGA')) then
        write(6,'(/, a)') " Hybrid DFT :: HSE06"
        write(6, '(a, F8.3)') ' Omega = ', hfx_parameter%omega 
       endif

      enddo

        write(6,'(a, D6.1)') 'Schwarz tolerance = ', hfx_parameter%eps_schwarz 
        write(6,'(a, D6.1)') 'pairlist tolerance = ', hfx_parameter%eps_pairlist
        write(6,'(a, D6.1)') 'Stored tolerance = ',  hfx_parameter%eps_stored
        write(6,'(a, D6.1)') 'Farfield tolerance =', hfx_parameter%eps_far
        write(6,'(a, L3)') 'On_the_fly = ', hfx_parameter%on_the_fly
        write(6,'(a, L3)') 'use_DM_trunc = ', hfx_parameter%DM_trunc
        write(6,'(a, L3)') 'Farfield = ', hfx_parameter%far_field
        write(6,'(a, L3)') 'Dynamic parallel = ', hfx_parameter%parallel
        write(6,'(a, I8)') 'Dynamic tasksize = ', hfx_parameter%fragsize
        write(*,'(A,9X,L12)') "HFX.UseFittedNAOs", hfx_parameter%is_fitted_nao
!        write(*,'(A,8X,L12)') "HFX.DumpFitData", hfx_parameter%dump_fit_data
        write(*,'(A,6X,I12)') "HFX.FitDataPoints", hfx_parameter%npts_fit
        write(*,'(A,14X,I5)') "HFX.MinimumNumberGaussians", hfx_parameter%min_num_gaus
        write(*,'(A,14X,I5)') "HFX.MaximumNumberGaussians", hfx_parameter%max_num_gaus
        write(*,'(A,14X,E12.3)') "HFX.SeparationExponents", hfx_parameter%threshold_exp_gaus
        write(*,'(A,14X,E12.3)') "HFX.ToleranceFit", hfx_parameter%tolerance_gaus
        write(*,'(A,14X,E12.3)') "HFX.GaussianEPS", hfx_parameter%gto_eps

      endif

      end subroutine  print_hfx_info

      end module hfx_types
