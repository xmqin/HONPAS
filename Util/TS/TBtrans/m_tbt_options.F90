! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code segment has been fully created by:
! Nick Papior, 2014

module m_tbt_options

  use precision, only : dp

  use units, only: Ang

  use m_ts_tdir, only: ts_tidx

  use m_ts_electype
  use m_ts_chem_pot

  use dictionary

  implicit none

  ! Common flags for parameters
  public
  save

  ! The standard name_prefix
#ifdef TBT_PHONON
  character(len=*), parameter :: name_prefix = 'PHT'
#else
  character(len=*), parameter :: name_prefix = 'TBT'
#endif

  ! The temperature
  real(dp) :: kT

  ! Electrodes and different chemical potentials
  integer :: N_Elec = 0
  type(Elec), allocatable, target :: Elecs(:)
  integer :: N_mu = 0
  type(ts_mu), allocatable, target :: mus(:)

  ! Whether we should stop right after having created
  ! the Green's function files
  logical :: stop_after_GS = .false.

  ! Dictionary to contain the data saving methods
  ! Each key corresponds to some data calculation
  ! algorithm.
  ! To check whether data should be calculated do:
  ! if ( 'DOS-Gf' .in. save_DATA ) then
  !   calculate DOS of Gf
  ! end fi
  type(dictionary_t) :: save_DATA

  ! Number of eigenchannels to calculate
  integer :: N_eigen = 0

  ! If the energy-contour is not perfectly divisable by the number of nodes then adjust
  integer :: BTD_method = 0 ! Optimization method for determining the best tri-diagonal matrix split
  ! 0  == We optimize for speed
  ! 1  == We optimize for memory

  ! A quantity describing the accuracy of the coordinates of the 
  ! electrodes.
  ! * Should only be edited by experienced users *
  real(dp) :: Elecs_xa_EPS = 1.e-4_dp

  ! Every 5% of the calculation progress it will print an estimation
  real :: percent_tracker = 5.

#ifdef NCDF_4
  ! Save file names for data files
  character(len=250) :: cdf_fname = ' '
  character(len=250) :: cdf_fname_sigma = ' '
  character(len=250) :: cdf_fname_proj = ' '
#endif


  ! List of private formats for printing information
  character(len=*), parameter, private :: f1 ='(''tbt: '',a,t53,''='',tr4,l1)'
  character(len=*), parameter, private :: f10='(''tbt: '',a,t53,''='',tr4,a)'
  character(len=*), parameter, private :: f11='(''tbt: '',a)'
  character(len=*), parameter, private :: f12='(''tbt: '',a,t53,''='',tr2,i0)'
  character(len=*), parameter, private :: f5 ='(''tbt: '',a,t53,''='',i5,a)'
  character(len=*), parameter, private :: f20='(''tbt: '',a,t53,''='',i0,'' -- '',i0)'
  character(len=*), parameter, private :: f6 ='(''tbt: '',a,t53,''='',f10.4,tr1,a)'
  character(len=*), parameter, private :: f7 ='(''tbt: '',a,t53,''='',f12.6,tr1,a)'
  character(len=*), parameter, private :: f8 ='(''tbt: '',a,t53,''='',f10.4)'
  character(len=*), parameter, private :: f9 ='(''tbt: '',a,t53,''='',tr1,e9.3)'
  character(len=*), parameter, private :: f15='(''tbt: '',a,t53,''='',2(tr1,i0,'' x''),'' '',i0)'


contains

  subroutine read_tbt_generic(na_u, lasto)

    use m_region, only: rgn_delete
    use fdf, only: fdf_defined
    use m_ts_method, only: ts_init_regions
    ! This array is never used used, so we delete it
    use m_ts_method, only: r_pvt

    use m_tbt_diag, only: init_diag

    ! The number of atoms
    integer, intent(in) :: na_u
    ! A summated list of last orbitals on atoms.
    integer, intent(in) :: lasto(0:na_u)

    ! Initialize the buffer regions
    if ( fdf_defined('TBT.Atoms.Buffer') ) then
       call ts_init_regions('TBT',na_u,lasto)
    else
       call ts_init_regions('TS',na_u,lasto)
    end if

    call rgn_delete(r_pvt)

    ! Initialize the diagonalization method.
    call init_diag( )

  end subroutine read_tbt_generic


  ! > Reads the chemical potentials as well as the applied
  ! Bias.
  ! The bias is an intricate part of the chemical potential why it
  ! is read in here.
  subroutine read_tbt_chem_pot( )

    use fdf, only : fdf_get
    use units, only: eV, Kelvin

    use m_ts_chem_pot, only : fdf_nmu, fdffake_mu, fdf_mu, name
    
    use m_tbt_hs, only: Volt, IsVolt

    implicit none

    ! *******************
    ! * LOCAL variables *
    ! *******************
    logical :: err
    integer :: i

    ! Read in the temperature
    kT = fdf_get('ElectronicTemperature',1.9e-3_dp,'Ry')
    kT = fdf_get('TS.ElectronicTemperature',kT,'Ry')
    kT = fdf_get('TBT.ElectronicTemperature',kT,'Ry')

    ! Read in the chemical potentials
    N_mu = fdf_nmu('TBT',kT,mus)
    if ( N_mu < 1 ) then
       N_mu = fdf_nmu('TS',kT,mus)
    end if
    err = .true.
    if ( N_mu < 1 ) then
       err = .false.
       if ( IsVolt ) then
          ! There is a bias: default
          ! to two chemical potentials with the
          ! applied bias
          N_mu = fdffake_mu(mus,kT,Volt)
          
       else
          ! There is no bias. Simply create
          ! one chemical potential.
          ! This will make the electrodes
          ! default to the one chemical potential.
          N_mu = 1
          allocate(mus(1))
          mus(1)%kT = kT
          mus(1)%ckT = ' '
          mus(1)%name = 'Fermi-level'
          mus(1)%mu = 0._dp
          mus(1)%cmu = '0. eV'
          mus(1)%ID = 1
          ! These are not used, but we populate
          ! them anyway
          mus(1)%N_poles = 1
          allocate(mus(1)%Eq_seg(1))
          mus(1)%Eq_seg(1) = '*NONE'
          
       end if
       
    end if
    do i = 1 , N_mu
       ! Default things that could be of importance
       if ( fdf_mu('TBT',mus(i),kT,Volt) ) then
          ! success
       else if ( fdf_mu('TS',mus(i),kT,Volt) ) then
          ! success
       else if ( err ) then
          ! only error out if it couldn't be found and forced
          ! created
          call die('Could not find chemical potential: ' &
               //trim(name(mus(i))))
       end if
    end do

#ifdef TBT_PHONON
    ! Phonon transport cannot define different chemical potentials
    ! Furthermore, they should be zero
    do i = 1 , N_mu
       if ( abs(mus(i)%mu) > 1.e-10 * eV ) then
          call die('Phonon transport does not define chemical &
               &potentials. I.e. you cannot lift the frequency spectra.')
       end if
    end do
#endif

  end subroutine read_tbt_chem_pot


  ! Reads all information regarding the electrodes, nothing more.
  subroutine read_tbt_elec( cell, na_u, xa, lasto)

    use fdf, only : fdf_get, fdf_obsolete, fdf_deprecated, leqi
    use parallel, only : IONode
    use intrinsic_missing, only : IDX_SPC_PROJ, EYE
    use intrinsic_missing, only : VEC_PROJ_SCA, VNORM

    use m_os, only : file_exist

    use files, only: slabel
    use units, only: eV

    use m_tbt_hs, only: spin_idx

    use m_ts_chem_pot, only : copy, chem_pot_add_Elec

    use m_ts_electype, only : fdf_nElec, fdf_Elec
    use m_ts_electype, only : Name, TotUsedOrbs, TotUsedAtoms
    use m_ts_electype, only : init_Elec_sim

    use m_ts_method, only : ts_init_electrodes, a_isBuffer

    implicit none

    ! *******************
    ! * INPUT variables *
    ! *******************
    real(dp), intent(in) :: cell(3,3)
    integer,  intent(in) :: na_u, lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)

    ! *******************
    ! * LOCAL variables *
    ! *******************
    integer :: i, j
    real(dp) :: rtmp
    logical :: err, bool

    if ( N_mu == 0 ) call die('read_tbt_elecs: error in programming')

    ! To determine the same coordinate nature of the electrodes
    Elecs_xa_EPS = fdf_get('TS.Elecs.Coord.Eps',0.001_dp*Ang, 'Bohr')
    Elecs_xa_EPS = fdf_get('TBT.Elecs.Coord.Eps',Elecs_xa_EPS,'Bohr')

    ! detect how many electrodes we have
    N_Elec = fdf_nElec('TBT', Elecs)
    if ( N_Elec < 1 ) N_Elec = fdf_nElec('TS', Elecs)
    if ( N_Elec < 1 ) then
       ! We initialize to 2 electrodes (Left/Right)
       N_Elec = 2
       allocate(Elecs(N_Elec))
       Elecs(1)%name = 'Left'
       Elecs(1)%ID = 1
       Elecs(2)%name = 'Right'
       Elecs(2)%ID = 2
       ! if they do-not exist, the user will be told
       if ( IONode ) then
          write(*,'(/,''tbt: *** '',a)') 'No electrode names were found, &
               &default Left/Right are expected'
       end if
    end if

    ! Setup default parameters for the electrodes
    ! first electrode is the "left"
    ! last electrode is the "right"
    ! the remaining electrodes have their chemical potential at 0
    ! Currently the transport direction for all electrodes is the default
    ! We should probably warn if +2 electrodes are used and t_dir is the
    ! same for all electrodes... Then the user needs to know what (s)he is doing...
    Elecs(:)%Bulk = fdf_get('TS.Elecs.Bulk',.true.) ! default everything to bulk electrodes
    Elecs(:)%Bulk = fdf_get('TBT.Elecs.Bulk',Elecs(1)%Bulk)

    rtmp = fdf_get('TS.Elecs.Eta',0.001_dp*eV,'Ry')
    rtmp = fdf_get('TBT.Elecs.Eta',rtmp,'Ry')
#ifdef TBT_PHONON
    ! eta value needs to be squared as it is phonon spectrum
    if ( rtmp > 0._dp ) rtmp = rtmp ** 2
#endif
    Elecs(:)%Eta = rtmp
    
    rtmp = fdf_get('TS.Elecs.Accuracy',1.e-13_dp*eV,'Ry')
    rtmp = fdf_get('TBT.Elecs.Accuracy',rtmp,'Ry')
    Elecs(:)%accu = rtmp

    ! whether all calculations should be performed
    ! "out-of-core" i.e. whether the GF files should be created or not
    ! In tbtrans this is now defaulted to in-core
    Elecs(:)%out_of_core = fdf_get('TBT.Elecs.Out-of-core',.false.)

    ! Whether we should try and re-use the surface Green function 
    ! files
    Elecs(:)%ReUseGF = fdf_get('TS.Elecs.GF.ReUse',.true.)
    Elecs(:)%ReUseGF = fdf_get('TBT.Elecs.GF.ReUse',Elecs(1)%ReUseGF)

    ! Will stop after creating the GF files.
    stop_after_GS = fdf_get('TBT.Elecs.GF.Only',.false.)

    do i = 1 , N_Elec

       ! If we only have 2 electrodes we take them 
       ! as though the atomic indices are the first and last
       ! respectively.
       if ( N_Elec == 2 ) then
          if ( i == 1 ) then
             err = fdf_Elec('TBT',slabel,Elecs(i),N_mu,mus,idx_a= 1, &
                  name_prefix = name_prefix)
             if ( .not. err ) &
                  err = fdf_Elec('TS',slabel,Elecs(i),N_mu,mus,idx_a= 1, &
                  name_prefix = name_prefix)
          else
             err = fdf_Elec('TBT',slabel,Elecs(i),N_mu,mus,idx_a=-1, &
                  name_prefix = name_prefix)
             if ( .not. err ) &
                  err = fdf_Elec('TS',slabel,Elecs(i),N_mu,mus,idx_a=-1, &
                  name_prefix = name_prefix)
          end if
       else
          ! Default things that could be of importance
          err = fdf_Elec('TBT',slabel,Elecs(i),N_mu,mus)
          if ( .not. err ) &
               err = fdf_Elec('TS',slabel,Elecs(i),N_mu,mus, &
               name_prefix = name_prefix)
       end if
       if ( .not. err ) then
          call die('Could not find electrode: '//trim(name(Elecs(i))))
       end if
       if ( Elecs(i)%idx_a < 0 ) &
            Elecs(i)%idx_a = na_u + Elecs(i)%idx_a + 1
       if ( Elecs(i)%idx_a < 1 .or. &
            na_u < Elecs(i)%idx_a ) then
          print *,Elecs(i)%idx_a,na_u
          call die("Electrode position does not exist")
       end if
       if ( N_Elec == 2 ) then
          ! Correct for buffer atoms, first electrode steps "up"
          ! second electrode steps "down"
          if ( i == 1 ) then
             j = Elecs(i)%idx_a
             do while ( a_isBuffer(j) )
                j = j + 1
             end do
             Elecs(i)%idx_a = j
          else
             j = Elecs(i)%idx_a + TotUsedAtoms(Elecs(i)) - 1
             do while ( a_isBuffer(j) )
                j = j - 1
             end do
             Elecs(i)%idx_a = j - TotUsedAtoms(Elecs(i)) + 1
          end if
       end if
       ! set the placement in orbitals
       Elecs(i)%idx_o = lasto(Elecs(i)%idx_a-1)+1

       ! Initialize electrode parameters
       call init_Elec_sim(Elecs(i),cell,na_u,xa)

    end do

    ! Initialize the electrode regions
    call ts_init_electrodes(na_u,lasto,N_Elec,Elecs)

    ! If many electrodes, no transport direction can be specified
    ! Hence we use this as an error-check (also for N_Elec == 1)
    if ( any(Elecs(:)%t_dir > 3) ) then
      ts_tidx = - N_Elec
    else
      
      if ( N_Elec /= 2 ) then
        ! Signals no specific unit-cell direction of transport
        ts_tidx = - N_Elec
      else

        ! Retrieve the indices of the unit-cell directions
        ! according to the electrode transport directions.
        ! We have already calculated the pivoting table for
        ! the electrodes
        i = Elecs(1)%pvt(Elecs(1)%t_dir)
        j = Elecs(2)%pvt(Elecs(2)%t_dir)

        bool = i == j

        ! For a single transport direction to be true,
        ! both the projections _has_ to be 1, exactly!
        rtmp = VEC_PROJ_SCA(cell(:,i), Elecs(1)%cell(:,Elecs(1)%t_dir))
        rtmp = rtmp / VNORM(Elecs(1)%cell(:,Elecs(1)%t_dir))
        bool = bool .and. abs(abs(rtmp) - 1._dp) < 1.e-5_dp
        rtmp = VEC_PROJ_SCA(cell(:,j), Elecs(2)%cell(:,Elecs(2)%t_dir))
        rtmp = rtmp / VNORM(Elecs(2)%cell(:,Elecs(2)%t_dir))
        bool = bool .and. abs(abs(rtmp) - 1._dp) < 1.e-5_dp

        if ( bool ) then

          ! The transport direction for the electrodes are the same...
          ! And fully encompassed! We have a single transport
          ! direction.
          ts_tidx = i

        else

          ! In case we have a skewed transport direction
          ! we have some restrictions...
          ts_tidx = -N_Elec

        end if
      end if
    end if

    ! Populate the electrodes in the chemical potential type
    do i = 1 , N_Elec
       err = .true.
       do j = 1 , N_mu
          if ( associated(Elecs(i)%mu,target=mus(j)) ) then
             call chem_pot_add_Elec(mus(j),i)
             err = .false.
             exit
          end if
       end do
       if ( err ) then
          call die('We could not attribute a chemical potential &
               &to electrode: '//trim(Elecs(i)%name))
       end if
    end do

    ! check that all electrodes and chemical potentials are paired in
    ! some way.
    if ( any(mus(:)%N_El == 0) ) then
       call die('A/Some chemical potential(s) has/have not been assigned any electrodes. &
            &All chemical potentials *MUST* be assigned an electrode')
    end if

    if ( na_u <= sum(TotUsedAtoms(Elecs)) ) then
       write(*,'(a)') 'Please stop this madness. What where you thinking?'
       call die('Electrodes occupy the entire device!!!')
    end if

  end subroutine read_tbt_elec


  subroutine read_tbt_after_Elec(nspin, cell, na_u, lasto, xa, no_u, kscell, kdispl)

    use fdf, only: fdf_get, leqi
    use parallel, only: IONode

    use m_ts_method, only: ts_method, TS_BTD
    use m_ts_method, only: ts_A_method, TS_BTD_A_COLUMN, TS_BTD_A_PROPAGATION

    use m_tbt_contour, only: read_contour_options

    use m_tbt_sigma_save, only: init_Sigma_options
    use m_tbt_dH, only: init_dH_options
    use m_tbt_dSE, only: init_dSE_options

    ! *******************
    ! * INPUT variables *
    ! *******************
    integer, intent(in) :: nspin
    real(dp), intent(in) :: cell(3,3)
    integer,  intent(in) :: na_u, lasto(0:na_u)
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: no_u
    integer,  intent(in) :: kscell(3,3)
    real(dp), intent(in) :: kdispl(3)

    integer :: i
    logical :: ltmp, Gamma3(3), only_T_Gf
    character(len=150) :: chars
    
    ! we must have read the electrodes first
    if ( N_Elec == 0 ) call die('read_tbt_options: Error in programming')

    percent_tracker = fdf_get('TBT.Progress',5.)
    percent_tracker = max(0., percent_tracker)
    
    ! Reading the Transiesta solution method
    chars = fdf_get('TBT.SolutionMethod','BTD')
    if ( leqi(chars,'BTD').or.leqi(chars,'tri') ) then
      ts_method = TS_BTD
    else
      call die('Unrecognized TBtrans solution method: '//trim(chars))
    end if

    chars = fdf_get('TS.BTD.Optimize','speed')
    chars = fdf_get('TBT.BTD.Optimize',trim(chars))
    if ( leqi(chars,'speed') ) then
      BTD_method = 0
    else if ( leqi(chars,'memory') ) then
      BTD_method = 1
    else
      call die('Could not determine flag TBT.BTD.Optimize, please &
          &see manual.')
    end if
    
    ! Read spectral calculation method for BTD method
    if ( N_Elec > 3 ) then
      chars = fdf_get('TS.BTD.Spectral','column')
    else
      chars = fdf_get('TS.BTD.Spectral','propagation')
    end if
    chars = fdf_get('TBT.BTD.Spectral',trim(chars))
    if ( leqi(chars,'column') ) then
      ts_A_method = TS_BTD_A_COLUMN
    else if ( leqi(chars,'propagation') ) then
      ts_A_method = TS_BTD_A_PROPAGATION
    else
      call die('TBT.BTD.Spectral option is not column or propagation. &
          &Please correct input.')
    end if

    ! initial
    only_T_Gf = .true.
    
    ! Whether we should assert and calculate
    ! all transmission amplitudes
    ltmp = fdf_get('TBT.T.Elecs.All',N_Elec == 1)
    ltmp = fdf_get('TBT.T.All',N_Elec == 1)
    if ( ltmp ) then
      save_DATA = save_DATA // ('T-all'.kv.1)
    end if
    
    N_eigen = fdf_get('TBT.T.Eig',0)
    if ( N_eigen > 0 ) then
      save_DATA = save_DATA // ('T-eig'.kv.N_eigen)
      only_T_Gf = .false.
    end if
    
    ! Should we calculate DOS of electrode bulk Green function
    ltmp = fdf_get('TBT.DOS.Elecs', .false. )
    if ( ltmp ) then
      save_DATA = save_DATA // ('DOS-Elecs'.kv.1)
    end if
    
    ltmp = fdf_get('TBT.T.Bulk', N_Elec == 1)
    if ( ltmp ) then
      ! when calculating the DOS for the electrode
      ! we also get the bulk transmission.
      save_DATA = save_DATA // ('DOS-Elecs'.kv.1)
    end if
    
    ! Should we calculate DOS of Green function
    ltmp = fdf_get('TBT.DOS.Gf', N_Elec == 1)
    if ( ltmp ) then
      save_DATA = save_DATA // ('DOS-Gf'.kv.1)
      only_T_Gf = .false.
    end if

    ! Should we calculate DOS of spectral function
    ltmp = fdf_get('TBT.DOS.A', N_Elec == 1)
    if ( ltmp ) then
      save_DATA = save_DATA // ('DOS-A'.kv.1)
      only_T_Gf = .false.
    end if

    ! Should we calculate DOS of all spectral functions
    ltmp = fdf_get('TBT.DOS.A.All', N_Elec == 1)
    if ( ltmp ) then
      save_DATA = save_DATA // ('DOS-A'.kv.1)
      save_DATA = save_DATA // ('DOS-A-all'.kv.1)
      only_T_Gf = .false.
    end if
    
    ! Should we calculate orbital current
    ltmp = fdf_get('TBT.Current.Orb', .false. )
    if ( ltmp ) then
      save_DATA = save_DATA // ('DOS-A'.kv.1)
      save_DATA = save_DATA // ('orb-current'.kv.1)
      only_T_Gf = .false.
    end if

    ! Options for density-matrix calculations
    ltmp = fdf_get('TBT.DM.Gf', .false.)
    if ( ltmp ) then
      save_DATA = save_DATA // ('DOS-Gf'.kv.1)
      save_DATA = save_DATA // ('DM-Gf'.kv.1)
    end if

    ltmp = fdf_get('TBT.DM.A', .false.)
    if ( ltmp ) then
      save_DATA = save_DATA // ('DOS-A'.kv.1)
      save_DATA = save_DATA // ('DM-A'.kv.1)
    end if
    
    
    ! Options for COOP and COHP curves.
    ! These are orbital (energy) populations that can be used to determine the
    ! bonding nature of the material.
    ltmp = fdf_get('TBT.COOP.Gf', .false.)
    if ( ltmp ) then
      save_DATA = save_DATA // ('DOS-Gf'.kv.1)
      save_DATA = save_DATA // ('COOP-Gf'.kv.1)
    end if

    ltmp = fdf_get('TBT.COOP.A', .false. )
    if ( ltmp ) then
      save_DATA = save_DATA // ('DOS-A'.kv.1)
      save_DATA = save_DATA // ('COOP-A'.kv.1)
    end if

    ltmp = fdf_get('TBT.COHP.Gf', .false. )
    if ( ltmp ) then
      save_DATA = save_DATA // ('DOS-Gf'.kv.1)
      save_DATA = save_DATA // ('COHP-Gf'.kv.1)
    end if

    ltmp = fdf_get('TBT.COHP.A', .false. )
    if ( ltmp ) then
      save_DATA = save_DATA // ('DOS-A'.kv.1)
      save_DATA = save_DATA // ('COHP-A'.kv.1)
    end if

    ltmp = fdf_get('TBT.T.Out',N_Elec == 1)
    if ( ltmp ) then
      save_DATA = save_DATA // ('T-sum-out'.kv.1)
    end if
    
    ! We cannot calculate the transmission for more than 3
    ! electrodes if using the diagonal
    if ( only_T_Gf .and. N_Elec > 3 ) then
      only_T_Gf = .false.
    end if

    if ( only_T_Gf ) then
      ! TODO, consider changing this to .true.
      only_T_Gf = fdf_get('TBT.T.Gf',.false.)
    end if
    if ( only_T_Gf ) then
      save_DATA = save_DATA // ('T-Gf'.kv.1)
      ts_A_method = TS_BTD_A_COLUMN

      ! We get sum-out for free, so save it
      save_DATA = save_DATA // ('T-sum-out'.kv.1)

      if ( N_Elec > 2 ) then
        ! When having more than one electrode
        ! we *must* calculate all transmissions
        ! to be able to get the actual transmissions
        ! using the diagonal Green function
        save_DATA = save_DATA // ('T-all'.kv.1)
      end if

    end if

    call init_Sigma_options( save_DATA )

    if ( IONode ) write(*,*) ! new-line
    ! The init_dH_options also checks for the
    ! consecutive sparse patterns.
    call init_dH_options( no_u )
    call init_dSE_options( )

    
    ! read in contour options
    call read_contour_options( N_Elec, Elecs, N_mu, mus )

    ! Check for Gamma in each direction
    do i = 1 , 3
      if ( kdispl(i) /= 0._dp ) then
        Gamma3(i) = .false.
      else if ( sum(kscell(i,:)) > 1 ) then
        ! Note it is the off-diagonal for this unit-cell
        ! direction
        Gamma3(i) = .false.
      else
        Gamma3(i) = .true.
      end if
    end do
    
    do i = 1 , N_Elec
      ! Initialize the electrode quantities for the stored values
      call check_Elec_sim(Elecs(i), nspin, cell, na_u, xa, &
          Elecs_xa_EPS, lasto, Gamma3)
    end do

  end subroutine read_tbt_after_Elec

  subroutine print_tbt_options(nspin)

    use units, only: Kelvin, eV
    use parallel, only: IONode
    use files, only: slabel

    use m_ts_method, only: ts_A_method, TS_BTD_A_COLUMN, TS_BTD_A_PROPAGATION

    use m_tbt_contour, only: print_contour_tbt_options, io_contour_tbt
    use m_tbt_contour, only: print_contour_tbt_block
    use m_tbt_save, only: print_save_options
    use m_tbt_diag, only: print_diag
    use m_tbt_dH, only: print_dH_options
    use m_tbt_dSE, only: print_dSE_options
    use m_tbt_sigma_save, only: print_Sigma_options
    use m_tbt_proj, only: print_proj_options
    use m_tbt_hs, only: Volt, IsVolt, spin_idx

    integer, intent(in) :: nspin
    
    integer :: i
    
    if ( .not. IONode ) return

    write(*,*)
    write(*,f11) repeat('*', 62)

    write(*,f6) 'Electronic temperature (reference)',kT/Kelvin,'K'
    if ( IsVolt ) then
       write(*,f6) 'Voltage', Volt/eV,'Volts'
    else
       write(*,f11) 'No applied bias'
    end if
    write(*,f1) 'Calculate transmission only using diag(Gf)',('T-Gf'.in.save_DATA)
    write(*,f1) 'Saving bulk transmission for electrodes',('DOS-Elecs'.in.save_DATA)
    write(*,f1) 'Saving DOS from bulk electrodes',('DOS-Elecs'.in.save_DATA)
    write(*,f1) 'Saving DOS from Green function',('DOS-Gf'.in.save_DATA)
    if ( 'DOS-A-all' .in. save_DATA ) then
       write(*,f1) 'Saving DOS from all spectral functions',.true.
    else
       write(*,f1) 'Saving DOS from spectral functions',('DOS-A' .in. save_DATA)
    end if
    write(*,f1) 'Saving bond currents (orb-orb)',('orb-current'.in.save_DATA)

    ! DM
    write(*,f1) 'Saving DM from Green function',('DM-Gf'.in.save_DATA)
    write(*,f1) 'Saving DM from spectral functions',('DM-A'.in.save_DATA)

    ! COOP/COHP curves
    write(*,f1) 'Saving COOP from Green function',('COOP-Gf'.in.save_DATA)
    write(*,f1) 'Saving COOP from spectral functions',('COOP-A'.in.save_DATA)
    write(*,f1) 'Saving COHP from Green function',('COHP-Gf'.in.save_DATA)
    write(*,f1) 'Saving COHP from spectral functions',('COHP-A'.in.save_DATA)

    write(*,f12) 'Calc. # transmission eigenvalues',N_eigen
    write(*,f1) 'Calc. T between all electrodes',('T-all'.in.save_DATA)
    write(*,f1) 'Calc. total T out of electrodes',('T-sum-out'.in.save_DATA)
    if ( nspin > 1 ) then
       if ( spin_idx == 0 ) then
          write(*,f11) 'Calculate spin UP and DOWN'
       else if ( spin_idx == 1 ) then
          write(*,f10) 'Calculate spin ','UP'
       else if ( spin_idx == 2 ) then
          write(*,f10) 'Calculate spin ','DOWN'
       else
          call die('Error in spin_idx')
       end if
    else
       write(*,f11) 'Non-polarized Hamiltonian'
    end if


    ! Algorithm choices...
    select case ( BTD_method )
    case ( 0 )
       write(*,f10)'BTD creation algorithm', 'speed'
    case ( 1 )
       write(*,f10)'BTD creation algorithm', 'memory'
    end select
    select case ( ts_A_method )
    case ( TS_BTD_A_PROPAGATION )
       write(*,f10)'BTD spectral function algorithm','propagation'
    case ( TS_BTD_A_COLUMN )
       write(*,f10)'BTD spectral function algorithm','column'
    case default
       call die('Error in setup BTD. A calc')
    end select

    
    call print_diag()
    call print_Sigma_options( save_DATA )
    call print_dH_options( )
    call print_dSE_options( )

    call print_save_options()
    call print_proj_options( save_DATA )

    write(*,f11)'          >> Electrodes << '
    do i = 1 , size(Elecs)
       call print_settings(Elecs(i),'tbt')
    end do
    
    call print_contour_tbt_options( 'TBT' )
    
    write(*,f11) repeat('*', 62)
    write(*,*)
    
    call io_contour_tbt(slabel)

    write(*,f11) repeat('<', 62)

    call print_contour_tbt_block( 'TBT' )

    write(*,f11) repeat('<', 62)

  end subroutine print_tbt_options

  subroutine print_tbt_warnings( Gamma )

    use fdf, only: fdf_get, fdf_defined
    use units, only: eV
    use parallel, only: IONode

    use m_tbt_save, only: save_parallel

    use m_tbt_dH, only: print_dH_warnings
    use m_tbt_hs, only: Volt

    ! Whether the user requests a Gamma calculation
    logical, intent(in) :: Gamma

    integer :: i
    logical :: ltmp, rem_DOS_Elecs, rem_T_Gf, has

    ! Removal of keys
    rem_DOS_Elecs = 'DOS-Elecs' .in. save_DATA
    if ( rem_DOS_Elecs ) then
      ! Currently the TBTGF files does not contain the DOS
      if ( all(Elecs(:)%out_of_core) ) &
          call delete(save_DATA,'DOS-Elecs')
      if ( 'Sigma-only' .in. save_DATA ) &
          call delete(save_DATA,'DOS-Elecs')
    end if
    
    rem_T_Gf = 'T-Gf' .in. save_DATA
    if ( N_Elec > 3 ) then
      call delete(save_DATA,'T-Gf')
    end if

    if ( .not. IONode ) return

    write(*,*) ! new-line
    write(*,'(3a)') repeat('*',24),' Begin: TBT CHECKS AND WARNINGS ',repeat('*',24)

    if ( N_eigen < 0 ) then
       call die('Number of transmission eigenvalues MUST be &
            &zero or positive.')
    end if

    if ( .not. Gamma ) then
       ltmp = .not. fdf_get('SpinSpiral',.false.)
       ltmp = fdf_get('TBT.Symmetry.TimeReversal',ltmp)
       if ( ('orb-current' .in.save_DATA) ) then
          if ( IONode .and. ltmp ) then
             write(*,'(a,/,a)') 'WARNING: k-averaging orbital currents with &
                  &time-reversal symmetry will not reproduce','the correct &
                  &orbital current. Set TBT.Symmetry.TimeReversal F'
          end if
       end if

       ! Also for COOP analysis
       has = 'COOP-Gf' .in. save_DATA
       has = has .or. ('COOP-A' .in. save_DATA)
       has = has .or. ('COHP-Gf' .in. save_DATA)
       has = has .or. ('COHP-A' .in. save_DATA)
       if ( IONode .and. ltmp .and. has ) then
          write(*,'(a,/,a)') 'WARNING: k-averaging COOP/COHP with &
               &time-reversal symmetry will not reproduce','the correct &
               &populations. Set TBT.Symmetry.TimeReversal F'
       end if
    end if

    ! CHECK THIS (we could allow it by only checking the difference...)
    if (  maxval(mus(:)%mu) - minval(mus(:)%mu) - abs(Volt) > 1.e-8_dp * eV ) then
       if ( IONode ) then
          write(*,'(a)') 'Chemical potentials [eV]:'
          do i = 1 , N_Elec
             write(*,'(a,f10.5,a)') trim(Name(Elecs(i)))//' at ',Elecs(i)%mu%mu/eV,' eV'
          end do
          write(*,'(a)') 'The difference must satisfy: "max(ChemPots)-min(ChemPots) - abs(Volt) < 1e-8 eV"'
          write(*,'(a,f10.5,a)') 'max(ChemPots) at ', maxval(mus(:)%mu)/eV,' eV'
          write(*,'(a,f10.5,a)') 'min(ChemPots) at ', minval(mus(:)%mu)/eV,' eV'
          write(*,'(a,f10.5,a)') '|V| at ', abs(Volt)/eV,' eV'
       end if
       call die('Chemical potentials are not consistent with the bias applied.')
    end if


    ! Check that the bias does not introduce a gating
    if ( any(abs(mus(:)%mu) - abs(Volt) > 1.e-9_dp ) ) then
       write(*,'(a)') 'Chemical potentials must lie in the range [-V;V] with the maximum &
            &difference being V'
       call die('Chemical potentials must not introduce consistent Ef shift to the system.')
    end if

    if ( rem_DOS_Elecs ) then
       if ( any(Elecs(:)%out_of_core) ) then
          write(*,'(a)')' Disabling electrode DOS calculation, only &
               &enabled for in-core self-energy calculations.'
       end if
       if ( 'Sigma-only' .in. save_DATA ) then
          write(*,'(a)')' Disabling electrode DOS calculation, only &
               &enabled when calculating transmission (not TBT.SelfEnergy.Only).'
       end if
    end if
    
    if ( rem_T_Gf .and. N_Elec > 3 ) then
       write(*,'(a)')' ** Disabling transport calculation using diagonal, &
            &not possible with N_elec > 3.'
    end if

    do i = 1, N_Elec
      if ( Elecs(i)%repeat .and. Elecs(i)%bloch%size() > 1 ) then
        write(*,'(a)') 'Electrode '//trim(Elecs(i)%name)//' is &
            &using Bloch unfolding using the repeat scheme! &
            &Please use the tiling scheme (it is orders of magnitudes faster!).'
      end if
    end do

#ifdef MPI
#ifdef NCDF_PARALLEL
    if ( .not. save_parallel ) then
       write(*,'(a)') ' ** Speed up the execution by utilizing parallel I/O'
       write(*,'(a)') '  > TBT.CDF.MPI true'
    end if
#endif
#endif

    call print_dH_warnings( save_DATA )

    write(*,'(3a,/)') repeat('*',24),' End: TBT CHECKS AND WARNINGS ',repeat('*',26)

  end subroutine print_tbt_warnings

end module m_tbt_options
