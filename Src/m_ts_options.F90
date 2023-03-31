! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

module m_ts_options

  use precision, only : dp

  use m_mixing, only: tMixer

  use m_ts_electype, only : Elec
  use m_ts_chem_pot, only : ts_mu
  use m_ts_tdir

  implicit none

  public
  save

  ! ###### SIESTA-options ######
  ! The following options override the siesta settings upon
  ! entering the transiesta SCF

  ! The tolerance
  real(dp) :: ts_Dtol ! = tolerance for density matrix
  real(dp) :: ts_Htol ! = tolerance for Hamiltonian
  logical :: ts_converge_dQ = .true. ! whether we should converge the charge
  real(dp) :: ts_dQtol ! = tolerance for charge in the device region (after SCF)
  integer :: ts_hist_keep = 0

  ! ###### end SIESTA-options #####

  ! Whether we should stop before transiesta begins...
  logical :: TS_siesta_stop = .false.

  !< Whether TranSiesta is allowed to start without the 0-bias calculation.
  logical :: TS_start_bias = .false.

  ! Controls to save the TSHS file
  logical :: TS_HS_save = .true.
  logical :: TS_DE_save = .false.
  ! whether we will use the bias-contour
  logical :: IsVolt = .false.
  ! maximum difference between chemical potentials
  real(dp) :: Volt = 0._dp
  ! The temperature for transiesta calculations
  real(dp) :: ts_kT
  ! Electrodes and different chemical potentials
  integer :: N_Elec = 0
  type(Elec), allocatable, target :: Elecs(:)
  integer :: N_mu = 0
  type(ts_mu), allocatable, target :: mus(:)

  integer :: TS_scf_mode = 0

  ! Flag to control whether we should update the forces (i.e. calculate energy-density matrix)
  logical :: Calc_Forces = .true.

  ! If the energy-contour is not perfectly divisable by the number of nodes then adjust
  integer :: BTD_method = 0 ! Optimization method for determining the best tri-diagonal matrix split
  ! 0  == We optimize for speed
  ! 1  == We optimize for memory

  ! File name for reading in the grid for the Hartree potential
  character(len=150) :: Hartree_fname = ' '

  ! A quantity describing the accuracy of the coordinates of the 
  ! electrodes.
  ! * Should only be edited by experienced users *
  real(dp) :: Elecs_xa_EPS = 1.e-4_dp

  ! The user can request to analyze the system, returning information about the 
  ! tri-diagonalization partition and the contour
  logical :: TS_Analyze = .false.

  ! List of private formats for printing information
  character(len=*), parameter, private :: f1 ='(''ts: '',a,t53,''='',tr4,l1)'
  character(len=*), parameter, private :: f10='(''ts: '',a,t53,''='',tr4,a)'
  character(len=*), parameter, private :: f11='(''ts: '',a)'
  character(len=*), parameter, private :: f5 ='(''ts: '',a,t53,''='',i5,a)'
  character(len=*), parameter, private :: f20='(''ts: '',a,t53,''='',i0,'' -- '',i0)'
  character(len=*), parameter, private :: f6 ='(''ts: '',a,t53,''='',f10.4,tr1,a)'
  character(len=*), parameter, private :: f7 ='(''ts: '',a,t53,''='',f12.6,tr1,a)'
  character(len=*), parameter, private :: f8 ='(''ts: '',a,t53,''='',f10.4)'
  character(len=*), parameter, private :: f9 ='(''ts: '',a,t53,''='',tr1,e9.3)'
  character(len=*), parameter, private :: f15='(''ts: '',a,t53,''='',2(tr1,i0,'' x''),'' '',i0)'

  ! set copy of SCF mixing
  type(tMixer), pointer :: ts_scf_mixs(:) => null()

contains


  ! We have separated the options routines to only deal 
  ! with their respective parts

  ! > Read generic options for transiesta which has
  ! nothing to do with the electrodes or the chemical potentials
  subroutine read_ts_generic( cell )

    use fdf, only : fdf_get, leqi
    use intrinsic_missing, only : VNORM

    use siesta_options, only : dDtol, dHtol

    use m_ts_global_vars, only: TSmode, onlyS
    use m_ts_method, only: TS_FULL, TS_BTD, TS_MUMPS, ts_method

    use m_ts_weight, only : read_ts_weight
    use ts_dq_m, only : ts_dq_read

#ifdef SIESTA__MUMPS
    use m_ts_mumps_init, only : read_ts_mumps
#endif

    use m_mixing, only: mixers_init
    use m_mixing_scf, only: scf_mixs

    ! Input variables
    real(dp), intent(in) :: cell(3,3)

    ! Local variables
    character(len=200) :: chars

    ! This has to be the first routine to be read
    if ( N_mu /= 0 ) call die('read_ts_generic: error in programming')
    if ( N_Elec /= 0 ) call die('read_ts_generic: error in programming')

    ! Read in general values that should be used in the electrode generation
    ! I.e. these FDF-parameters are used for diagon runs with transiesta
#ifdef TRANSIESTA
    TS_HS_save = fdf_get('TS.HS.Save',.true.)
    TS_DE_save = fdf_get('TS.DE.Save',.true.)
#else
    TS_HS_save = fdf_get('TS.HS.Save',.false.)
    TS_DE_save = fdf_get('TS.DE.Save',.false.)
#endif
    onlyS = fdf_get('TS.onlyS',.false.)
    onlyS = fdf_get('TS.S.Save',onlyS)

    ! Immediately return from transiesta when this occurs
    ! no settings from the intrinsic transiesta routines
    ! are needed.
    if ( onlyS .or. .not. TSmode ) return

    ! When running TSmode we FORCE TS.HS.Save and TS.DE.Save
    TS_HS_save = .true.
    TS_DE_save = .true.

    ! Force the run of a biased TranSiesta run
    ! from a pristine siesta calculation.
    TS_start_bias = fdf_get('TS.Voltage.FromSiesta', .false.)

    ! Read in the transiesta SCF mixing options
    call mixers_init('TS.SCF', ts_scf_mixs )
    if ( .not. associated(ts_scf_mixs) ) then
       ts_scf_mixs => scf_mixs
    end if
    
    ! Read in the mixing for the transiesta cycles
    ts_Dtol = fdf_get('TS.SCF.DM.Tolerance',dDTol)
    ts_Htol = fdf_get('TS.SCF.H.Tolerance',dHTol)
    ts_hist_keep = fdf_get('TS.SCF.Mix.History.Keep',0)

    ! Stop after siesta has converged
    TS_siesta_stop = fdf_get('TS.SIESTA.Only',.false.)

    ! Reading the Transiesta solution method
    chars = fdf_get('TS.SolutionMethod','BTD')
    if ( leqi(chars,'full') ) then
       ts_method = TS_FULL
    else if ( leqi(chars,'BTD') .or. leqi(chars,'tri') ) then
       ts_method = TS_BTD
#ifdef SIESTA__MUMPS
    else if ( leqi(chars,'mumps') ) then
       ts_method = TS_MUMPS
#endif
    else
       call die('Unrecognized TranSiesta solution method: '//trim(chars))
    end if

    ! currently this does not work
    chars = fdf_get('SCF.Initialize','diagon')
    chars = fdf_get('TS.SCF.Initialize',chars)
    if ( leqi(chars,'diagon') ) then
       TS_scf_mode = 0
       chars = 'none'
    else if ( leqi(chars,'transiesta') ) then
       TS_scf_mode = 1
       chars = 'init'
    end if

    chars = fdf_get('TS.BTD.Optimize','speed')
    if ( leqi(chars,'speed') ) then
       BTD_method = 0
    else if ( leqi(chars,'memory') ) then
       BTD_method = 1
    else
       call die('Could not determine flag TS.BTD.Optimize, please &
            &see manual.')
    end if

    ! Determine whether the user wishes to only do an analyzation
    TS_Analyze = fdf_get('TS.Analyze',.false.)

    ! Read charge-correction methods
    call ts_dq_read( )

#ifdef SIESTA__MUMPS
    call read_ts_mumps( )
#endif

    call read_ts_weight( )

    ! whether to calculate the forces or not (default calculate everything)
    Calc_Forces = fdf_get('TS.Forces',.true.)

  end subroutine read_ts_generic

  
  ! > Reads the chemical potentials as well as the applied
  ! Bias.
  ! The bias is an intricate part of the chemical potential why it
  ! is read in here.
  subroutine read_ts_chem_pot( )

    use fdf, only : fdf_get
    use units, only: eV, Kelvin

    use siesta_options, only : kT => Temp

    use m_ts_global_vars, only: TSmode, onlyS
    
    use m_ts_chem_pot, only : fdf_nmu, fdffake_mu, fdf_mu, name

    implicit none
    
! *******************
! * LOCAL variables *
! *******************
    logical :: err
    integer :: i, j
    real(dp) :: rtmp

    if ( onlyS .or. .not. TSmode ) return

    ts_kT = fdf_get('TS.ElectronicTemperature',kT,'Ry')
    if ( ts_kT / Kelvin < 10._dp ) then
       call die('TranSiesta electronic temperature *must* &
            &be larger than 10 kT')
    end if
    
    ! The sign can not be chosen from this (several mu, where to define it)
    Volt   = fdf_get('TS.Voltage',0._dp,'Ry') 
    ! Voltage situation is above 0.01 mV
    IsVolt = abs(Volt) > 0.00001_dp * eV

    ! Read in the chemical potentials
    N_mu = fdf_nmu('TS',ts_kT,mus)
    err = .true.
    if ( N_mu < 1 ) then
       err = .false.
       N_mu = fdffake_mu(mus,ts_kT,Volt)
    end if
    do i = 1 , N_mu
       ! Default things that could be of importance
       if ( err .and. .not. fdf_mu('TS',mus(i),ts_kT,Volt) ) then
          call die('Could not find chemical potential: ' &
               //trim(name(mus(i))))
       end if
       
       ! We do not allow the electronic temperature
       ! to be below 10 kT
       if ( mus(i)%kT / Kelvin < 10._dp ) then
          call die('TranSiesta electronic temperature *must* &
               &be larger than 10 kT')
       end if
       
    end do

    ! We consider 10 Kelvin to be the minimum allowed
    ! temperature difference of the leads.
    rtmp = 10._dp * Kelvin
    do j = 1 , N_mu - 1
       do i = j + 1 , N_mu
          if ( abs(mus(i)%kT - mus(j)%kT) > rtmp ) then
             ! If there exists a temperature gradient
             ! we are in non-equilibrium, hence we need the 
             ! bias-setup no matter V!
             IsVolt = .true.
             exit
          end if
       end do
    end do

  end subroutine read_ts_chem_pot


  ! Reads all information regarding the electrodes, nothing more.
  subroutine read_ts_elec( cell, na_u, xa, lasto)

    use fdf, only : fdf_get, fdf_obsolete, fdf_deprecated, leqi
    use parallel, only : IONode
    use intrinsic_missing, only : IDX_SPC_PROJ, EYE, VNORM
    use intrinsic_missing, only : VEC_PROJ_SCA

    use m_os, only : file_exist

    use files, only: slabel
    use units, only: eV, Ang

    use m_ts_global_vars, only: TSmode, onlyS
    
    use m_ts_chem_pot, only : copy, chem_pot_add_Elec

    use m_ts_electype, only : fdf_nElec, fdf_Elec
    use m_ts_electype, only : Name, TotUsedOrbs, TotUsedAtoms
    use m_ts_electype, only : init_Elec_sim

    use m_ts_method, only : ts_init_electrodes, a_isBuffer
    use m_cite, only : add_citation

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
    character(len=200) :: chars, c
    type(ts_mu) :: tmp_mu

    if ( onlyS .or. .not. TSmode ) return

    if ( N_mu == 0 ) call die('read_ts_elecs: error in programming')

    ! the title of the green's functions are now non-generic
    call fdf_obsolete('TS.GFTitle')

    ! The regular options for describing the electrodes can not be 
    ! used anymore...
    call fdf_obsolete('TS.HSFileLeft')
    call fdf_obsolete('TS.GFFileLeft')
    call fdf_obsolete('TS.NumUsedAtomsLeft')
    call fdf_obsolete('TS.ReplicateA1Left')
    call fdf_obsolete('TS.ReplicateA2Left')
    call fdf_obsolete('TS.HSFileRight')
    call fdf_obsolete('TS.GFFileRight')
    call fdf_obsolete('TS.NumUsedAtomsRight')
    call fdf_obsolete('TS.ReplicateA1Right')
    call fdf_obsolete('TS.ReplicateA2Right')

    ! notice that this does not have the same meaning... 
    call fdf_deprecated('TS.UpdateDMCROnly','TS.Elecs.DM.Update')
    call fdf_deprecated('TS.UseBulk','TS.Elecs.Bulk')

    ! whether or not the electrodes should be re-instantiated
    call fdf_deprecated('TS.CalcGF','TS.Elecs.GF.ReUse')
    call fdf_deprecated('TS.ReUseGF','TS.Elecs.GF.ReUse')

    ! To determine the same coordinate nature of the electrodes
    Elecs_xa_EPS = fdf_get('TS.Elecs.Coord.Eps',0.001_dp*Ang, 'Bohr')

    ! detect how many electrodes we have
    N_Elec = fdf_nElec('TS',Elecs)
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
        c = '(''transiesta: *** '',a,/)'
        write(*,c)'No electrode names were found, default Left/Right are expected'
      end if
    end if
    
    ! If only one electrode you are not allowed to move the Fermi-level
    ! of the electrode. That should be done by other means (i.e. use NetCharge)
    if ( N_Elec == 1 ) then
      ! Notice that below the chemical potential gets corrected
      ! EVEN if the user supplied a bias.
      if ( IsVolt .and. IONode ) then
        c = '(''transiesta: *** '',a)'
        write(*,c) 'Single electrode calculations does not allow shifting the chemical potential.'
        write(*,c) 'You should do that by changing the states filled in the system.'
        write(*,c) 'Consult the manual on how to do this.'
        call die('Please set the chemical potential to zero for your one electrode')
      end if
    end if
    
    ! Setup default parameters for the electrodes
    ! first electrode is the "left"
    ! last electrode is the "right"
    ! the remaining electrodes have their chemical potential at 0
    ! We should probably warn if +2 electrodes are used and t_dir is the
    ! same for all electrodes... Then the user needs to know what (s)he is doing...
    ! Accuracy required for self-energy convergence
    Elecs(:)%accu = fdf_get('TS.Elecs.Accuracy',1.e-13_dp*eV,'Ry')
    Elecs(:)%Eta  = fdf_get('TS.Elecs.Eta',0.001_dp*eV,'Ry')
    Elecs(:)%Bulk = fdf_get('TS.Elecs.Bulk',.true.) ! default everything to bulk electrodes
    if ( Elecs(1)%Bulk ) then
      ! Default is cross-terms if we use bulk electrodes
      chars = 'cross-terms'
    else
      ! For non-bulk systems, we default to all
      chars = 'all'
    end if
    c = fdf_get('TS.Elecs.DM.Update',trim(chars))
    if ( leqi(c,'none') ) then
      Elecs(:)%DM_update = 0
    else if ( leqi(c,'cross-terms') .or. &
        leqi(c,'cross-term') ) then
      Elecs(:)%DM_update = 1
    else if ( leqi(c,'all') ) then
      Elecs(:)%DM_update = 2
    else
      call die('TS.Elecs.DM.Update [cross-terms,none,all]: &
          &unrecognized option: '//trim(c))
    end if

    ! Whether we should always set the DM to bulk
    ! values (by reading in from electrode DM)
    if ( TS_scf_mode == 1 .and. .not. IsVolt ) then
      chars = 'bulk'
    else
      chars = 'diagon'
    end if
    chars = fdf_get('TS.Elecs.DM.Init',trim(chars))
    if ( leqi(chars,'diagon') ) then
      Elecs(:)%DM_init = 0
    else if ( leqi(chars,'bulk') ) then
      Elecs(:)%DM_init = 1
    else if ( leqi(chars,'force-bulk') ) then
      Elecs(:)%DM_init = 2
    else
      call die('TS.Elecs.DM.Init unknown value [diagon,bulk,force-bulk]')
    end if
    if ( IsVolt ) then
      if ( Elecs(1)%DM_init == 1 .and. IONode) then
        if ( IONode ) then
          c = '(''transiesta: *** '',a,/)'
          write(*,c)'Will default to not read in electrode DM, only applicable for V = 0 calculations'
        end if
      end if
      Elecs(:)%DM_init = 0
    end if

    ! Whether we should try and re-use the surface Green function 
    ! files
    Elecs(:)%ReUseGF = fdf_get('TS.Elecs.GF.ReUse',.true.)

    ! whether all calculations should be performed
    ! "out-of-core" i.e. whether the GF files should be created or not
    Elecs(:)%out_of_core = fdf_get('TS.Elecs.Out-of-core',.true.)

    do i = 1 , N_Elec

       ! If we only have 2 electrodes we take them 
       ! as though the atomic indices are the first and last
       ! respectively.
       if ( N_Elec == 2 ) then
          if ( i == 1 ) then
             err = fdf_Elec('TS',slabel,Elecs(i),N_mu,mus,idx_a= 1)
          else if ( i == 2 ) then
             err = fdf_Elec('TS',slabel,Elecs(i),N_mu,mus,idx_a=-1)
          end if
       else
          ! Default things that could be of importance
          err = fdf_Elec('TS',slabel,Elecs(i),N_mu,mus)
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

       ! Initialize the electrode quantities for the
       ! stored values
       call init_Elec_sim(Elecs(i), cell, na_u, xa)

    end do

    ! If many electrodes, no transport direction can be specified
    ! Hence we use this as an error-check (also for N_Elec == 1)
    if ( any(Elecs(:)%t_dir > 3) ) then
      ts_tidx = - N_Elec

      ! We add the real-space self-energy article
      if ( IONode ) then
        call add_citation("10.1103/PhysRevB.100.195417")
      end if

    else
      select case ( N_Elec )
      case ( 1 )
        ! The easy case
        ! We simple need to figure out if the electrode
        ! has its transport direction aligned with the
        ! lattice vectors

        i = Elecs(1)%pvt(Elecs(1)%t_dir)

        ! For a single transport direction to be true,
        ! both the projections _has_ to be 1, exactly!
        rtmp = VEC_PROJ_SCA(cell(:,i), Elecs(1)%cell(:,Elecs(1)%t_dir))
        rtmp = rtmp / VNORM(Elecs(1)%cell(:,Elecs(1)%t_dir))
        bool = abs(abs(rtmp) - 1._dp) < 1.e-5_dp

        if ( bool ) then

          ! The transport direction for the electrodes are the same...
          ! And fully encompassed! We have a single transport
          ! direction.
          ts_tidx = i

        else

          ! In case we have a skewed transport direction
          ! we have some restrictions...
          ts_tidx = - N_Elec

        end if

      case ( 2 )

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
          ts_tidx = - N_Elec

        end if

      case default

        ! N_Elec > 2
        ! Here we always have these settings
        ts_tidx = - N_Elec

      end select
    end if

    ! The user can selectively decide how the bias
    ! is applied.
    ! For N-terminal calculations we advice the user
    ! to use a Poisson solution they add.
    if ( ts_tidx > 0 ) then
       ! We have a single unified semi-inifinite direction
       chars = fdf_get('TS.Poisson','ramp')
    else
       chars = fdf_get('TS.Poisson','elec-box')
    end if

#ifdef NCDF_4
    if ( file_exist(chars, Bcast = .true.) ) then
       
       Hartree_fname = trim(chars)
       ts_tidx = 0
       
    else
#endif
       Hartree_fname = ' '
       if ( leqi(chars,'ramp') .or. leqi(chars, 'ramp-cell') ) then
          if ( ts_tidx <= 0 ) then
             call die('TS.Poisson cannot be ramp for &
                  &anything but 2-electrodes with aligned transport direction.')
          end if
       else if ( leqi(chars,'elec-box') ) then
          ts_tidx = - N_Elec
       else
#ifdef NCDF_4
          call die('Error in specifying how the Hartree potential &
               &should be placed. [ramp|elec-box|NetCDF-file]')
#else
          call die('Error in specifying how the Hartree potential &
               &should be placed. [ramp|elec-box]')
#endif
       end if
#ifdef NCDF_4
    end if
#endif

    ! In case we are doing equilibrium, fix the chemical potential to the first
    if ( .not. IsVolt ) then

       ! force it to be zero... can be necessary if considering single electrode
       ! calculations (assures V == 0)
       Volt = 0._dp

       ! copy over electrode...
       call copy(mus(1),tmp_mu)
       
       ! Deallocate all strings
       do i = 1 , N_mu
          deallocate(mus(i)%Eq_seg)
       end do
       deallocate(mus)

       ! create the first chemical potential again
       N_mu = 1
       allocate(mus(1))
       call copy(tmp_mu,mus(1))
       deallocate(tmp_mu%Eq_seg)

       ! Firmly assure the chemical potential to be zero
       mus(1)%mu = 0._dp
       mus(1)%ID = 1
       
       ! Assign all electrodes to the same chemical potential
       do i = 1 , N_Elec
          Elecs(i)%mu => mus(1)
       end do

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

    if ( na_u <= sum(TotUsedAtoms(Elecs)) ) then
      write(*,'(a)') 'Please stop this madness. What where you thinking?'
      write(*,*) na_u, sum(TotUsedAtoms(Elecs))
      call die('Electrodes occupy the entire device!!!')
    end if
    
    ! Initialize the electrode regions
    call ts_init_electrodes(na_u,lasto,N_Elec,Elecs)

  end subroutine read_ts_elec


  ! This routine reads options for transiesta after
  ! having read in the electrodes
  ! It allows one to do things between reading options
  subroutine read_ts_after_Elec( cell, nspin, na_u, xa, lasto, &
       ts_kscell, ts_kdispl)

    use fdf, only : fdf_get, leqi
    use intrinsic_missing, only : VNORM, IDX_SPC_PROJ, EYE
    use atomlist, only : qa
    use siesta_options, only : charnet

    use m_ts_global_vars, only: TSmode, onlyS

    use m_ts_method, only: TS_BTD_A_PROPAGATION, TS_BTD_A_COLUMN
    use m_ts_method, only: ts_A_method, a_isBuffer

    use m_ts_electype, only: check_Elec_sim
    
    use m_ts_contour, only: read_contour_options
    use m_ts_contour_eq, only: N_Eq, Eq_c
    
    use m_ts_weight, only : read_ts_weight
    use ts_dq_m, only : ts_dq_read
    
    use m_ts_hartree, only: read_ts_hartree_options

#ifdef SIESTA__MUMPS
    use m_ts_mumps_init, only : read_ts_mumps
#endif
    
    ! Input variables
    real(dp), intent(in) :: cell(3,3)
    integer, intent(in) :: nspin, na_u
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: lasto(0:na_u)
    integer, intent(in) :: ts_kscell(3,3)
    real(dp), intent(in) :: ts_kdispl(3)

    ! Local variables
    character(len=50) :: chars
    integer :: i
    real(dp) :: dev_q
    logical :: Gamma3(3)

    if ( onlyS .or. .not. TSmode ) return

    if ( N_Elec == 0 ) call die('read_ts_after_Elecs: error in programming')

    ! Read spectral calculation method for BTD method
    if ( N_Elec > 3 ) then
       chars = fdf_get('TS.BTD.Spectral','column')
    else
       chars = fdf_get('TS.BTD.Spectral','propagation')
    end if
    if ( leqi(chars,'column') ) then
       ts_A_method = TS_BTD_A_COLUMN
    else if ( leqi(chars,'propagation') ) then
       ts_A_method = TS_BTD_A_PROPAGATION
    else
       call die('TS.BTD.Spectral option is not column or propagation. &
            &Please correct input.')
    end if
    
    ! Read in options again, at this point we have
    ! the correct ts_tidx
    call read_ts_hartree_options()
    
    ! read in contour options
    call read_contour_options( N_Elec, Elecs, N_mu, mus, ts_kT, IsVolt, Volt )

    ! Quickly update the forces update if continued fraction
    ! is used, one need an extra Green function evaluation to
    ! get the forces correct (iR and -R)
    do i = 1 , N_Eq
       if ( leqi(Eq_c(i)%c_io%part,'cont-frac') ) then
          ! the forces are not updated, dispite the user requests
          calc_forces = .false.
       end if
    end do

    ! Check for Gamma in each direction
    do i = 1 , 3
       if ( ts_kdispl(i) /= 0._dp ) then
          Gamma3(i) = .false.
       else if ( sum(ts_kscell(i,:)) > 1 ) then
          ! Note it is the off-diagonal for this direction
          Gamma3(i) = .false.
       else
          Gamma3(i) = .true.
       end if
    end do
          
    do i = 1 , N_Elec
      ! Initialize the electrode quantities for the
      ! stored values
      call check_Elec_sim(Elecs(i), nspin, cell, na_u, xa, &
          Elecs_xa_EPS, lasto, Gamma3, ts_kscell, ts_kdispl)
    end do

    ! Now we know which atoms are buffer's and electrodes
    ! This allows us to decide the tolerance for convergence of
    ! the charges.
    ! By default we allow a difference up to q(device) / 1000
    dev_q = - charnet
    do i = 1, na_u
      if ( .not. a_isBuffer(i) ) then
        ! count this atom
        dev_q = dev_q + qa(i)
      end if
    end do
    ts_dQtol = fdf_get('TS.SCF.dQ.Tolerance',dev_q / 1000._dp)
    ts_converge_dQ = fdf_get('TS.SCF.dQ.Converge', .true.)

  end subroutine read_ts_after_Elec

  
  subroutine print_ts_options( cell )

    use fdf, only: fdf_get, leqi
    use parallel, only: IOnode

    use intrinsic_missing, only : IDX_SPC_PROJ, EYE
    use units, only: eV, Kelvin

    use m_mixing, only: mixers_print
    use m_mixing_scf, only: scf_mixs

    use m_ts_electype, only: print_settings

    use m_ts_global_vars, only: TSmode, onlyS

    use m_ts_contour, only: print_contour_options

    use m_ts_method, only: TS_FULL, TS_BTD, TS_MUMPS, ts_method, na_Buf
    use m_ts_method, only: TS_BTD_A_COLUMN, TS_BTD_A_PROPAGATION
    use m_ts_method, only: ts_A_method

    use ts_dq_m, only: TS_DQ_METHOD, TS_DQ_METHOD_BUFFER, TS_DQ_METHOD_FERMI
    use ts_dq_m, only: TS_DQ_FACTOR, TS_DQ_FERMI_TOLERANCE
    use ts_dq_m, only: TS_DQ_FERMI_MAX, TS_DQ_FERMI_ETA

    use m_ts_weight, only: TS_W_METHOD, TS_W_CORRELATED
    use m_ts_weight, only: TS_W_ORB_ORB, TS_W_TR_ATOM_ATOM, TS_W_SUM_ATOM_ATOM
    use m_ts_weight, only: TS_W_TR_ATOM_ORB, TS_W_SUM_ATOM_ORB, TS_W_MEAN
    use m_ts_weight, only: TS_W_K_METHOD
    use m_ts_weight, only: TS_W_K_CORRELATED, TS_W_K_UNCORRELATED

#ifdef SIESTA__MUMPS
    use m_ts_mumps_init, only: MUMPS_mem, MUMPS_ordering, MUMPS_block
#endif

    use m_ts_hartree, only: TS_HA_PLANES, TS_HA_frac, TS_HA_offset

    implicit none

    ! The unit-cell
    real(dp), intent(in) :: cell(3,3)
    
! *******************
! * LOCAL variables *
! *******************
    character(len=200) :: chars
    logical :: ltmp
    real(dp) :: tmp33(3,3)
    integer :: i, tdir

    if ( .not. IONode ) return

    write(*,*)
    write(*,f11) repeat('*', 62)

    if ( onlyS .or. .not. TSmode ) then
       
       write(*,f1) 'Save H and S matrices', TS_HS_save
       write(*,f1) 'Save DM and EDM matrices', TS_DE_save
       write(*,f1) 'Only save the overlap matrix S', onlyS

       write(*,f11) repeat('*', 62)
       write(*,*)

       return
    end if

    write(*,f6) 'Voltage', Volt/eV,'Volts'
    write(*,f1) 'Save H and S matrices', TS_HS_save
    if ( TS_Analyze ) then
       write(*,f11)'Will analyze bandwidth of LCAO sparse matrix and quit'
    end if
    write(*,f7) 'Electronic temperature',ts_kT/Kelvin,'K'
    if ( ts_tidx < 1 ) then
       write(*,f11) 'Transport individually selected for electrodes'
    else
       write(chars,'(a,i0)') 'A',ts_tidx
       write(*,f10) 'Transport along unit-cell vector',trim(chars)
       ! Calculate Cartesian transport direction
       call eye(3,tmp33)
       tdir = IDX_SPC_PROJ(tmp33,cell(:,ts_tidx),mag=.true.)
       select case ( tdir )
       case ( 1 )
          write(*,f10) 'Transport along Cartesian vector','X'
       case ( 2 )
          write(*,f10) 'Transport along Cartesian vector','Y'
       case ( 3 )
          write(*,f10) 'Transport along Cartesian vector','Z'
       end select
     end if
     chars = ' '
     if ( TS_HA_PLANES(1, 1) ) chars = trim(chars) // '-A'
     if ( TS_HA_PLANES(2, 1) ) chars = trim(chars) // '+A'
     if ( TS_HA_PLANES(1, 2) ) chars = trim(chars) // '-B'
     if ( TS_HA_PLANES(2, 2) ) chars = trim(chars) // '+B'
     if ( TS_HA_PLANES(1, 3) ) chars = trim(chars) // '-C'
     if ( TS_HA_PLANES(2, 3) ) chars = trim(chars) // '+C'
     write(*,f10) 'Fixing Hartree potential at cell boundary', trim(chars)
    write(*,f8) 'Fix Hartree potential fraction', TS_HA_frac
    write(*,f7) 'Hartree potential offset', TS_HA_offset/eV, 'eV'

    if ( ts_method == TS_FULL ) then
       write(*,f10)'Solution method', 'Full inverse'
    else if ( ts_method == TS_BTD ) then
       write(*,f10)'Solution method', 'BTD'
       chars = fdf_get('TS.BTD.Pivot','atom+'//trim(Elecs(1)%name))
       write(*,f10)'BTD pivoting method', trim(chars)
       if ( BTD_method == 0 ) then
          chars = 'speed'
       else if  ( BTD_method == 1 ) then
          chars = 'memory'
       end if
       write(*,f10)'BTD creation algorithm', trim(chars)
       if ( IsVolt ) then
          select case ( ts_A_method )
          case ( TS_BTD_A_PROPAGATION )
             write(*,f10)'BTD spectral function algorithm','propagation'
          case ( TS_BTD_A_COLUMN )
             write(*,f10)'BTD spectral function algorithm','column'
          case default
             call die('Error in setup BTD. A calc')
          end select
       end if
#ifdef SIESTA__MUMPS
    else if ( ts_method == TS_MUMPS ) then
       write(*,f10)'Solution method', 'MUMPS'
       write(*,f5)'MUMPS extra memory', MUMPS_mem,'%'
       write(*,f5)'MUMPS blocking factor', MUMPS_block,''
       select case ( MUMPS_ordering ) 
       case ( 7 )
          write(*,f10)'MUMPS ordering', 'auto'
       case ( 6 )
          write(*,f10)'MUMPS ordering', 'QAMD'
       case ( 5 )
          write(*,f10)'MUMPS ordering', 'METIS'
       case ( 4 )
          write(*,f10)'MUMPS ordering', 'PORD'
       case ( 3 )
          write(*,f10)'MUMPS ordering', 'SCOTCH'
       case ( 2 )
          write(*,f10)'MUMPS ordering', 'AMF'
       case ( 0 )
          write(*,f10)'MUMPS ordering', 'AMD'
       end select
#endif
    end if
    write(*,f9) 'SCF DM tolerance',ts_Dtol
    write(*,f7) 'SCF Hamiltonian tolerance',ts_Htol/eV, 'eV'
    write(*,f1) 'SCF converge charge',ts_converge_dQ
    write(*,f9) 'SCF charge tolerance',ts_dQtol

    select case ( TS_scf_mode )
    case ( 0 )
       write(*,f10) 'Initialize DM by','diagon'
    case ( 1 )
       write(*,f10) 'Initialize DM by','transiesta'
    end select
    if ( IsVolt ) then
       if ( len_trim(Hartree_fname) > 0 ) then
          write(*,f10) 'User supplied Hartree potential', &
               trim(Hartree_fname)
       else
          if ( ts_tidx > 0 ) then
             write(*,f11) 'Hartree potential ramp across entire cell'
          else
             write(*,f11) 'Hartree potential will be placed in electrode box'
          end if
       end if

       chars = 'Non-equilibrium contour weight method'
       select case ( TS_W_METHOD )
       case ( TS_W_ORB_ORB )
          write(*,f10) trim(chars),'orb-orb'
       case ( TS_W_CORRELATED + TS_W_TR_ATOM_ATOM )
          write(*,f10) trim(chars),'Correlated Tr[atom]-Tr[atom]'
       case ( TS_W_TR_ATOM_ATOM )
          write(*,f10) trim(chars),'Uncorrelated Tr[atom]-Tr[atom]'
       case ( TS_W_CORRELATED + TS_W_TR_ATOM_ORB )
          write(*,f10) trim(chars),'Correlated Tr[atom]-orb'
       case ( TS_W_TR_ATOM_ORB )
          write(*,f10) trim(chars),'Uncorrelated Tr[atom]-orb'
       case ( TS_W_CORRELATED + TS_W_SUM_ATOM_ATOM )
          write(*,f10) trim(chars),'Correlated Sum[atom]-Sum[atom]'
       case ( TS_W_SUM_ATOM_ATOM )
          write(*,f10) trim(chars),'Uncorrelated Sum[atom]-Sum[atom]'
       case ( TS_W_CORRELATED + TS_W_SUM_ATOM_ORB )
          write(*,f10) trim(chars),'Correlated Sum[atom]-orb'
       case ( TS_W_SUM_ATOM_ORB )
          write(*,f10) trim(chars),'Uncorrelated Sum[atom]-orb'
       case ( TS_W_MEAN )
          write(*,f10) trim(chars),'Algebraic mean'
       case default
          ! This is an easy place for cathing mistakes
          call die('Error in code, weighting method unrecognized.')
       end select
       chars = 'Non-equilibrium contour weight k-method'
       select case ( TS_W_K_METHOD ) 
       case ( TS_W_K_CORRELATED )
          write(*,f10) trim(chars),'Correlated k-points'
       case ( TS_W_K_UNCORRELATED )
          write(*,f10) trim(chars),'Uncorrelated k-points'
       end select
    end if
    if ( .not. Calc_Forces ) then
       write(*,f11) '*** TranSiesta will NOT update forces ***'
    end if

    if ( TS_DQ_METHOD == 0 ) then
       write(*,f11)'Will not correct charge fluctuations'
    else if ( TS_DQ_METHOD == TS_DQ_METHOD_BUFFER ) then ! Correct in buffer
       if ( 0 < na_Buf ) then
          write(*,f10)'Charge correction','buffer'
       else
          call die('Charge correction can not happen in buffer as no buffer &
               &atoms exist.')
       end if
       write(*,f8)'Charge correction factor',TS_DQ_FACTOR
    else if ( TS_DQ_METHOD == TS_DQ_METHOD_FERMI ) then ! Correct fermi-lever
       write(*,f10)'Charge correction','Fermi-level'
       write(*,f8)'Charge correction dQ tolerance',TS_DQ_FERMI_TOLERANCE
       write(*,f7)'Fermi-level extrapolation eta value ',TS_DQ_FERMI_ETA/eV, 'eV'
       write(*,f8)'Charge correction factor',TS_DQ_FACTOR
       write(*,f7)'Max change in Fermi-level allowed', &
            TS_DQ_FERMI_MAX / eV,'eV'
    end if

    ! Print mixing options
    if ( associated(ts_scf_mixs, target=scf_mixs) ) then
       write(*,f11)'TS.SCF mixing options same as SCF'
    else
       call mixers_print('TS.SCF', ts_scf_mixs)
    end if

    write(*,f11)'          >> Electrodes << '
    ltmp = ts_tidx < 1 .and. IsVolt
    do i = 1 , size(Elecs)
       call print_settings(Elecs(i), 'ts', &
            box = ltmp)
    end do

    ! Print the contour information
    call print_contour_options( 'TS' , IsVolt )

    write(*,f11) repeat('*', 62)
    write(*,*)

  end subroutine print_ts_options


  subroutine print_ts_warnings( Gamma, cell, na_u, xa, Nmove )

    use parallel, only: IONode, Nodes
    use intrinsic_missing, only : VNORM, VEC_PROJ_SCA

    use m_os, only: file_exist

    use units, only: Kelvin, eV, Ang
    use siesta_options, only: FixSpin

    use m_ts_global_vars, only: TSmode, onlyS
    use m_ts_chem_pot, only : Name, Eq_segs
    use m_ts_electype, only : TotUsedAtoms, Name, Elec_frac

    use m_ts_method, only : a_isElec, a_isBuffer
    use m_ts_method, only : ts_A_method, TS_BTD_A_COLUMN

    use m_ts_contour_eq, only: N_Eq_E
    use m_ts_contour_neq, only: contour_neq_warnings

    use m_ts_hartree, only: TS_HA_frac

    ! Input variables
    logical, intent(in) :: Gamma
    real(dp), intent(in) :: cell(3,3)
    integer, intent(in) :: na_u
    real(dp), intent(in) :: xa(3,na_u)
    integer, intent(in) :: Nmove

    ! Local variables
    integer :: i, j, iEl, idx, idx1, idx2, itmp3(3)
    real(dp) :: rtmp, tmp3(3), tmp33(3,3), bdir(2)
    real(dp) :: p(3)
    logical :: err, warn, ltmp

    if ( .not. IONode ) return
    
    warn = .false.
    err = .false.

    write(*,'(3a)') repeat('*',24),' Begin: TS CHECKS AND WARNINGS ',repeat('*',24)
    if ( FixSpin ) then
      if ( TSmode ) then
        write(*,'(a)') 'Fixed spin for transiesta calculations is not implemented!'
        call die('Fixing spin is not possible in TranSiesta')
      end if
      if ( TS_HS_save ) then
        write(*,'(a)') 'Fixed spin aligns the Fermi-levels in the output TSHS to spin-UP!'
      end if
      if ( TS_DE_save ) then
        write(*,'(a)') 'Fixed spin aligns the Fermi-levels in the output TSDE to spin-UP!'
      end if
    end if

    if ( TS_HA_frac /= 1._dp ) then
      write(*,'(a)') 'Fraction of Hartree potential is NOT 1.'
      warn = .true.
    end if
    if ( TS_HA_frac < 0._dp .or. 1._dp < TS_HA_frac ) then
      write(*,'(a)') 'Fraction of Hartree potential is below 0.'
      write(*,'(a)') '  MUST be in range [0;1]'
      call die('Vha fraction erronously set.')
    end if
    
    ! Return if not a transiesta calculation
    if ( onlyS .or. .not. TSmode ) then
       write(*,'(3a,/)') repeat('*',24),' End: TS CHECKS AND WARNINGS ',repeat('*',26)
       return
    end if

    if ( ts_tidx < 1 .and. len_trim(Hartree_fname) == 0 .and. IsVolt ) then
       write(*,'(a)') 'Hartree potiental correction is the box solution &
            &which is not advised. Please supply your own Poisson solution.'
    end if

    if ( ts_A_method == TS_BTD_A_COLUMN ) then
       write(*,'(a)') 'Memory usage can be reduced by setting:'
       write(*,'(a)') '   TS.BTD.Spectral propagation'
    end if


    ! Check that all chemical potentials are really different
    ! their energy difference has to be below 0.1 meV and the
    ! temperature difference has to be below 10 K.
    rtmp = 10._dp * Kelvin
    do i = 1 , N_mu - 1
       do j = i + 1 , N_mu
          if ( abs(mus(i)%mu - mus(j)%mu) < 0.0001_dp * eV .and. &
               abs(mus(i)%kT - mus(j)%kT) <= rtmp ) then
             write(*,'(a)') 'Two chemical potentials: '//trim(name(mus(i)))//' and ' &
                  //trim(name(mus(j)))//' are the same, in bias calculations this &
                  &is not allowed.'
             err = .true.
          end if
       end do
    end do

    ! Check that all chemical potentials are in use
    if ( any(mus(:)%N_El == 0) ) then
       write(*,'(a)') 'A/Some chemical potential(s) have not been assigned any electrodes. &
            &All chemical potentials *MUST* be assigned an electrode'
       err = .true.
    end if

    ! check that all have at least 2 contour points on the equilibrium contour
    ! the 3rd is the fictive pole segment
    if ( .not. all(Eq_segs(mus(:)) > 0) ) then
       write(*,'(a)') 'All chemical potentials does not have at least &
            &1 equilibrium contours'
       err = .true.
    end if
    if ( .not. all(Eq_segs(mus(:)) /= 2) ) then
       write(*,'(a)') 'No chemical potential can have only two &
            &equilibrium contours'
       write(*,'(a)')'Either of these:'
       write(*,'(a)')'   1. continued fraction'
       write(*,'(a)')'  or'
       write(*,'(a)')'   1. Circle contour'
       write(*,'(a)')'   2. Tail contour'
       write(*,'(a)')'   3. Residuals (poles)'
       err = .true.
    end if

    ! we need to check that they indeed do not overlap
    do i = 1 , N_Elec
       idx1 = Elecs(i)%idx_a
       idx2 = idx1 + TotUsedAtoms(Elecs(i)) - 1
       ! we need to check every electrode,
       ! specifically because if one of the electrodes is fully located
       ! inside the other and we check the "small" one 
       ltmp = .false.
       do j = 1 , N_Elec
          if ( i == j ) cycle
          idx = Elecs(j)%idx_a
          if ( (idx <= idx1 .and. &
               idx1 < idx + TotUsedAtoms(Elecs(j))) ) then
             ltmp = .true.
          end if
          if ( (idx <= idx2 .and. &
               idx2 < idx + TotUsedAtoms(Elecs(j))) ) then
             ltmp = .true.
          end if
          if ( ltmp ) then
             write(*,'(a)') 'Electrode: '//trim(Elecs(i)%name)
             write(*,'(a,i0,a,i0)') 'Positions: ',idx1,' -- ',idx2 
             idx1 = Elecs(j)%idx_a
             idx2 = idx1 + TotUsedAtoms(Elecs(j)) - 1
             write(*,'(a)') 'Electrode: '//trim(Elecs(j)%name)
             write(*,'(a,i0,a,i0)') 'Positions: ',idx1,' -- ',idx2 
             write(*,'(a)') 'Overlapping electrodes is not physical, please correct.'
             err = .true.
          end if
       end do

       ! Warn if using non-bulk electrodes and one does not update everything!
       if ( .not. Elecs(i)%Bulk ) then
          select case ( Elecs(i)%DM_update )
          case ( 0 )
             write(*,'(a)') 'Electrode: '//trim(Elecs(i)%name)
             write(*,'(a)') '  is using a non-bulk electrode Hamiltonian and does not'
             write(*,'(a)') '  update the electrode region or the cross-term DM!'
             write(*,'(a)') '  Please consider setting "DM-update all"'
          case ( 1 )
             write(*,'(a)') 'Electrode: '//trim(Elecs(i)%name)
             write(*,'(a)') '  is using a non-bulk electrode Hamiltonian and does not'
             write(*,'(a)') '  update the electrode region DM!'
             write(*,'(a)') '  Please consider setting "DM-update all"'
          end select
       end if

    end do
    
    ! CHECK THIS (we could allow it by only checking the difference...)
    if (  maxval(mus(:)%mu) - minval(mus(:)%mu) - abs(Volt) > 1.e-4_dp * eV ) then
       write(*,'(a)') 'Chemical potentials [eV]:'
       do i = 1 , N_Elec
          write(*,'(a,f10.5,a)') trim(Elecs(i)%name)//' at ',Elecs(i)%mu%mu/eV,' eV'
       end do
       write(*,'(a)') 'The difference must satisfy: "max(ChemPots)-min(ChemPots) - abs(Volt) < 1e-4 eV"'
       write(*,'(a,f10.5,a)') 'max(ChemPots) at ', maxval(mus(:)%mu)/eV,' eV'
       write(*,'(a,f10.5,a)') 'min(ChemPots) at ', minval(mus(:)%mu)/eV,' eV'
       write(*,'(a,f10.5,a)') '|V| at ', abs(Volt)/eV,' eV'
       write(*,'(a)') 'Chemical potentials are not consistent with the bias applied.'
       err = .true.
    end if

    ! Check that the bias does not introduce a gating
    if ( any(abs(mus(:)%mu) - abs(Volt) > 1.e-9_dp ) ) then
       write(*,'(a)') 'Chemical potentials must lie in the range [-V;V] with the maximum &
            &difference being V'
       write(*,'(a)') 'Chemical potentials must not introduce consistent Ef shift to the system.'
       err = .true.
    end if


    ! Check that we can actually start directly in transiesta
    if ( TS_scf_mode == 1 ) then ! TS-start
       if ( .not. all(Elecs(:)%DM_update >= 1) ) then
          write(*,'(a)')'WARNING: Responsibility is now on your side'
          write(*,'(a)')'WARNING: Requesting immediate start, yet we &
               &do not update cross-terms.'
          warn = .true.
       end if
    end if
    
    ! Calculate the number of optimal contour points
    i = mod(N_Eq_E(), Nodes) ! get remaining part of equilibrium contour
    if ( i /= 0 ) then
       i = Nodes - i
       write(*,'(a)')'Without loosing performance you can increase &
            &the equilibrium integration precision.'
       write(*,'(a,i0,a)')'You can add ',i,' more energy points in the &
            &equilibrium contours, for FREE!'
       if ( i/N_mu > 0 ) then
          write(*,'(a,i0,a)')'This is ',i/N_mu, &
               ' more energy points per chemical potential.'
       end if
    end if
    
    call contour_nEq_warnings()
    
    if ( .not. Calc_Forces ) then
       write(*,f11) '***       TranSiesta will NOT update forces       ***'
       write(*,f11) '*** ALL FORCES AFTER TRANSIESTA HAS RUN ARE WRONG ***'
       if ( Nmove > 0 ) then
          write(*,'(a)')'Relaxation with TranSiesta *REQUIRES* an update of &
               &the energy density matrix. Will continue at your request.'
          err = .true.
       end if
    end if

    if ( Nmove > 0 .and. .not. all(Elecs(:)%DM_update > 0) ) then
       write(*,'(a)') 'TranSiesta relaxation is only allowed if you also &
            &update, at least, the cross terms, please set: &
            &TS.Elecs.DM.Update [cross-terms|all]'
       err = .true.
    end if

    ! A transiesta calculation requires that all atoms
    ! are within the unit-cell.
    ! However, for N_Elec > 2 calculations this
    ! is not a requirement as the potential is placed
    ! irrespective of the actual box.
    
    ltmp = .false.
    call reclat(cell,tmp33,0)

    do i = 1 , na_u

       ! If we use a boxed or user-supplied potential
       ! profile, we need not check the coordinates.
       ! It is the users responsibility for user-defined
       ! potential grids.
       ! The boxed one takes care of grid-points in the periodic picture
       ! And if the grid-plane does not cut the actual grid, it will
       ! error out.
       if ( ts_tidx < 1 ) cycle

       ! Buffer atoms can be where-ever
       if ( a_isBuffer(i) ) cycle
       
       ! Check the index
       do j = 1 , 3
          itmp3(j) = floor( dot_product(xa(:,i),tmp33(:,j)) )
       end do

       ! Only check the transport direction
       ! Note that we have projected onto the unit-cell
       ! vector, hence tidx and not tdir
       if ( itmp3(ts_tidx) /= 0 ) then
          write(*,'(i0,''('',i0,''='',i0,'')'',tr1)',advance='no')&
               i,ts_tidx,itmp3(ts_tidx)
          ltmp = .true.
       end if
       
    end do
    if ( ltmp ) then
       write(*,'(/a)')'atom(cell-dir=neighbour-cell)'
       write(*,'(a)') '*** Device atomic coordinates are not inside unit-cell.'
       write(*,'(a)') '*** This is a requirement for bias calculations'
       write(*,'(a)') '    as the Poisson equation cannot be correctly handled'
       write(*,'(a)') '    due to inconsistencies with the grid and atomic coordinates'
       if ( IsVolt ) then
          write(*,'(a)') 'Please move device atoms inside the device region in &
               &the transport direction.'
          err = .true.
       else
          write(*,'(a)') '*** Will continue, but will die when running V /= 0 ***'
          warn = .true.
       end if
    end if

    if ( ts_tidx < 1 ) then

       write(*,'(a)') '*** TranSiesta semi-infinite directions are individual ***'
       write(*,'(a)') '*** It is heavily advised to have any electrodes with no &
            &periodicity'
       write(*,'(a)') '    in the transverse directions be located as far from any &
            &cell-boundaries'
       write(*,'(a)') '    as possible. This has to do with the electrostatic potential &
            &correction. ***'
       if ( IsVolt ) then
          write(*,'(a)') '*** Please ensure electrode unit-cells are as confined as possible.'
          write(*,'(a)') '    I.e. do not add superfluous vacuum if not needed in the &
               &electrode calculation.'
          write(*,'(a)') '    The initial guess for the potential profile is heavily influenced'
          write(*,'(a)') '    by the electrode unit-cell sizes! ***'
       end if

    else

       ! The transport direction is well-defined
       ! and the Hartree potential is fixed at the bottom of the
       ! unit-cell of the A[ts_tidx] direction.
       ! We will let the user know if any atoms co-incide with
       ! the plane as that might hurt convergence a little.

       if ( N_Elec <= 2 ) then
          
          ! Get the electrode fraction of the position
          if ( N_Elec == 1 ) then
             
             call Elec_frac(Elecs(1),cell,na_u,xa,ts_tidx, fmin = bdir(1))
             call Elec_frac(Elecs(1),cell,na_u,xa,ts_tidx, fmax = bdir(2))
             iEl = 1
             if ( bdir(1) < bdir(2) ) then
                i = 1
             else
                i = 2
             end if
             
          else
             
             ! Get the electrode fraction of the position
             call Elec_frac(Elecs(1),cell,na_u,xa,ts_tidx, fmin = bdir(1))
             call Elec_frac(Elecs(2),cell,na_u,xa,ts_tidx, fmin = bdir(2))

             ! Determine the electrode closest to the
             ! lower boundary
             if ( bdir(1) < bdir(2) ) then
                ! The first electrode is closest
                i = 1
                iEl = 1
                call Elec_frac(Elecs(2),cell,na_u,xa,ts_tidx, fmax = bdir(2))
             else
                ! The second electrode is closest
                i = 2
                iEl = 2
                call Elec_frac(Elecs(1),cell,na_u,xa,ts_tidx, fmax = bdir(1))
             end if
             
          end if

          ! Get the fraction
          tmp3 = cell(:,ts_tidx) / VNORM(cell(:,ts_tidx))

          ! Check lower limit
          if ( bdir(i) < 0._dp ) then

             ! Tell the user to shift the entire structure
             ! by the offset to origo + 1/2 a bond-length

             ! Origo offset:
             p = - cell(:,ts_tidx) * bdir(i)
             p = p + tmp3 * Elecs(iEl)%dINF_layer * 0.5_dp
             ! Bond-length
             write(*,'(a,f9.5,a)') 'Electrode: '//trim(Elecs(iEl)%name)//' lies &
                  &outside the unit-cell.'
             write(*,'(a)')'Please shift the entire structure using the &
                  &following recipe:'
             write(*,'(a)') 'If you already have AtomicCoordinatesFormat, add these'
             write(*,'(tr1,a)') 'AtomicCoordinatesFormat Ang'
             write(*,'(tr1,a)') '%block AtomicCoordinatesOrigin'
             write(*,'(tr1,3(tr2,f12.4))') p / Ang
             write(*,'(tr1,a)') '%endblock AtomicCoordinatesOrigin'
             err = .true.

          end if

          ! Check upper limit
          if ( i == 1 ) then
             i = 2
          else
             i = 1
          end if
          if ( N_Elec == 2 ) iEl = i
          if ( 1._dp < bdir(i) ) then

             bdir(i) = bdir(i) - 1._dp

             ! Tell the user to shift the entire structure
             ! by the offset to origo + 1/2 a bond-length

             ! Origo offset:
             p = - cell(:,ts_tidx) * bdir(i)
             p = p - tmp3 * Elecs(iEl)%dINF_layer * 0.5_dp
             ! Bond-length
             write(*,'(a,f9.5,a)') 'Electrode: '//trim(Elecs(iEl)%name)//' lies &
                  &outside the unit-cell.'
             write(*,'(a)')'Please shift the entire structure using the &
                  &following recipe:'
             write(*,'(a)') 'If you already have AtomicCoordinatesFormat, add these'
             write(*,'(tr1,a)') 'AtomicCoordinatesFormat Ang'
             write(*,'(tr1,a)') '%block AtomicCoordinatesOrigin'
             write(*,'(tr1,3(tr2,f12.4))') p / Ang
             write(*,'(tr1,a)') '%endblock AtomicCoordinatesOrigin'
             err = .true.

          end if

       else

          call die("ts_options: ts_tidx < 0 with N_Elec > 2")
          
       end if

    end if
    
    ! If the user has requested to initialize using transiesta
    ! and the user does not utilize the bulk DM, they should be
    ! warned
    if ( TS_scf_mode == 1 .and. any(Elecs(:)%DM_init == 0) ) then
       write(*,'(a)') 'You are not initializing the electrode DM/EDM. &
            &This may result in very wrong electrostatic potentials close to &
            &the electrode/device boundary region.'
       if ( IsVolt ) then
         write(*,'(a)') '    This warning is only applicable for V == 0 calculations!'
         write(*,'(a)') '    I give this warning because it is not clear how your V = 0 calcualtion was done.'
       end if
       warn = .true.
    end if

    ! warn the user about suspicous work regarding the electrodes
    do i = 1 , N_Elec

       idx1 = Elecs(i)%idx_a
       idx2 = idx1 + TotUsedAtoms(Elecs(i)) - 1

       if ( .not. Elecs(i)%Bulk ) then
          write(*,'(a)') 'Electrode '//trim(Elecs(i)%name)//' will &
               &not use bulk Hamiltonian.'
          warn = .true.
       end if

       if ( Elecs(i)%DM_update == 0 ) then
          write(*,'(a)') 'Electrode '//trim(Elecs(i)%name)//' will &
               &not update cross-terms or local region.'
          warn = .true.
       end if
       
       if ( .not. Elecs(i)%Bulk .and. Elecs(i)%DM_update /= 2 ) then
          write(*,'(3a)') 'Electrode ',trim(Elecs(i)%name), &
               ' has non-bulk Hamiltonian and does not update all'
          warn = .true.
       end if

       if ( .not. Elecs(i)%kcell_check ) then
          write(*,'(a)') 'Electrode '//trim(Elecs(i)%name)//' will &
               &not check the k-grid sampling vs. system k-grid &
               &sampling. Ensure appropriate sampling.'
       end if
       if ( Elecs(i)%V_frac_CT /= 0._dp ) then
          write(*,'(a)') 'Electrode '//trim(Elecs(i)%name)//' will &
               &shift coupling Hamiltonian with a shift in energy &
               &corresponding to the applied bias. Be careful here.'
          warn = .true.
       end if
       if ( Elecs(i)%delta_Ef /= 0._dp ) then
          write(*,'(a)') 'Electrode '//trim(Elecs(i)%name)//' will &
               &shift the electronic structure manually. Be careful here.'
          warn = .true.
       end if

       if ( Elecs(i)%repeat .and. Elecs(i)%bloch%size() > 1 ) then
         write(*,'(a)') 'Electrode '//trim(Elecs(i)%name)//' is &
             &using Bloch unfolding using the repeat scheme! &
             &Please use the tiling scheme (it is orders of magnitudes faster!).'
         warn = .true.
       end if

       ! if any buffer atoms exist, we should suggest to the user
       ! to use TS.Elec.<elec> [DM-update cross-terms|all]
       ! in case any buffer atoms are too close
       ltmp = .false.
       do j = 1 , na_u
          ! skip non-buffer atoms
          if ( .not. a_isBuffer(j) ) cycle
          do idx = idx1 , idx2
             ! Proximity of 4 Ang enables this check
             ltmp = VNORM(xa(:,idx)-xa(:,j)) < 4._dp * Ang
             if ( ltmp ) exit
          end do
          if ( ltmp ) exit
       end do
       if ( ltmp .and. Elecs(i)%DM_update == 0 ) then
          ! some buffer atoms are close to this electrode
          ! Advice to use dm_update
          write(*,'(a,/,a)') 'Electrode '//trim(Elecs(i)%name)//' is &
               &likely terminated by buffer atoms. It is HIGHLY recommended to add this:', &
               '  TS.Elec.'//trim(Elecs(i)%name)//' DM-update [cross-terms|all]'
          warn = .true.
       end if

       ! In case DM_bulk is requested we assert that the file exists
       ltmp = file_exist(Elecs(i)%DEfile)
       if ( Elecs(i)%DM_init > 0 .and. .not. ltmp ) then
          write(*,'(a,/,a)') 'Electrode '//trim(Elecs(i)%name)//' TSDE &
               &file cannot be located in: '//trim(Elecs(i)%DEfile)//'.', &
               '  Please add TS.DE.Save T to the electrode calculation or &
               &specify the exact file position using ''TSDE-file'' in the&
               & TS.Elec block.'
          err = .true.
       end if

    end do

    if ( N_Elec /= 2 .and. any(Elecs(:)%DM_update == 0) ) then
       write(*,'(a,/,a)') 'Consider updating more elements when doing &
            &N-electrode calculations. The charge conservation typically &
            &increases.','  TS.Elecs.DM.Update [cross-terms|all]'
       warn = .true.
    end if

    ! Check that the pivoting table is unique
    do iEl = 1, N_Elec
       if ( sum(Elecs(iEl)%pvt) /= 6 .or. count(Elecs(iEl)%pvt==2) /= 1 ) then
          write(*,'(a,/,a)') 'The pivoting table for the electrode unit-cell, &
               &onto the simulation unit-cell is not unique: '//trim(Elecs(iEl)%name), &
               '  Please check your electrode and device cell parameters!'
          write(*,'(a)') '  Combining this with electric fields or dipole-corrections is NOT advised!'
          warn = .true.
       end if
    end do

    write(*,'(3a,/)') repeat('*',24),' End: TS CHECKS AND WARNINGS ',repeat('*',26)

    if ( warn ) then
       ! Print BIG warning sign

       write(*,'(tr18,a)') repeat('*',40)
       write(*,'(tr19,a)') 'TRANSIESTA REPORTED IMPORTANT WARNINGS'
       write(*,'(tr18,a)') repeat('*',40)

    end if

    if ( err ) then
       write(*,'(/tr18,a)') repeat('*',30)
       write(*,'(tr19,a)') 'TRANSIESTA REPORTED ERRORS'
       write(*,'(tr18,a)') repeat('*',30)

       call die('One or more errors have occured doing TranSiesta &
            &initialization, check the output')
    end if
    
  end subroutine print_ts_warnings


  subroutine print_ts_blocks( na_u, xa ) 

    use parallel, only : IONode
    use files, only : slabel

    use m_ts_global_vars, only: TSmode, onlyS
    use m_ts_chem_pot, only: print_mus_block
    use m_ts_contour, only: print_contour_block, io_contour


    ! Input variables
    integer, intent(in) :: na_u
    real(dp), intent(in) :: xa(3,na_u)

    ! Local variables
    integer :: i

    if ( .not. IONode ) return

    if ( onlyS .or. .not. TSmode ) return

    write(*,'(/,a,/)') '>>> TranSiesta block information for FDF-file START <<<'
    
    call print_mus_block( 'TS' , N_mu , mus)
    
    call print_contour_block( 'TS' , IsVolt )
    
    write(*,'(/,a,/)') '>>> TranSiesta block information for FDF-file END <<<'

    ! write out the contour
    call io_contour(IsVolt, mus, slabel)

  end subroutine print_ts_blocks

  subroutine val_swap(v1,v2)
    real(dp), intent(inout) :: v1, v2
    real(dp) :: tmp
    tmp = v1
    v1  = v2
    v2  = tmp
  end subroutine val_swap
  
end module m_ts_options
