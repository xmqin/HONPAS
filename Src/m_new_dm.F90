!>  Prepares a starting density matrix for a new geometry iteration

!>     This DM can be:
!>         1. Synthesized directly from atomic occupations (not idempotent)
!>         2. Read from file
!>         3. Re-used (with possible extrapolation) from previous geometry step(s).

!>     In cases 2 and 3, the structure of the read or extrapolated DM 
!>     is automatically adapted to the current sparsity pattern.

!>     Special cases:
!>            Harris: The matrix is always initialized
!>            Force calculation: The DM should be written to disk
!>                               at the time of the "no displacement"
!>                               calculation and read from file at
!>                               every subsequent step.

!>              In all other cases (including "server operation"), the
!>              default is to allow DM re-use (with possible extrapolation)
!>              from previous geometry steps.
!>              The fdf variables 'DM.AllowReuse' and 'DM.AllowExtrapolation'
!>              (mapped to 'allow_dm_reuse' and 'allow_dm_extrapolation', and
!>              both 'true' by default) can be used to change this behavior.
!>
!>              There is no re-use of the DM for "Forces", and "Phonon"
!>              dynamics types (i.e., the DM is re-initialized)
!>
!>              For "CG" calculations, the default is not to extrapolate the
!>              DM (unless requested by setting 'DM.AllowExtrapolation' to
!>              "true"). The previous step's DM is reused.

module m_new_dm

  use sys, only: die
  use parallel, only: IONode
  use precision, only: dp
  use alloc, only: re_alloc, de_alloc

  implicit none

  private
  
  public :: new_dm
  public :: get_allowed_history_depth

contains

  subroutine new_DM(SC_changed, DM_history, DM_2D, EDM_2D)

    use siesta_options
    use siesta_geom,      only: ucell, xa, na_u, isc_off
    use siesta_geom,      only: nsc, nsc_old

    use sparse_matrices,  only: sparse_pattern, block_dist
    use atomlist,         only: Datm, iaorb, lasto, no_u
    use m_steps,          only: istp
    use m_spin,   only: spin
    use m_handle_sparse, only : bulk_expand

    use class_dSpData2D
    use class_Sparsity
    use class_Pair_Geometry_dSpData2D
    use class_Fstack_Pair_Geometry_dSpData2D

    use fdf, only: fdf_get, fdf_defined
    use units, only : eV
    use m_ts_global_vars,only: TSrun
    use m_ts_electype, only : copy_DM
    use m_ts_options, only : TS_analyze
    use m_ts_options, only : N_Elec, Elecs
    use m_ts_method
    use m_energies, only: Ef

    logical, intent(in) :: SC_changed ! Has auxiliary supercell changed?
    type(Fstack_Pair_Geometry_dSpData2D), intent(inout) :: DM_history
    type(dSpData2D), intent(inout) :: DM_2D, EDM_2D

    ! Local variables
    logical :: DM_init ! Initialize density matrix from file or from atomic density
    logical :: read_DM
    integer :: DM_in_history, n_depth

    ! Variable to signal how the initialization was performed.
    ! The following table lists the (current) possibilities:
    !  0: DM is filled from atomic information (no information present)
    !  1: .DM file is read, Siesta continuation
    !  2: An extrapolation of previous geometries.
    !     The DM's from the previous coordinates are kept in memory
    !  3: .TSDE file is read, TranSiesta continuation
    ! If the value is negative it means that one cannot mix in the
    ! first step, unless SCF.Mix.First.Force is true!
    integer :: init_method

    real(dp), pointer :: DM(:,:), EDM(:,:)
    integer :: iElec
    integer :: na_a, na_dev
    integer, allocatable :: allowed_a(:)
    logical :: set_Ef
    real(dp) :: old_Ef, diff_Ef

    if ( IONode ) then
       write(*,"(a,i5)") "new_DM -- step: ", istp
    end if

    ! In principle we allow the re-use of the DM (i.e, we do not initialize it)
    ! Initialization is either:
    !  1) read of DM
    !  2) atomic initialization of DM (possibly user-defined spin-configuration)
    DM_init = .false.
    
    ! As per defaults
    read_DM = UseSaveDM

    ! Except if there are explicit instructions
    if ( .not. allow_DM_reuse ) then
       DM_init = .true.
       
       ! Do not allow to re-use the dm
       if (IONode) then
          write(*,"(a)") "DM re-use not allowed. Resetting DM at every geometry step"
          if ( read_DM ) then
             write(*,"(a)") "DM.UseSaveDM  overridden!!"
          end if
       end if
       read_DM = .false.
       
    end if

    ! or using Harris...
    if ( harrisfun ) DM_init = .true.

    ! or we are in the first step, or performing force-constant calculations
    DM_in_history = n_items(DM_history)

    ! Force initialization when
    !  1) no DM are accummulated in the history
    !  2) when force-constants runs are made
    !     In this case reading the DM is equivalent to
    !     reading the un-displaced DM which should be
    !     closer to SCF than the previously displaced one.
    if ( DM_in_history == 0 ) then
       DM_init = .true.
    else if ( UseSaveDM ) then
       if ( idyn == 6 .and. writeDM )  & ! Force Constants
            DM_init = .true.
       if ( idyn == 7 .or. idyn == 9 ) & ! Phonon series (writedm??)
            DM_init = .true.
    end if

#ifdef TO_BE_REMOVED
    ! ... or if the auxiliary cell has changed
    ! (in this case we have to avoid reading back saved copy from file)
    if ( SC_changed ) then
       
       if ( initDMaux ) then
          ! Reset DM history, thus we MUST
          ! also re-initialize DM
          DM_init = .true.
          read_DM = .false.
          if ( IONode ) then
             write(*,"(a)") "DM history reset as auxiliary supercell changed."
          end if
          call get_allowed_history_depth(n_depth)
          call new(DM_history, n_depth, "(reset DM history stack)")
       else
          if ( IONode ) then
             write(*,"(a)") "** Warning: DM history NOT reset upon supercell change"
             write(*,"(a)") "** Warning: since 'ReinitialiseDM' is set to .false."
          end if
       end if
     end if
#endif

    if ( DM_init ) then
       
       if ( IOnode ) then
          write(*,"(a)") "Initializing Density Matrix..."
       end if

       call init_DM(spin, na_u, no_u, nsc, lasto, iaorb, &
            Datm, &
            block_dist, sparse_pattern, &
            DM_2D, EDM_2D, &
            read_DM, &
            init_anti_ferro, &
            init_method)

    else

       if (IOnode) then
          write(*,"(a)") "Re-using DM from previous geometries..."
       endif

       ! Extrapolation or simple re-structuring
       if ( IONode ) then
          write(*,'(a,i0)') "Number of DMs in history: ", DM_in_history
       end if
       call extrapolate_DM_with_coords(DM_history, na_u, xa(:,1:na_u), &
            sparse_pattern, DM_2D)
       if ( IONode ) then
          write(*,'(a)') "New DM after history re-use:"
          call print_type(DM_2D)
       end if

       ! A negative value specifies that one cannot mix in the initial step
       init_method = -2

    end if

    if ( init_method == 0 ) then

       ! In case we have initialized from atomic fillings
       ! we allow the user to supply different files for
       ! starting the calculation in a state of bulk 
       ! calculations
       call bulk_expand(na_u,xa,lasto,ucell,nsc,isc_off,DM_2D)

       ! Note that normalization is performed later

    end if

    ! In case we will immediately start a transiesta run
    if ( TSrun .and. DM_Init .and. (.not. TS_analyze ) ) then
       ! In transiesta we can always initialize
       ! the density matrix with the bulk-values
       ! so as to "fix" the density in the leads

       ! If the Fermi-level has not been
       ! set, we initialize it to the mean of the
       ! electrode chemical potentials
       if ( any(Elecs(:)%DM_init > 0) ) then

          set_Ef = abs(Ef) < 0.00001_dp .and. &
               (.not. fdf_defined('TS.Fermi.Initial') ) 
          if ( IONode ) then
             write(*,'(/,a)') 'transiesta: Will read in bulk &
                  &density matrices for electrodes'
             if ( set_Ef ) then
                write(*,'(a)') &
                     'transiesta: Will average Fermi-levels of electrodes'
             end if
          end if

          if ( init_method == 0 ) then
            ! We are starting from atomic-filled orbitals
            ! We are now allowed to overwrite everything (even buffer!)
            na_a = na_u
            allocate(allowed_a(na_a))
            do iElec = 1 , na_a
              allowed_a(iElec) = iElec
            end do
          else
            ! We will only overwrite elements in the electrodes
            ! Not in buffer
            na_a = 0
            do iElec = 1 , na_u
              if ( a_isElec(iElec) ) then
                if ( Elecs(atom_type(iElec))%DM_init > 0 ) then
                  na_a = na_a + 1
                end if
              end if
            end do
            allocate(allowed_a(na_a))
            na_a = 0 
            do iElec = 1 , na_u
              if ( a_isElec(iElec) ) then
                if ( Elecs(atom_type(iElec))%DM_init > 0 ) then
                  na_a = na_a + 1
                  allowed_a(na_a) = iElec
                end if
              end if
            end do
          end if

          do iElec = 1 , N_Elec

            ! We shift the mean by one fraction of the electrode
            if ( set_Ef ) then
              Ef = Ef + Elecs(iElec)%Ef / N_Elec
            end if

            if ( Elecs(iElec)%DM_init == 0 ) cycle
            
            if ( IONode ) then
              write(*,'(/,2a)') 'transiesta: Reading in electrode TSDE for ', &
                  trim(Elecs(iElec)%Name)
            end if
            
            ! Copy over the DM in the lead
            ! Notice that the EDM matrix that is copied over
            ! will be equivalent at Ef == 0
            call copy_DM(Elecs(iElec),na_u,xa,lasto,nsc,isc_off, &
                ucell, DM_2D, EDM_2D, na_a, allowed_a)
            
          end do

          ! Clean-up
          deallocate(allowed_a)

       end if
        
       ! Initialize the Fermi-level
       diff_Ef = Ef

       if ( fdf_defined('TS.Fermi.Initial') ) then
          ! Write out some information regarding
          ! how the Ef is set

          old_Ef = Ef
          Ef = fdf_get('TS.Fermi.Initial',Ef,'Ry')
          if ( abs(init_method) == 2 ) then
             ! As the fermi-level has been read in from a previous
             ! calculation (TSDE), the EDM should only be shifted by the difference
             diff_Ef = Ef - old_Ef
          else
             ! Setting diff_Ef has no meaning
             diff_Ef = Ef
          end if

          if ( IONode ) then
             write(*,*) ! new-line
             if ( abs(init_method) < 2 ) then
                write(*,'(a,f9.5,a)')'transiesta: Setting the Fermi-level to: ', &
                     Ef / eV,' eV'
             else if ( abs(init_method) == 2 ) then
                write(*,'(a,2(f10.6,a))')'transiesta: Changing Fermi-level from -> to: ', &
                     old_Ef / eV,' -> ',Ef / eV, ' eV'
             end if
          end if

       end if

       ! The electrode EDM is aligned at Ef == 0
       ! We need to align the energy matrix
       iElec = nnzs(sparse_pattern) * spin%EDM
       DM => val(DM_2D)
       EDM => val(EDM_2D)
       call daxpy(iElec,diff_Ef,DM(1,1),1,EDM(1,1),1)

    end if


    ! Determine how the mixing of the first step should be performed
    !
    ! The idea is that one will allow mixing on the first SCF step
    ! only for atomically filled orbitals. This will in most cases
    ! be a good initial choice, while in certain systems (spin-orbit)
    ! it may not be so good.
    ! Subsequently for any MD steps beyond the initial step, we will not
    ! allow mixing the first SCF step if the DM has been extrapolated
    ! or re-used.
    ! If allowing mixing on the first SCF step it has been known that
    ! it may lead to "wrong" convergence due to memory effects in cases of few
    ! iterations. For instance when performing FC runs (externally driven)
    ! this would lead to non-degenerate transverse eigenmodes for simple molecules
    if ( init_method == 0 ) then ! atomicly filled data
      ! allow mix_scf_first
    else if ( mix_scf_first_force ) then
      ! user requested a mix in the first step
      mix_scf_first = .true.
      if ( IONode .and. init_method < 0 ) then
        ! Warn the user about this
        write(*,"(a)") "Mixing first scf iteration is *forced* although the sparsity pattern has changed..."
        write(*,"(a)") "Mixing non-compatible matrices might result in problems."
        write(*,"(a)") "Please use a linear mixer for the first few steps to remove history effects."
      end if
    else if ( mix_scf_first .and. init_method < 0 ) then
       ! Do not allow mixing first SCF step since we are reusing a DM
       ! from another geometry.
       mix_scf_first = .false.
    end if
    
  end subroutine new_DM


  ! Tries to initialize the DM by:
  !  1) reading a DM/TSDE(transiesta) file
  !  2) fall-back to atomic initialization using
  !     a possibly user-defined spin-configuration.
  subroutine init_DM(spin, na_u, no_u, nsc, lasto, iaorb, &
       DM_atom, &
       dit, sp, DM_2D, EDM_2D, &
       read_DM, &
       anti_ferro, init_method)
    
    use t_spin, only : tSpin

    use class_OrbitalDistribution
    use class_Sparsity
    use class_dSpData2D
    use class_Fstack_dData1D, only: reset
    use m_ts_global_vars,only: ts_method_init, TSinit, TSrun, TSmode
    use m_ts_options,   only : TS_scf_mode, ts_hist_keep
    use m_ts_options,   only : val_swap, ts_scf_mixs
    use m_ts_options,   only : ts_Dtol, ts_Htol
    use m_ts_options,   only : IsVolt, TS_start_bias
    use siesta_options, only : dDtol, dHtol

    use m_mixing, only: mixers_history_init
    use m_mixing_scf, only: scf_mixs, scf_mix

    ! ********* INPUT ***************************************************
    ! type(tSpin) spin              : spin configuration for this system
    ! integer na_u                  : number of atoms in unit-cell
    ! integer no_u                  : number of orbitals in unit-cell
    ! integer nsc(3)                : number of supercells along each direction
    ! integer lasto(0:na_u)         : last orbital of each atom
    ! integer iaorb(no_u)           : the atomic index of the corresponding orbital
    ! real(dp) DM_atom(no_u)        : atomic density based on atomic configuration
    ! type(OrbitalDistribution) dit : the distribution used for the orbitals
    ! type(Sparsity) sp             : sparsity pattern of DM
    ! type(dSpData2D) DM_2D         : the density matrix 
    ! logical read_DM               : true if init_DM should try and read the DM from disk
    ! logical anti_ferro            : whether the DM should be anti-ferro magnetic or not
    ! integer init_method           : returns method it has read the data with
    !                                   0 == atomic filling (possibly user-defined
    !                                   1 == .DM read
    !                                   2 == .TSDE read
    ! *******************************************************************

    ! The spin-configuration that is used to determine the spin-order.
    type(tSpin), intent(in) :: spin
    ! Number of atoms in the unit-cell
    integer, intent(in) :: na_u
    ! Number of orbitals in the unit-cell
    integer, intent(in) :: no_u
    ! Number of supercells along each direction
    integer, intent(in) :: nsc(3)
    ! The last orbital on each atom
    integer, intent(in) :: lasto(0:na_u)
    ! The atom containing orbital "io"
    integer, intent(in) :: iaorb(no_u)
    ! DM filled with atomic charges based on atom and orbitals
    real(dp), intent(in) :: DM_atom(no_u)
    ! Parallel distribution of DM
    type(OrbitalDistribution), intent(in) :: dit
    ! Sparse pattern for DM
    type(Sparsity), intent(inout) :: sp
    ! The DM, these will be initialiazed upon return by
    ! the atomic filling
    type(dSpData2D), intent(inout) :: DM_2D
    ! The EDM, these will be initialiazed upon return by
    ! the atomic filling
    type(dSpData2D), intent(inout) :: EDM_2D
    ! Whether we should try and read from DM/TSDE file
    logical, intent(in) :: read_DM
    ! Setting for the magnetic setup of DM
    ! If .true. the spin-configuration will be anti-ferro-magnetic
    logical, intent(in) :: anti_ferro
    integer, intent(out) :: init_method

    ! *** Local variables
    integer :: i

    ! Signal that we are using atomic fillings
    init_method = 0

    ! Try to read DM from disk if wanted
    if ( read_DM ) then
       
       ! Try and read the DM from the files
       call init_DM_file(spin, no_u, nsc, &
            dit, sp, DM_2D, EDM_2D, &
            init_method)

    end if

    if ( init_method == 0 ) then
       ! The reading of DM did not succeed...
       ! We will initialize with atomic charges
       ! This is guarenteed to succeed.
       call init_DM_atomic(spin, na_u, no_u, lasto, iaorb, &
            DM_atom, &
            dit, sp, DM_2D, &
            anti_ferro)

       ! The init_method now signals that it is atomic filling
       
    end if

    
    if ( TS_scf_mode == 1 .and. TSinit ) then

      ! if the user requests to start the transiesta SCF immediately.
      call ts_method_init( .true. )

    else 

      ! Print-out whether transiesta is starting, or siesta is starting
      call ts_method_init( abs(init_method) == 2 )

    end if

    if ( TSmode .and. IsVolt .and. abs(init_method) /= 2 ) then
      if ( IONode ) then
        write(*,'(a,/)') 'ts: We highly recommend you to perform 0-bias calculations before *any* bias calculations'
      end if
      if ( .not. TS_start_bias ) then
        call die('ts: You have to calculate the 0 V and re-use the TSDE from &
            &that calculation.')
      end if
    end if

    if ( TSrun ) then

      ! Correct the convergence parameters in transiesta
      call val_swap(dDtol,ts_Dtol)
      call val_swap(dHtol,ts_Htol)

      ! From now on, a new mixing cycle starts,
      ! Check in mixer.F for new mixing schemes.
      if ( associated(ts_scf_mixs, target=scf_mixs) ) then
        do i = 1 , size(scf_mix%stack)
          call reset(scf_mix%stack(i), -ts_hist_keep)
        end do
      else
        call mixers_history_init(scf_mixs)
      end if
      ! Transfer scf_mixing to the transiesta mixing routine
      scf_mix => ts_scf_mixs(1)

    end if

  end subroutine init_DM

  !> Routine for reading the DM from a file.  
  !> This is a simple read-inset routine which reads
  !> a DM/TSDE file, and inserts the quantities
  !> into the the resulting DM (and/or EDM).
  subroutine init_DM_file(spin, no_u, nsc, &
       dit, sp, DM_2D, EDM_2D, &
       init_method)


    ! If the readed DM file has a different number of spin-components,
    ! this routine will easily extrapolate the quantities:
    !
    !   none => polarized/non-collinear/spin-orbit
    !       DM(:,1) = DM_read(:,1) / 2
    !       DM(:,2) = DM_read(:,1) / 2
    !
    !   polarized => none
    !       DM(:,1) = DM_read(:,1) + DM_read(:,2)
    !   polarized => non-collinear/spin-orbit
    !       DM(:,1) = DM_read(:,1)
    !       DM(:,2) = DM_read(:,2)
    !
    !   non-collinear => none
    !       DM(:,1) = DM_read(:,1) + DM_read(:,2)
    !   non-collinear => polarized
    !       DM(:,1) = DM_read(:,1)
    !       DM(:,2) = DM_read(:,2)
    !   non-collinear => spin-orbit
    !       DM(:,1) = DM_read(:,1)
    !       DM(:,2) = DM_read(:,2)
    !       DM(:,3) = DM_read(:,3)
    !       DM(:,4) = DM_read(:,4)
    !
    !   spin-orbit => none
    !       DM(:,1) = DM_read(:,1) + DM_read(:,2)
    !   spin-orbit => polarized
    !       DM(:,1) = DM_read(:,1)
    !       DM(:,2) = DM_read(:,2)
    !   spin-orbit => non-collinear
    !       DM(:,1) = DM_read(:,1)
    !       DM(:,2) = DM_read(:,2)
    !       DM(:,3) = DM_read(:,3)
    !       DM(:,4) = DM_read(:,4)

    ! Data-structures
    use class_Sparsity
    use class_dSpData2D
    use class_OrbitalDistribution
    use t_spin, only: tSpin

    use fdf, only: fdf_get
    use files, only : slabel
    use m_iodm, only : read_DM
    use m_ts_iodm, only: read_TS_DM
    use m_energies, only: Ef  ! Transiesta uses the EF obtained in a initial SIESTA run
    ! to place the electrodes and scattering region energy
    ! levels at the appropriate relative position, so it is
    ! stored in the TSDE file.
    use m_ts_global_vars,only: TSmode

    use m_handle_sparse, only: correct_supercell_SpD
    use m_handle_sparse, only: unfold_noauxiliary_supercell_SpD
    use m_handle_sparse, only: fold_auxiliary_supercell_SpD
    
    ! ********* INPUT ***************************************************
    ! type(tSpin) spin              : spin configuration for this system
    ! integer no_u                  : number of orbitals in unit-cell
    ! integer nsc(3)                : number of supercells along each direction
    ! type(OrbitalDistribution) dit : the distribution used for the orbitals
    ! type(Sparsity) sp             : sparsity pattern of DM
    ! type(dSpData2D) DM_2D         : the density matrix 
    ! type(dSpData2D) DM_2D         : the energy density matrix 
    ! integer init_method           : returned value to determine the
    !                                 method used to read the DM/EDM
    !                                   0 == not read
    !                                   1 == .DM read
    !                                   2 == .TSDE read
    ! *******************************************************************

    !> The spin-configuration that is used to determine the spin-order.
    type(tSpin), intent(in) :: spin
    !> Number of orbitals in the unit-cell
    integer, intent(in) :: no_u
    !> Number of supercells along each direction
    integer, intent(in) :: nsc(3)
    !> Parallel distribution of DM/EDM
    type(OrbitalDistribution), intent(in) :: dit
    !> Sparse pattern for DM/EDM
    type(Sparsity), intent(inout) :: sp
    !> The DM and EDM, these will be initialiazed upon return
    !> if the routine could read the files
    type(dSpData2D), intent(inout) :: DM_2D, EDM_2D

    !> To signal the method by which we have read DM/EDM
    integer, intent(out) :: init_method

    
    ! *** Local variables:
    logical :: corrected_nsc
    logical :: DM_found
    logical :: TSDE_found
    ! The file we should read
    character(len=256) :: fname
    ! Number of spin-components read from the DM/TSDE file
    integer :: nspin_read
    integer :: nsc_read(3)

    ! The currently read stuff
    type(dSpData2D) :: DM_read
    type(Sparsity), pointer :: sp_read
    type(dSpData2D) :: EDM_read

    ! Signal the file has not been found.
    init_method = 0
    
    if ( IONode ) write(*,*) ! New line

    DM_found = .false.
    TSDE_found = .false.
    
    if ( TSmode ) then
       if ( IONode ) write(*,'(a)',advance='no') &
            'Attempting to read DM, EDM from TSDE file... '

       ! Retrieve the name of the initialization file.
       fname = fdf_get('File.TSDE.Init',trim(slabel)//'.TSDE')

       ! Try and read the file
       call read_ts_dm(trim(fname), dit, nsc_read, DM_read, EDM_read, Ef, TSDE_found)

       if ( TSDE_found ) then
          ! Signal we have read TSDE
          init_method = 2

          DM_found = .true.

       else if ( IONode ) then
          write(*,'(a)') 'Failed...'
       end if

    end if

    if ( .not. DM_found ) then
       if ( IONode ) write(*,'(a)',advance='no') &
            'Attempting to read DM from file... '

       ! Retrieve the name of the initialization file.
       fname = fdf_get('File.DM.Init',trim(slabel)//'.DM')

       call read_DM(trim(fname), dit, nsc_read, DM_read, DM_found)

       if ( DM_found ) then
         ! Signal that the DM file has been found
         init_method = 1
         
       end if
       
    end if

    ! If DM is found, we check that it indeed is conformant with the number
    ! of orbitals
    if ( DM_found ) then

       if ( nrows_g(DM_read) /= nrows_g(sp) ) then
          if ( IONode ) then
             write(*,'(a)') 'Failed...'
             write(*,"(2(a,i0,/),a)") &
                  "WARNING: Wrong number of orbitals in DM file: ",nrows_g(DM_read), &
                  "WARNING: Expected number of orbitals in DM file: ",nrows_g(sp), &
                  "WARNING: Falling back to alternate initialization of DM."
          end if
          
          DM_found = .false.
          TSDE_found = .false.
          
       end if

    else if ( IONode ) then
       
       ! The DM/TSDE was not found... Signal it has not been found
       write(*,'(a)') 'Failed...'
       
    end if


    ! Density matrix size checks
    if ( DM_found ) then

      ! Specify the *correct* init_method
      ! In cases where the sparsity pattern changes, we assume this has to do with
      ! MD simulations and thus a geometry change.
      ! By default we do not allow mixing the first step when re-using a DM from another
      ! geometry. This is because of idempotency.
      sp_read => spar(DM_read)
      if ( .not. equivalent(sp, sp_read) ) then
        init_method = -init_method
      end if

      corrected_nsc = .false.
      if ( nsc_read(1) /= 0 .and. any(nsc /= nsc_read) ) then

        ! There are three cases:
        if ( all(nsc_read == 1) .and. fdf_get('DM.Init.Unfold', .true.) ) then

          ! 1. The read DM was created from a Gamma-only calculation and thus nsc == 1, always.
          !    In this case we know that the DM elements in the Gamma calculation (io,jo)
          !    are replicated in all image cells (io, jo + i_s * no_u) where i_s is the image cell
          !    offfset, since there are no k-points and thus no phases to worry about.
          !
          !    In very special circumstances the user may avoid this behavior by
          !    by setting 'DM.Init.Unfold false' 
          
          call unfold_noauxiliary_supercell_SpD(sp, DM_read)
          if ( TSDE_found ) then
            call unfold_noauxiliary_supercell_SpD(sp, EDM_read)
          end if

        else if ( all(nsc == 1) ) then

          ! 2. The read DM was created from a calculation with a non-trivial auxiliary cell, 
          !    and the current calculation is Gamma-only. In this case it is necessary to fold back
          !    all supercell DM entries.

          call fold_auxiliary_supercell_SpD(sp, DM_read)
          if ( TSDE_found ) then
            call fold_auxiliary_supercell_SpD(sp, EDM_read)
          end if
            
        else

          ! 3. The read DM has a different supercell size. In this case we only copy those
          !    elements that we know exists in the call below.
         
          ! Correct the supercell information
          ! Even for EDM this will work because correct_supercell_SpD
          ! changes the sparse pattern in-place and EDM and DM have
          ! a shared sp
          call correct_supercell_SpD(nsc_read, DM_read, nsc)
          corrected_nsc = .true.
          
        end if

      end if

      nspin_read = size(DM_read, 2)

      if ( IONode ) then
        if ( spin%DM == nspin_read ) then
          write(*,'(a)') 'Succeeded...'
        else if ( spin%DM < nspin_read ) then
          write(*,'(a)') 'Succeeded by reducing the number of spin-components...'
        else
          write(*,'(a)') 'Succeeded by increasing the number of spin-components...'
        end if
      end if
      
      if ( IONode ) then
        write(*,'(a)') "DM from file:"
        call print_type(DM_read)
      end if

      call restruct_Data(spin%DM, DM_read, DM_2D, .not. corrected_nsc)
      if ( IONode ) then
        write(*,'(a)') "DM to be used:"
        call print_type(DM_2D)
      end if
      if ( TSDE_found ) then
        call restruct_Data(spin%EDM, EDM_read, EDM_2D, .false.)
      end if

    end if

    ! Put deletes here to avoid complicating the logic
    call delete(DM_read)
    call delete(EDM_read)

    if ( .not. DM_found ) init_method = 0

  contains

    !> Driver to fix both the spin dimension and the sparsity pattern
    !> of a DM object, which typically has been read from file.
    !>
    !> Note that only the cases in which one of the spin dimensions is 1
    !> are treated.
    
    subroutine restruct_Data(nspin, in_2D, out_2D, show_warning)
      
      use m_restruct_SpData2D, only: restruct_dSpData2D
      
      !> Current spin dimension in the program
      integer, intent(in) :: nspin
      !> Input (DM) bud, to be mined for info
      type(dSpData2D), intent(inout) :: in_2D
      !> Output (DM) bud, created
      type(dSpData2D), intent(inout) :: out_2D
      !> Whether to show sanity-check warnings 
      logical, intent(in) :: show_warning
      
      integer :: nspin_read, i
      real(dp), pointer :: A2(:,:)

      nspin_read = size(in_2D, 2)
      
      if ( nspin == 1 .and. nspin /= nspin_read ) then
         ! This SCF has 1 spin-component.
         
         ! The read DM has at least 2!
         ! Thus we sum the spinors to form the non-polarized version
         A2 => val(in_2D)
!$OMP parallel do default(shared), private(i)
         do i = 1 , size(A2, 1)
            A2(i,1) = A2(i,1) + A2(i,2)
         end do
!$OMP end parallel do

      end if

      !> Calls [[restruct_dSpData2D]] to re-structure the sparsity
      !> data to match the output DM, with maximum spin%DM number of
      !> spin-components. It returns a new out_2D bud.
      !> @note
      !> The current sparsity pattern `sp` is known here by host association from
      !> the parent routine, which is a bit confusing. 
      !> `out_2D` in the called routine. Here it is the `out` sparsity, corresponding
      !> to the current target sparsity in the program.
      !> @endnote
      call restruct_dSpData2D(in_2D, sp, out_2D, nspin, show_warning=show_warning)

      if ( nspin_read == 1 .and. nspin /= nspin_read ) then
         ! This SCF has more than 2 spin-components.
         ! The read DM has 1.
         ! Thus we divide the spinors to form the polarized case.
         A2 => val(out_2D)
!$OMP parallel do default(shared), private(i)
         do i = 1 , size(A2, 1)
            A2(i,1) = A2(i,1) * 0.5_dp
            A2(i,2) = A2(i,1)
         end do
!$OMP end parallel do
         
      end if
       
    end subroutine restruct_Data

  end subroutine init_DM_file

  subroutine init_DM_atomic(spin, na_u, no_u, lasto, iaorb, &
       DM_atom, &
       dit, sp, DM_2D, &
       anti_ferro)
    
    ! The DM is generated assuming atomic charging
    ! (filling up atomic orbitals). The DM originated that way is
    ! not a good DM due to overlaps, but the SCF cycling corrects
    ! that for the next cycle.

    ! Spin polarized calculations starting from atoms:
    !    Default: All atoms with maximum polarization compatible with
    !             atomic configuration. In Ferromagnetic ordering (up).
    !    If DM.InitSpinAF is true, as default but in Antiferro order:
    !             even atoms have spin down, odd up.
    !    If fdf %block DM.InitSpin is present it overwrites previous
    !         schemes: magnetic moments are explicitly given for some atoms.
    !         Atoms not mentioned in the block are initialized non polarized.

    ! Written by E. Artacho. December 1997. Taken from the original piece
    ! of siesta.f written by P. Ordejon.
    ! Non-collinear spin added by J.M.Soler, May 1998.
    ! Updated by N.R.Papior, December, 2016
    ! ********* INPUT ***************************************************
    ! type(tSpin) spin              : spin configuration for this system
    ! integer na_u                  : number of atoms in unit-cell
    ! integer no_u                  : number of orbitals in unit-cell
    ! integer lasto(0:na_u)         : last orbital of each atom
    ! integer iaorb(no_u)           : the atomic index of the corresponding orbital
    ! real(dp) DM_atom(no_u)        : atomic density based on atomic configuration
    ! type(OrbitalDistribution) dit : the distribution used for the orbitals
    ! type(Sparsity) sp             : sparsity pattern of DM
    ! type(dSpData2D) DM_2D         : the density matrix 
    ! logical anti_ferro            : whether the DM should be anti-ferro magnetic or not
    ! *******************************************************************

    use class_Sparsity
    use class_OrbitalDistribution
    use class_dSpData2D
    use class_OrbitalDistribution
    use class_dData2D
    use t_spin, only: tSpin

    use fdf
    use parsing
    use sparse_matrices, only: S  ! Note direct import of (associated now) pointer

#ifdef MPI
    use mpi_siesta
#endif

    ! The spin-configuration that is used to determine the spin-order.
    type(tSpin), intent(in) :: spin
    ! Number of atoms in the unit-cell
    integer, intent(in) :: na_u
    ! Number of orbitals in the unit-cell
    integer, intent(in) :: no_u
    ! The last orbital on each atom
    integer, intent(in) :: lasto(0:na_u)
    ! The atom containing orbital "io"
    integer, intent(in) :: iaorb(no_u)
    ! DM filled with atomic charges based on atom and orbitals
    real(dp), intent(in) :: DM_atom(no_u)
    ! Parallel distribution of DM
    type(OrbitalDistribution), intent(in) :: dit
    ! Sparse pattern for DM
    type(Sparsity), intent(inout) :: sp
    ! The DM, these will be initialiazed upon return by
    ! the atomic filling
    type(dSpData2D), intent(inout) :: DM_2D
    ! Setting for the magnetic setup of DM
    ! If .true. the spin-configuration will be anti-ferro-magnetic
    logical, intent(in) :: anti_ferro

    ! *** Local variables
    type(block_fdf) :: bfdf
    logical :: init_spin_block
    type(dData2D) :: arr_2D
    ! The pointers for the sparse pattern
    integer :: no_l, nnz
    integer :: io, gio, i, ind, jo
    integer, pointer :: ncol(:), ptr(:), col(:)
    real(dp), pointer :: DM(:,:)
    
    ! First we need to figure out the options
    if ( spin%DM > 1 ) then
       init_spin_block = fdf_block('DM.InitSpin', bfdf)
    else
       ! Do not read the initspin block in case of
       ! non-polarized calculation
       init_spin_block = .false.
    end if

    ! Retrieve pointers to data
    call attach(sp, nrows=no_l, nnzs=nnz, &
         n_col=ncol, list_ptr=ptr, list_col=col)

    ! First we create the DM
    call newdData2D(arr_2D, nnz, spin%DM, "DM")
    call newdSpData2D(sp, arr_2D, dit, DM_2D, &
         "DM initialized from atoms")
    call delete(arr_2D)
    DM => val(DM_2D)

    if ( .not. init_spin_block ) then
       call init_atomic()
    else
       call init_user()
    end if

    ! We have initialized with atomic information. Correct in case we
    ! are using such a small cell that there are direct interactions
    ! of orbitals with their own images, and we insist on using the
    ! Gamma-point only. Otherwise S(diagonal) is always 1.0 and the
    ! simple atomic-orbital superposition works as intended.

      do io = 1, no_l
         ! Retrieve global orbital index
         gio = index_local_to_global(dit, io)
         do i = 1 , ncol(io)
            ind = ptr(io) + i
            jo = col(ind)
            ! Check for diagonal element
            if ( gio == jo ) then
               DM(ind,:) = DM(ind,:) / S(ind)
            endif
         end do
      end do

  contains

    subroutine init_atomic()
      
      integer :: io, gio, i, ind, jo, i1, i2
      ! This subroutine initializes the DM based on
      ! the anti_ferro
      
!$OMP parallel default(shared), private(i1,i2)

      ! Initialize to 0
      do i2 = 1 , spin%DM
!$OMP do
         do i1 = 1 , nnz
            DM(i1,i2) = 0._dp
         end do
!$OMP end do nowait
      end do
       
!$OMP barrier
      
      ! Automatic, for non magnetic (nspin=1) or for Ferro or Antiferro
!$OMP single
      do io = 1, no_l
         
         ! Retrieve global orbital index
         gio = index_local_to_global(dit, io)
         
!$OMP task firstprivate(io,gio), private(i,ind,jo,i1,i2)
         do i = 1 , ncol(io)
            ind = ptr(io) + i
            jo = col(ind)
            ! Check for diagonal element
            if ( gio /= jo ) cycle
            
            if ( spin%DM == 1 ) then
               
               ! No spin polarization
               DM(ind,1) = DM_atom(gio)
               
            else
               
               ! Spin polarization
               i1 = 1
               i2 = 2
               
               ! Ferro or antiferro
               if ( anti_ferro .and. &
                    mod(iaorb(gio),2) == 0 ) then
                  i1 = 2
                  i2 = 1
               end if
               
               DM(ind,i1) = min( DM_atom(gio), 1._dp )
               DM(ind,i2) = DM_atom(gio) - DM(ind,i1)
               
            end if

         end do
!$OMP end task
      end do
!$OMP end single nowait
      
!$OMP end parallel

      if ( IONode ) then
         write(*,'(a)') "DM filled with atomic data:"
         call print_type(DM_2D)
      end if

    end subroutine init_atomic

    subroutine init_user()

      use units, only: deg

      type(parsed_line), pointer :: pline

      integer, pointer :: atom_dat(:) => null()
      integer, pointer :: atom
      real(dp), pointer :: spin_dat(:,:) => null()
      real(dp), pointer :: phi, theta, spin_val

      real(dp), parameter :: epsilon = 1.e-8_dp
      ! To grab +/-
      character(len=1) :: updo
      character(len=80) :: msg
      logical :: bad_syntax, non_col

      integer :: ni, nn, nr, na
      integer :: i, io, jo, gio, ia, ind, is
      real(dp) :: cosph, sinph, costh, sinth
      real(dp) :: qio, rate, spin_at, spio


      ! Define defaults
      non_col = .false.
      bad_syntax = .false.

      ! Allocate options
      call re_alloc( atom_dat, 1, na_u, 'atom_dat', 'init_dm_user')
      call re_alloc( spin_dat, 1, na_u, 1, 3, 'spin_dat', 'init_dm_user' )
      
      ! Read the data from the block and then we populate DM
      na = 0
      do while( fdf_bline(bfdf,pline) )
         
         ! Read number of names, integers and reals on this line
         ni = fdf_bnintegers(pline)
         nn = fdf_bnnames(pline)
         nr = fdf_bnreals(pline)

         ! Simple check for bad syntax
         if ( ni /= 1 ) then
            bad_syntax = .true.
            exit
         end if

         na = na + 1
         if ( na > na_u ) then
            call die('There can only be one initial spin-component per atom.')
         end if

         ! Point to the data
         atom     => atom_dat(na)
         spin_val => spin_dat(na,1)
         phi      => spin_dat(na,2)
         theta    => spin_dat(na,3)
         
         ! Default values for the not always supplied data
         phi   = 0._dp
         theta = 0._dp

         ! Read atom involved
         atom = fdf_bintegers(pline,1)
         
         ! Check that the atom is correct
         ! Allow for negative indices...
         if ( atom < 0 ) atom = atom + na_u + 1

         ! Now check that it is within range...
         if ( atom < 1 .or. na_u < atom ) then
            bad_syntax = .true.
            exit
         end if


         ! Now it depends on the type of input...
         if ( nn == 0 ) then
            ! There are no names in this one...
            
            ! Read value of spin
            if ( nr == 3 ) then
               spin_val = fdf_breals(pline,1)
               theta    = fdf_breals(pline,2) * deg
               phi      = fdf_breals(pline,3) * deg
            else if ( nr == 1 ) then
               spin_val = fdf_breals(pline,1)
            else
               ! Print bad-syntax error and stop
               bad_syntax = .true.
               exit
            end if
            
         else if ( nn == 1 ) then
            ! Read spin as + or - (maximun value)
            updo = fdf_bnames(pline,1)
            if ( updo .eq. '+' ) then
               spin_val =  100._dp
            else if ( updo .eq. '-' ) then
               spin_val = -100._dp
            else
               bad_syntax = .true.
               exit
            end if
            
            if ( nr == 2 ) then
               theta = fdf_breals(pline,1) * deg
               phi   = fdf_breals(pline,2) * deg
            else if ( nr /= 0 ) then
               bad_syntax = .true.
               exit
            end if
         else
            ! Print bad-syntax error and stop
            bad_syntax = .TRUE.
            exit 
         end if

         if ( abs(theta) > 1.d-12 ) non_col = .true.
         if ( abs(phi) > 1.d-12 ) non_col = .true.
         
      end do

      if ( bad_syntax ) then
         write(msg,'(a,i4)') &
              'initDM: ERROR: bad syntax in DM.InitSpin, line', na
         call die(trim(msg))
      end if

      if ( non_col .and. spin%DM < 4 ) then
         ! The user has requested non-collinear spin.
         ! This is not applicable to this simulation
         
          if ( IONode ) then
             write(*,'(/,2a)') 'initDM: WARNING: noncollinear spins ',  &
                  'in DM.InitSpin not used because nspin < 4'
          end if
          non_col = .false.
          
       end if

       ! Now we may fill DM

!$OMP parallel default(shared), private(i,ind,is,jo,gio,ia,io)

       ! Initialize to 0
       do is = 1 , spin%DM
!$OMP do
          do i = 1 , nnz
             DM(i,is) = 0._dp
          end do
!$OMP end do nowait
       end do

!$OMP barrier
       
       ! Initialize the paramagnetic case
!$OMP single
       do io = 1 , no_l
          ! Get global index
          gio = index_local_to_global(dit, io)

          ! Get atomic index
          ia = iaorb(gio)
          
!$OMP task firstprivate(io,gio)
          do i = 1, ncol(io)
             ind = ptr(io) + i
             jo = col(ind)
             if ( gio == jo ) then
                DM(ind, 1) = 0.5_dp * DM_atom(gio)
                DM(ind, 2) = DM(ind, 1)
             end if
          end do
!$OMP end task

       end do
!$OMP end single

       ! Now correct for the user input 
!$OMP single
       do ia = 1, na
          
          ! Get data
          atom     => atom_dat(ia)
          spin_val => spin_dat(ia,1)
          phi      => spin_dat(ia,2)
          theta    => spin_dat(ia,3)

          ! Find maximum atomic moment that the atoms involved can carry
          spin_at = 0._dp
          do gio = lasto(atom-1) + 1, lasto(atom)
             spin_at = spin_at + min( DM_atom(gio), 2._dp - DM_atom(gio) )
          end do

         
          if ( spin_at < epsilon ) then
             ! This is the simple case.
             ! The atom cannot be spin-polarized and thus
             ! we may easily set the contribution
             
             if ( IONode ) then
                write(*,'(a,i0,a)')  &
                     'initDM: WARNING: atom ', atom, &
                     ' has a closed-shell and cannot be polarized'
             end if

             ! We have already initialized for the
             ! paramagnetic case (no polarization)

          else

             ! If given spin is larger than possible, make it to max atomic
             if ( abs(spin_val) > spin_at .and. &
                  abs(spin_val) > epsilon ) &
                  spin_val = spin_at * spin_val / abs(spin_val)
             
             ! Initialize orbitals with same rate as atom
             rate = spin_val / spin_at
             do gio = lasto(atom-1) + 1, lasto(atom)
                io = index_global_to_local(dit, gio)
             
                ! Immediately skip if not on this node... 
                if ( io < 1 ) cycle

!$OMP task firstprivate(gio,io,rate), &
!$OMP& private(qio,spio,costh,sinth,cosph,sinph)
                ! The orbital charge per spin
                qio = DM_atom(gio) * 0.5_dp
                spio = rate * min( qio, 1._dp - qio )
                do i = 1 , ncol(io)
                   ind = ptr(io) + i
                   jo = col(ind)
                   
                   ! Immediately skip if not diagonal term..
                   if ( gio /= jo ) cycle

                   if ( non_col ) then
                      
                      ! Store non-collinear-spin density matrix as
                      !   ispin=1 => D11
                      !   ispin=2 => D22
                      !   ispin=3 => Real(D12)
                      !   ispin=4 => Imag(D12)
                      costh = cos(theta)
                      sinth = sin(theta)
                      cosph = cos(phi)
                      sinph = sin(phi)
                      
                      DM(ind,1) = qio + spio * costh
                      DM(ind,2) = qio - spio * costh
                      DM(ind,3) =       spio * sinth * cosph
                      DM(ind,4) =       spio * sinth * sinph
                      if ( spin%SO ) then ! spin-orbit coupling
                         DM(ind,5) = 0._dp
                         DM(ind,6) = 0._dp
                         DM(ind,7) = DM(ind,3)
                         DM(ind,8) = DM(ind,4)
                      end if
                      
                   else
                      
                      DM(ind,1) = qio + spio
                      DM(ind,2) = qio - spio
                      
                   end if
                   
                end do
!$OMP end task
             end do
             
          end if
       enddo
!$OMP end single nowait

!$OMP end parallel
       
      call de_alloc(atom_dat, 'atom_dat', 'init_dm_user')
      call de_alloc(spin_dat, 'spin_dat', 'init_dm_user')

      if ( IONode ) then
         write(*,'(a)') "DM filled with atomic data (user-defined):"
         call print_type(DM_2D)
      end if

    end subroutine init_user

  end subroutine init_DM_atomic

  
  subroutine extrapolate_dm_with_coords(DM_history,na_u,xa,sparse_pattern,DMnew)
    use class_Sparsity
    use class_dData2D
    use class_OrbitalDistribution
    use class_dSpData2D
    use class_Geometry
    use class_Pair_Geometry_dSpData2D
    use class_Fstack_Pair_Geometry_dSpData2D

    use m_restruct_SpData2D, only: restruct_dSpData2D
    use fdf, only: fdf_get

    type(Fstack_Pair_Geometry_dSpData2D), intent(in) :: DM_history
    integer, intent(in)                             :: na_u 
    real(dp), intent(in)                            :: xa(:,:)
    type(Sparsity), intent(in)                      :: sparse_pattern
    type(dSpData2D), intent(inout)                   :: DMnew

    integer :: n, i, nspin, nnzs_out
    real(dp), allocatable   :: c(:)
    real(dp), allocatable   :: xan(:,:,:), dummy_cell(:,:,:)
    type(Geometry), pointer :: geom 
    type(dSpData2D), pointer :: dm
    type(OrbitalDistribution), pointer    :: orb_dist
    type(Pair_Geometry_dSpData2D), pointer :: pair

    type(dSpData2D)       :: DMtmp
    type(dData2D)       :: a_out

    real(dp), dimension(:,:) , pointer  :: a, ai, xp

    n = n_items(DM_history)
    allocate(c(n))

    allocate(xan(3,na_u,n),dummy_cell(3,3,n))

    do i = 1, n
       pair => get_pointer(DM_history,i)
       call firstp(pair,geom)
       xp => coords(geom)
       xan(:,:,i) = xp(:,:)
       dummy_cell(:,:,i) = 1.0_dp
    enddo
    if (fdf_get("UseDIISforDMExtrapolation",.true.)) then
       ! Cast Jose Soler's idea into a DIIS framework

       if (fdf_get("UseSVDExperimental",.false.)) then
          ! Attempt to use the "alternate" KF method with
          ! first differences.  It does not work well yet
          call extrapolate_diis_svd_new(na_u,n,dummy_cell,  &
               xan,dummy_cell(:,:,1),xa,c)
       else if (fdf_get("UseSVD",.true.)) then
          ! Straightforward SVD
          call extrapolate_diis_svd(na_u,n,dummy_cell,  &
               xan,dummy_cell(:,:,1),xa,c)
       else
          call extrapolate_diis(na_u,n,dummy_cell, &
               xan,dummy_cell(:,:,1),xa,c)
       endif
    else  
       ! Use Jose Soler's original method
       call extrapolate(na_u,n,dummy_cell,xan,dummy_cell(:,:,1),xa,c)
    endif
    if (ionode) then
       print *, "DM extrapolation coefficients: "
       do i = 1, n
          print "(i0,f10.5)", i, c(i)
       enddo
    endif

    pair => get_pointer(DM_history,1)
    call secondp(pair,dm)

    ! We assume that all DMs in the history stack have the same orbital distribution...
    orb_dist => dist(dm)
    a => val(dm)
    nspin = size(a,dim=2)
    nnzs_out = nnzs(sparse_pattern)

    ! Scratch array to accumulate the elements
    call newdData2D(a_out,nnzs_out, nspin,name="(temp array for extrapolation)")
    a => val(a_out)
    a(:,:) = 0.0_dp

    do i = 1, n
       pair => get_pointer(DM_history,i)
       call secondp(pair,dm)
       !           if (.not. associated(orb_dist,dist(dm))) then
       !              call die("Different orbital distributions in DM history stack")
       !           endif
       call restruct_dSpData2D(dm,sparse_pattern,DMtmp)
       ai => val(DMtmp)
       a(:,:) = a(:,:) + c(i) * ai(:,:)
    enddo

    call newdSpData2D(sparse_pattern,a_out,orb_dist, &
         DMnew,name="SpM extrapolated using coords")
    call delete(a_out)
    call delete(DMtmp)
    deallocate(xan,c,dummy_cell)

  end subroutine extrapolate_dm_with_coords

  subroutine get_allowed_history_depth(n)
    ! 
    ! Encapsulates the logic of DM extrapolation and re-use
    ! (work in progress)
    use fdf, only: fdf_get
    use siesta_options, only: DM_history_depth, harrisfun, idyn

    integer, intent(out)         :: n

    n = DM_history_depth     ! As set by default or by the user

    if (harrisfun) then
       n = 0
       if (ionode) print "(a)", &
            "DM_history_depth set to zero for 'Harris' run"
       return
    else if (.not. fdf_get("DM.AllowReuse",.true.)) then
       n = 0
       if (ionode) print "(a)",   &
            "DM_history_depth set to zero since no re-use is allowed"
       return
    else if (.not. fdf_get("DM.AllowExtrapolation",.true.)) then
       n = 1
       if (ionode) print "(a)", &
            "DM_history_depth set to one since no extrapolation is allowed"
       return
    else if (idyn .eq. 0)  then   ! Geometry relaxation
       if (fdf_get("DM.AllowExtrapolation",.false.)) then
          if (ionode) print "(a,i0)", "Requested Extrapolation for geometry relaxation. DM_history_depth: ", n
       else	
          n = 1
          if (ionode) print "(a)", &
               "DM_history_depth set to one: no extrapolation allowed by default for geometry relaxation"
       endif
       return
    else if ((idyn .eq. 6) .OR. (idyn .eq. 7))  then   ! Forces or Phonon series
       n = 1
       if (ionode) print "(a)", &
            "DM_history_depth set to one for 'Forces' run"
       return
    endif
  end subroutine get_allowed_history_depth


  ! *******************************************************************
  ! SUBROUTINE extrapolate( na, n, cell, xa, cell0, x0, c )
  ! *******************************************************************
  ! Finds optimal coefficients for an approximate expasion of the form
  ! D_0 = sum_i c_i D_i, where D_i is the density matrix in the i'th
  ! previous iteration, or any other function of the atomic positions. 
  ! In practice, given points x_0 and x_i, i=1,...,n, it finds the
  ! coefficients c_i such that sum_i c_i*x_i, is closest to x_0, with
  ! sum_i c_i = 1. Unit cell vectors are included for completeness of
  ! the geometry specification only. Though not used in this version,
  ! this routine can be used if cell vectors change (see algorithms).
  ! Written by J.M.Soler. Feb.2010.
  ! ************************* INPUT ***********************************
  !  integer  na          ! Number of atoms
  !  integer  n           ! Number of previous atomic positions
  !  real(dp) cell(3,3,n) ! n previous unit cell vectors
  !  real(dp) xa(3,na,n)  ! n previous atomic positions (most recent first)
  !  real(dp) cell0(3,3)  ! Present unit cell vectors
  !  real(dp) x0(3,na)    ! Present atomic positions
  ! ************************* OUTPUT **********************************
  !  real(dp) c(n)        ! Expansion coefficients
  ! ************************ UNITS ************************************
  ! Unit of distance is arbitrary, but must be the same in all arguments
  ! ********* BEHAVIOUR ***********************************************
  ! - Returns without any action if n<1 or na<1
  ! - Stops with an error message if matrix inversion fails
  ! ********* DEPENDENCIES ********************************************
  ! Routines called: 
  !   inver     : Matrix inversion
  ! Modules used:
  !   precision : defines parameter 'dp' (double precision real kind)
  !   sys       : provides the stopping subroutine 'die'
  ! ********* ALGORITHMS **********************************************
  ! Using a linear approximation for D(x), imposing that sum(c)=1, and
  ! assuming that <dD/dx_i*dD/dx_j>=0, <(dD/dx_i)**2>=<(dD/dx_j)**2>
  ! (i.e. assuming no knowledge on the parcial derivatives), we find
  ! (D(x0)-D(sum_i c_i*x_i))**2 = const * (x0-sum_i c_i*x_i)**2
  ! Therefore, given points x_0 and x_i, i=1,...,n, we look for the
  ! coefficients c_i such that sum_i c_i*x_i, is closest to x_0, with
  ! sum_i c_i = 1. It is straighforward to show that this reduces to 
  ! solving S*c=s0, where S_ij=x_i*x_j and s0_i=x_i*x0, under the
  ! constraint sum(c)=1. To impose it, we rewrite the expansion as
  ! xmean + sum_i c_i*(x_i-xmean), where xmean=(1/n)*sum_i x_i.
  ! Since the vectors (x_i-xmean) are linearly dependent, the matrix
  ! S_ij=(x_i-xmean)*(x_j-xmean) is singular. Therefore, we substitute
  ! the last row of S and s0 by 1, thus imposing sum(c)=1.
  ! Unit cell vectors are not used in this version, but this routine
  ! can be safely used even if cell vectors change, since D depends
  ! on the interatomic distances, and it can be shown that
  ! sum_ia,ja (xa(:,ia,i)-xa(:,ja,i))*(xa(:,ia,j)-xa(:,ja,j))
  !   = 2*na*sum_ia (xa(:,ia,i)-xmean(:,ia))*(xa(:,ia,j)-xmean(:,ia))
  !   = 2*na*S_ij
  ! This implies that approximating the present ineratomic distances,
  ! as an expansion of previous ones, is equivalent to approximating
  ! the atomic positions.
  ! *******************************************************************

  SUBROUTINE extrapolate( na, n, cell, xa, cell0, x0, c )

    ! Used procedures and parameters
    USE sys,       only: message           ! Termination routine
    USE precision, only: dp            ! Double precision real kind
    use parallel,  only: Node

    ! Passed arguments
    implicit none
    integer, intent(in) :: na          ! Number of atoms
    integer, intent(in) :: n           ! Number of previous atomic positions
    real(dp),intent(in) :: cell(3,3,n) ! n previous unit cell vectors
    real(dp),intent(in) :: xa(3,na,n)  ! n previous atomic positions
    real(dp),intent(in) :: cell0(3,3)  ! Present unit cell vectors
    real(dp),intent(in) :: x0(3,na)    ! Present atomic positions
    real(dp),intent(out):: c(n)        ! Expansion coefficients

    ! Internal variables and arrays
    integer :: i, ierr, j, m, ix, ia
    real(dp):: s(n,n), s0(n), si(n,n), xmean(3,na)

    ! Trap special cases
    if (na<1 .or. n<1) then
       return
    else if (n==1) then
       c(1) = 1
       return
    end if

    ! Find average of previous positions
    do ia=1,na
       do ix=1,3
          xmean(ix,ia) = sum(xa(ix,ia,1:n)) / n
       enddo
    enddo

    ! The above is equivalent to
    !   xmean = sum(xa,dim=3)/n

    ! Find matrix s of dot products. Subtract xmean to place origin within the
    ! hyperplane of x vectors

    do j = 1,n
       s0(j) = sum( (x0-xmean) * (xa(:,:,j)-xmean) )
       do i = j,n
          s(i,j) = sum( (xa(:,:,i)-xmean) * (xa(:,:,j)-xmean) )
          s(j,i) = s(i,j)
       end do
    end do

    ! Find the largest number of (first) m linearly independent vectors xa-xmean.
    ! Notice that m<n because we have subtracted xmean. Optimally, we should not
    ! restrict ourselves to the first (most recent) m vectors, but that would
    ! complicate the code too much.
    do m = n-1,0,-1
       if (m==0) exit                   ! Do not try to invert a 0x0 'matrix'
       call inver( s, si, m, n, ierr )
       if (ierr==0) exit                ! Not too elegant, but it will do.
    end do

    if (node==0) print *, " # of linearly independent vectors: ", m

    ! Trap the case in which the first two xa vectors are equal
    if (m==0) then  ! Just use most recent point only
       c(n) = 1
       c(1:n-1) = 0
       return
    end if

    ! Set one more row for equation sum(c)=1.
    m = m+1
    s0(m) = 1
    s(m,1:m) = 1

    ! Invert the full equations matrix
    call inver( s, si, m, n, ierr )
    if (ierr/=0) then
       c(:) = 0.0_dp
       c(n) = 1.0_dp
       call message('WARNING','extrapolate: matrix inversion failed')
       call message('WARNING','extrapolate: using last item in history')
       return
    endif

    ! Find expansion coefficients
    ! This is wrong regarding the order...
    c(1:m) = matmul( si(1:m,1:m), s0(1:m) )
    c(m+1:n) = 0

  END SUBROUTINE extrapolate

  ! *******************************************************************
  ! SUBROUTINE extrapolate_diis( na, n, cell, xa, cell0, x0, c )
  ! *******************************************************************
  ! Finds optimal coefficients for an approximate expasion of the form
  ! D_0 = sum_i c_i D_i, where D_i is the density matrix in the i'th
  ! previous iteration, or any other function of the atomic positions. 
  ! In practice, given points x_0 and x_i, i=1,...,n, it finds the
  ! coefficients c_i such that sum_i c_i*x_i, is closest to x_0, with
  ! sum_i c_i = 1. Unit cell vectors are included for completeness of
  ! the geometry specification only. Though not used in this version,
  ! this routine can be used if cell vectors change (see algorithms).
  ! Original version Written by J.M.Soler. Feb.2010.
  ! Couched in DIIS language by A. Garcia, Nov. 2012.
  ! ************************* INPUT ***********************************
  !  integer  na          ! Number of atoms
  !  integer  n           ! Number of previous atomic positions
  !  real(dp) cell(3,3,n) ! n previous unit cell vectors
  !  real(dp) xa(3,na,n)  ! n previous atomic positions (most recent first)
  !  real(dp) cell0(3,3)  ! Present unit cell vectors
  !  real(dp) x0(3,na)    ! Present atomic positions
  ! ************************* OUTPUT **********************************
  !  real(dp) c(n)        ! Expansion coefficients
  ! ************************ UNITS ************************************
  ! Unit of distance is arbitrary, but must be the same in all arguments
  ! ********* BEHAVIOUR ***********************************************
  ! - Returns without any action if n<1 or na<1
  ! - Stops with an error message if matrix inversion fails
  ! ********* DEPENDENCIES ********************************************
  ! Routines called: 
  !   inver     : Matrix inversion
  ! Modules used:
  !   precision : defines parameter 'dp' (double precision real kind)
  !   sys       : provides the stopping subroutine 'die'
  ! ********* ALGORITHMS **********************************************

  SUBROUTINE extrapolate_diis( na, n, cell, xa, cell0, x0, c )

    ! Used procedures and parameters
    USE sys,       only: message, die
    USE precision, only: dp            ! Double precision real kind
    use parallel,  only: Node

    ! Passed arguments
    implicit none
    integer, intent(in) :: na          ! Number of atoms
    integer, intent(in) :: n           ! Number of previous atomic positions
    real(dp),intent(in) :: cell(3,3,n) ! n previous unit cell vectors
    real(dp),intent(inout) :: xa(3,na,n)  ! n previous atomic positions
    real(dp),intent(in) :: cell0(3,3)  ! Present unit cell vectors
    real(dp),intent(in) :: x0(3,na)    ! Present atomic positions
    real(dp),intent(out):: c(n)        ! Expansion coefficients

    ! Internal variables and arrays
    integer :: i, info, j, m, ix, ia
    real(dp), allocatable, dimension(:,:) :: s, si

    allocate (s(n+1,n+1), si(n+1,n+1))

    ! Trap special cases
    if (na<1 .or. n<1) then
       return
    else if (n==1) then
       c(1) = 1
       return
    end if

    ! Find residuals with respect to x0 
    do i = 1, n
       do ia=1,na
          do ix=1,3
             xa(ix,ia,i) = xa(ix,ia,i) - x0(ix,ia)
          enddo
       enddo
    enddo

    ! Find matrix s of dot products.

    do j = 1,n
       do i = j,n
          s(i,j) = sum( xa(:,:,i) * xa(:,:,j) )
          s(j,i) = s(i,j)
       end do
       ! Now extend the matrix with ones in an extra column
       ! and row ...
       s(j,n+1)=1.0_dp   ! This value is really arbitrary
       s(n+1,j)=1.0_dp   ! This represents the Sum_{Ci} = 1 constraint
    end do
    s(n+1,n+1) = 0.0_dp

    if (Node == 0) then
       call print_mat(s,n+1)
    endif

    ! Find rank of the xi*xj matrix
    do m = n,0,-1
       if (m==0) exit                   ! Do not try to invert a 0x0 'matrix'
       call inver( s, si, m, n+1, info )
       if (info==0) exit                ! Not too elegant, but it will do.
    end do
    if (node==0) print *, " Rank of DIIS matrix: ", m

    if (m<n) then
       info = -1
    else
       ! Invert the full matrix
       call inver( s, si, n+1, n+1, info)
    endif
    c(:) = 0.0_dp

    ! If inver was successful, get coefficients for DIIS/Pulay mixing
    ! (Last column of the inverse matrix, corresponding to solving a
    ! linear system with (0,0,0,...,0,1) in the right-hand side)
    if (info .eq. 0) then
       do i=1,n
          c(i)=si(i,n+1)
       enddo
    else
       ! Otherwise, use only last step
       if (Node == 0) then
          write(*,"(a,i5)")  &
               "Warning: unstable inversion in DIIS - use last item only"
       endif
       c(n) = 1.0_dp
    endif

    deallocate(s,si)
  end SUBROUTINE extrapolate_diis

  SUBROUTINE extrapolate_diis_svd( na, n, cell, xa, cell0, x0, c )

    ! Used procedures and parameters
    USE sys,       only: message, die
    USE precision, only: dp            ! Double precision real kind
    use parallel,  only: Node
    use m_svd,     only: solve_with_svd

    ! Passed arguments
    implicit none
    integer, intent(in) :: na          ! Number of atoms
    integer, intent(in) :: n           ! Number of previous atomic positions
    real(dp),intent(in) :: cell(3,3,n) ! n previous unit cell vectors
    real(dp),intent(inout) :: xa(3,na,n)  ! n previous atomic positions
    real(dp),intent(in) :: cell0(3,3)  ! Present unit cell vectors
    real(dp),intent(in) :: x0(3,na)    ! Present atomic positions
    real(dp),intent(out):: c(n)        ! Expansion coefficients

    ! Internal variables and arrays
    integer :: i, info, j, m, ix, ia, rank
    real(dp), allocatable, dimension(:,:) :: s, si
    real(dp), allocatable, dimension(:)   :: rhs, sigma, beta

    allocate (s(n+1,n+1), si(n+1,n+1), sigma(n+1), rhs(n+1))
    allocate (beta(n+1))  ! full set of coefficients

    ! Trap special cases
    if (na<1 .or. n<1) then
       return
    else if (n==1) then
       c(1) = 1
       return
    end if

    ! Find residuals with respect to x0 
    do i = 1, n
       do ia=1,na
          do ix=1,3
             xa(ix,ia,i) = xa(ix,ia,i) - x0(ix,ia)
          enddo
       enddo
    enddo

    ! Find matrix s of dot products.

    do j = 1,n
       do i = j,n
          s(i,j) = sum( xa(:,:,i) * xa(:,:,j) )
          s(j,i) = s(i,j)
       end do
       ! Now extend the matrix with ones in an extra column
       ! and row ...
       s(j,n+1)=1.0_dp   ! This value is really arbitrary
       s(n+1,j)=1.0_dp   ! This represents the Sum_{Ci} = 1 constraint
    end do
    s(n+1,n+1) = 0.0_dp

    if (Node == 0) then
       ! call print_mat(s,n+1)
    endif

    ! Find rank of the xi*xj matrix
    do m = n,0,-1
       if (m==0) exit                   ! Do not try to invert a 0x0 'matrix'
       call inver( s, si, m, n+1, info )
       if (info==0) exit                ! Not too elegant, but it will do.
    end do
    if (node==0) print *, " Estimated Rank of xi*xj matrix: ", m

    rhs(1:n) = 0.0_dp
    rhs(n+1) = 1.0_dp
    !
    call solve_with_svd(s,rhs,beta,info,sigma=sigma,rank_out=rank)
    !
    if (node==0) then
       print *, 'matrix rank: ', rank
       print "(a,/,6g12.5)", 'Singular values: ', sigma(:)
    endif
    if (info == 0) then
       c(1:n) = beta(1:n)
    else
       if (node==0) then
          print *, 'The SVD algorithm failed to converge'
          print *, 'Re-using last item only...'
          print *, "info: ", info
       endif
       c(:) = 0.0_dp
       c(n) = 1.0_dp
    endif

    deallocate(s,si,sigma,rhs,beta)
  end SUBROUTINE extrapolate_diis_svd

  SUBROUTINE extrapolate_diis_svd_new( na, n, cell, xa, cell0, x0, c )

    ! Used procedures and parameters
    USE sys,       only: message, die
    USE precision, only: dp            ! Double precision real kind
    use parallel,  only: Node
    use m_svd,     only: solve_with_svd

    ! Passed arguments
    implicit none
    integer, intent(in) :: na          ! Number of atoms
    integer, intent(in) :: n           ! Number of previous atomic positions
    real(dp),intent(in) :: cell(3,3,n) ! n previous unit cell vectors
    real(dp),intent(inout) :: xa(3,na,n)  ! n previous atomic positions
    real(dp),intent(in) :: cell0(3,3)  ! Present unit cell vectors
    real(dp),intent(in) :: x0(3,na)    ! Present atomic positions
    real(dp),intent(out):: c(n)        ! Expansion coefficients

    ! Internal variables and arrays
    integer :: i, info, j, ix, ia, rank
    real(dp), allocatable, dimension(:,:) :: s
    real(dp), allocatable, dimension(:)   :: rhs, sigma
    real(dp), allocatable, dimension(:)   :: beta ! intermediate coeffs

    allocate (s(n-1,n-1), rhs(n-1), beta(n-1), sigma(n-1))

    ! Trap special cases
    if (na<1 .or. n<1) then
       return
    else if (n==1) then
       c(1) = 1
       return
    end if

    ! Find residual of Nth term with respect to x0 
    do ia=1,na
       do ix=1,3
          xa(ix,ia,n) = xa(ix,ia,n) - x0(ix,ia)
       enddo
    enddo

    ! Find first differences
    ! As in alternate method in KF
    do j=1,n-1
       do ia=1,na
          do ix=1,3
             xa(ix,ia,j) = xa(ix,ia,j+1) - xa(ix,ia,j)
          enddo
       enddo
    enddo

    ! Find matrix s of dot products.

    do j = 1,n-1
       do i = j,n-1
          s(i,j) = sum( xa(:,:,i) * xa(:,:,j) )
          s(j,i) = s(i,j)
       end do
       ! And right-hand side
       rhs(j) = -sum( xa(:,:,j) * xa(:,:,n) )
    end do

    if (Node == 0) then
       ! call print_mat(s,n+1)
    endif

    call solve_with_svd(s,rhs,beta,info,rank_out=rank,sigma=sigma)
    if (node==0) then
       print *, 'matrix rank: ', rank
       print "(a,/,6g12.5)", 'Singular values: ', sigma(:)
       print "(a,/,6g12.5)", 'Beta coeffs: ', beta(:)
    endif
    if (info /= 0) then
       if (node==0) then
          print *, 'The SVD algorithm failed to converge'
          print *, 'Re-using last item only...'
          print *, "info: ", info
       endif
       c(:) = 0.0_dp
       c(n) = 1.0_dp
    endif

    ! Re-shuffle coefficients
    c(n) = 1.0_dp
    c(1) = -beta(1)
    do j = 2, n-1
       c(j) = beta(j-1)-beta(j)
    enddo

    deallocate(s,beta,rhs,sigma)
  end SUBROUTINE extrapolate_diis_svd_new

  subroutine print_mat(a,n)
    integer, intent(in)  :: n
    real(dp), intent(in) :: a(n,n)
    integer :: i, j

    print *, "mat:"
    do i = 1, n
       print "(6g15.7)", (a(i,j),j=1,n)
    enddo
    print *, "------"
  end subroutine print_mat

end module m_new_dm
