! Also the mixing container
module m_mixing_scf

  use class_Fstack_dData1D
  use m_mixing, only: tMixer

  implicit none

  private
  save

  type(tMixer), pointer :: scf_mixs(:) => null()
  type(tMixer), pointer :: scf_mix => null()


  ! Default mixing, no discrepancy between spin-components
  integer, parameter :: MIX_SPIN_ALL = 1
  ! Only use spinor components for mixing
  integer, parameter :: MIX_SPIN_SPINOR = 2
  ! Only use spin-sum for mixing (implicit on spinor)
  integer, parameter :: MIX_SPIN_SUM = 3
  ! Use both spin-sum and spin-difference density for mixing (implicit on spinor)
  integer, parameter :: MIX_SPIN_SUM_DIFF = 4
  ! It makes little sense to only mix difference as for spin-polarised
  ! calculations with no difference it will converge immediately

  ! How the spin mixing algorthim is chosen
  integer :: mix_spin = MIX_SPIN_ALL 

  public :: scf_mixs, scf_mix

  public :: mix_spin
  public :: MIX_SPIN_ALL, MIX_SPIN_SPINOR, MIX_SPIN_SUM, MIX_SPIN_SUM_DIFF

  public :: mixers_scf_init
  public :: mixers_scf_print, mixers_scf_print_block
  public :: mixers_scf_history_init
  public :: mixers_scf_reset
  
  public :: mixing_scf_converged
  
contains

  subroutine mixers_scf_init( nspin, Comm )

    use fdf
    use precision, only: dp
#ifdef MPI
    use mpi_siesta, only: MPI_Comm_World
#endif
    use m_mixing, only: mixers_reset, mixers_init
    use m_mixing, only: mix_method, mix_method_variant
    use m_mixing, only: mixer_init
    use m_mixing, only: mixers_history_init

    ! The number of spin-components
    integer, intent(in) :: nspin
    
    ! The communicator used for the mixer
    integer, intent(in), optional :: Comm

    ! Block constructs
    type(block_fdf) :: bfdf

    ! Get number of history steps
    integer :: n_hist, n_kick, n_restart, n_save
    real(dp) :: w, w_kick
    integer :: n_lin_after
    real(dp) :: w_lin_after
    logical :: lin_after

    ! number of history steps saved
    type(tMixer), pointer :: m
    integer :: nm, im, im2, tmp
    logical :: is_broyden
    character(len=70) :: method, variant, opt
    
    ! If the mixers are denoted by a block, then
    ! the entire logic *MUST* be defined in the blocks
    
    opt = fdf_get('SCF.Mix.Spin','all')
    if      ( leqi(opt, 'all') ) then
       mix_spin = MIX_SPIN_ALL
    else if ( leqi(opt, 'spinor') ) then
       mix_spin = MIX_SPIN_SPINOR
    else if ( leqi(opt, 'sum') ) then
       mix_spin = MIX_SPIN_SUM
    else if ( leqi(opt, 'sum+diff') ) then
       mix_spin = MIX_SPIN_SUM_DIFF
    else
       call die("Unknown option given for SCF.Mix.Spin &
            &all|spinor|sum|sum+diff")
    end if
    ! If there is only one spinor we should mix all...
    if ( nspin == 1 ) mix_spin = MIX_SPIN_ALL

    ! Initialize to ensure debug stuff read
    call mixers_init('SCF', scf_mixs, Comm = Comm )

    ! Check for existance of the SCF.Mix block
    if ( associated(scf_mixs) ) then

       if ( size(scf_mixs) > 0 ) then
          return
       end if

       ! Something has gone wrong...
       ! The user has supplied a block, but
       ! haven't added any content to the block...
       ! However, we fall-back to the default mechanism
       
    end if

    ! ensure nullification
    call mixers_reset(scf_mixs)


    ! >>>*** FIRST ***<<<
    ! Read in compatibility options
    
    ! Figure out if we are dealing with
    ! Broyden or Pulay
    n_hist = fdf_get('DM.NumberPulay', 2)
    tmp = fdf_get('DM.NumberBroyden', 0)
    is_broyden = tmp > 0
    if ( is_broyden ) then
       n_hist = tmp
    end if

    ! Define default mixing weight (used for 
    ! Pulay, Broyden and linear mixing)
    w      = fdf_get('DM.MixingWeight', 0.25_dp)

    ! Default kick-options
    n_kick = fdf_get('DM.NumberKick', 0)
    w_kick = fdf_get('DM.KickMixingWeight', 0.5_dp)

    lin_after = fdf_get('SCF.LinearMixingAfterPulay', .false.)
    w_lin_after = fdf_get('SCF.MixingWeightAfterPulay', w)
    ! >>>*** END ***<<<


    ! Read options in new format

    ! Get history length
    n_hist = fdf_get('SCF.Mixer.History', n_hist)

    ! update mixing weight and kick mixing weight
    w      = fdf_get('SCF.Mixer.Weight',w)
    n_kick = fdf_get('SCF.Mixer.Kick',n_kick)
    w_kick = fdf_get('SCF.Mixer.Kick.Weight',w_kick)

    ! Restart after this number of iterations
    n_restart = fdf_get('SCF.Mixer.Restart', 0)
    n_save    = fdf_get('SCF.Mixer.Restart.Save', 1)
    ! negative savings are not allowed
    n_save = max(0, n_save)

    ! Get the variant of the mixing method
    if ( is_broyden ) then
       method = 'Broyden'
    else if ( n_hist > 0 ) then
       method = 'Pulay'
    else
       method = 'Linear'
    end if
    method = fdf_get('SCF.Mixer.Method', trim(method))
    variant = fdf_get('SCF.Mixer.Variant', 'original')

    
    ! Determine whether linear mixing should be
    ! performed after the "advanced" mixing
    n_lin_after = fdf_get('SCF.Mixer.Linear.After', -1)
    w_lin_after = fdf_get('SCF.Mixer.Linear.After.Weight', w_lin_after)
    
    ! Determine total number of mixers
    nm = 1
    if ( n_lin_after >= 0 .or. lin_after ) nm = nm + 1
    if ( n_kick > 0 ) nm = nm + 1

    ! Initiailaze all mixers
    allocate(scf_mixs(nm))
    scf_mixs(:)%w = w
    scf_mixs(:)%n_hist = n_hist
    scf_mixs(:)%restart = n_restart
    scf_mixs(:)%restart_save = n_save
    
    ! 1. Current mixing index
    im = 1
    ! Store the advanced mixer index (for references to
    ! later mixers)
    im2 = im
    m => scf_mixs(im)
    m%name = method
    m%m = mix_method(method)
    m%v = mix_method_variant(m%m, variant)

    ! 2. Setup the linear mixing after the actual mixing
    if ( n_lin_after > 0 .or. lin_after ) then
       im = im + 1
       m => scf_mixs(im)
       ! Signal to switch to this mixer after
       ! convergence
       scf_mixs(im2)%next_conv => m
       m%name = 'Linear-After'
       m%m = mix_method('linear')
       m%w = w_lin_after
       m%n_itt = n_lin_after
       ! jump back to previous after having run a
       ! few iterations
       m%next => scf_mixs(im2)

    end if
    
    ! In case we have a kick, apply the kick here
    ! This overrides the "linear.after" option
    if ( n_kick > 0 ) then
       
       im = im + 1
       m => scf_mixs(im)
       m%name = 'Linear-Kick'
       m%n_itt = 1
       m%n_hist = 0
       m%m = mix_method('linear')
       m%w = w_kick
       m%next => scf_mixs(im2)

       ! set the default mixer to kick
       scf_mixs(im2)%n_itt = n_kick - 1
       scf_mixs(im2)%next => m
       scf_mixs(im2)%restart = n_kick - 1

    end if

    ! Correct the input
    do im = 1 , nm
       call mixer_init( scf_mixs(im) )
    end do

    ! Initialize the allocation of each mixer
    call mixers_history_init( scf_mixs )
    
#ifdef MPI
    if ( present(Comm) ) then
       scf_mixs(:)%Comm = Comm
    else
       scf_mixs(:)%Comm = MPI_Comm_World
    end if
#endif
    
  end subroutine mixers_scf_init

  subroutine mixers_scf_print( nspin )
    use parallel, only: IONode
    use m_mixing, only: mixers_print
    integer, intent(in) :: nspin

    ! Print mixing options
    call mixers_print( 'SCF' , scf_mixs )
    
    if ( IONode .and. nspin > 1 ) then
       select case ( mix_spin )
       case ( MIX_SPIN_ALL )
          write(*, '(a,t50,a)') 'mix.SCF: Spin-component mixing','all'
       case ( MIX_SPIN_SPINOR )
          write(*, '(a,t50,a)') 'mix.SCF: Spin-component mixing','spinor'
          if ( nspin <= 2 ) then
             call die("SCF.Mixer.Spin spinor option only valid for &
                  &non-collinear and spin-orbit calculations")
          end if
       case ( MIX_SPIN_SUM )
          write(*, '(a,t50,a)') 'mix.SCF: Spin-component mixing','sum'
       case ( MIX_SPIN_SUM_DIFF ) 
          write(*, '(a,t50,a)') 'mix.SCF: Spin-component mixing','sum and diff'
       end select
    end if
    
  end subroutine mixers_scf_print

  subroutine mixers_scf_print_block( )
    use m_mixing, only: mixers_print_block

    ! Print mixing options
    call mixers_print_block( 'SCF' , scf_mixs )
    
  end subroutine mixers_scf_print_block

  subroutine mixing_scf_converged( SCFconverged )

    use parallel, only: IONode

    logical, intent(inout) :: SCFconverged
    integer :: i

    ! Return if no convergence
    if ( .not. SCFconverged ) return

    if ( associated(scf_mix%next_conv) ) then

       ! this means that we skip to the 
       ! following algorithm
       scf_mix => scf_mix%next_conv
       SCFconverged = .false.

       if ( allocated(scf_mix%stack) ) then

          do i = 1 , size(scf_mix%stack)

             ! delete all but one history
             ! This should be fine
             call reset(scf_mix%stack(i), -1)

          end do

       end if

       if ( IONode ) then
         write(*,'(/,2a)') ':!: SCF cycle continuation mixer: ', &
              trim(scf_mix%name)
       end if

    end if

  end subroutine mixing_scf_converged


  subroutine mixers_scf_reset()
    use m_mixing, only: mixers_reset

    nullify(scf_mix)
    call mixers_reset( scf_mixs )

  end subroutine mixers_scf_reset


  subroutine mixers_scf_history_init( )
    use m_mixing, only: mixers_history_init
    
    call mixers_history_init( scf_mixs )
    scf_mix => scf_mixs(1)
    
  end subroutine mixers_scf_history_init

end module m_mixing_scf
