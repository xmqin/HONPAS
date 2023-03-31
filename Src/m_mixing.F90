
! Module for all mixing methods in a standard way

! This module implements mixing of the Pulay and Broyden
! type.

! The Pulay method is implemented in the fast calculation
! setup and in the stable method.
! The stable method is executed if the inversion fails.
!  - Stable: G.Kresse and J.Furthmuller, Comp. Mat. Sci. 6, 15, 1996
!  - gr (guarenteed-reduction) : http://arxiv.org/pdf/cond-mat/0005521.pdf

! All implemented methods employ a restart with variable
! history saving.

module m_mixing
  
  use precision, only: dp

#ifdef MPI
  ! MPI stuff
  use mpi_siesta
#endif

  ! Intrinsic classes for retaining history
  use class_dData1D
  use class_Fstack_dData1D

  implicit none
  
  private

  save

  integer, parameter :: MIX_LINEAR = 1
  integer, parameter :: MIX_PULAY = 2
  integer, parameter :: MIX_BROYDEN = 3
  integer, parameter :: MIX_FIRE = 4

  ! Action tokens (binary: 0, 1, 2, 4, 8, ...!)
  integer, parameter :: ACTION_MIX = 0
  integer, parameter :: ACTION_RESTART = 1
  integer, parameter :: ACTION_NEXT = 2

  type tMixer

     ! Name of mixer
     character(len=24) :: name

     ! The different saved variables per iteration
     ! and their respective stacks
     type(Fstack_dData1D), allocatable :: stack(:)
     
     ! The method of the mixer
     integer :: m = MIX_PULAY
     
     ! In case the mixing method has a variant
     ! this denote the variant
     ! This value is thus specific for each method
     integer :: v = 0

     ! The currently reached iteration
     integer :: cur_itt = 0, start_itt = 0

     ! Different mixers may have different histories
     integer :: n_hist = 2
     
     ! Number of iterations using this mixer
     ! There are a couple of signals here
     !  == 0 :
     !     only use this mixer until convergence
     !   > 0 :
     !     after having runned n_itt step to "next"
     integer :: n_itt = 0

     ! When mod(cur_itt,restart_itt) == 0 the history will
     ! be _reset_
     integer :: restart = 0
     integer :: restart_save = 0

     ! This is an action token specifying the current
     ! action
     integer :: action = ACTION_MIX

     ! The next mixing method following this method
     type(tMixer), pointer :: next => null()

     ! The next mixing method following this method
     ! Only used if mixing method achieved convergence
     ! using this method
     type(tMixer), pointer :: next_conv => null()

     ! ** Parameters specific for the method:

     ! The mixing parameter used for this mixer
     real(dp) :: w = 0._dp

     ! linear array of real variables used specifically
     ! for this mixing type
     real(dp), pointer :: rv(:) => null()
     integer, pointer :: iv(:) => null()

#ifdef MPI
     ! In case we have MPI the mixing scheme
     ! can implement a reduction scheme.
     ! This can be MPI_Comm_Self to not employ any
     ! reductions
     integer :: Comm = MPI_Comm_Self
#endif

  end type tMixer

  
  ! Indices for special constanst
  integer, parameter :: I_PREVIOUS_RES = 0
  integer, parameter :: I_P_RESTART = -1
  integer, parameter :: I_P_NEXT = -2
  ! This index should always be the lowest index
  ! This is used to allocate the correct bounds for the
  ! additional array of information
  integer, parameter :: I_SVD_COND = -3


  ! Debug mixing runs
  logical :: debug_mix = .false.
  ! In case of parallel mixing this also contains the node number
  character(len=20) :: debug_msg = 'mix:'
  
  public :: tMixer

  ! Routines are divided in three sections

  ! 1. Routines used to construct the mixers
  !    Routines used to print information regarding
  !    the mixers
  public :: mixers_init, mixer_init
  public :: mixers_print, mixers_print_block
  public :: mixers_history_init
  public :: mixers_reset

  ! 2. Public functions for retrieving information
  !    from external routines
  public :: mix_method, mix_method_variant

  ! 3. Actual mixing methods
  public :: mixing

  public :: MIX_LINEAR, MIX_FIRE, MIX_PULAY, MIX_BROYDEN

  interface mixing
     module procedure mixing_1d, mixing_2d
  end interface mixing

contains

  !> Initialize a set of mixers by reading in fdf information.
  !! @param[in] prefix the fdf-label prefixes
  !! @param[pointer] mixers the mixers that are to be initialized
  !! @param[in] Comm @opt optional MPI-communicator
  subroutine mixers_init( prefix, mixers, Comm )

    use parallel, only: IONode, Node
    use fdf
    
    ! FDF-prefix for searching keywords
    character(len=*), intent(in) :: prefix
    ! The array of mixers (has to be nullified upon entry)
    type(tMixer), pointer :: mixers(:)
    integer, intent(in), optional :: Comm

    ! Block constructs
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline

    ! number of history steps saved
    integer :: n_hist, n_restart, n_save
    real(dp) :: w

    integer :: nm, im, im2
    character(len=10) :: lp
    character(len=70) :: method, variant

    ! Default mixing options...
    if ( fdf_get('Mixer.Debug',.false.) ) then
       debug_mix = IONode
       debug_msg = 'mix:'
    end if
    if ( fdf_get('Mixer.Debug.MPI',.false.) ) then
       debug_mix = .true.
       write(debug_msg,'(a,i0,a)') 'mix (',Node,'):'
    end if

    lp = trim(prefix)//'.Mixer'

    ! ensure nullification
    call mixers_reset(mixers)

    ! Return immediately if the user hasn't defined
    ! an fdf-block for the mixing options...
    if ( .not. fdf_block(trim(lp)//'s', bfdf) ) return

    ! update mixing weight and kick mixing weight
    w      = fdf_get(trim(lp)//'.Weight',0.1_dp)
    ! Get history length
    n_hist = fdf_get(trim(lp)//'.History',6)
    ! Restart after this number of iterations
    n_restart = fdf_get(trim(lp)//'.Restart',0)
    n_save    = fdf_get(trim(lp)//'.Restart.Save',1)
    ! negative savings are not allowed
    n_save = max(0, n_save)



    ! Read in the options regarding the mixing options
    nm = 0
    do while ( fdf_bline(bfdf,pline) ) 
       if ( fdf_bnnames(pline) == 0 ) cycle
       nm = nm + 1
    end do
    if ( nm == 0 ) then
       call die('mixing: No mixing schemes selected. &
            &Please at least add one mixer.')
    end if

    
    ! Allocate all denoted mixers...
    allocate(mixers(nm))
    mixers(:)%w = w
    mixers(:)%n_hist = n_hist
    mixers(:)%restart = n_restart
    mixers(:)%restart_save = n_save

    
    ! Rewind to grab names.
    call fdf_brewind(bfdf)
    nm = 0
    do while ( fdf_bline(bfdf,pline) ) 
       if ( fdf_bnnames(pline) == 0 ) cycle
       
       nm = nm + 1
       mixers(nm)%name = fdf_bnames(pline,1)
       
    end do
    
    ! Now read all mixers for this segment and their options
    do im = 1 , nm
       
       call read_block( mixers(im) )
       
    end do
    
    ! Create history stack and associate correct 
    ! stack pointers
    call mixers_history_init(mixers)

#ifdef MPI
    if ( present(Comm) ) then
       mixers(:)%Comm = Comm
    else
       mixers(:)%Comm = MPI_Comm_World
    end if
#endif
    
  contains

    subroutine read_block( m )
      type(tMixer), intent(inout), target :: m

      character(len=64) :: opt
      
      ! create block string
      opt = trim(lp)//'.'//trim(m%name)

      if ( .not. fdf_block(opt,bfdf) ) then
         call die('Block: '//trim(opt)//' does not exist!')
      end if

      ! Default to the pulay method...
      ! This enables NOT writing this in the block
      method = 'pulay'
      variant = ' '
      
      ! read method
      do while ( fdf_bline(bfdf,pline) )
         if ( fdf_bnnames(pline) == 0 ) cycle
         
         opt = fdf_bnames(pline,1)
         
         if ( leqi(opt,'method') ) then
            
            method = fdf_bnames(pline,2)
            
         else if ( leqi(opt,'variant') ) then
            
            variant = fdf_bnames(pline,2)
            
         end if
         
      end do
      
      ! Retrieve the method and the variant
      m%m = mix_method(method)
      m%v = mix_method_variant(m%m, variant)

      ! Define separate defaults which are
      ! not part of the default input options
      select case ( m%m )
      case ( MIX_LINEAR )
         m%n_hist = 0
      end select

      call fdf_brewind(bfdf)

      ! read options
      do while ( fdf_bline(bfdf,pline) )
         if ( fdf_bnnames(pline) == 0 ) cycle
         
         opt = fdf_bnames(pline,1)
         
         if ( leqi(opt,'iterations') &
              .or. leqi(opt,'itt') ) then

            m%n_itt = fdf_bintegers(pline,1)
            
         else if ( leqi(opt,'history') ) then
            
            m%n_hist = fdf_bintegers(pline,1)

         else if ( leqi(opt,'weight') .or. leqi(opt, 'w') ) then
            
            m%w = fdf_breals(pline,1)
            
         else if ( leqi(opt,'restart') ) then
            
            m%restart = fdf_bintegers(pline,1)

         else if ( leqi(opt,'restart.save') ) then
            
            m%restart_save = fdf_bintegers(pline,1)
            m%restart_save = max(0,m%restart_save)

         end if
         
      end do

      ! Initialize the mixer by setting the correct
      ! standard options and allocate space in the mixers...
      call mixer_init( m )

      ! Read the options for this mixer
      call fdf_brewind(bfdf)
      
      ! read options
      do while ( fdf_bline(bfdf,pline) )
         if ( fdf_bnnames(pline) == 0 ) cycle
         
         opt = fdf_bnames(pline,1)
         
         if ( leqi(opt,'next') ) then
            
            nullify(m%next)
            
            opt = fdf_bnames(pline,2)
            do im2 = 1 , nm
               if ( leqi(opt,mixers(im2)%name) ) then
                  m%next => mixers(im2)
                  exit
               end if
            end do
            
            if ( .not. associated(m%next) ) then
               call die('mixing: Could not find next mixer. &
                    &Ensure all mixers exist and their names.')
            end if

            if ( associated(m%next, target=m) ) then
               call die('mixing: Next *must* not be it-self. &
                    &Please change accordingly.')
            end if
            
         else if ( leqi(opt,'next.conv') ) then
            
            nullify(m%next_conv)
            
            opt = fdf_bnames(pline,2)
            do im2 = 1 , nm
               if ( leqi(opt,mixers(im2)%name) ) then
                  m%next_conv => mixers(im2)
                  exit
               end if
            end do
            
            if ( .not. associated(m%next_conv) ) then
               call die('mixing: Could not find next convergence mixer. &
                    &Ensure all mixers exist and their names.')
            end if

            if ( associated(m%next_conv,target=m) ) then
               call die('mixing: next.conv *must* not be it-self. &
                    &Please change accordingly.')
            end if

         end if

      end do

      ! Ensure that if a next have not been specified
      ! it will continue indefinitely.
      if ( .not. associated(m%next) ) then
         m%n_itt = 0
      end if

      
      ! Read the options for this mixer
      call fdf_brewind(bfdf)
      
      ! read options
      do while ( fdf_bline(bfdf,pline) )
         ! skip lines without associated content
         if ( fdf_bnnames(pline) == 0 ) cycle
         
         opt = fdf_bnames(pline,1)

         ! Do options so that a pulay option may refer to
         ! the actual names of the constants
         
         if ( m%m == MIX_PULAY ) then
           if ( leqi(opt,'weight.linear') &
               .or. leqi(opt,'w.linear') ) then
             
             m%rv(1) = fdf_breals(pline,1)
             if ( m%rv(1) <= 0._dp .or. 1._dp < m%rv(1) ) then
               call die("m_mixing: Mixing weight should be 0 < weight <= 1")
             end if
             
           else if ( leqi(opt,'svd.cond') ) then
             m%rv(I_SVD_COND) = fdf_bvalues(pline,1)
           end if
           
         else if ( m%m == MIX_BROYDEN ) then
           
           ! The linear mixing weight
           if ( leqi(opt,'weight.linear') &
               .or. leqi(opt,'w.linear') ) then
             
             m%rv(1) = fdf_breals(pline,1)
             if ( m%rv(1) <= 0._dp .or. 1._dp < m%rv(1) ) then
               call die("m_mixing: Linear mixing weight should be 0 < weight <= 1")
             end if

           else if ( leqi(opt,'weight.prime') &
               .or. leqi(opt,'w.prime') ) then
             
             m%rv(2) = fdf_breals(pline,1)
             if ( m%rv(2) < 0._dp .or. 1._dp < m%rv(2) ) then
               call die("m_mixing: Prime weight should be 0 <= weight <= 1")
             end if

           else if ( leqi(opt,'svd.cond') ) then
             
             m%rv(I_SVD_COND) = fdf_bvalues(pline,1)
             
           end if
           
         end if
         
         ! Generic options for all advanced methods...
         if ( leqi(opt,'next.p') ) then

            ! Only allow stepping to the next when
            ! having a next associated
            if ( associated(m%next) ) then
               m%rv(I_P_NEXT) = fdf_bvalues(pline,1)
            end if

         else if ( leqi(opt,'restart.p') ) then

            m%rv(I_P_RESTART) = fdf_bvalues(pline,1)

         end if

      end do

    end subroutine read_block

  end subroutine mixers_init


  !> Initialize a single mixer depending on the preset
  !! options. Useful for external correct setup.
  !!
  !! @param[inout] mix mixer to be initialized
  subroutine mixer_init( mix )
    type(tMixer), intent(inout) :: mix
    integer :: n

    if ( mix%w <= 0._dp .or. 1._dp < mix%w ) then
      call die("m_mixing: Mixing weight should be: 0 < weight <= 1")
    end if

    ! Correct amount of history in the mixing.
    if ( 0 < mix%restart .and. &
         mix%restart < mix%n_hist ) then
       ! This is if we restart this scheme,
       ! then it does not make sense to have a history
       ! greater than the restart count
       mix%n_hist = mix%restart
    end if
    if ( 0 < mix%n_itt .and. &
         mix%n_itt < mix%n_hist ) then
       ! If this only runs for n_itt itterations,
       ! it makes no sense to have a history greater
       ! than this.
       mix%n_hist = mix%n_itt
    end if

    select case ( mix%m )
    case ( MIX_LINEAR )
       
       allocate(mix%rv(I_SVD_COND:0))
       ! Kill any history settings that do not apply to the
       ! linear mixer.
       mix%restart = 0
       mix%restart_save = 0
              
    case ( MIX_PULAY )
       
       allocate(mix%rv(I_SVD_COND:1))
       mix%rv(1) = mix%w
       ! We allocate the double residual (n_hist-1)
       mix%n_hist = max(2, mix%n_hist)
       if ( mix%v == 1 .or. mix%v == 3 ) then
          
         ! The GR method requires an even number
         ! of restart steps
         ! And then we ensure the history to be aligned
         ! with a restart (restart has precedence)
         if ( mix%restart /= 0 ) then
           mix%restart = mix%restart + mod(mix%restart, 2)
         end if

       end if
       
    case ( MIX_BROYDEN )
       
       ! allocate temporary array
       mix%n_hist = max(2, mix%n_hist)
       n = 2 + mix%n_hist
       allocate(mix%rv(I_SVD_COND:n))
       mix%rv(1:n) = mix%w

    end select

    if ( mix%restart < 0 ) then
       call die('mixing: restart count must be positive')
    end if
    
    mix%restart_save = min(mix%n_hist - 1, mix%restart_save)
    mix%restart_save = max(0, mix%restart_save)

    ! This is the restart parameter
    ! I.e. if |f_k / f - 1| < rp
    ! only works for positive rp
    mix%rv(I_PREVIOUS_RES) = huge(1._dp)
    mix%rv(I_P_RESTART) = -1._dp
    mix%rv(I_P_NEXT) = -1._dp
    mix%rv(I_SVD_COND) = 1.e-8_dp
    
  end subroutine mixer_init


  !> Initialize all history for the mixers
  !!
  !! Routine for clearing all history and setting up the
  !! arrays so that they may be used subsequently.
  !!
  !! @param[inout] mixers the mixers to be initialized
  subroutine mixers_history_init( mixers )
    type(tMixer), intent(inout), target :: mixers(:)

    type(tMixer), pointer :: m
    integer :: im, is, ns
    logical :: is_GR

    do im = 1 , size(mixers)
       m => mixers(im)
       if ( debug_mix .and. current_itt(m) >= 1 ) then
          write(*,'(a,a)') trim(debug_msg), &
               ' resetting history of all mixers'
          exit
       end if
    end do

    ! Clean up all arrays and reference counted
    ! objects
    do im = 1 , size(mixers)

       m => mixers(im)

       ! reset history track
       m%start_itt = 0
       m%cur_itt = 0

       ! do not try and de-allocate something not
       ! allocated
       if ( allocated(m%stack) ) then

          ns = size(m%stack)
          do is = 1 , ns
             call delete(m%stack(is))
          end do
          
          ! clean-up
          deallocate(m%stack)
          
       end if
       
       ! Re-populate 
       select case ( m%m )
       case ( MIX_LINEAR )

          ! do nothing

       case ( MIX_PULAY )

          is_GR = (m%v == 1) .or. (m%v == 3)

          if ( .not. is_GR ) then
             allocate(m%stack(3))
          else
             allocate(m%stack(2))
          end if

          ! These arrays contains these informations
          !   s1 = m%stack(1)
          !   s2 = m%stack(2)
          !   s3 = m%stack(3)
          ! Here <> is input function, x[in], and
          ! <>' is the corresponding output, x[out].
          ! First iteration:
          !   s1 = { 1' - 1 }
          !   s3 = { 1' }
          ! Second iteration
          !   s2 = { 2' - 2 - (1' - 1) }
          !   s1 = { 2 - 1 , 2' - 2 }
          !   s3 = { 2' }
          ! Third iteration
          !   s2 = { 2' - 2 - (1' - 1) , 3' - 3 - (2' - 2) }
          !   s1 = { 2 - 1 , 3 - 2, 3' - 3 }
          !   s3 = { 3' }
          ! and so on

          ! allocate x[i+1] - x[i]
          call new(m%stack(1), m%n_hist)
          ! allocate F[i+1] - F[i]
          call new(m%stack(2), m%n_hist-1)
          
          if ( .not. is_GR ) then
             call new(m%stack(3), 1)
          end if

       case ( MIX_BROYDEN )

          ! Same as original Pulay
          allocate(m%stack(3))
          call new(m%stack(1), m%n_hist)
          call new(m%stack(2), m%n_hist-1)
          call new(m%stack(3), 1)

       end select

    end do

  end subroutine mixers_history_init


  !> Reset the mixers, i.e. clean _everything_
  !!
  !! Also deallocates (and nullifies) the input array!
  !!
  !! @param[inout] mixers array of mixers to be cleaned
  subroutine mixers_reset( mixers )
    type(tMixer), pointer :: mixers(:)
    
    type(tMixer), pointer :: m

    integer :: im, is, ns

    if ( .not. associated(mixers) ) return

    do im = 1 , size(mixers)
       m => mixers(im)
       
       if ( allocated(m%stack) ) then
          ns = size(m%stack)
          do is = 1 , ns
             call delete(m%stack(is))
          end do
          deallocate(m%stack)
       end if
       
       if ( associated(m%rv) ) then
          deallocate(m%rv)
          nullify(m%rv)
       end if
       
       if ( associated(m%iv) ) then
          deallocate(m%iv)
          nullify(m%iv)
       end if
       
    end do

    deallocate(mixers)
    nullify(mixers)

  end subroutine mixers_reset
  

  !> Print (to std-out) information regarding the mixers
  !!
  !! @param[in] prefix the prefix (fdf) for the mixers
  !! @param[in] mixers array of mixers allocated
  subroutine mixers_print( prefix, mixers )

    use parallel, only: IONode
    
    character(len=*), intent(in) :: prefix
    type(tMixer), intent(in), target :: mixers(:)

    type(tMixer), pointer :: m
    character(len=64) :: fmt

    logical :: bool
    integer :: i

    if ( .not. IONode ) return

    fmt = 'mix.'//trim(prefix)//':'

    if ( debug_mix ) then
       write(*,'(2a,t50,''= '',l)') trim(fmt), &
            ' Debug messages',debug_mix
    end if

    ! Print out options for all mixers
    do i = 1 , size(mixers)

       m => mixers(i)

       select case ( m%m )
       case ( MIX_LINEAR )
          
          write(*,'(2a,t50,''= '',a)') trim(fmt), &
               ' Linear mixing',trim(m%name)
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Mixing weight',m%w

          if ( m%n_hist > 0 .and. (&
               associated(m%next) &
               .or. associated(m%next_conv)) ) then
             write(*,'(2a,t50,''= '',i0)') trim(fmt), &
                  '    Carried history steps',m%n_hist
          end if

       case ( MIX_PULAY )
          
          write(*,'(2a,t50,''= '',a)') trim(fmt), &
               ' Pulay mixing',trim(m%name)

          select case ( m%v )
          case ( 0 )
             write(*,'(2a,t50,''= '',a)') trim(fmt), &
                  '    Variant','stable'
          case ( 1 ) 
             write(*,'(2a,t50,''= '',a)') trim(fmt), &
                  '    Variant','GR'
          case ( 2 )
             write(*,'(2a,t50,''= '',a)') trim(fmt), &
                  '    Variant','stable-SVD'
          case ( 3 )
             write(*,'(2a,t50,''= '',a)') trim(fmt), &
                  '    Variant','GR-SVD'
          end select

          write(*,'(2a,t50,''= '',i0)') trim(fmt), &
               '    History steps',m%n_hist
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Linear mixing weight',m%rv(1)
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Mixing weight',m%w
          write(*,'(2a,t50,''= '',e10.4)') trim(fmt), &
               '    SVD condition',m%rv(I_SVD_COND)
          if ( m%rv(I_P_NEXT) > 0._dp ) then
             write(*,'(2a,t50,''= '',f6.4)') trim(fmt), &
                  '    Step mixer parameter',m%rv(I_P_NEXT)
          end if
          bool = .false.
          if ( m%rv(I_P_RESTART) > 0._dp ) then
             write(*,'(2a,t50,''= '',f6.4)') trim(fmt), &
                  '    Restart parameter',m%rv(I_P_RESTART)
             bool = .true.
          end if
          if ( m%restart > 0 ) then
             write(*,'(2a,t50,''= '',i0)') trim(fmt), &
                  '    Restart steps',m%restart
             bool = .true.
          end if
          if ( bool ) then
             write(*,'(2a,t50,''= '',i0)') trim(fmt), &
                  '    Restart save steps',m%restart_save
          end if

       case ( MIX_BROYDEN )
          
          write(*,'(2a,t50,''= '',a)') trim(fmt), &
               ' Broyden mixing',trim(m%name)

          !write(*,'(2a,t50,''= '',a)') trim(fmt), &
          !     '    Variant','original'

          write(*,'(2a,t50,''= '',i0)') trim(fmt), &
               '    History steps',m%n_hist
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Linear mixing weight',m%rv(1)
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Inverse Jacobian weight',m%w
          write(*,'(2a,t50,''= '',f12.6)') trim(fmt), &
               '    Weight prime',m%rv(2)
          write(*,'(2a,t50,''= '',e10.4)') trim(fmt), &
               '    SVD condition',m%rv(I_SVD_COND)
          if ( m%rv(I_P_NEXT) > 0._dp ) then
             write(*,'(2a,t50,''= '',f6.4)') trim(fmt), &
                  '    Step mixer parameter',m%rv(I_P_NEXT)
          end if
          bool = .false.
          if ( m%rv(I_P_RESTART) > 0._dp ) then
             write(*,'(2a,t50,''= '',f6.4)') trim(fmt), &
                  '    Restart parameter',m%rv(I_P_RESTART)
             bool = .true.
          end if
          if ( m%restart > 0 ) then
             write(*,'(2a,t50,''= '',i0)') trim(fmt), &
                  '    Restart steps',m%restart
             bool = .true.
          end if
          if ( bool ) then
             write(*,'(2a,t50,''= '',i0)') trim(fmt), &
                  '    Restart save steps',m%restart_save
          end if
          
       case ( MIX_FIRE )
          
          write(*,'(2a,t50,''= '',a)') trim(fmt), &
               ' Fire mixing',trim(m%name)

       end select

       if ( m%n_itt > 0 ) then
          write(*,'(2a,t50,''= '',i0)') trim(fmt), &
               '    Number of mixing iterations',m%n_itt
          if ( associated(m%next) ) then
             write(*,'(2a,t50,''= '',a)') trim(fmt), &
                  '    Following mixer',trim(m%next%name)
          else
             call die('Something went wrong, if the mixer does not go &
                  &indefinitely it should have a following method.')
          end if
       end if
       
       if ( associated(m%next_conv) ) then
          write(*,'(2a,t50,''= '',a)') trim(fmt), &
               '    Following mixer upon convergence',trim(m%next_conv%name)
       end if
          
    end do

  end subroutine mixers_print


  !> Print (to std-out) the fdf-blocks that recreate the mixer settings
  !!
  !! @param[in] prefix the fdf-prefix for reading the blocks
  !! @param[in] mixers array of mixers that should be printed
  !!    their fdf-blocks
  subroutine mixers_print_block( prefix, mixers )
    use parallel, only: IONode
    
    character(len=*), intent(in) :: prefix
    type(tMixer), intent(in), target :: mixers(:)

    type(tMixer), pointer :: m

    logical :: bool
    integer :: i

    if ( .not. IONode ) return

    ! Write block of input
    write(*,'(/3a)')'%block ',trim(prefix), '.Mixers'
    do i = 1 , size(mixers)
       m => mixers(i)
       write(*,'(t3,a)') trim(m%name)
    end do
    write(*,'(3a)')'%endblock ',trim(prefix), '.Mixers'


    ! Print out options for all mixers
    do i = 1 , size(mixers)

       m => mixers(i)

       ! Write out this block
       write(*,'(/4a)')'%block ',trim(prefix),'.Mixer.',trim(m%name)

       write(*,'(t3,a)') '# Mixing method'

       ! Write out method
       select case ( m%m )
       case ( MIX_LINEAR )
          
          write(*,'(t2,2(tr1,a))') 'method','linear'

       case ( MIX_PULAY )

          write(*,'(t2,2(tr1,a))') 'method','pulay'
          select case ( m%v )
          case ( 0 )
             write(*,'(t2,2(tr1,a))') 'variant','stable'
          case ( 1 )
             write(*,'(t2,2(tr1,a))') 'variant','GR'
          case ( 2 )
             write(*,'(t2,2(tr1,a))') 'variant','stable+SVD'
          case ( 3 )
             write(*,'(t2,2(tr1,a))') 'variant','GR+SVD'
          end select

       case ( MIX_BROYDEN )

          write(*,'(t2,2(tr1,a))') 'method','broyden'

          ! currently no variants exists
          
       end select


       
       ! remark
       write(*,'(/,t3,a)') '# Mixing options'

       ! Weight
       ! For Broyden this is the inverse Jacobian
       write(*,'(t3,a,f6.4)') 'weight ', m%w
       select case ( m%m )
       case ( MIX_PULAY )
         write(*,'(t3,a,f6.4)') 'weight.linear ', m%rv(1)
       case ( MIX_BROYDEN )
         write(*,'(t3,a,f6.4)') 'weight.linear ', m%rv(1)
         write(*,'(t3,a,f6.4)') 'weight.prime ', m%rv(2)
       end select
          
       if ( m%n_hist > 0 ) then
          write(*,'(t3,a,i0)') 'history ', m%n_hist
       end if

       bool = .false.
       if ( m%restart > 0 ) then
          write(*,'(t3,a,i0)') 'restart ', m%restart
          bool = .true.
       end if
       select case ( m%m )
       case ( MIX_PULAY, MIX_BROYDEN )
          if ( m%rv(I_P_RESTART) > 0._dp ) then
             write(*,'(t3,a,e10.5)')'restart.p ', m%rv(I_P_RESTART)
             bool = .true.
          end if
       end select
       if ( bool ) then
          write(*,'(t3,a,i0)')'restart.save ', m%restart_save
       end if



       ! remark

       bool = .false.
       if ( m%n_itt > 0 ) then
          write(*,'(/,t3,a)') '# Continuation options'
          write(*,'(t3,a,i0)')'iterations ', m%n_itt
          bool = .true.
       end if
       select case ( m%m )
       case ( MIX_PULAY , MIX_BROYDEN )
          if ( m%rv(I_P_NEXT) > 0._dp ) then
             if ( .not. bool ) &
                  write(*,'(/,t3,a)') '# Continuation options'
             write(*,'(t3,a,f6.4)')'next.p ', m%rv(I_P_NEXT)
             bool = .true.
          end if
       end select
       if ( bool .and. associated(m%next) ) then
          write(*,'(t2,2(tr1,a))') 'next', trim(m%next%name)
       else if ( bool ) then
          call die('Something went wrong, if the mixer does not go &
               &indefinitely it should have a following method.')
       end if
       
       if ( associated(m%next_conv) ) then
          if ( .not. bool ) &
               write(*,'(/,t3,a)') '# Continuation options'
          write(*,'(t2,2(tr1,a))') 'next.conv', trim(m%next_conv%name)
       end if


       write(*,'(4a)')'%endblock ',trim(prefix),'.Mixer.',trim(m%name)
       
    end do

    write(*,*) ! new-line
    
  end subroutine mixers_print_block


  !> Return the integer specification of the mixing type
  !!
  !! @param[in] str the character representation of the mixing type
  !! @return the integer corresponding to the mixing type
  function mix_method(str) result(m)
    use fdf, only: leqi
    character(len=*), intent(in) :: str
    integer :: m

    if ( leqi(str,'linear') ) then
       m = MIX_LINEAR
    else if ( leqi(str,'pulay') .or. &
         leqi(str,'diis') .or. &
         leqi(str, 'anderson') ) then
       m = MIX_PULAY
    else if ( leqi(str,'broyden') ) then
       m = MIX_BROYDEN
    else if ( leqi(str,'fire') ) then
       m = MIX_FIRE
       call die('mixing: FIRE currently not supported.')
    else
       call die('mixing: Unknown mixing variant.')
    end if

  end function mix_method


  !> Return the variant of the mixing method
  !!
  !! @param[in] m the integer type of the mixing method
  !! @param[in] str the character specification of the mixing method variant
  !! @return the variant of the mixing method
  function mix_method_variant(m, str) result(v)
    use fdf, only: leqi
    integer, intent(in) :: m
    character(len=*), intent(in) :: str
    integer :: v

    v = 0
    select case ( m )
    case ( MIX_LINEAR )
       ! no variants
    case ( MIX_PULAY ) 
       v = 0
       ! We do not implement tho non-stable version
       ! There is no need to have an inferior Pulay mixer...
       if ( leqi(str,'original') .or. &
            leqi(str,'kresse') .or. leqi(str,'stable') ) then
          ! stable version, will nearly always succeed on inversion
          v = 0
       else if ( leqi(str,'original+svd') .or. &
            leqi(str,'kresse+svd') .or. leqi(str,'stable+svd') ) then
          ! stable version, will nearly always succeed on inversion
          v = 2
       else if ( leqi(str,'gr') .or. &
            leqi(str,'guarenteed-reduction') .or. &
            leqi(str,'bowler-gillan') ) then
          ! Guarenteed reduction version
          v = 1
       else if ( leqi(str,'gr+svd') .or. &
            leqi(str,'guarenteed-reduction+svd') .or. &
            leqi(str,'bowler-gillan+svd') ) then
          ! Guarenteed reduction version
          v = 3
       end if
    case ( MIX_BROYDEN )
       ! Currently only one variant
       v = 0
    case ( MIX_FIRE ) 
       ! no variants
    end select

  end function mix_method_variant




  ! The basic mixing procedure is this:
  ! 1. Initialize the mixing algorithm
  !    This will typically mean that one needs to
  !    push input and output matrices
  ! 2. Calculate the mixing coefficients
  ! 3. Use coefficients to calculate the
  !    next (optimized) guess
  ! 4. Finalize the mixing method
  !
  ! Having the routines split up in this manner
  ! allows one to skip step 2 and use coefficients
  ! from another set of input/output to retrieve the
  ! mixing coefficients.
  ! Say we may retrieve mixing coefficients from
  ! the Hamiltonian, but use them for the density-matrix

  !> Initialize the mixing algorithm
  !!
  !! @param[pointer] mix the mixing method
  !! @param[in] n size of the arrays to be used in the algorithm
  !! @param[in] xin array of the input variables
  !! @param[in] xout array of the output variables
  subroutine mixing_init( mix, n, xin, F)
    
    ! The current mixing method 
    type(tMixer), pointer :: mix
    integer, intent(in) :: n
    ! In/out of the function
    real(dp), intent(in) :: xin(n), F(n)

    real(dp), pointer :: res(:), rres(:)

    integer :: i, ns
    real(dp) :: dnorm, dtmp
    logical :: p_next, p_restart

    ! Initialize action for mixer
    mix%action = ACTION_MIX
    
    ! Step iterator (so first mixing has cur_itt == 1)
    mix%cur_itt = mix%cur_itt + 1

    ! If we are going to skip to next, we signal it
    ! before entering
    if ( mix%n_itt > 0 .and. &
         mix%n_itt <= current_itt(mix) ) then
       mix%action = IOR(mix%action, ACTION_NEXT)
    end if


    ! Check whether the residual norm is below a certain
    ! criteria
    p_next = mix%rv(I_P_NEXT) > 0._dp
    p_restart = mix%rv(I_P_RESTART) > 0._dp

    ! Check whether a parameter next/restart is required
    if ( p_restart .or. p_next ) then

       ! Calculate norm: ||f_k||
       dnorm = norm(n, F, F)
#ifdef MPI
       dtmp = dnorm
       call MPI_AllReduce(dtmp,dnorm,1, &
            MPI_double_precision, MPI_Sum, &
            mix%Comm, i)
#endif

       ! Calculate the relative difference
       dtmp = abs(dnorm / mix%rv(I_PREVIOUS_RES) - 1._dp)

       ! We first check for next, that has precedence
       if ( p_next ) then

          if ( dtmp < mix%rv(I_P_NEXT) ) then
             ! Signal stepping mixer
             mix%action = IOR(mix%action, ACTION_NEXT)
          end if
          
          if ( debug_mix .and. current_itt(mix) > 1 ) &
               write(*,'(a,2(a,e8.3))') trim(debug_msg), &
               ' | ||f_k|| - ||f_k-1|| |/||f_k-1|| < np  :  ', &
               dtmp, ' < ', mix%rv(I_P_NEXT)
          
       end if
       
       if ( p_restart ) then

          if ( dtmp < mix%rv(I_P_RESTART) ) then
             ! Signal restart
             mix%action = IOR(mix%action, ACTION_RESTART)
          end if

          if ( debug_mix .and. current_itt(mix) > 1 ) &
               write(*,'(a,2(a,e8.3))') trim(debug_msg), &
               ' | ||f_k|| - ||f_k-1|| |/||f_k-1|| < rp  :  ', &
               dtmp, ' < ', mix%rv(I_P_RESTART)

       end if
       
       ! Store the new residual norm
       mix%rv(I_PREVIOUS_RES) = dnorm
       
    end if


    ! Push information to the stack
    select case ( mix%m )
    case ( MIX_LINEAR )
       if ( debug_mix ) &
            write(*,'(2a)') trim(debug_msg),' linear'
       call init_linear()
    case ( MIX_PULAY )
       if ( debug_mix ) then
          select case ( mix%v )
          case ( 0 )
            write(*,'(2a)') trim(debug_msg),' Pulay'
          case ( 1 )
            write(*,'(2a)') trim(debug_msg),' Pulay, GR'
          case ( 2 )
            write(*,'(2a)') trim(debug_msg),' Pulay-SVD'
          case ( 3 )
            write(*,'(2a)') trim(debug_msg),' Pulay-SVD, GR'
          end select
       end if
       call init_pulay()
    case ( MIX_BROYDEN )
      if ( debug_mix ) &
           write(*,'(2a)') trim(debug_msg),' Broyden'
       call init_broyden()
    end select

  contains

    subroutine init_linear()
      
      ! information for this depends on the
      ! following method
      call fake_history_from_linear(mix%next)
      call fake_history_from_linear(mix%next_conv)

    end subroutine init_linear

    subroutine init_pulay()
      logical :: GR_linear

      select case ( mix%v )
      case ( 0 , 2 ) ! Stable Pulay

         ! Add the residual to the stack
         call push_F(mix%stack(1), n, F)
         
         ns = n_items(mix%stack(1))
         
         ! Add the residuals of the residuals if applicable
         if ( ns > 1 ) then
            
            ! Create F[i+1] - F[i]
            call push_diff(mix%stack(2), mix%stack(1))

            ! Update the residual to reflect the input residual
            res => getstackval(mix, 1, ns-1)
            rres => getstackval(mix, 3)
         
!$OMP parallel do default(shared), private(i)
            do i = 1 , n
               res(i) = res(i) - rres(i) + xin(i)
            end do
!$OMP end parallel do
            
         end if

      case ( 1 , 3 )

         ! Whether this is the linear cycle...
         GR_linear = mod(current_itt(mix), 2) == 1

         ! Add the residual to the stack
         call push_F(mix%stack(1), n, F, mix%rv(1))

         ns = n_items(mix%stack(1))

         if ( GR_linear .and. current_itt(mix) > 1 .and. &
              ns > 1 ) then
            
            res => getstackval(mix, 1)
            rres => getstackval(mix, 2)
!$OMP parallel do default(shared), private(i)
            do i = 1 , n
               rres(i) = rres(i) + res(i)
            end do
!$OMP end parallel do

         else if ( ns > 1 .and. .not. GR_linear ) then

            ! now we can calculate RRes[i]
            call push_diff(mix%stack(2),mix%stack(1))

         end if

      end select
      
    end subroutine init_pulay

    subroutine init_broyden()

      ! Add the residual to the stack
      call push_F(mix%stack(1), n, F)

      ns = n_items(mix%stack(1))

      ! Add the residuals of the residuals if applicable
      if ( ns > 1 ) then

         ! Create F[i+1] - F[i]
         call push_diff(mix%stack(2), mix%stack(1))

         ! Update the residual to reflect the input residual
         res => getstackval(mix, 1, ns-1)
         rres => getstackval(mix, 3)

!$OMP parallel do default(shared), private(i)
         do i = 1 , n
            res(i) = res(i) - rres(i) + xin(i)
         end do
!$OMP end parallel do

      else

         ! Store F[x_in] (used to create the input residual)
         call push_stack_data(mix%stack(3), n)
         res => getstackval(mix, 3)
!$OMP parallel do default(shared), private(i)
         do i = 1 , n
            res(i) = xin(i) + F(i)
         end do
!$OMP end parallel do 

      end if

    end subroutine init_broyden
    
    subroutine fake_history_from_linear(next)
      type(tMixer), pointer :: next
      real(dp), pointer :: t1(:), t2(:)
      integer :: ns, nh, i, nhl

      if ( .not. associated(next) ) return

      ! Reduce to # history of linear
      nhl = mix%n_hist
      ! if the number of fake-history steps saved is
      ! zero we immediately return.
      ! Only if mix%n_hist > 0 will the below
      ! occur.
      if ( nhl == 0 ) return

      ! Check for the type of following method
      select case ( next%m )
      case ( MIX_PULAY )

         ! Here it depends on the variant
         select case ( next%v ) 
         case ( 0 , 2 ) ! stable pulay mixing

            ! Add the residual to the stack
            call push_F(next%stack(1), n, F)

            ns = n_items(next%stack(1))
            
            ! Add the residuals of the residuals if applicable
            if ( ns > 1 ) then

               ! Create F[i+1] - F[i]
               call push_diff(next%stack(2), next%stack(1))

               ! Update the residual to reflect the input residual
               t1 => getstackval(next, 1, ns-1)
               t2 => getstackval(next, 3)

!$OMP parallel do default(shared), private(i)
               do i = 1 , n
                  t1(i) = t1(i) - t2(i) + xin(i)
                  t2(i) = xin(i) + F(i)
               end do
!$OMP end parallel do

            else

               call push_stack_data(next%stack(3), n)
               t1 => getstackval(next, 3)
!$OMP parallel do default(shared), private(i)
               do i = 1 , n
                  t1(i) = xin(i) + F(i)
               end do
!$OMP end parallel do
               
            end if

            ! Clean up the data...
            if ( ns >= nhl ) then
               call reset(next%stack(1), 1)
               call reset(next%stack(2), 1)
               ns = ns - 1
            end if

            nh = max_size(next%stack(1))

            if ( debug_mix ) &
                 write(*,'(a,2(a,i0))') trim(debug_msg), &
                 ' next%n_hist = ', ns, ' / ', nh
            
         end select
         
      case ( MIX_BROYDEN )

         ! Add the residual to the stack
         call push_F(next%stack(1), n, F)

         ns = n_items(next%stack(1))

         ! Add the residuals of the residuals if applicable
         if ( ns >= 2 ) then

            ! Create F[i+1] - F[i]
            call push_diff(next%stack(2), next%stack(1))

            ! Update the residual to reflect the input residual
            t1 => getstackval(next, 1, ns-1)
            t2 => getstackval(next, 3)

!$OMP parallel do default(shared), private(i)
            do i = 1 , n
               t1(i) = t1(i) - t2(i) + xin(i)
               t2(i) = xin(i) + F(i)
            end do
!$OMP end parallel do

         else

            call push_stack_data(next%stack(3), n)
            t1 => getstackval(next, 3)
!$OMP parallel do default(shared), private(i)
            do i = 1 , n
               t1(i) = xin(i) + F(i)
            end do
!$OMP end parallel do

         end if

         ! Clean up the data...
         if ( ns >= nhl ) then
            call reset(next%stack(1), 1)
            call reset(next%stack(2), 1)
            ns = ns - 1
         end if

         nh = max_size(next%stack(1))

         if ( debug_mix ) &
              write(*,'(a,2(a,i0))') trim(debug_msg), &
              ' next%n_hist = ', ns, ' / ', nh
         
      end select

    end subroutine fake_history_from_linear
    
  end subroutine mixing_init

  
  !> Function to retrieve the number of coefficients
  !! calculated in this iteration.
  !! This is so external routines can query the size
  !! of the arrays used.
  !!
  !! @param[in] mix the used mixer
  function mixing_ncoeff( mix ) result(n)
    type(tMixer), intent(in) :: mix
    integer :: n

    n = 0

    select case ( mix%m )
    case ( MIX_PULAY )
       n = n_items(mix%stack(2))
    case ( MIX_BROYDEN )
       n = n_items(mix%stack(2))
    end select

  end function mixing_ncoeff


  !> Calculate the mixing coefficients for the
  !! current mixer
  !!
  !! @param[in] mix the current mixer
  !! @param[in] n the number of elements used to calculate
  !!           the coefficients
  !! @param[in] xin the input value
  !! @param[in] F xout - xin, (residual)
  !! @param[out] coeff the coefficients
  subroutine mixing_coeff( mix, n, xin, F, coeff)

    use parallel, only: IONode

    type(tMixer), intent(inout) :: mix
    integer, intent(in) :: n
    real(dp), intent(in) :: xin(n), F(n)
    real(dp), intent(out) :: coeff(:)
    integer :: ncoeff

    ncoeff = size(coeff)

    if ( ncoeff < mixing_ncoeff(mix) ) then
       write(*,'(a)') 'mix: Error in calculating coefficients'
       ! Do not allow this...
       return
    end if

    select case ( mix%m )
    case ( MIX_LINEAR )
       call linear_coeff()
    case ( MIX_PULAY )
       call pulay_coeff()
    case ( MIX_BROYDEN )
       call broyden_coeff()
    end select

  contains

    subroutine linear_coeff()
      integer :: i
      do i = 1 , ncoeff
         coeff(i) = 0._dp
      end do
    end subroutine linear_coeff

    subroutine pulay_coeff()
      integer :: ns, nh, nmax
      integer :: i, j, info
      logical :: lreturn

      ! Calculation quantities
      real(dp) :: dnorm, G
      real(dp), pointer :: res(:), rres(:), rres1(:), rres2(:)
      real(dp), allocatable :: A(:,:), Ainv(:,:)

      ns = n_items(mix%stack(1))
      nmax = max_size(mix%stack(1))
      nh = n_items(mix%stack(2))

      lreturn = .false.
      
      ! Easy check for initial step...
      select case ( mix%v )
      case ( 0 , 2 ) ! Stable Pulay

        lreturn = ns == 1

      case ( 1 , 3 ) ! Guaranteed Pulay

        lreturn = mod(current_itt(mix), 2) == 1

      end select

      ! In case we return we are actually doing
      ! linear mixing
      if ( lreturn ) return

      
      ! Print out number of currently used history steps
      if ( debug_mix ) &
          write(*,'(a,2(a,i0))') trim(debug_msg), &
          ' n_hist = ',ns, ' / ',nmax

      ! Allocate arrays for calculating the
      ! coefficients
      allocate(A(nh,nh), Ainv(nh,nh))

      
      ! Calculate A_ij coefficients for inversion
      do i = 1 , nh
        
        ! Get RRes[i] array
        rres1 => getstackval(mix, 2, i)

        do j = 1 , i - 1

          ! Get RRes[j] array
          rres2 => getstackval(mix, 2, j)

          ! A(i,j) = A(j,i) = norm(RRes[i],RRes[j])
          A(i,j) = norm(n, rres1, rres2)
          A(j,i) = A(i,j)

        end do

        ! Diagonal
        A(i,i) = norm(n, rres1, rres1)

      end do

#ifdef MPI
      ! Global operations, but only for the non-extended entries
      call MPI_AllReduce(A(1,1),Ainv(1,1),nh*nh, &
          MPI_double_precision, MPI_Sum, &
          mix%Comm, i)
      ! copy over reduced arrays
      A(:,:) = Ainv(:,:)
#endif

      ! Get inverse of matrix
      select case ( mix%v )
      case ( 0 , 1 )

        call inverse(nh, A, Ainv, info)

        if ( info /= 0 ) then

          ! We will first try the SVD routine
          call svd(nh, A, Ainv, mix%rv(I_SVD_COND), j, info)

          ! only inform if we should not use SVD per default
          if ( IONode ) &
              write(*,'(2a,i0,a,i0)') trim(debug_msg), &
              ' Pulay -- inversion failed, > SVD [rank/size] ', &
              j, ' / ', nh

        end if

      case ( 2 , 3 )

        ! We forcefully use the SVD routine
        call svd(nh, A, Ainv, mix%rv(I_SVD_COND), j, info)

      end select

      
      ! NOTE, although mix%stack(1) contains
      ! the x[i] - x[i-1], the tip of the stack
      ! contains F[i]!
      ! res == F[i]
      res => getstackval(mix, 1)

      ! Initialize the coefficients
      do i = 1 , nh
        coeff(i) = 0._dp
      end do

      if ( info == 0 ) then
        
        ! Calculate the coefficients on all processors
        do j = 1 , nh

          ! res  == F[i]
          ! rres == F[j+1] - F[j]
          rres => getstackval(mix, 2, j)
          dnorm = norm(n, rres, res)

          do i = 1 , nh
            coeff(i) = coeff(i) - Ainv(i,j) * dnorm
          end do

        end do

#ifdef MPI
        ! Reduce the coefficients
        call MPI_AllReduce(coeff(1),A(1,1),nh, &
            MPI_double_precision, MPI_Sum, &
            mix%Comm, i)
        do i = 1 , nh
          coeff(i) = A(i,1)
        end do
#endif

      else

        info = 0

        ! reset to linear mixing
        write(*,'(2a)') trim(debug_msg), &
            ' Pulay -- inversion failed, SVD failed, > linear'

      end if

      ! Clean up memory
      deallocate(A, Ainv)
            
    end subroutine pulay_coeff

    subroutine broyden_coeff()
      integer :: ns, nh, nmax
      integer :: i, j, k, info

      ! Calculation quantities
      real(dp) :: dnorm, dtmp
      real(dp), pointer :: w(:), res(:), rres(:), rres1(:), rres2(:)
      real(dp), allocatable :: c(:), A(:,:), Ainv(:,:)

      ns = n_items(mix%stack(1))
      nmax = max_size(mix%stack(1))
      nh = n_items(mix%stack(2))

      ! Easy check for initial step...
      if ( ns == 1 ) then

         ! reset
         coeff = 0._dp

         return
         
      end if

      
      ! Print out number of currently used history steps
      if ( debug_mix ) &
           write(*,'(a,2(a,i0))') trim(debug_msg), &
           ' n_hist = ',ns, ' / ',nmax


      ! This is the modified Broyden algorithm...

      ! Retrieve the previous weights
      w => mix%rv(3:2+nh)
      select case ( mix%v )
      case ( 2 ) ! Unity Broyden
         
         w(nh) = 1._dp
         
      case ( 1 ) ! RMS Broyden
         
         dnorm = norm(n, F, F)
#ifdef MPI
         call MPI_AllReduce(dnorm, dtmp, 1, &
              MPI_Double_Precision, MPI_Max, &
              mix%Comm, i)
         dnorm = dtmp
#endif
         w(:) = 1._dp / sqrt(dnorm)
         if ( debug_mix ) &
              write(*,'(2(a,e10.4))') &
              trim(debug_msg)//' weight = ',w(1), &
              ' , norm = ',dnorm

      case ( 0 ) ! Varying weight
         
         dnorm = 0._dp
!$OMP parallel do default(shared), private(i), &
!$OMP& reduction(max:dnorm)
         do i = 1 , n
            dnorm = max( dnorm, abs(F(i)) )
         end do
!$OMP end parallel do
#ifdef MPI
         call MPI_AllReduce(dnorm, dtmp, 1, &
              MPI_Double_Precision, MPI_Max, &
              mix%Comm, i)
         dnorm = dtmp
#endif
         ! Problay 0.2 should be changed to user-defined
         w(nh) = exp( 1._dp / (dnorm + 0.2_dp) )
         if ( debug_mix ) &
              write(*,'(2a,1000(tr1,e10.4))') &
              trim(debug_msg),' weights = ',w(1:nh)

      end select

      ! Allocate arrays used
      allocate(c(nh))
      allocate(A(nh,nh), Ainv(nh,nh))

      
      !  < RRes[i] | Res[n] >
      do i = 1 , nh
         rres => getstackval(mix, 2, i)
         c(i) = norm(n, rres, F)
      end do
#ifdef MPI
      call MPI_AllReduce(c(1), A(1,1), nh, &
           MPI_Double_Precision, MPI_Sum, &
           mix%Comm, i)
      do i = 1 , nh
         c(i) = A(i,1)
      end do
#endif

      ! Create A_ij coefficients for inversion
      do i = 1 , nh
         
         ! Get RRes[i] array
         rres1 => getstackval(mix, 2, i)
         
         do j = 1 , i -1 

            ! Get RRes[j] array
            rres2 => getstackval(mix, 2, j)

            ! A(i,j) = A(j,i) = dot_product(RRes[i],RRes[j])
            A(i,j) = w(i) * w(j) * norm(n, rres1, rres2)
            A(j,i) = A(i,j)

         end do

         ! Do the diagonal term
         A(i,i) = w(i) * w(i) * norm(n, rres1, rres1)
         
      end do

#ifdef MPI
      call MPI_AllReduce(A(1,1), Ainv(1,1), nh*nh, &
           MPI_double_precision, MPI_Sum, &
           mix%Comm, i)
      A(:,:) = Ainv(:,:)
#endif

      ! Add the diagonal term
      ! This should also prevent it from being
      ! singular (unless mix%rv(2) == 0)
      do i = 1 , nh
         A(i,i) = mix%rv(2) ** 2 + A(i,i)
      end do

      ! Calculate the inverse
      call inverse(nh, A, Ainv, info)

      if ( info /= 0 ) then
                  
         ! We will first try the SVD routine
         call svd(nh, A, Ainv, mix%rv(I_SVD_COND), j, info)

         ! only inform if we should not use SVD per default
         if ( IONode ) &
             write(*,'(2a,i0,a,i0)') trim(debug_msg), &
             ' Broyden -- inversion failed, > SVD [rank/size] ', &
             j, ' / ', nh

      end if

      do i = 1 , nh
         coeff(i) = 0._dp
      end do

      if ( info == 0 ) then

         ! Calculate the coefficients...
         do i = 1 , nh

            do j = 1 , nh
               ! Ainv should be symmetric (A is)
               coeff(i) = coeff(i) + w(j) * c(j) * Ainv(j,i)
            end do

            ! Calculate correct weight...
            coeff(i) = - w(i) * coeff(i)
            
         end do

      else

         ! reset to linear mixing
         write(*,'(2a)') trim(debug_msg), &
              ' Broyden -- inversion failed, SVD failed, > linear'

      end if

      deallocate(A, Ainv)
      
    end subroutine broyden_coeff

  end subroutine mixing_coeff



  !> Calculate the guess for the next iteration
  !!
  !! Note this gets passed the coefficients. Hence,
  !! they may be calculated from another set of history
  !! steps.
  !! This may be useful in certain situations.
  !!
  !! @param[in] mix the current mixer
  !! @param[in] n the number of elements used to calculate
  !!           the coefficients
  !! @param[in] xin the input value
  !! @param[in] F the xin residual
  !! @param[out] xnext the input for the following iteration
  !! @param[in] coeff the coefficients
  subroutine mixing_calc_next( mix, n, xin, F, xnext, coeff )
    ! The current mixing method 
    type(tMixer), pointer :: mix
    integer, intent(in) :: n
    real(dp), intent(in) :: xin(n)
    real(dp), intent(in) :: F(n)
    real(dp), intent(out) :: xnext(n)
    real(dp), intent(in) :: coeff(:)

    select case ( mix%m )
    case ( MIX_LINEAR )
       call mixing_linear()
    case ( MIX_PULAY )
       call mixing_pulay()
    case ( MIX_BROYDEN )
       call mixing_broyden()
    end select

  contains

    subroutine mixing_linear()
      integer :: i
      real(dp) :: w
      w = mix%w

      if ( debug_mix ) write(*,'(2a,e10.4)') &
           trim(debug_msg),' alpha = ', w
      
!$OMP parallel do default(shared), private(i)
      do i = 1 , n
         xnext(i) = xin(i) + w * F(i)
      end do
!$OMP end parallel do
      
    end subroutine mixing_linear

    subroutine mixing_pulay()
      integer :: ns, nh
      integer :: i, j
      logical :: lreturn
      real(dp) :: G
      real(dp), pointer :: res(:), rres(:)

      ns = n_items(mix%stack(1))
      nh = size(coeff)
      if ( nh /= n_items(mix%stack(2)) ) then
         write(*,'(a)') 'mix: Error in mixing of Pulay'
         xnext = 0._dp
         return
      end if

      ! Easy check for initial step...
      select case ( mix%v )
      case ( 0 , 2 ) ! Stable Pulay
         
         lreturn = ns == 1
         
         if ( lreturn .and. debug_mix ) &
              write(*,'(2a,e10.4)') trim(debug_msg), &
              ' Pulay (initial), alpha = ', mix%rv(1)
         
      case ( 1 , 3 ) ! Guaranteed Pulay
         
         lreturn = mod(current_itt(mix), 2) == 1

         if ( lreturn .and. debug_mix ) &
              write(*,'(2a,e10.4)') trim(debug_msg), &
              ' Direct mixing, alpha = ', mix%rv(1)

      end select

      ! In case we return we are actually doing
      ! linear mixing
      if ( lreturn ) then

         ! We are doing a linear mixing
!$OMP parallel do default(shared), private(i)
         do i = 1 , n
            xnext(i) = xin(i) + F(i) * mix%rv(1)
         end do
!$OMP end parallel do

         return

      end if

      ! Get the linear mixing term...
      G = mix%w

      ! if debugging print out the different variables
      if ( debug_mix ) then
         write(*,'(2a,f10.6,a,e10.4,a,100(tr1,e10.4))') &
              trim(debug_msg),&
              ' G = ',G,', sum(alpha) = ',sum(coeff), &
              ', alpha = ',coeff
      end if


!$OMP parallel default(shared), private(i, j)

      !  x^opt[i+1] = x[i] + G F[i]
!$OMP do
      do i = 1 , n
         xnext(i) = xin(i) + G * F(i)
      end do
!$OMP end do
      
      do j = 1 , nh
      
         !  res == x[j] - x[j-1]
         ! rres == F[j+1] - F[j]
!$OMP single
        res => getstackval(mix, 1, j)
        rres => getstackval(mix, 2, j)
!$OMP end single ! implicit barrier

         !  x^opt[i+1] = x^opt[i+1] +
         !       alpha_j ( x[j] - x[j-1] + G (F[j+1] - F[j]) )
!$OMP do
         do i = 1 , n
            xnext(i) = xnext(i) + coeff(j) * ( res(i) + G * rres(i) )
         end do
!$OMP end do

      end do

!$OMP end parallel

    end subroutine mixing_pulay

    subroutine mixing_broyden()
      integer :: ns, nh
      integer :: i, j
      real(dp) :: G
      real(dp), pointer :: res(:), rres(:)

      ns = n_items(mix%stack(1))
      nh = size(coeff)
      if ( nh /= n_items(mix%stack(2)) ) then
         write(*,'(a)') 'mix: Error in mixing of Broyden'
         xnext = 0._dp
         return
      end if

      if ( ns == 1 ) then
         if ( debug_mix ) &
              write(*,'(2a,e10.4)') trim(debug_msg), &
              ' Broyden (initial), alpha = ', mix%rv(1)

         ! We are doing a linear mixing
!$OMP parallel do default(shared), private(i)
         do i = 1 , n
            xnext(i) = xin(i) + F(i) * mix%rv(1)
         end do
!$OMP end parallel do

         return
         
      end if

      ! Get the inverse Jacobian term...
      G = mix%w

      ! if debugging print out the different variables
      if ( debug_mix ) then
         write(*,'(2a,f10.6,a,e10.4,a,100(tr1,e10.4))') &
              trim(debug_msg), ' G = ', G, &
              ', sum(coeff) = ',sum(coeff), &
              ', coeff = ',coeff
      end if


!$OMP parallel default(shared), private(i, j)

      !  x^opt[i+1] = x[i] + G F[i]
!$OMP do
      do i = 1 , n
         xnext(i) = xin(i) + G * F(i)
      end do
!$OMP end do
      
      do j = 1 , nh
      
         !  res == x[j] - x[j-1]
         ! rres == F[j+1] - F[j]
!$OMP single
         res => getstackval(mix, 1, j)
         rres => getstackval(mix, 2, j)
!$OMP end single ! implicit barrier

         !  x^opt[i+1] = x^opt[i+1] +
         !       alpha_j ( x[j] - x[j-1] + G (F[j+1] - F[j]) )
!$OMP do
         do i = 1 , n
            xnext(i) = xnext(i) + coeff(j) * ( res(i) + G * rres(i) )
         end do
!$OMP end do

      end do

!$OMP end parallel

    end subroutine mixing_broyden

  end subroutine mixing_calc_next


  !> Finalize the mixing algorithm
  !!
  !! @param[inout] mix mixer to be finalized
  !! @param[in] n size of the input arrays
  !! @param[in] xin the input for this iteration
  !! @param[in] F the residual for this iteration
  !! @param[in] xnext the optimized input for the next iteration
  subroutine mixing_finalize(mix, n, xin, F, xnext)
    use parallel, only: IONode
    
    type(tMixer), pointer :: mix
    integer, intent(in) :: n
    real(dp), intent(in) :: xin(n), F(n), xnext(n)
    integer :: rsave, is

    select case ( mix%m )
    case ( MIX_LINEAR )
       call fin_linear()
    case ( MIX_PULAY )
       call fin_pulay()
    case ( MIX_BROYDEN )
       call fin_broyden()
    end select


    ! Fix the action to finalize it..
    if ( mix%restart > 0 ) then
      if ( mod(current_itt(mix),mix%restart) == 0 ) then
        mix%action = IOR(mix%action, ACTION_RESTART)
      end if
    end if

    
    ! Check the actual finalization...
    ! First check whether we should restart history
    if ( IAND(mix%action, ACTION_RESTART) == ACTION_RESTART ) then

       ! The user has requested to restart the
       ! mixing scheme now
       rsave = mix%restart_save

       select case ( mix%m )
       case ( MIX_PULAY )

          if ( IONode ) then
             write(*,'(a)')'mix: Pulay -- resetting history'
          end if
          
          if ( rsave == 0 ) then
            do is = 1, size(mix%stack)
              call reset(mix%stack(is))
            end do
          else
            call reset(mix%stack(1),-rsave)
            call reset(mix%stack(2),-rsave+1)
          end if
          
       case ( MIX_BROYDEN )

          if ( IONode ) then
             write(*,'(a)')'mix: Broyden -- resetting history'
          end if

          if ( rsave == 0 ) then
             call reset(mix%stack(1))
             call reset(mix%stack(2))
             call reset(mix%stack(3))
          else
             call reset(mix%stack(1),-rsave)
             call reset(mix%stack(2),-rsave+1)
          end if
          
       end select

       if ( allocated(mix%stack) ) then
          if ( debug_mix ) &
               write(*,'(a,a,i0)') trim(debug_msg), &
               ' saved hist = ',n_items(mix%stack(1))
       end if
       
    end if

    
    ! check whether we should change the mixer
    if ( IAND(mix%action, ACTION_NEXT) == ACTION_NEXT ) then
       
       call mixing_step( mix )
       
    end if
    
  contains

    subroutine fin_linear()

      ! do nothing...

    end subroutine fin_linear

    subroutine fin_pulay()

      integer :: ns
      integer :: i
      logical :: GR_linear
      real(dp), pointer :: res(:), rres(:)

      ns = n_items(mix%stack(1))
      
      select case ( mix%v )
      case ( 0 , 2 ) ! stable Pulay

         if ( n_items(mix%stack(3)) == 0 ) then
            call push_stack_data(mix%stack(3), n)
         end if

         res => getstackval(mix, 3)

!$OMP parallel do default(shared), private(i)
         do i = 1 , n
            ! Input:
            !  res == x[i-1] + F[i-1]
            res(i) = xin(i) + F(i)
            ! Output:
            !  res == x[i] + F[i]
         end do
!$OMP end parallel do

      case ( 1 , 3 ) ! GR Pulay
         
         GR_linear = mod(current_itt(mix), 2) == 1
         
         if ( n_items(mix%stack(2)) > 0 .and. &
              .not. GR_linear ) then
            
            res => getstackval(mix, 1)
            rres => getstackval(mix, 2)
            
!$OMP parallel do default(shared), private(i)
            do i = 1 , n
               ! Input:
               !  rres == F[i] - F[i-1]
               rres(i) = rres(i) - res(i)
               ! Output:
               !  rres == - F[i-1]
            end do
!$OMP end parallel do

            call pop(mix%stack(1))
            
            ! Note that this is Res[i-1] = (F^i-1_out - F^i-1_in)
            res => getstackval(mix,1)
!$OMP parallel do default(shared), private(i)
            do i = 1 , n
               res(i) = res(i) - xin(i) + xnext(i)
            end do
!$OMP end parallel do

         end if
         
      end select
      
    end subroutine fin_pulay

    subroutine fin_broyden()

      integer :: ns, nh
      integer :: i
      real(dp), pointer :: res(:), rres(:)

      ns = current_itt(mix)
      nh = n_items(mix%stack(2))

      if ( ns >= 2 .and. n_items(mix%stack(3)) > 0 ) then
         
         ! Update the residual to reflect the input residual
         res => getstackval(mix, 3)
         
!$OMP parallel do default(shared), private(i)
         do i = 1 , n
            ! Input:
            !  res == x[i-1] + F[i-1]
            res(i) = xin(i) + F(i)
            ! Output:
            !  res == x[i] + F[i]
         end do
!$OMP end parallel do

      end if

      ! Update weights (if necessary)
      if ( nh > 1 ) then
         do i = 3 , nh + 1
            mix%rv(i) = mix%rv(i+1)
         end do
      end if
      
    end subroutine fin_broyden

  end subroutine mixing_finalize
  

  ! Perform the actual mixing...
  subroutine mixing_1d( mix, n, xin, F, xnext, nsub)
    ! The current mixing method 
    type(tMixer), pointer :: mix
    ! The current step in the SCF and size of arrays
    integer, intent(in) :: n
    ! x1 == Input function,
    ! F1 == Residual from x1
    real(dp), intent(in) :: xin(n), F(n)
    ! x2 == Next input function
    real(dp), intent(inout) :: xnext(n)
    ! Number of elements used for calculating the mixing
    ! coefficients
    integer, intent(in), optional :: nsub

    ! Coefficients
    integer :: ncoeff
    real(dp), allocatable :: coeff(:)

    call mixing_init( mix, n, xin, F )

    ncoeff = mixing_ncoeff( mix )
    allocate(coeff(ncoeff))

    ! Calculate coefficients
    if ( present(nsub) ) then
       call mixing_coeff( mix, nsub, xin, F, coeff)
    else
       call mixing_coeff( mix, n, xin, F, coeff)
    end if

    ! Calculate the following output
    call mixing_calc_next( mix, n, xin, F, xnext, coeff)

    ! Coefficients are not needed anymore...
    deallocate(coeff)

    ! Finalize the mixer
    call mixing_finalize( mix, n, xin, F, xnext)

  end subroutine mixing_1d
  

  subroutine mixing_2d( mix, n1, n2, xin, F, xnext, nsub)
    
    type(tMixer), pointer :: mix
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: xin(n1,n2), F(n1,n2)
    real(dp), intent(inout) :: xnext(n1,n2)
    integer, intent(in), optional :: nsub

    ! Simple wrapper for 1D
    if ( present(nsub) ) then
       call mixing_1d( mix, n1*n2 , xin(1,1), F(1,1), xnext(1,1) ,&
            nsub = n1 * nsub)
    else
       call mixing_1d( mix, n1*n2 , xin(1,1), F(1,1), xnext(1,1))
    end if
    
  end subroutine mixing_2d




  ! Step the mixing object and ensure that
  ! the old history is either copied over or freed
  subroutine mixing_step( mix )

    use parallel, only: IONode 

    type(tMixer), pointer :: mix

    type(tMixer), pointer :: next => null()

    type(dData1D), pointer :: d1D
    integer :: i, is, n, init_itt
    logical :: reset_stack, copy_stack

    ! First try and
    next => mix%next
    if ( associated(next) ) then
      
      ! Whether or not the two methods are allowed
      ! to share history
      copy_stack = mix%m == next%m
      ! Only Pulay original and Broyden methods are compatible.
      ! The Pulay-GR variant does not have the same data stored.
      select case ( mix%m )
      case ( MIX_PULAY )
        select case ( next%m )
        case ( MIX_PULAY )
          copy_stack = mod(mix%v, 2) == mod(next%v, 2)
        case ( MIX_BROYDEN )
          copy_stack = mod(mix%v, 2) == 0
        end select
      case ( MIX_BROYDEN )
        select case ( next%m )
        case ( MIX_PULAY )
          copy_stack = mod(mix%v, 2) == 0
          ! otherwise both are broyden
        end select
      end select

      copy_stack = copy_stack .and. allocated(mix%stack)

      ! If the two methods are similar
      if ( copy_stack ) then

        ! They are similar, copy over the history stack
        do is = 1 , size(mix%stack)

          ! Get maximum size of the current stack,
          n = n_items(mix%stack(is))

          ! Note that this will automatically take care of
          ! wrap-arounds and delete the unneccesry elements
          do i = 1 , n
            d1D => get_pointer(mix%stack(is),i)
            call push(next%stack(is),d1D)
          end do

          ! nullify
          nullify(d1D)

        end do

      end if

    end if

    reset_stack = .true.
    if ( associated(next) ) then
      if ( associated(next%next, mix) .and. &
          next%n_itt > 0 ) then
        ! if this is a circular mixing routine
        ! we should not reset the history...
        reset_stack = .false.
      end if
    end if

    if ( reset_stack ) then
      select case ( mix%m )
      case ( MIX_PULAY , MIX_BROYDEN )

        n = size(mix%stack)
        do is = 1 , n
          call reset(mix%stack(is))
        end do

      end select
    end if

    if ( associated(next) ) then
      init_itt = 0
      ! Set-up the next mixer
      select case ( next%m )
      case ( MIX_PULAY , MIX_BROYDEN )
        init_itt = n_items(next%stack(1))
      end select
      next%start_itt = init_itt
      next%cur_itt = init_itt
      if ( IONode ) then
        if ( copy_stack ) then
          write(*,'(3a)') trim(debug_msg),' switching mixer (with history) --> ', &
              trim(next%name)
        else
          write(*,'(3a)') trim(debug_msg),' switching mixer --> ', &
              trim(next%name)
        end if
      end if
      mix => mix%next
    end if

  end subroutine mixing_step



  ! Calculate the inverse of a matrix
  subroutine inverse(n, A, B, info )

    integer, intent(in) :: n
    real(dp), intent(in)  :: A(n,n)
    real(dp), intent(out) :: B(n,n)
    integer, intent(out) :: info

    integer :: i, j
    ! Local arrays
    real(dp) :: pm(n,n), work(n*4), err
    ! Relative tolerance dependent on the magnitude
    ! For now we retain the old tolerance
    real(dp), parameter :: etol = 1.e-4_dp
    integer :: ipiv(n)

    ! initialize info
    info = 0

    ! simple check and fast return
    if ( n == 1 ) then

       B(1,1) = 1._dp / A(1,1)

       return

    end if

    call lapack_inv()
    if ( info /= 0 ) call simple_inv()

  contains

    subroutine lapack_inv()

      B(:, :) = A(:, :)

      call dgetrf(n,n,B,n,ipiv,info)
      if ( info /= 0 ) return
      
      call dgetri(n,B,n,ipiv,work,n*4,info)
      if ( info /= 0 ) return

      ! This sets info appropriately
      call check_inv()
      
    end subroutine lapack_inv
    
    subroutine simple_inv()
      
      real(dp) :: x
      integer :: k
      
      ! Copy over A
      B(:, :) = A(:, :)
      
      do i = 1 , n
         if ( B(i,i) == 0._dp ) then
            info = -n
            return
         end if
         
         x = 1._dp / B(i,i)
         B(i,i) = 1._dp
         do j = 1 , n
            B(j,i) = B(j,i) * x
         end do

         do k = 1 , n
            if ( (k-i) /= 0 ) then
               x = B(i,k)
               B(i,k) = 0._dp
               do j = 1 , n
                  B(j,k) = B(j,k) - B(j,i) * x
               end do
            end if
         end do
      end do

      ! This sets info appropriately
      call check_inv()

    end subroutine simple_inv

    subroutine check_inv()
      
      ! Check correcteness
      pm = matmul(A,B)
      
      do j = 1 , n 
         do i = 1 , n
            if ( i == j ) then
               err = pm(i,j) - 1._dp
            else
               err = pm(i,j)
            end if

            ! This is pretty strict tolerance!
            if ( abs(err) > etol ) then
               ! Signal failure in inversion
               info = - n - 1
               return
            end if
            
         end do
      end do
      
    end subroutine check_inv
    
  end subroutine inverse


  ! Calculate the svd of a matrix
  ! With   ||(Ax - b)|| <<
  subroutine svd(n, A, B, cond, rank, info)

    integer, intent(in) :: n
    real(dp), intent(in)  :: A(n,n)
    real(dp), intent(out) :: B(n,n)
    real(dp), intent(in) :: cond
    integer, intent(out) :: rank, info

    ! Local arrays
    integer :: i
    character(len=64) :: fmt
    real(dp) :: AA(n,n), S(n), work(n*5)

    if ( n == 1 ) then
      ! Simple scheme
      B(1,1) = 1._dp / A(1,1)
      return
    end if

    ! Copy A matrix
    AA(:, :) = A(:, :)
    ! setup pseudo inverse solution for minimizing
    ! constraints
    B(:, :) = 0._dp
    do i = 1 , n
       B(i,i) = 1._dp
    end do
    call dgelss(n,n,n,AA,n,B,n,S,cond,rank,work,n*5,info)
    
    ! if debugging print out the different variables
    if ( debug_mix ) then
      ! also mark the rank
      
      if ( rank == n ) then
        ! complete rank
        write(*,'(2a,100(tr1,e10.4))') &
            trim(debug_msg),' SVD singular = ',S
      else
        ! this prints the location of the SVD rank, if not full
        write(fmt,'(i0,2a)') rank, '(tr1,e10.4),'' >'',100(tr1,e10.4)'
        write(*,'(2a,'//trim(fmt)//')') &
            trim(debug_msg),' SVD singular = ',S
      end if
    end if

    ! One could potentially search a smaller part of the iterations
    ! for a complete satisfactory SVD solution.
    ! However, typically the singular values are because the
    ! two latest iterations which means that one will end up with
    ! the linear mixing scheme with only a coefficient for the latest
    ! iteration.
    ! This, is per see not a bad choice since it forces a linear mixing
    ! in the end. However, the SVD solution seems better if one tries to
    ! really push the tolerances.

    ! One could do this (which requires the routine to be recursive)
    ! if ( rank < n ) then
    !   call svd(n-1,A(1:n-1,1:n-1), B(1:n-1, 1:n-1), cond, rank, info)
    !   ! zero out everything else
    !   B(n,:) = 0._dp
    !   B(:,n) = 0._dp
    ! end if

  end subroutine svd


  ! *************************************************
  ! *                Helper routines                *
  ! *                    LOCAL                      *
  ! *************************************************
  
  ! Returns the value array from the stack(:)
  ! Returns this array:
  !    mix%stack(sidx)(hidx) ! defaults to the last item
  function getstackval(mix,sidx,hidx) result(d1)
    type(tMixer), intent(in) :: mix
    integer, intent(in) :: sidx
    integer, intent(in), optional :: hidx

    real(dp), pointer :: d1(:)

    type(dData1D), pointer :: dD1
    
    if ( present(hidx) ) then
       dD1 => get_pointer(mix%stack(sidx),hidx)
    else
       dD1 => get_pointer(mix%stack(sidx), &
            n_items(mix%stack(sidx)))
    end if
    
    d1 => val(dD1)
    
  end function getstackval


  ! Returns true if the following 
  ! "advanced" mixer is 'method'
  function is_next(mix,method,next) result(bool)
    type(tMixer), intent(in), target :: mix
    integer, intent(in) :: method
    type(tMixer), pointer, optional :: next

    logical :: bool

    type(tMixer), pointer :: m

    bool = .false.
    m => mix%next
    do while ( associated(m) )

       if ( m%m == MIX_LINEAR ) then
          m => m%next
       else if ( m%m == method ) then
          bool = .true.
          exit
       else
          ! Quit if it does not do anything
          exit
       end if

       ! this will prevent cyclic combinations
       if ( associated(m,mix) ) exit

    end do

    if ( present(next) ) then
       next => m
    end if

  end function is_next


  !> Get current iteration count
  !!
  !! This is abstracted because the initial iteration
  !! and the current iteration may be uniquely defined.
  function current_itt(mix) result(itt)
    type(tMixer), intent(in) :: mix
    integer :: itt

    itt = mix%cur_itt - mix%start_itt

  end function current_itt

  
  ! Stack handling routines
  function stack_check(stack,n) result(check)
    type(Fstack_dData1D), intent(inout) :: stack
    integer, intent(in) :: n
    logical :: check

    ! Local arrays
    type(dData1D), pointer :: dD1

    if ( n_items(stack) == 0 ) then
       check = .true.
    else

       ! Check that the stack stored arrays are
       ! of same size...
       
       dD1 => get_pointer(stack,1)
       check = n == size(dD1)

    end if

  end function stack_check

  subroutine push_stack_data(s_F,n)
    type(Fstack_dData1D), intent(inout) :: s_F
    integer, intent(in) :: n

    type(dData1D) :: dD1

    if ( .not. stack_check(s_F,n) ) then
       call die('mixing: history has changed size...')
    end if

    call newdData1D(dD1, n, '(F)')
    
    ! Push the data to the stack
    call push(s_F,dD1)
    
    ! Delete double reference
    call delete(dD1)
       
  end subroutine push_stack_data

  subroutine push_F(s_F,n,F,fact)
    type(Fstack_dData1D), intent(inout) :: s_F
    integer, intent(in) :: n
    real(dp), intent(in) :: F(n)
    real(dp), intent(in), optional :: fact

    type(dData1D) :: dD1
    real(dp), pointer :: sF(:)
    integer :: in, ns
    integer :: i

    if ( .not. stack_check(s_F,n) ) then
       call die('mixing: history has changed size...')
    end if

    in = n_items(s_F)
    ns = max_size(s_F)

    if ( in == ns ) then
       
       ! we have to cycle the storage
       call get(s_F,1,dD1)

    else

       call newdData1D(dD1, n, '(F)')

    end if

    sF => val(dD1)

    if ( present(fact) ) then
!$OMP parallel do default(shared), private(i)
       do i = 1 , n
          sF(i) = F(i) * fact
       end do
!$OMP end parallel do
    else
       call dcopy(n, F, 1, sF, 1)
    end if

    ! Push the data to the stack
    call push(s_F,dD1)

    ! Delete double reference
    call delete(dD1)

  end subroutine push_F

  subroutine update_F(s_F,n,F)
    type(Fstack_dData1D), intent(inout) :: s_F
    integer, intent(in) :: n
    real(dp), intent(in) :: F(n)

    type(dData1D), pointer :: dD1
    real(dp), pointer :: FF(:)
    integer :: in

    if ( .not. stack_check(s_F,n) ) then
       call die('mixing: history has changed size...')
    end if

    in = n_items(s_F)

    if ( in == 0 ) then
       
       ! We need to add it as it does not exist
       call push_F(s_F,n,F)

    else

       ! we have an entry, update the latest
       dD1 => get_pointer(s_F,in)

       FF => val(dD1)

       call dcopy(n, F, 1, FF, 1)

    end if

  end subroutine update_F

  subroutine push_diff(s_rres,s_res,alpha)
    type(Fstack_dData1D), intent(inout) :: s_rres
    type(Fstack_dData1D), intent(in) :: s_res
    real(dp), intent(in), optional :: alpha

    type(dData1D) :: dD1
    type(dData1D), pointer :: pD1
    real(dp), pointer :: res1(:), res2(:), rres(:)
    integer :: in, ns, i, n

    if ( n_items(s_res) < 2 ) then
       call die('mixing: Residual residuals cannot be calculated, &
            &inferior residual size.')
    end if

    in = n_items(s_res)

    ! First get the value of in
    pD1 => get_pointer(s_res,in-1)
    res1 => val(pD1)
    ! get the value of in
    pD1 => get_pointer(s_res,in)
    res2 => val(pD1)

    in = n_items(s_rres)
    ns = max_size(s_rres)

    if ( in == ns ) then
       
       ! we have to cycle the storage
       call get(s_rres,1,dD1)

    else

       call newdData1D(dD1, size(res1), '(res)')

    end if

    ! Get the residual of the residual
    rres => val(dD1)
    n = size(rres)

    if ( present(alpha) ) then
!$OMP parallel do default(shared), private(i)
       do i = 1 , n
          rres(i) = (res2(i) - res1(i)) * alpha
       end do
!$OMP end parallel do
    else
!$OMP parallel do default(shared), private(i)
       do i = 1 , n
          rres(i) = res2(i) - res1(i)
       end do
!$OMP end parallel do
    end if
    
    ! Push the data to the stack
    call push(s_rres,dD1)

    ! Delete double reference
    call delete(dD1)

  end subroutine push_diff


  !> Calculate the norm of two arrays
  function norm(n, x1, x2)
    integer, intent(in) :: n
    real(dp), intent(in) :: x1(n), x2(n)
    real(dp) :: norm

    ! Currently we use an external routine
    integer :: i

    ! Calculate dot product
    
    norm = 0._dp
!$OMP parallel do default(shared), private(i) &
!$OMP& reduction(+:norm)
    do i = 1 , n
       norm = norm + x1(i) * x2(i)
    end do
!$OMP end parallel do
    
  end function norm


end module m_mixing
  
