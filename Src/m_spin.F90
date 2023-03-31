! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module t_spin

  implicit none

  !> Type containing a simulations spin configuration
  !!
  !! Thus this type contains _all_ relevant information
  !! regarding the spin-configuration.
  type tSpin

     !> Dimensionality of the Hamiltonian
     integer :: H = 1
     !> Dimensionality of the density matrix
     integer :: DM = 1
     !> Dimensionality of the energy density matrix
     integer :: EDM = 1
     !> Dimensionality of the grid operations
     integer :: Grid = 1
     !> Dimensionality of the diagonal spin-components
     integer :: spinor = 1

     ! Storing the type of calculation
     
     !> Whether the simulation has no spin
     logical :: none = .true.
     !> Collinear spin
     logical :: Col = .false.
     !> Non-colinear spin
     logical :: NCol = .false.
     !> Spin-orbit coupling
     logical :: SO = .false.

     ! Perhaps one could argue that one may
     ! associate a symmetry to the spin which
     ! then denotes whether the spin-configuration
     ! is assumed time-reversal symmetric...
     ! Hmm... 

  end type tSpin

end module t_spin

module m_spin
  use precision, only: dp

  use t_spin, only: tSpin

  implicit none
  
  private

  !> Spin configuration for SIESTA
  type(tSpin), public, save :: spin

  ! Use plain integers instead of pointers, to avoid problems
  ! in the PEXSI-only nodes, which might not call the spin_init
  ! routine. The values are copied at the end of that routine.
  
  ! Create short-hands for the spin-configuration
  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%Grid
    integer, save, public, pointer :: nspin => null() ! (Grid)
  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%spinor
    integer, save, public, pointer :: spinor_dim => null() ! (spinor)
  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%H, spin%DM
    integer, save, public, pointer :: h_spin_dim => null() ! (H and DM)
  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%EDM
    integer, save, public, pointer :: e_spin_dim => null() ! (EDM)
  

  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%none
  logical, save, public, pointer :: NoMagn ! (none)
  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%Col
  logical, save, public, pointer :: SPpol ! (Col)
  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%NCol
  logical, save, public, pointer :: NonCol ! (NCol)
  ! DO NOT USE THIS VARIABLE, USE -> type(tSpin) :: spin%SO
  logical, save, public, pointer :: SpOrb ! (SO)

  ! TODO : this is unrelated to the spin-configuration...
  ! Consider moving this to some other place...
  logical, save, public :: TrSym   = .true.

  ! Whether we are performing spiral arrangement of spins
  logical, save, public :: Spiral  = .false.
  ! Pitch wave vector for spiral configuration
  real(dp), save, public :: qSpiral(3) = 0._dp

  public :: init_spin

  public :: print_spin_options
  public :: init_spiral

  public :: fname_spin

contains
  
  subroutine init_spin(default_nspin)
    
    use sys, only: die
    use fdf, only : fdf_get, leqi, fdf_deprecated
    use alloc, only: re_alloc

    !< Externally option for a specific spin configuration for the Hamiltonian.
    !! Upon entry this value is used to setup the default settings, but
    !! fdf-flags still has precedence.
    !! Upon exit, it contains the requested spin-number for the Hamiltonian
    integer, intent(inout), optional :: default_nspin

    character(len=32) :: opt, opt_old

    ! Create pointer assignments...
    call int_pointer(spinor_dim, spin%spinor)
    call int_pointer(nspin     , spin%grid)
    call int_pointer(h_spin_dim, spin%H)
    call int_pointer(e_spin_dim, spin%EDM)

    ! Create pointer assignments...
    call log_pointer(NoMagn, spin%none)
    call log_pointer(SPpol , spin%Col)
    call log_pointer(NonCol, spin%NCol)
    call log_pointer(SpOrb , spin%SO)

    ! Time reversal symmetry
    TrSym  = .true.

    ! All components of the 'spin' variable
    ! is initially correct...
    spin%none = .false.
    spin%Col = .false.
    spin%NCol = .false.
    spin%SO = .false.
    
    ! Read in old flags (discouraged)
    if ( present(default_nspin) ) then
      select case ( default_nspin )
      case ( 1 )
        ! default to non-polarized (the default)
      case ( 2 )
        spin%Col = .true.
      case ( 4 )
        spin%NCol = .true.
      case ( 8 )
        spin%SO = .true.
      end select
    end if
    spin%Col  = fdf_get('SpinPolarized', spin%Col)
    spin%NCol = fdf_get('NonCollinearSpin', spin%NCol)
    spin%SO   = fdf_get('SpinOrbit', spin%SO)

    ! Announce the deprecated flags (if used)...
    call fdf_deprecated('SpinPolarized','Spin')
    call fdf_deprecated('NonCollinearSpin','Spin')
    call fdf_deprecated('SpinOrbit','Spin')

    ! Set default option from "old" options
    if ( spin%SO ) then
       opt_old = 'spin-orbit'
    else if ( spin%NCol ) then
       opt_old = 'non-collinear'
    else if ( spin%Col ) then
       opt_old = 'polarized'
    else
       opt_old = 'none'
    end if

    ! In order to enable text input (and obsolete the
    ! 4 different options we use a single value now)
    opt = fdf_get('Spin', opt_old)
    
    if ( leqi(opt, 'none') .or. &
         leqi(opt, 'non-polarized') .or. &
         leqi(opt, 'non-polarised') .or. &
         leqi(opt, 'NP') .or. leqi(opt,'N-P') ) then

       spin%none = .true.
       call warn_and_set_to_false(Spin%Col)
       call warn_and_set_to_false(Spin%NCol)
       call warn_and_set_to_false(Spin%SO)
       
    else if ( leqi(opt, 'polarized') .or. &
         leqi(opt, 'collinear') .or. leqi(opt, 'colinear') .or. &
         leqi(opt, 'polarised') .or. leqi(opt, 'P') ) then
       
       spin%Col = .true.
       call warn_and_set_to_false(Spin%none)
       call warn_and_set_to_false(Spin%NCol)
       call warn_and_set_to_false(Spin%SO)
       
    else if ( leqi(opt, 'non-collinear') .or. leqi(opt, 'non-colinear') .or. &
         leqi(opt, 'NC') .or. leqi(opt, 'N-C') ) then
       
       spin%NCol = .true.
       call warn_and_set_to_false(Spin%none)
       call warn_and_set_to_false(Spin%Col)
       call warn_and_set_to_false(Spin%SO)
       
    else if ( leqi(opt, 'spin-orbit') .or. leqi(opt, 'S-O') .or. &
         leqi(opt, 'SOC') .or. leqi(opt, 'SO') ) then
       call warn_and_set_to_false(Spin%none)
       call warn_and_set_to_false(Spin%Col)
       call warn_and_set_to_false(Spin%NCol)
       
       spin%SO = .true.
       call warn_and_set_to_false(Spin%none)
       call warn_and_set_to_false(Spin%Col)
       call warn_and_set_to_false(Spin%NCol)
       
    else
       write(*,*) 'Unknown spin flag: ', trim(opt)
       call die('Spin: unknown flag, please assert the correct input.')
    end if

    ! Note that, in what follows,
    !   spinor_dim = min(h_spin_dim,2)
    !   e_spin_dim = min(h_spin_dim,4)
    !   nspin      = min(h_spin_dim,4)  ! Relevant for dhscf, local DM
    !      should probably be called nspin_grid
    !
    ! so everything can be determined if h_spin_dim is known.
    ! It is tempting to go back to the old 'nspin' overloading,
    ! making 'nspin' mean again 'h_spin_dim'.
    ! But this has to be done carefully, as some routines expect
    ! an argument 'nspin' which really means 'spinor_dim' (like diagon),
    ! and others (such as dhscf) expect 'nspin' to mean 'nspin_grid'.

    if ( spin%SO ) then
       ! Spin-orbit case

       ! Dimensions
       spin%H      = 8
       spin%DM     = 8
       spin%EDM    = 4
       spin%Grid   = 4
       spin%spinor = 2

       ! Flags
       spin%none = .false.
       spin%Col  = .false.
       spin%NCol = .false.
       spin%SO   = .true.

       ! should be moved...
       TRSym      = .false.

    else if ( spin%NCol ) then
       ! Non-collinear case

       ! Dimensions
       spin%H      = 4
       spin%DM     = 4
       spin%EDM    = 4
       spin%Grid   = 4
       spin%spinor = 2

       ! Flags
       spin%none = .false.
       spin%Col  = .false.
       spin%NCol = .true.
       spin%SO   = .false.

       ! should be moved...
       TRSym      = .false.

    else if ( spin%Col ) then
       ! Collinear case

       ! Dimensions
       spin%H      = 2
       spin%DM     = 2
       spin%EDM    = 2
       spin%Grid   = 2
       spin%spinor = 2

       ! Flags
       spin%none = .false.
       spin%Col  = .true.
       spin%NCol = .false.
       spin%SO   = .false.

       ! should be moved...
       TRSym      = .true.

    else if ( spin%none ) then
       ! No spin configuration...

       ! Dimensions
       spin%H      = 1
       spin%DM     = 1
       spin%EDM    = 1
       spin%Grid   = 1
       spin%spinor = 1

       ! Flags
       spin%none = .true.
       spin%Col  = .false.
       spin%NCol = .false.
       spin%SO   = .false.

       ! should be moved...
       TRSym      = .true.

    end if

    ! Get true time reversal symmetry
    TRSym = fdf_get('TimeReversalSymmetry',TRSym)

    if ( present(default_nspin) ) then
      default_nspin = spin%H
    end if

  contains

    subroutine int_pointer(from, to)
      integer, pointer :: from
      integer, intent(in), target :: to

      from => to

    end subroutine int_pointer

    subroutine log_pointer(from, to)
      logical, pointer :: from
      logical, intent(in), target :: to

      from => to

    end subroutine log_pointer
    
    subroutine warn_and_set_to_false(var)
      logical, intent(inout) :: var
      if (var) then
         call message("WARNING","Deprecated spin keyword overridden by new-style 'Spin' input")
         call message("WARNING","Option from deprecated keyword: "//trim(opt_old))
         call message("WARNING","Option from 'Spin' keyword input: "//trim(opt))
         var = .false.
      endif
    end subroutine warn_and_set_to_false

  end subroutine init_spin


  ! Print out spin-configuration options
  subroutine print_spin_options( )
    use parallel, only: IONode
    use m_cite, only : add_citation

    character(len=32) :: opt

    if ( .not. IONode ) return

    if ( spin%SO ) then
       opt = 'spin-orbit'
       call add_citation("10.1088/0953-8984/19/19/489001")
    else if ( spin%NCol ) then
       opt = 'non-collinear'
    else if ( spin%Col ) then
       opt = 'collinear'
    else 
       opt = 'none'
    end if

    write(*,'(a,t53,''= '',a)') 'redata: Spin configuration',trim(opt)
    write(*,'(a,t53,''= '',i0)')'redata: Number of spin components',spin%H
    write(*,'(a,t53,''= '',l1)')'redata: Time-Reversal Symmetry',TRSym
    write(*,'(a,t53,''= '',l1)')'redata: Spin spiral', Spiral
    if ( Spiral ) then
      write(*,'(a,t53,''='',3(tr1,e12.6))') &
          'redata: Spin spiral pitch wave vector', qSpiral
      if ( .not. spin%NCol ) then
        write(*,'(a)') 'redata: WARNING: spin spiral requires non-collinear spin'
        call die("Spin spiral requires non-collinear spin")
      end if
    end if

    if ( spin%SO ) then
       write(*,'(a)') repeat('#',60)
       call message("WARNING","This spin-orbit implementation uses a local approximation.")
       call message("WARNING","From 4.2 and onwards the full non-local &
           &approximation is implemented")
       call message("WARNING","You are strongly advised to use >=4.2 versions!")
       write(*,'(a)') repeat('#',60)
    end if

  end subroutine print_spin_options
  

  subroutine init_spiral( ucell )
    use fdf, only: fdf_get
    use fdf, only: block_fdf, parsed_line
    use fdf, only: fdf_block, fdf_bline, fdf_bclose
    use fdf, only: fdf_bvalues
    use m_get_kpoints_scale, only: get_kpoints_scale

    ! Unit cell lattice vectors
    real(dp), intent(in) :: ucell(3,3)

    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline

    ! Reciprocal cell vectors
    real(dp) :: rcell(3,3)
    integer :: ierr

    Spiral = fdf_block('Spin.Spiral', bfdf)

    if ( .not. Spiral ) return
    ! We cannot use TRS for spirals, regardless of user!
    TRSym = .false.

    ! Retrieve scale
    call get_kpoints_scale('Spin.Spiral.Scale', rcell, ierr)

    if ( ierr /= 0 ) then
      call die('init_spiral: ERROR in Spin.Spiral.Scale, could &
          &not find scale for spiral wave vector.')
    end if

    if (.not. fdf_bline(bfdf,pline)) &
        call die('init_spiral: ERROR in Spin.Spiral block, could &
        &not find spiral wave vector.')

    ! Read pitch wave-vector
    qSpiral(1) = fdf_bvalues(pline,1)
    qSpiral(2) = fdf_bvalues(pline,2)
    qSpiral(3) = fdf_bvalues(pline,3)
    qSpiral(:) = matmul(rcell, qSpiral)

    call fdf_bclose(bfdf)

  end subroutine init_spiral

  function fname_spin(nspin,ispin,slabel,suffix,basename) result(fname)
    integer, intent(in) :: nspin, ispin
    character(len=*), intent(in), optional :: slabel, suffix, basename
    character(len=200) :: fname
    
    if ( present(basename) ) then
       if ( nspin == 1 ) then
          fname = trim(basename)
       else
          if ( ispin == 1 ) fname = trim(basename)//'_UP'
          if ( ispin == 2 ) fname = trim(basename)//'_DN'
       end if
    else
       if ( .not. &
            ( present(slabel) .and. present(suffix) ) ) &
            call die('Error in filename input')
       if ( nspin == 1 ) then
          fname = trim(slabel)//'.'//trim(suffix)
       else
          if ( ispin == 1 ) fname = trim(slabel)//'.'//trim(suffix)//'_UP'
          if ( ispin == 2 ) fname = trim(slabel)//'.'//trim(suffix)//'_DN'
       end if
    end if
    
  end function fname_spin


end module m_spin
