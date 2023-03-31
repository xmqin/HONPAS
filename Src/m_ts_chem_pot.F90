!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com


module m_ts_chem_pot

  use precision, only : dp

  use units, only : Pi, eV, Kelvin

  use m_ts_io_ctype, only : C_N_NAME_LEN
  
  implicit none

  integer,  public, parameter :: NAME_MU_LEN = 32

  real(dp), public, parameter :: mu_same = 1.e-8_dp
  integer, private, parameter :: def_poles = 8

  type :: ts_mu
     ! name of the chemical potential
     character(len=NAME_MU_LEN) :: name = ' '
     ! ID
     integer  :: ID = 0
     ! Number of poles for the equilibrium contour
     integer  :: N_poles = 0
     ! the chemical potential
     real(dp) :: mu = 0._dp
     character(len=NAME_MU_LEN) :: cmu
     ! the temperature associated with this chemical potential
     real(dp) :: kT = 0._dp
     character(len=NAME_MU_LEN) :: ckT
     ! number of electrodes having this chemical potential
     integer  :: N_El = 0
     ! array of electrode indices (conforming with the Elecs-array)
     integer, pointer :: el(:) => null()
     ! We must have a container which determines the contour segments
     ! that gets attached to the chemical potential
     character(len=C_N_NAME_LEN), allocatable :: Eq_seg(:)
  end type ts_mu
  public :: ts_mu

  interface hasC
     module procedure hasCio
     module procedure hasCeq
  end interface hasC
  public :: hasC

  interface hasEl
     module procedure hasEl_i
  end interface hasEl
  public :: hasEl

  interface name
     module procedure name_
  end interface name
  public :: name

  interface copy
     module procedure copy_
  end interface copy
  public :: copy

  public :: Eq_segs

  public :: chem_pot_add_Elec

  public :: fdf_nmu, fdffake_mu, fdf_mu

  public :: print_mus_block

  private

contains

  function fdf_nmu(prefix,kT,this_n) result(n)
    use fdf

    character(len=*), intent(in) :: prefix
    real(dp), intent(in) :: kT
    type(ts_mu), intent(inout), allocatable :: this_n(:)
    integer :: n

    ! prepare to read in the data...
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    integer :: i

    real(dp) :: E_pole
    logical :: found

    n = 0
    found = fdf_block(trim(prefix)//'.ChemPots',bfdf)
    if ( .not. found ) return

    ! first count the number of electrodes
    n = 0
    do while ( fdf_bline(bfdf,pline) )
       if ( fdf_bnnames(pline) == 0 ) cycle
       n = n + 1 
    end do

    allocate(this_n(n))
    this_n(:)%N_poles = fdf_get('TS.Contours.Eq.Pole.N',def_poles)
    E_pole = fdf_get('TS.Contours.Eq.Pole',1.5_dp*eV,'Ry')
    if ( E_pole > 0._dp ) then
       call E2Npoles(E_pole,kT,i)
       this_n(:)%N_poles = i
    end if

    ! rewind to read again
    call fdf_brewind(bfdf)

    n = 0
    do while ( fdf_bline(bfdf,pline) )
       if ( fdf_bnnames(pline) == 0 ) cycle
       n = n + 1 
       ! Attach the name
       this_n(n)%Name = trim(fdf_bnames(pline,1))
       ! Save the ID of the chemical potential
       this_n(n)%ID = n
       if ( n > 1 ) then
          ! Check that no name is the same
          do i = 1 , n - 1 
             if ( leqi(name(this_n(i)),name(this_n(n))) ) then
                call die('Chemical potential names must not be the same')
             end if
          end do
       end if
    end do

  end function fdf_nmu

  function fdffake_mu(this_n,kT,Volt) result(n)
    use fdf

    type(ts_mu), intent(inout), allocatable :: this_n(:)
    ! SIESTA temperature
    real(dp), intent(in) :: kT, Volt
    real(dp) :: E_pole
    integer :: n, n_pole
    ! In case the user whishes to utilise the standard 
    ! transiesta setup we fake the chemical potentials
    n = 2
    allocate(this_n(n))

    ! Read in number of poles
    this_n(:)%N_poles = fdf_get('TS.Contours.Eq.Pole.N',def_poles)
    E_pole = fdf_get('TS.Contours.Eq.Pole',1.5_dp*eV,'Ry')
    ! If the energy is larger than zero, the user requests
    ! number of poles
    if ( E_pole > 0._dp ) then
       call E2Npoles(E_pole,kT,n_pole)
       this_n(:)%N_poles = n_pole
    end if

    ! Set the temperature
    this_n(:)%kT = kT
    this_n(:)%ckT = ' '

    ! We star-mark the defaultet contours...
    ! this will let us construct them readily...

    ! *** NOTE this is hard-coded together with the ts_contour_eq

    ! Create the left chemical potential...
    this_n(1)%name = 'Left'
    this_n(1)%mu   = Volt * 0.5_dp
    this_n(1)%cmu  = 'V/2'
    allocate(this_n(1)%Eq_seg(3)) ! one fake for poles
    this_n(1)%Eq_seg(1) = '*c-left'
    this_n(1)%Eq_seg(2) = '*t-left'
    this_n(1)%ID = 1

    this_n(2)%name = 'Right'
    this_n(2)%mu   = - Volt * 0.5_dp
    this_n(2)%cmu  = '-V/2'
    allocate(this_n(2)%Eq_seg(3)) ! one fake for poles
    this_n(2)%Eq_seg(1) = '*c-right'
    this_n(2)%Eq_seg(2) = '*t-right'
    this_n(2)%ID = 2

  end function fdffake_mu
    
  function fdf_mu(prefix,this, kT, Volt) result(found)
    use fdf
    use m_ts_io_ctype, only: pline_E_parse

    character(len=*), intent(in) :: prefix
    type(ts_mu), intent(inout) :: this
    ! The SIESTA Electronic temperature, NOT the chemical potential
    real(dp), intent(in) :: kT
    real(dp), intent(in) :: Volt
    logical :: found

    ! prepare to read in the data...
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline => null()
    logical :: info(2), bool_pole(2)
    real(dp) :: E_pole
    character(len=256) :: ln

    found = fdf_block(trim(prefix)//'.ChemPot.'//trim(Name(this)),bfdf)
    if ( .not. found ) return

    info(:) = .false.
#ifdef TBTRANS
    ! TBtrans does not need the equilbrium contour information
    info(2) = .true.
#ifdef TBT_PHONON
    ! PHtrans does not need to specify the chemical potential
    info(1) = .true.
#endif
#endif

    bool_pole = .false.
    
    ! Initialize the temperature for this chemical potential
    this%kT  = kT
    this%ckT = ' '

    do while ( fdf_bline(bfdf,pline) )
       if ( fdf_bnnames(pline) == 0 ) cycle

       ln = trim(fdf_bnames(pline,1))
       
       ! We select the input
       if ( leqi(ln,'chemical-shift') .or. &
            leqi(ln,'mu') ) then
          if ( fdf_bnvalues(pline) < 1 .and. &
               fdf_bnnames(pline) < 2 ) call die('Chemical-shift not supplied')

          ! utilize the io_ctype E-parser to grap the chemical potential
          ! we force the energies to be before the last assignment
          call pline_E_parse(pline,1,ln, &
               val = this%mu, V = Volt, kT = kT, before=3)
          this%cmu = trim(ln)

          info(1) = .true.

       else if ( leqi(ln,'contour.eq') ) then
          ! we automatically make room for one pole contour
          call read_contour_names('Equilibrium',this%Eq_seg,fakes=1)
          info(2) = .true.

       else if ( leqi(ln,'temp') .or. leqi(ln,'kT') .or. &
            leqi(ln,'Electronic.Temperature') .or. &
            leqi(ln,'ElectronicTemperature') ) then

          call pline_E_parse(pline,1,ln, &
               val = this%kT, kT = kT, before=3)
          this%ckT = trim(ln)

       else if ( leqi(ln,'contour.eq.pole.n') ) then

          if ( fdf_bnintegers(pline) < 1 ) &
               call die('You have not specified a number for &
               &number of poles.')

          this%N_poles = fdf_bintegers(pline,1)
          bool_pole(1) = .true.

       else if ( leqi(ln,'contour.eq.pole') ) then

          call pline_E_parse(pline,1,ln, &
               val = E_pole, kT = kT, before=3)

          bool_pole(2) = .true.

       else
          
          call die('Unrecognized option "'//trim(ln)//'" &
               &for chemical potential: '//trim(this%name))

       end if

    end do

    ! Check whether the Eq.Pole has been set,
    if ( bool_pole(2) ) then
       ! Update the number of poles
       call E2Npoles(E_pole, this%kT, this%N_poles)
    else if ( .not. bool_pole(1) ) then
       E_pole = fdf_get('TS.Contours.Eq.Pole',1.5_dp*eV,'Ry')
       if ( E_pole > 0._dp ) then
          call E2Npoles(E_pole,this%kT,this%N_poles)
       end if
    end if

    if ( .not. info(2) ) then
       ! The contour.eq hasn't been defined
       ! we assume the user want to use the continued
       ! fraction method
       
       ! we allow this
       allocate(this%Eq_seg(1))
       ! The equilibrium here
       this%Eq_seg(1) = 'cont-frac-'//trim(this%name)
       
       if ( bool_pole(2) ) then
          ! Update the number of poles (only up to 70% as
          ! they become sparse
          this%N_poles = int(E_pole * 0.7_dp / Pi / this%kT)
       else if ( .not. bool_pole(1) ) then
          ! Another default is 60 * Pi * kT ~ 60 poles
          ! scale for the sparse top
          E_pole = Pi * 60._dp * this%kT * 0.7_dp
          E_pole = fdf_get('TS.Contours.Eq.Pole',E_pole,'Ry')
          if ( E_pole > 0._dp ) then
             this%N_poles = int(E_pole / Pi / this%kT)
          end if
       end if

#ifndef TRANSIESTA_OVERRIDE
       if ( this%N_poles < 20 ) then
          call die('The continued fraction method requires at least 20 poles.')
       end if
#endif
       info(2) = .true.
    end if

    if ( .not. all(info) ) then
       write(*,*)'You need to supply at least:'
       write(*,*)' - chemical-shift'
       write(*,*)' - contour.eq'
       call die('You have not supplied all chemical potential information')
    end if


    if ( this%N_poles < 1 ) then
       print '(a)','Electrode: '//trim(this%name)
       print '(a,i0)',' Number of poles: ',this%N_poles
       print '(a,f10.2,a)',' Temperature: ', this%kT/Kelvin,' K'
       call die('Number of poles must be larger than or equal to 1')
    end if

  contains
    
    subroutine read_contour_names(name,con,fakes)
      character(len=*), intent(in) :: name
      character(len=C_N_NAME_LEN), intent(inout), allocatable :: con(:)
      integer, intent(in), optional :: fakes
      integer :: i, empty

      character(len=256) :: ln

      if ( allocated(con) ) deallocate(con)

      ! we need to read in the equilibrium contour
      ! skip to "begin"
      if ( .not. fdf_bline(bfdf,pline) ) &
           call die("Chemical potential block ended prematurely.")

      ! read in the begin ... end block
      ln = fdf_bnames(pline,1)
      if ( .not. leqi(ln,"begin") ) &
           call die(trim(name)//" contour errorneously formatted. &
           &First line *must* be begin!")

      ! Count lines
      i = 0
      empty = 0
      do 
         if ( .not. fdf_bline(bfdf,pline) ) &
              call die("Chemical potential block ended prematurely.")
         if ( fdf_bnnames(pline) < 1 ) then
            empty = empty + 1
         else
            ln = fdf_bnames(pline,1)
            if ( leqi(ln,'end') ) exit
            i = i + 1
         end if
      end do
      
      ! allocate names
      if ( present(fakes) ) then
         allocate(con(i+fakes))
         empty = empty - fakes
      else
         allocate(con(i))
      end if
      con = ' '
      do i = 0 , size(con) + empty
         if ( .not. fdf_bbackspace(bfdf) ) &
              call die("Backspacing too much...")
      end do
      i = 0
      do 
         if ( .not. fdf_bline(bfdf,pline) ) &
              call die("Chemical potential block ended prematurely.")
         if ( fdf_bnnames(pline) < 1 ) cycle
         ln = fdf_bnames(pline,1)
         if ( leqi(ln,'end') ) exit
         i = i + 1
         if ( len_trim(ln) > C_N_NAME_LEN ) then
            call die('Contour name: '//trim(ln)//' is too long, please use a &
                 &shorter name.')
         end if
         con(i) = trim(ln)
      end do
      
    end subroutine read_contour_names
    
  end function fdf_mu

  elemental function Eq_segs(this) result(count)
    type(ts_mu), intent(in) :: this
    integer :: count
    count = size(this%Eq_seg)
  end function Eq_segs

  subroutine chem_pot_add_Elec(this,iEl) 
    type(ts_mu), intent(inout) :: this
    integer, intent(in) :: iEl
    integer, pointer :: tmp(:), tmp2(:)

    nullify(tmp)
    allocate(tmp(this%N_El+1))

    if ( this%N_El == 0 ) then
       this%N_El = 1
       this%el => tmp
       this%el(1) = iEl
    else if ( all(this%el /= iEl) ) then
       ! copy over
       tmp(1:this%N_El) = this%el(:)
       this%N_El = this%N_El + 1
       tmp(this%N_El) = iEl
       ! clean up
       tmp2 => this%el
       this%el => tmp
       deallocate(tmp2)
    end if
  end subroutine chem_pot_add_Elec

  elemental function hasEl_i(this,iEl) result(has)
    type(ts_mu), intent(in) :: this
    integer, intent(in) :: iEl
    logical :: has
    has = this%N_El > 0
    if ( has ) has = any(this%el == iEl)
  end function hasEl_i

  elemental function hasCio(this,c_io) result(has)
    use m_ts_io_ctype, only : ts_c_io
    type(ts_mu), intent(in) :: this
    type(ts_c_io), intent(in) :: c_io 
    logical :: has
    integer :: i
    do i = 1 , Eq_segs(this)
       has = this%Eq_seg(i) .eq. c_io%name
       if ( has ) return
    end do
  end function hasCio

  elemental function hasCeq(this,c) result(has)
    use m_ts_cctype, only : ts_cw
    type(ts_mu), intent(in) :: this
    type(ts_cw), intent(in) :: c 
    logical :: has
    has = hasCio(this,c%c_io)
  end function hasCeq

  elemental function name_(this) result(name)
    type(ts_mu), intent(in) :: this
    character(len=NAME_MU_LEN) :: name
    name = this%name
  end function name_

  subroutine print_mus_block(prefix,N_mu,mus)
    use parallel, only : IONode
    character(len=*), intent(in) :: prefix
    integer, intent(in) :: N_mu
    type(ts_mu), intent(in) :: mus(N_mu)
    integer :: i

    if (IONode) then
       write(*,'(2a)') '%block ',trim(prefix)//'.ChemPots'
       do i = 1 , N_mu
          write(*,'(t3,a)') trim(mus(i)%name)
       end do
       write(*,'(2a,/)') '%endblock ',trim(prefix)//'.ChemPots'
    end if

    do i = 1 , N_mu
       call print_mu_block(prefix,mus(i))
    end do
    
  end subroutine print_mus_block

  subroutine print_mu_block(prefix,this)
    use parallel, only : IONode
    use fdf
    character(len=*), intent(in) :: prefix
    type(ts_mu), intent(in) :: this
    real(dp) :: E_pole
    integer :: i, def_pole
    character(len=64) :: ln

    if ( .not. IONode ) return

    def_pole = fdf_get('TS.Contours.Eq.Pole.N',def_poles)
    E_pole = fdf_get('TS.Contours.Eq.Pole',1.5_dp*eV,'Ry')
    if ( E_pole > 0._dp ) then
       call E2Npoles(E_pole,this%kT, def_pole)
    end if
    
    ! Start by writing out the block beginning
    write(*,'(a,a)') '%block '//trim(prefix)//'.ChemPot.',trim(this%name)
    write(*,'(t3,a,tr2,a)') 'mu',trim(this%cmu)
    write(*,'(t3,a)') 'contour.eq'
    write(*,'(t4,a)') 'begin'
    do i = 1 , Eq_segs(this) - 1 ! no pole-allowed
       ln = this%Eq_seg(i)
       if ( ln(1:1) == '*' ) then
          write(*,'(t5,a)') trim(ln(2:))
       else
          write(*,'(t5,a)') trim(ln)
       end if
    end do
    write(*,'(t4,a)') 'end'
    if ( def_pole /= this%N_poles ) then
       write(*,'(t3,a,tr2,i0)') 'contour.eq.pole.n',this%N_poles
    end if
    if ( len_trim(this%ckT) > 0 ) then
       write(*,'(t3,2a)') 'Temp ',this%ckT
    end if
    write(*,'(a,a)') '%endblock '//trim(prefix)//'.ChemPot.',trim(this%name)
    
  end subroutine print_mu_block

  subroutine E2Npoles(E_pole,kT,n_pole)
    real(dp), intent(in) :: kT, E_pole
    integer, intent(out) :: n_pole
    real(dp) :: tmp

    ! Calculate number of \pi kT
    tmp = E_pole / ( Pi * kT )
    ! Using ceiling should also force it to at least
    ! 1 pole.
    n_pole = ceiling( tmp / 2._dp )
    
  end subroutine E2Npoles

  subroutine copy_(this,copy)
    type(ts_mu), intent(inout) :: this, copy
    copy%name    = this%name
    copy%ID      = this%ID
    copy%N_poles = this%N_poles
    copy%mu      = this%mu
    copy%cmu     = this%cmu
    copy%kT      = this%kT
    copy%ckT     = this%ckT
    copy%N_El    = this%N_El
    allocate(copy%Eq_seg(size(this%Eq_seg)))
    copy%Eq_seg(:) = this%Eq_seg(:)
  end subroutine copy_

end module m_ts_chem_pot
  
