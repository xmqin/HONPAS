! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---

! This code has been fully implemented by:
! Nick Papior, 2014
!
! Please attribute the original author in case of dublication
module m_tbt_contour

  use precision, only : dp
  use fdf, only : leqi, fdf_convfac

! Use the type associated with the contour
! Maybe they should be collected to this module.
! However, I like this partition.
  use m_ts_electype

  use m_ts_chem_pot
  use m_ts_cctype
  use m_ts_io_contour
  use m_ts_io_ctype

  implicit none

  ! Contour path
  integer, save, public :: N_tbt
  type(ts_c_io),  pointer, save, public :: tbt_io(:) => null()
  type(ts_cw), pointer, save, public :: tbt_c(:) => null()

  ! The contour specific variables
  real(dp), save, public :: tbt_Eta

  ! this is heavily linked with the CONTOUR_EQ from m_ts_contour_eq
  integer, parameter, public :: CONTOUR_TBT = 5

  public :: read_contour_options
  public :: print_contour_tbt_options
  public :: print_contour_tbt_block
  public :: io_contour_tbt
  public :: N_tbt_E
  public :: tbt_E
  public :: c2weight

  private

contains

  subroutine read_contour_options(N_Elec,Elecs,N_mu,mus)

    use parallel, only : Node
    use fdf
    use m_ts_electype

    ! only to gain access to the chemical shifts
    integer, intent(in) :: N_Elec
    type(Elec), intent(in), target :: Elecs(N_Elec)
    integer, intent(in) :: N_mu
    type(ts_mu), intent(in), target :: mus(N_mu)
    
    real(dp) :: Volt, tmp, max_kT, Ry2eV
    character(len=C_N_NAME_LEN) :: ctmp
#ifdef TBT_PHONON
    integer :: i
#endif

    Ry2eV = fdf_convfac('Ry', 'eV')

    ! broadening
    tbt_Eta = fdf_get('TBT.Contours.Eta',0._dp,'Ry')
    if ( tbt_Eta < 0._dp .and. Node == 0 ) then
      call die('tbtrans: error cannot use the advanced Green function')
       write(*,'(a)')'*** NOTICE ***'
       write(*,'(a)')'tbtrans will use the advanced Green function'
    end if
    max_kT = maxval(mus(:)%kT)

    ! Get applied bias
    tmp = fdf_get('TS.Voltage',0._dp,'Ry')
    Volt = fdf_get('TBT.Voltage',tmp,'Ry')
    if ( abs(tmp-Volt) > 1.e-5_dp .and. Node == 0 ) then
       write(*,'(a)') '*** WARNING: transiesta and tbtrans bias &
            &not equivalent!'
       write(*,'(a)') '*** WARNING: Be sure to use an interpolation scheme!'
    end if

    ! Bias-window setup
    call my_setup(N_tbt,tbt_c,tbt_io,max_kT)

    ! We only allow the user to either use the old input format, or the new
    ! per-electrode input

    if ( N_tbt < 1 ) then

       ! We create the default version
       ! *** NOTE requires an even splitting of the bias
       N_tbt = 1
       nullify(tbt_io,tbt_c)
       allocate(tbt_io(N_tbt),tbt_c(N_tbt))
       tbt_c(1)%c_io => tbt_io(1)
       ctmp = 'line'
       if ( ts_exists_contour_block('TBT',' ',ctmp) ) then
          call ts_read_contour_block('TBT',' ',ctmp,tbt_io(1),max_kT,Volt)
       else
          tbt_io(1)%part = 'line'
          tbt_io(1)%name = ctmp
#ifdef TBT_PHONON
          tbt_io(1)%ca = '0. eV'
          tbt_io(1)%a  = 0._dp
          tbt_io(1)%cb = '0.5 eV'
          tbt_io(1)%b  =  0.5_dp / Ry2eV
          tbt_io(1)%cd = '0.0025 eV'
          tbt_io(1)%d = 0.0025_dp / Ry2eV
#else
          tbt_io(1)%ca = '-2. eV'
          tbt_io(1)%a  = - 2._dp / Ry2eV
          tbt_io(1)%cb = '2. eV'
          tbt_io(1)%b  =  2._dp / Ry2eV
          tbt_io(1)%cd = '0.01 eV'
          tbt_io(1)%d = 0.01_dp / Ry2eV
#endif
          tbt_io(1)%method = 'mid-rule'
       end if
       call ts_fix_contour(tbt_io(1))

       ! setup the contour
       allocate(tbt_c(1)%c(tbt_c(1)%c_io%N),tbt_c(1)%w(tbt_c(1)%c_io%N,1))
       call setup_tbt_contour(tbt_c(1), tbt_Eta)

    end if

#ifdef TBT_PHONON
    ! The half-width is squared as the phonon energy
    ! Important to do it here since the Eta should be propagated
    ! in units of Ry (and not Ry ** 2). All contours should be
    ! created with units Ry.
    tbt_Eta = tbt_Eta ** 2
    
    do i = 1 , N_tbt
       if ( tbt_io(i)%a < 0._dp ) then
          call die('Phonon transport is only defined for positive &
               &frequencies. Please correct your input.')
       end if
    end do
#endif

    ! TODO correct empty cycles, i.e. if two line contours are neighbours 
    ! then we have overlying energy points...

  contains 

    subroutine my_setup(N_tbt,tbt_c,tbt_io,max_kT)
      integer, intent(inout) :: N_tbt
      type(ts_cw), pointer :: tbt_c(:)
      type(ts_c_io), pointer :: tbt_io(:)
      real(dp), intent(in) :: max_kT

      ! Local variables
      integer :: i, j
      character(len=C_N_NAME_LEN), allocatable :: tmp(:)
      logical :: connected

      N_tbt = fdf_nc_iotype('TBT',' ')
      if ( N_tbt < 1 ) return

      allocate(tmp(N_tbt))

      tmp(1) = fdf_name_c_iotype('TBT',' ',1)
      do i = 2 , N_tbt
         tmp(i) = fdf_name_c_iotype('TBT',' ',i)
         do j = 1 , i - 1
            if ( leqi(tmp(j),tmp(i)) ) then
               call neq_die('You cannot have two names from the window &
                    &to be the same...')
            end if
         end do
      end do
      
      ! allocate all required objects
      nullify(tbt_io,tbt_c)
      allocate(tbt_io(N_tbt),tbt_c(N_tbt))
      
      do i = 1 , N_tbt

         ! assign pointer
         tbt_c(i)%c_io => tbt_io(i)
         ! read in the contour
         call ts_read_contour_block('TBT',' ',tmp(i),tbt_io(i), max_kT, Volt)

         ! transport type contour
         tbt_c(i)%c_io%type = 'tran'

      end do
      deallocate(tmp)

      do i = 1 , N_tbt - 1
         if ( i == 1 ) then
            call ts_fix_contour(tbt_io(i), next=tbt_io(i+1) , &
                 connected=connected )
         else
            call ts_fix_contour(tbt_io(i), &
                 prev=tbt_io(i-1), next=tbt_io(i+1), &
                 connected=connected )
         end if
         if ( Node == 0 .and. .not. connected ) then
            write(*,'(a)') 'tbt: *** Contours are not connected ***'
         end if
            
      end do
      if ( N_tbt > 1 ) then
         call ts_fix_contour(tbt_io(N_tbt), prev = tbt_io(N_tbt-1) )
      else
         call ts_fix_contour(tbt_io(N_tbt))
      end if

      ! setup the contour
      do i = 1 , N_tbt
         if ( tbt_c(i)%c_io%N < 1 ) then
            write(*,*) 'Contour: '//trim(tbt_c(i)%c_io%Name)//' has 0 points.'
            write(*,*) 'Please ensure at least 1 point in each contour...'
            call die('Contour number of points, invalid')
         end if

         ! allocate contour
         allocate(tbt_c(i)%c(tbt_c(i)%c_io%N),tbt_c(i)%w(tbt_c(i)%c_io%N,1))
         call setup_tbt_contour(tbt_c(i), tbt_Eta)

      end do

    end subroutine my_setup

  end subroutine read_contour_options
                                         

  ! This routine assures that we have setup all the 
  ! equilibrium contours for the passed electrode
  subroutine setup_tbt_contour(c, Eta)
    type(ts_cw), intent(inout) :: c
    real(dp), intent(in) :: Eta

    if ( leqi(c%c_io%part,'user') ) then
       
      call contour_file(c,Eta)
      
    else if ( leqi(c%c_io%part,'line') ) then
      
      call contour_line(c,Eta)
      
    else
      
      call neq_die('Unrecognized contour type for &
          &tbtrans, MUST be a line part.')
      
    end if
    
  end subroutine setup_tbt_contour

  subroutine contour_line(c,Eta)
    use m_integrate
    use m_gauss_quad
    type(ts_cw), intent(inout) :: c
    real(dp), intent(in) :: Eta

    ! local variables
    character(len=c_N) :: tmpC
    real(dp) :: a,b, tmp
    real(dp), allocatable :: ce(:), cw(:)

    if ( .not. leqi(c%c_io%part,'line') ) &
         call die('Contour is not a line')

    if ( c%c_io%N < 1 ) then
       call die('Contour: '//trim(c%c_io%Name)//' has &
            &an erroneous number of points (<1).')
    end if

    ! get bounds
    a = c%c_io%a
    b = c%c_io%b

    allocate(ce(c%c_io%N))
    allocate(cw(c%c_io%N))

    select case ( method(c%c_io%method) )
    case ( CC_MID )
       
       call Mid_Rule(c%c_io%N,ce,cw,a,b)
       
    case ( CC_SIMP_MIX )
       
       call Simpson_38_3_rule(c%c_io%N,ce,cw,a,b)
       
    case ( CC_BOOLE_MIX )
       
       call Booles_Simpson_38_3_rule(c%c_io%N,ce,cw,a,b)       
       
    case ( CC_G_LEGENDRE ) 
       
       call Gauss_Legendre_Rec(c%c_io%N,0,a,b,ce,cw)
       
    case ( CC_TANH_SINH ) 

       ! we should also gain an option for this
       if ( c_io_has_opt(c%c_io,'precision') ) then
          tmpC = c_io_get_opt(c%c_io,'precision')
          read(tmpC,'(g20.10)') tmp
       else
          tmp = 2.e-2_dp * abs(b-a) / real(c%c_io%N,dp)
          write(tmpC,'(g20.10)') tmp
          call c_io_add_opt(c%c_io,'precision',tmpC)
       end if

       call TanhSinh_Exact(c%c_io%N,ce,cw,a,b, p=tmp)

    case ( CC_USER )

      call contour_file(c,Eta)

      deallocate(ce, cw)
      return
      
    case default

       call die('Could not determine the line-integral')

    end select

    c%c = cmplx(ce,Eta, dp)
    c%w(:,1) = cmplx(cw,0._dp, dp)

    deallocate(ce,cw)
    
  end subroutine contour_line

  ! This routine will read the contour points from a file
  subroutine contour_file(c,Eta)
    use m_io, only: io_assign, io_close
    use fdf, only: fdf_convfac

    type(ts_cw), intent(inout) :: c
    ! The lifting into the complex plane
    real(dp), intent(in) :: Eta

    integer :: iu, iostat, ne
    logical :: exist
    complex(dp) :: E , W
    real(dp) :: rE, iE, rW, iW, conv
    character(len=512) :: file, line
    character(len=16) :: unit
    
    ! The contour type contains the file name in:
    !  c%c_io%cN (weirdly enough)
    file = c%c_io%cN

    call io_assign(iu)
    
    ! Open the contour file
    inquire(file=trim(file), exist=exist)
    if ( .not. exist ) then
      call die('The file: '//trim(file)//' could not be found to read contour points!')
    end if
    
    ! Open the file
    open(iu, file=trim(file), form='formatted', status='old')
    
    ne = 0
    ! The default unit is eV.
    ! On every line an optional unit-specificer may be used to specify the
    ! subsequent lines units (until another unit is specified)
    conv = fdf_convfac('eV', 'Ry')
    do
      ! Now we have the line
      read(iu, '(a)', iostat=iostat) line
      if ( iostat /= 0 ) exit
      if ( len_trim(line) == 0 ) cycle
      line = trim(adjustl(line))
      if ( line(1:1) == '#' ) cycle

      ! We have a line with energy and weight
      ne = ne + 1
      ! There are three optional ways of reading this
      ! 1.  ReE ImE, ReW ImW [unit]
      read(line, *, iostat=iostat) rE, iE, rW, iW, unit
      if ( iostat == 0 ) then
        conv = fdf_convfac(unit, 'Ry')
      else
        read(line, *, iostat=iostat) rE, iE, rW, iW
      end if
      if ( iostat == 0 ) then
        c%c(ne) = cmplx(rE,iE, dp) * conv
        c%w(ne,1) = cmplx(rW,iW, dp) * conv
        cycle
      end if

      ! 2.  ReE ImE, ReW [unit]
      iW = 0._dp
      read(line, *, iostat=iostat) rE, iE, rW, unit
      if ( iostat == 0 ) then
        conv = fdf_convfac(unit, 'Ry')
      else
        read(line, *, iostat=iostat) rE, iE, rW
      end if
      if ( iostat == 0 ) then
        c%c(ne) = cmplx(rE,iE, dp) * conv
        c%w(ne,1) = cmplx(rW,iW, dp) * conv
        cycle
      end if

      ! 3.  ReE , ReW [unit]
      iE = Eta
      iW = 0._dp
      read(line, *, iostat=iostat) rE, rW, unit
      if ( iostat == 0 ) then
        conv = fdf_convfac(unit, 'Ry')
      else
        read(line, *, iostat=iostat) rE, rW
      end if
      if ( iostat == 0 ) then
        c%c(ne) = cmplx(rE * conv,iE, dp)
        c%w(ne,1) = cmplx(rW,iW, dp) * conv
        cycle
      end if

      call die('Contour file: '//trim(file)//' is not formatted correctly. &
          &Please read the documentation!')

    end do

    call io_close(iu)

    if ( c%c_io%N /= ne ) then
      call die('Error in reading the contour points from file: '//trim(file))
    end if
    
  end subroutine contour_file

  function TBT_E(id,step) result(c)
    integer, intent(in) :: id
    integer, intent(in), optional :: step
    type(ts_c_idx) :: c ! the configuration of the energy-segment
    integer :: lstep, i, PN
    lstep = 1
    if ( present(step) ) lstep = step
    PN = N_tbt_E()
    if ( id <= PN ) then
       c = get_c(id)
       return
    end if
    c = get_c(-1)
    i = MOD(PN,lstep)
    if ( i /= 0 ) PN = PN + lstep - i
    if ( id <= PN ) then
       c%exist = .true.
       c%fake  = .true.
    end if
  end function TBT_E

  function get_c(id) result(c)
    integer, intent(in) :: id
    type(ts_c_idx) :: c
    integer :: i,j,iE
    c%exist = .false.
    c%fake  = .false.
    c%e     = cmplx(0._dp,0._dp, dp)
    c%idx   = 0
    if ( id < 1 ) return

    iE = 0
    do j = 1 , N_tbt ! number of contours
       if ( iE + tbt_c(j)%c_io%N < id ) then
          iE = iE + tbt_c(j)%c_io%N
          cycle
       end if
       i = id - iE
       if ( i <= tbt_c(j)%c_io%N ) then
          c%exist = .true.
          c%e     = tbt_c(j)%c(i)
          c%idx(1) = 2 ! designates the non-equilibrium contours
          c%idx(2) = j ! designates the index of the non-equilibrium contour
          c%idx(3) = i ! is the index of the non-equilibrium contour
          return
       end if
    end do

  end function get_c

  subroutine c2weight(c,W)
    type(ts_c_idx), intent(in) :: c
    real(dp), intent(out) :: W
    ! local variables
    type(ts_cw), pointer :: cw

    cw => tbt_c(c%idx(2))
    W = cw%w(c%idx(3),1)

  end subroutine c2weight
  
  function N_TBT_E() result(N)
    integer :: N, i
    N = 0
    do i = 1 , N_tbt
       N = N + size(tbt_c(i)%c)
    end do
  end function N_TBT_E

  subroutine print_contour_tbt_block(prefix)
    use parallel, only : IONode
    character(len=*), intent(in) :: prefix

    integer :: i

    if ( IONode ) then
       write(*,'(2a)') '%block ',trim(prefix)//'.Contours'
       do i = 1 , N_tbt
          write(*,'(tr4,a)') trim(tbt_io(i)%name)
       end do
       write(*,'(2a,/)') '%endblock ',trim(prefix)//'.Contours'
    end if

    do i = 1 , N_tbt
       call ts_print_contour_block(trim(prefix)//'.Contour.',tbt_io(i))
    end do

  end subroutine print_contour_tbt_block


  subroutine print_contour_tbt_options(prefix)

    use parallel, only : IONode
    use m_ts_io_contour

    character(len=*), intent(in) :: prefix
    character(len=200) :: chars
    integer :: i
    type(ts_c_opt_ll), pointer :: opt
    real(dp) :: Ry2eV

    if ( .not. IONode ) return

    Ry2eV = fdf_convfac('Ry', 'eV')
    
    write(*,opt_n) '             >> TBtrans contour << '
#ifdef TBT_PHONON
    write(*,opt_g_u) 'Device Green function imaginary Eta', &
         tbt_Eta*Ry2eV**2,'eV**2'
#else
    write(*,opt_g_u) 'Device Green function imaginary Eta',tbt_Eta*Ry2eV,'eV'
#endif
    do i = 1 , N_tbt
       chars = '  '//trim(tbt_io(i)%part)
       write(*,opt_c) 'Contour name',trim(prefix)// &
            '.Contour.'//trim(tbt_io(i)%name)
       call write_e(trim(chars)//' contour E_min',tbt_io(i)%a)
       call write_e(trim(chars)//' contour E_max',tbt_io(i)%b)
       write(*,opt_int) trim(chars)//' contour points',tbt_io(i)%N
       write(*,opt_c) trim(chars)//' contour method', &
            trim(longmethod2str(tbt_io(i)))
       opt => tbt_io(i)%opt
       do while ( associated(opt) )
          write(*,opt_c) '   Option for contour method', trim(opt%opt)
          opt => opt%next
       end do
    end do

  end subroutine print_contour_tbt_options

  subroutine io_contour_tbt(slabel)
    use parallel, only : IONode
    use m_tbt_save, only : name_save
    character(len=*), intent(in) :: slabel

! *********************
! * LOCAL variables   *
! *********************
    integer :: iu, i
    type(ts_c_idx) :: cidx
    character(len=200) :: fname

    if ( .not. IONode ) return

    ! Get save name
    call name_save(1,1,fname,'CC')

    call io_assign( iu )
    open( iu, file=trim(fname), status='unknown' )
    write(iu,'(a)') '# Contour path for the transport part'
    write(iu,'(a,a24,3(tr1,a25))') '#','Re(c) [eV]','Im(c) [eV]','w [eV]'

    cidx%idx(1) = CONTOUR_TBT

    do i = 1 , N_tbt
       
       cidx%idx(2) = i
       call io_contour_c(iu,cidx)
       
    end do

    call io_close( iu )

  end subroutine io_contour_tbt

  subroutine io_contour_c(iu,cidx)
    use fdf, only: fdf_convfac
    use m_ts_aux, only : nf
    integer, intent(in) :: iu
    type(ts_c_idx), intent(inout) :: cidx
    type(ts_cw), pointer :: c
    integer :: i
    real(dp) :: Ry2eV
    
    Ry2eV = fdf_convfac('Ry', 'eV')

    if ( cidx%idx(1) == CONTOUR_TBT ) then
       c => tbt_c(cidx%idx(2))
    else
       call die('io_contour_c: Error in code')
    end if

    do i = 1 , size(c%c)
       write(iu,'(3(e25.17,tr1))') c%c(i) * Ry2eV, real(c%w(i,1),dp) * Ry2eV
    end do

  end subroutine io_contour_c

  subroutine neq_die(msg)
    character(len=*), intent(in) :: msg
    
    write(*,*) 'Killing... printing out so-far gathered information'
    call print_contour_tbt_options('TBT')
    call print_contour_tbt_block('TBT')
    call die(msg)

  end subroutine neq_die
    
end module m_tbt_contour
