!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2015, nickpapior@gmail.com

module m_ts_contour_neq

  use precision, only : dp

! Use the type associated with the contour
! Maybe they should be collected to this module.
! However, I like this partition.
  use m_ts_electype

  use m_ts_chem_pot
  use m_ts_cctype
  use m_ts_io_contour
  use m_ts_io_ctype
  use m_ts_aux

  implicit none

  ! The non-equilibrium density integration are attributed discussions with
  ! Antti-Pekka Jauho. 

  ! For further attributions see the original paper by Brandbyge, et. al, 2002: DOI: 10.1103/PhysRevB.65.165401
  
  ! The non-equilibrium contour IO-segments
  integer, save, public :: N_nEq
  type(ts_c_io), pointer, save, public :: nEq_io(:) => null()
  type(ts_cw)  , pointer, save, public :: nEq_c(:) => null()

  ! Instead of defining several separate contours for
  ! each part, we define the bias window contour and
  ! define a cut-off radius in units of kT to determine
  ! the range of the Fermi-functions.
  ! It is not allowed to be less than 5 kT
  real(dp), save, public :: cutoff_kT = 5._dp

  ! The contour specific variables
  real(dp), save, public :: nEq_Eta = 0._dp

  ! A table type for containing the ID's of different
  ! segments
  type :: ts_nEq_ID
     ! This corresponds to this correction:
     !   \Delta = G * \Gamma_El * G^\dagger * \
     !        ( nF(E - El_\mu) - nF(E - \mu) )
     ! And hence this correction belongs to the \mu equilibrium
     ! contour.

     ! The index of the equilibrium chemical
     ! potential that this contour ID has correction
     ! contributions for
     type(ts_mu), pointer :: mu => null()
     ! The index for the electrode attached to this
     ! integration quantity
     type(Elec), pointer :: El => null()
     ! The index associated with the non-equilibrium 
     ! array index as calculated in the tri|full|mumps[gk] routines
     integer :: ID = 0

     ! Min and Max energy allowed for this ID to contain
     ! an energy point
     real(dp) :: E(2)

  end type ts_nEq_ID
  integer, save, public :: N_nEq_id = 0
  type(ts_nEq_id), allocatable, save, public :: nEq_ID(:)

  ! this is heavily linked with the CONTOUR_EQ from m_ts_contour_eq
  integer, parameter, public :: CONTOUR_NEQ = 2

  public :: read_contour_neq_options
  public :: print_contour_neq_options
  public :: print_contour_neq_block
  public :: contour_nEq_warnings
  public :: io_contour_neq
  public :: N_nEq_E
  public :: nEq_E
  public :: has_cE_nEq
  public :: c2weight_neq
  public :: ID2mu

  private

contains

  subroutine read_contour_neq_options(N_Elec,Elecs,N_mu,mus,kT,Volt)

    use fdf
    use m_ts_electype

    ! only to gain access to the chemical shifts
    integer, intent(in) :: N_Elec
    type(Elec), intent(in), target :: Elecs(N_Elec)
    integer, intent(in) :: N_mu
    type(ts_mu), intent(in), target :: mus(N_mu)
    ! This temperature is the SIESTA electronic temperature
    real(dp), intent(in) :: kT, Volt
    
    integer :: i, j, k, N, cur_mu
    real(dp) :: min_E, max_E, rtmp, rtmp2, Ry2eV
    logical :: err

    ! Notify the user of obsolete keys.
    call fdf_obsolete('TS.biasContour.Eta')
    Ry2eV = fdf_convfac('Ry', 'eV')

    ! check that we in-fact have a bias calculation
    if ( N_mu < 2 ) then
       call die('Something has gone wrong. We can only find one chemical potential')
    end if

    ! broadening defaults to the electrodes Eta values and their
    ! propagation. However, the user can denote an eta value in the 
    ! device region as well.
    nEq_Eta = fdf_get('TS.Contours.nEq.Eta',0._dp,'Ry')
    if ( nEq_Eta < 0._dp ) call die('ERROR: nEq_Eta < 0, we do not allow &
        &for using the advanced Green function, please correct.')
    if ( nEq_Eta == 0._dp ) then
      do i = 1, N_Elec
        if ( Elecs(i)%Eta <= 0._dp ) then
          call die('ERROR: A non-equilibrium contour requires a positive finite imaginary &
              &number to ensure we do not invert an eigenvalue energy, please correct eta values!')
        end if
      end do
    end if

    ! Temperature cut-off for the tail integrals
    cutoff_kT = fdf_get('TS.Contours.nEq.Fermi.Cutoff',5._dp)
    if ( cutoff_kT < 2.5_dp ) then
       call die('The cutoff must be larger than or equal to 5 as the &
            &Fermi-function still has a weight at Ef + 2.5 kT.')
    end if

    ! We only allow the user to either use the old input format, or the new
    ! per-electrode input

    ! Bias-window setup
    call my_setup('nEq',N_nEq,nEq_c,nEq_io)
    if ( N_nEq < 1 ) then

      ! Find the minimum integration energy and the maximum
      j = 1
      k = 1
      min_E = mus(j)%mu - mus(j)%kT * cutoff_kT
      max_E = mus(k)%mu + mus(k)%kT * cutoff_kT
      do i = 1, N_mu
        if ( mus(i)%mu - mus(i)%kT * cutoff_kT < min_E ) then
          min_E = mus(i)%mu - mus(i)%kT * cutoff_kT
          j = i
        end if
        if ( max_E < mus(i)%mu + mus(i)%kT * cutoff_kT ) then
          max_E = mus(i)%mu + mus(i)%kT * cutoff_kT
          k = i
        end if
      end do
      
      ! We create the default version
      N_nEq = 1
      nullify(nEq_io,nEq_c)
      allocate(nEq_io(N_nEq),nEq_c(N_nEq))
      
      ! Create the integrals
      nEq_io(1)%part = 'line'
      nEq_io(1)%name = 'line'
      write(nEq_io(1)%ca,'(a,g12.5,a)') &
          '- ',cutoff_kT * mus(j)%kT / kT ,' kT - |V|/2'
      ! mus(j)%mu is negative
      nEq_io(1)%a  = - cutoff_kT * mus(j)%kT + mus(j)%mu
      write(nEq_io(1)%cb,'(a,g12.5,a)') &
          '|V|/2 + ',cutoff_kT * mus(k)%kT / kT,' kT'
      nEq_io(1)%b  = mus(k)%mu + cutoff_kT * mus(k)%kT
      nEq_io(1)%cd = '0.01 eV'
      nEq_io(1)%d = 0.01_dp / Ry2eV
      nEq_io(1)%method = 'mid-rule'

      do i = 1 , N_nEq
        nEq_c(i)%c_io => nEq_io(i)

        ! Fix the contours checking for their consistency
        call ts_fix_contour(nEq_io(i))

        ! setup the contour
        allocate(nEq_c(i)%c(nEq_c(i)%c_io%N),nEq_c(i)%w(nEq_c(i)%c_io%N,1))

        call setup_nEq_contour(nEq_c(i), nEq_Eta)

      end do

    end if

    ! At this point we have created all contours
    ! We need to check that the cut-off is obeyed.
    ! We check this by finding the max and min energy.
    min_E = huge(1._dp)
    max_E = -huge(1._dp)
    do i = 1 , N_nEq
       if ( nEq_c(i)%c_io%a < min_E ) then
          min_E = nEq_c(i)%c_io%a
       end if
       if ( max_E < nEq_c(i)%c_io%b ) then
          max_E = nEq_c(i)%c_io%b
       end if
    end do

    ! Check that no chemical potential lies more than cut-off from it
    err = .false.
    do i = 1 , N_mu
       rtmp = mus(i)%mu - cutoff_kT * mus(i)%kT
       if ( min_E - rtmp > 0.000001_dp ) then
          err = .true.
          write(*,'(a)')'transiesta: Lowest energy and chemical potential &
               &'//trim(mus(i)%name)//' does not conform with Fermi-function &
               &cut-off [eV].'
          write(*,'(a,g12.5)')'Specified cut-off: ',rtmp * Ry2eV
          write(*,'(a,g12.5)')'Minimum energy in bias contour: ',min_E * Ry2eV
          write(*,'(a)')'Assert that all contours have at least 2 points to &
               &contain both end-points.'
       end if
       rtmp = mus(i)%mu + cutoff_kT * mus(i)%kT
       if ( rtmp - max_E > 0.000001_dp ) then
          err = .true.
          write(*,'(a)')'transiesta: Highest energy and chemical potential &
               &'//trim(mus(i)%name)//' does not conform with Fermi-function &
               &cut-off [eV].'
          write(*,'(a,g12.5)')'Specified cut-off: ',rtmp * Ry2eV
          write(*,'(a,g12.5)')'Maximum energy in bias contour: ',max_E * Ry2eV
          write(*,'(a)')'Assert that all contours have at least 2 points to &
               &contain both end-points.'
       end if
    end do
    if ( err ) &
         call neq_die('The energy grid is not fulfilling the requirements &
         &please see the output.')

    ! 2 chemical potentials => 1 segment
    ! 3 chemical potentials => 3 segments
    ! 4 chemical potentials => 6 segments, etc.
    !N_nEq_segs = 0 ! count
    !do j = N_mu - 1 , 1 , -1
    !   N_nEq_segs = N_nEq_segs + j
    !end do

    ! count the number of non-equilibrium segments concerning the
    ! electrodes density matrix update
    ! 2 electrodes, 2 mu => 2
    ! 3 electrodes, 2 mu => 3
    ! 4 electrodes, 2 mu => 4
    ! 5 electrodes, 2 mu => 5 , etc.
    ! 3 electrodes, 3 mu => 6
    ! 6 electrodes, 2 mu => 6
    ! 7 electrodes, 2 mu => 7
    ! 4 electrodes, 3 mu => 8
    ! 5 electrodes, 3 mu => 10
    ! 4 electrodes, 4 mu => 12
    ! etc.
    ! For each chemical potential we need the contribution
    ! from each other electrode with different chemical potential
    N_nEq_ID = 0
    do i = 1 , N_mu
       N_nEq_ID = N_nEq_ID + N_Elec - mus(i)%N_El
    end do
    if ( N_nEq_ID < 2 ) then
       call neq_die('Could not find enough non-equilibrium segments. &
            &Please correct/could be a bug...')
    end if

    allocate(nEq_ID(N_nEq_ID))
    ! Create all ID's
    cur_mu = 1
    N = 0 ! current reached ID
    do j = 1 , N_mu - 1
       do i = j + 1 , N_mu 
          do k = 1 , mus(i)%N_El 
             N = N + 1 
             nEq_ID(N)%ID = N
             nEq_ID(N)%mu => mus(j)
             nEq_ID(N)%El => Elecs(mus(i)%El(k))
          end do
          do k = 1 , mus(j)%N_El 
             N = N + 1 
             nEq_ID(N)%ID = N
             nEq_ID(N)%mu => mus(i)
             nEq_ID(N)%El => Elecs(mus(j)%El(k))
          end do
       end do
    end do
    if ( N /= N_nEq_ID ) then
       call die('Could not find all available non-equilibrium &
            &IDs, bug in code')
    end if

    ! Update the min/max energies allowed
    ! to contain any contour elements for each ID
    ! We add an extra 0.0000001 eV to account for floating
    ! point accurracy
    do N = 1 , N_nEq_ID
       ! Figure out the actual minimum/maximum 
       ! by taking into account the cut-off Fermi function value
       rtmp  = nEq_ID(N)%mu%mu    - cutoff_kT * nEq_ID(N)%mu%kT
       rtmp2 = nEq_ID(N)%El%mu%mu - cutoff_kT * nEq_ID(N)%El%mu%kT
       min_E = min(rtmp,rtmp2)
       rtmp  = nEq_ID(N)%mu%mu    + cutoff_kT * nEq_ID(N)%mu%kT
       rtmp2 = nEq_ID(N)%El%mu%mu + cutoff_kT * nEq_ID(N)%El%mu%kT
       max_E = max(rtmp,rtmp2)
       nEq_ID(N)%E(1) = min_E - 0.0000001_dp / Ry2eV
       nEq_ID(N)%E(2) = max_E + 0.0000001_dp / Ry2eV
    end do

  contains 

    subroutine my_setup(suffix,N_nEq,nEq_c,nEq_io)
      character(len=*), intent(in) :: suffix
      integer, intent(inout) :: N_nEq
      type(ts_cw), pointer :: nEq_c(:)
      type(ts_c_io), pointer :: nEq_io(:)

      ! Local variables
      logical :: connected
      integer :: i, j
      character(len=C_N_NAME_LEN), allocatable :: tmp(:)

      N_nEq = fdf_nc_iotype('TS',suffix)
      if ( N_nEq < 1 ) return

      allocate(tmp(N_nEq))

      tmp(1) = fdf_name_c_iotype('TS',suffix,1)
      do i = 2 , N_nEq
         tmp(i) = fdf_name_c_iotype('TS',suffix,i)
         do j = 1 , i - 1
            if ( leqi(tmp(j),tmp(i)) ) then
               call neq_die('You cannot have two names from the bias-window &
                    &to be the same...')
            end if
         end do
      end do
      
      ! allocate all required objects
      nullify(nEq_io,nEq_c)
      allocate(nEq_io(N_nEq),nEq_c(N_nEq))
      
      do i = 1 , N_nEq

         ! assign pointer
         nEq_c(i)%c_io => nEq_io(i)
         ! read in the contour
         call ts_read_contour_block('TS',suffix,tmp(i),nEq_io(i), kT, Volt)

         ! Assign non-equilibrium contour
         nEq_c(i)%c_io%type = 'neq'
         
      end do
      deallocate(tmp)

      do i = 1 , N_nEq - 1
         if ( i == 1 ) then
            call ts_fix_contour(nEq_io(i), next=nEq_io(i+1), &
                 connected = connected )
         else
            call ts_fix_contour(nEq_io(i), &
                 prev=nEq_io(i-1), next=nEq_io(i+1), &
                 connected = connected )
         end if
         if ( .not. connected ) then
            call die('Contour: '//trim(nEq_io(i)%name)// &
                 ' and '//trim(nEq_io(i+1)%name)//' are not connected.')
         end if
      end do
      ! The fix_contour checks and corrects the neighbouring 
      ! contours.
      if ( N_nEq > 1 ) then
         call ts_fix_contour(nEq_io(N_nEq),prev=nEq_io(N_nEq-1))
      else
         call ts_fix_contour(nEq_io(N_nEq))
      end if

      ! setup the contour
      do i = 1 , N_nEq
         if ( nEq_c(i)%c_io%N < 1 ) then
            write(*,*) 'Contour: '//trim(nEq_c(i)%c_io%Name)//' has 0 points.'
            write(*,*) 'Please ensure at least 1 point in each contour...'
            call die('Contour number of points, invalid')
         end if

         ! allocate contour
         allocate(nEq_c(i)%c(nEq_c(i)%c_io%N),nEq_c(i)%w(nEq_c(i)%c_io%N,1))
         call setup_nEq_contour(nEq_c(i), nEq_Eta)
      end do

      if ( nEq_c(1)%c_io%a > nEq_c(N_nEq)%c_io%b ) then
         call neq_die('The non-equilibrium contours must be in increasing &
              &energy. Even if your bias is negative. Please correct.')
      end if

    end subroutine my_setup

  end subroutine read_contour_neq_options

  subroutine contour_nEq_warnings( ) 

    use parallel, only : IONode, Nodes

    integer :: i, imu, not_calc, N
    complex(dp) :: W, ZW
    type(ts_c_idx) :: cE
    logical :: err

    ! Print out a warning if there are any points that does not
    ! contribute to any chemical potentials
    not_calc = 0
    do i = 1 , N_nEq_E()

       ! Get the current energy point
       cE = nEq_E(i)
       
       err = .false.
       do N = 1 , N_nEq_ID
          call c2weight_nEq(cE,N,1._dp,W,imu,ZW)
          if ( imu > 0 ) err = .true.
          if ( err ) exit
       end do
       if ( .not. err ) not_calc = not_calc + 1
    end do
    if ( IONode .and. not_calc > 0 ) then
       write(*,'(a,i0,a)') &
            '*** You have ',not_calc,' unused non-equilibrium contour &
            &points which degrades performance considerably.'
       write(*,'(a)') &
            '    Consider correcting your nEq contours.'
    end if

    if ( IONode .and. not_calc == 0 ) then
       i = mod(N_nEq_E(), Nodes) ! get remaining part of equilibrium contour
       if ( IONode .and. i /= 0 ) then
          i = Nodes - i
          write(*,'(a)')'Without loosing performance you can increase &
               &the non-equilibrium integration precision.'
          write(*,'(a,i0,a)')'You can add ',i,' more energy points in the &
               &non-equilibrium contours, for FREE!'
       end if
    end if

  end subroutine contour_nEq_warnings

  ! This routine assures that we have setup all the 
  ! equilibrium contours for the passed electrode
  subroutine setup_nEq_contour(c, Eta)
    use fdf, only: leqi
    type(ts_cw), intent(inout) :: c
    real(dp), intent(in) :: Eta

    if ( leqi(c%c_io%part,'user') ) then
       
      call contour_file(c,Eta)

    else if ( leqi(c%c_io%part,'line') ) then

      call contour_line(c,Eta)

    else if ( leqi(c%c_io%part,'tail') ) then

      c%c_io%part = 'line'

      call contour_line(c,Eta)

      c%c_io%part = 'tail'

    else

      call neq_die('Unrecognized contour type for the &
          &non-equilibrium part.')
       
    end if
    
  end subroutine setup_nEq_contour

  function has_cE_nEq(cE,iEl,ID) result(has)
    type(ts_c_idx), intent(in) :: cE
    integer, intent(in) :: iEl, ID
    logical :: has
    real(dp) :: E
    has = .false.
    if ( cE%fake ) return
    if ( cE%idx(1) /= CONTOUR_NEQ ) return
    
    has = nEq_ID(ID)%El%ID == iEl

    if ( has ) then

       ! Retrieve the real energy of the contour index
       E = real(nEq_c(cE%idx(2))%c(cE%idx(3)),dp)

       ! We need to check that the energy is within the cutoff
       ! range of kT from either chemical potentials
       if ( E < nEq_ID(ID)%E(1) ) has = .false.
       if ( nEq_ID(ID)%E(2) < E ) has = .false.

    end if

  end function has_cE_nEq

  subroutine c2weight_nEq(c,ID,k,W,imu,ZW)
    use m_ts_aux, only : nf
    type(ts_c_idx), intent(in) :: c
    integer, intent(in) :: ID ! the non-equilibrium contour index
    real(dp), intent(in) :: k ! generic weight
    integer, intent(out) :: imu ! returns the equilibrium idx for EDM
    complex(dp), intent(out) :: W, ZW ! the weight returned
    ! local variables
    real(dp) :: E
    type(ts_cw), pointer :: cw

    imu = 0
    W  = 0._dp
    ZW = 0._dp

    ! Get the contour points and weights
    cw => nEq_c(c%idx(2))

    ! Retrieve the real energy of the contour index
    E = real(cw%c(c%idx(3)),dp)

    ! We need to check that the energy is within the cutoff
    ! range of kT from either chemical potentials
    if ( E < nEq_ID(ID)%E(1) ) return
    if ( nEq_ID(ID)%E(2) < E ) return

    ! Update the equilibrium contribution ID
    imu = nEq_ID(ID)%mu%ID

    ! Notice that we multiply with -i to move Gf.Gamma.Gf^\dagger
    ! to the negative imaginary part
    ! This means that we do not need to create different
    ! routines for the equilibrium density and the non-equilibrium
    ! density

    ! nf function is: nF(E-E1) - nF(E-E2) IMPORTANT
    W = k * real(cw%w(c%idx(3),1),dp) * &
         nf(E, &
         nEq_ID(ID)%El%mu%mu, nEq_ID(ID)%El%mu%kT, &
         nEq_ID(ID)%mu%mu, nEq_ID(ID)%mu%kT )

    ZW = E * W

  end subroutine c2weight_nEq

  function ID2mu(ID) result(imu)
    integer, intent(in) :: ID
    integer :: imu
    imu = nEq_ID(ID)%mu%ID
  end function ID2mu

  subroutine contour_line(c,Eta)
    use fdf, only: leqi
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
            &an errorneous number of points.')
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

       call contour_file(c, Eta)

       ! Immediately return as the user has specified *everything*
       deallocate(ce, cw)
       return
       
    case default

       call die('Could not determine the line-integral')

    end select

    c%c(:) = cmplx(ce, Eta, dp)
    c%w(:,1) = cmplx(cw, 0._dp, dp)

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

  function nEq_E(id,step) result(c)
    integer, intent(in) :: id
    integer, intent(in), optional :: step
    type(ts_c_idx) :: c ! the configuration of the energy-segment
    integer :: lstep, i, PN
    lstep = 1
    if ( present(step) ) lstep = step
    PN = N_nEq_E()
    c = get_c(id)
    if ( id <= PN ) return
    i = MOD(PN,lstep)
    if ( i /= 0 ) PN = PN + lstep - i
    if ( id <= PN ) then
       c%exist = .true.
       c%fake  = .true.
    end if
  end function nEq_E

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
    do j = 1 , N_nEq ! number of contours
       if ( iE + nEq_c(j)%c_io%N < id ) then
          iE = iE + nEq_c(j)%c_io%N
          cycle
       end if
       i = id - iE
       if ( i <= nEq_c(j)%c_io%N ) then
          c%exist = .true.
          c%e     = nEq_c(j)%c(i)
          c%idx(1) = CONTOUR_NEQ ! designates the non-equilibrium contours
          c%idx(2) = j ! designates the index of the non-equilibrium contour
          c%idx(3) = i ! is the index of the non-equilibrium contour
          return
       end if
    end do

  end function get_c

  function N_nEq_E() result(N)
    integer :: N, i
    N = 0
    do i = 1 , N_nEq
       N = N + size(nEq_c(i)%c)
    end do
  end function N_nEq_E

  subroutine print_contour_neq_block(prefix)
    use parallel, only : IONode
    character(len=*), intent(in) :: prefix

    integer :: i

    if ( IONode ) then
       write(*,'(2a)') '%block ',trim(prefix)//'.Contours.nEq'
       do i = 1 , N_nEq
          write(*,'(tr4,a)') trim(nEq_io(i)%name)
       end do
       write(*,'(2a,/)') '%endblock ',trim(prefix)//'.Contours.nEq'
    end if

    do i = 1 , N_nEq
       call ts_print_contour_block(trim(prefix)//'.Contour.nEq.',nEq_io(i))
    end do

  end subroutine print_contour_neq_block


  subroutine print_contour_neq_options(prefix)

    use parallel, only : IONode
    use fdf, only : fdf_convfac
    use m_ts_io_contour

    character(len=*), intent(in) :: prefix
    character(len=200) :: chars
    integer :: i
    type(ts_c_opt_ll), pointer :: opt
    real(dp) :: Ry2eV

    if ( .not. IONode ) return

    Ry2eV = fdf_convfac('Ry', 'eV')
    
    write(*,opt_n) '        >> non-Equilibrium contour << '
    write(*,opt_g_u) 'Device Green function imaginary Eta',nEq_Eta*Ry2eV,'eV'
    write(*,opt_g_u) 'Fermi-function cut-off',cutoff_kT,'kT'
    do i = 1 , N_nEq
       chars = '  '//trim(nEq_io(i)%part)
       write(*,opt_c) 'Contour name',trim(prefix)// &
            '.Contour.nEq.'//trim(neq_io(i)%name)
       call write_e(trim(chars)//' contour E_min',neq_io(i)%a)
       call write_e(trim(chars)//' contour E_max',neq_io(i)%b)
       write(*,opt_int) trim(chars)//' contour points',neq_io(i)%N
       write(*,opt_c) trim(chars)//' contour method', &
            trim(longmethod2str(neq_io(i)))
       opt => neq_io(i)%opt
       do while ( associated(opt) )
          write(*,opt_c) '   Option for contour method', trim(opt%opt)
          opt => opt%next
       end do
    end do
    
  end subroutine print_contour_neq_options

  subroutine io_contour_nEq(slabel,suffix)
    use parallel, only : IONode
    character(len=*), intent(in) :: slabel
    character(len=*), intent(in), optional :: suffix

! *********************
! * LOCAL variables   *
! *********************
    integer :: i

    if ( .not. IONode ) return

    do i = 1 , N_nEq_ID
       call io_contour_nEq_ID(nEq_ID(i),slabel,suffix)
    end do

  end subroutine io_contour_nEq


  subroutine io_contour_nEq_ID(nEq_ID,slabel,suffix)
    use parallel, only : IONode
    use units, only : Kelvin
    use fdf, only: fdf_convfac
    type(ts_nEq_ID), intent(in) :: nEq_ID
    character(len=*), intent(in) :: slabel
    character(len=*), intent(in), optional :: suffix

! *********************
! * LOCAL variables   *
! *********************
    character(len=256) :: fname
    integer        :: i, unit
    type(ts_c_idx) :: c
    real(dp) :: m1, m2, Ry2eV
    character(len=32) :: cm1, cm2, kT1, kT2
    complex(dp) :: W, ZW, Z
    integer :: imu
    
    if ( .not. IONode ) return
    Ry2eV = fdf_convfac('Ry', 'eV')

    ! format of file is:
    ! <slabel>.TSCCNEQ-<ID>-<correction chemical potential>_<from electrode>
    ! The ID is important since it corresponds to the "internal" order of chemical potentials.
    if ( present(suffix) ) then
      write(fname,'(a,".",a,"-",a,"_",a)') trim(slabel), trim(suffix), trim(nEq_ID%mu%name), trim(nEq_ID%El%name)
    else
      write(fname,'(a,".",a,"-",i0,"-",a,"_",a)') trim(slabel), 'TSCCNEQ', nEq_ID%ID, &
          trim(nEq_ID%mu%name), trim(nEq_ID%El%name)
    end if

    call io_assign( unit )
    open( unit, file=fname, status='unknown' )
    write(unit,'(a)') '# Contour path for the non-equilibrium part'
    write(unit,'(2a)')'# Correction for chemical potential ',trim(nEq_ID%mu%name)
    write(unit,'(2a)')'# Scattering states coming from ',trim(nEq_ID%El%name)
    
    m1 = nEq_ID%El%mu%mu * Ry2eV
    if ( m1 < 0._dp ) then
       write(cm1,'(a,g10.4)')'+ ',-m1
    else
       write(cm1,'(a,g10.4)')'- ',m1
    end if
    write(kT1,'(g10.4)') nEq_ID%El%mu%kT / Kelvin
    m2 = nEq_ID%mu%mu * Ry2eV
    if ( m2 < 0._dp ) then
       write(cm2,'(a,g10.4)')'+ ',-m2
    else
       write(cm2,'(a,g10.4)')'- ',m2
    end if
    write(kT2,'(g10.4)') nEq_ID%mu%kT / Kelvin
    write(unit,'(a,/,9a)') '# Window Fermi function: ', &
         '#   nF(E ', trim(cm1),' eV, ',trim(kT1),' K) &
         &- nF(E ',trim(cm2),' eV, ',trim(kT2),' K)'
    write(unit,'(a,a24,2(tr1,a25))') '#','Re(c) [eV]','Im(c) [eV]','w [eV]'
    
    do i = 1 , N_nEq_E()

       ! Retrieve energy-point
       c = nEq_E(i,1)

       ! Check if it exists on this ID
       call c2weight_nEq(c,nEq_ID%ID,1._dp,W,imu,ZW)
       if ( imu < 1 ) cycle
       
       Z = nEq_c(c%idx(2))%c(c%idx(3))
       ! The imaginary part reflects the imaginary part in the
       ! electrode self-energy.
       ! Since the energy-points are duplicated for different
       ! non-equilibrium parts.
       ! If the electrode eta value is negative it refers to the device
       ! eta value.
       if ( nEq_ID%El%Eta > 0._dp ) then
         Z = cmplx(real(Z, dp), nEq_ID%El%Eta, dp)
       end if

       write(unit,'(3(e25.17,tr1))') Z * Ry2eV, real(W, dp) * Ry2eV
       
    end do

    call io_close( unit )

  end subroutine io_contour_nEq_ID

  subroutine neq_die(msg)
    character(len=*), intent(in) :: msg
    
    write(*,*) 'Killing... printing out so-far gathered information'
    call print_contour_neq_options('TS')
    call print_contour_neq_block('TS')
    call die(msg)

  end subroutine neq_die
    
end module m_ts_contour_neq
