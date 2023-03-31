!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2013, nickpapior@gmail.com
!

module m_ts_contour_eq

  use precision, only : dp

! Use the type associated with the contour
! Maybe they should be collected to this module.
! However, I like this partition.
  use m_ts_electype

  use m_ts_cctype
  use m_ts_chem_pot
  use m_ts_io_contour
  use m_ts_io_ctype
  use m_ts_aux

  implicit none

  ! equilibrium contour IO-segments
  integer, save, public :: N_Eq
  type(ts_c_io), pointer, save, public :: Eq_io(:) => null()
  type(ts_cw)  , pointer, save, public :: Eq_c(:) => null()

  ! The contours for the equilibrium density are attributed a fruitful discussion with
  ! Hans Skriver. Previously the routine names reflected his contribution.
  ! However, for programming clarity we have employed a different naming scheme.
  ! In order to retain the contributions it is encouraged to keep this sentiment for his
  ! contributions.

  ! For further attributions see the original paper by Brandbyge, et. al, 2002: DOI: 10.1103/PhysRevB.65.165401

  ! this is heavily linked with CONTOUR_NEQ from m_ts_contour_neq
  integer, parameter, public :: CONTOUR_EQ = 1

  public :: read_contour_eq_options
  public :: print_contour_eq_options
  public :: print_contour_eq_block
  public :: io_contour_eq
  public :: N_Eq_E
  public :: get_c_io_index
  public :: Eq_linear_index
  public :: Eq_E
  public :: c2energy
  public :: c2weight_eq
  public :: ID2idx

  interface ID2idx
     module procedure ID2idx_cw
     module procedure ID2idx_cidx
  end interface

  private

contains

  subroutine read_contour_eq_options(N_mu, mus, Volt)

    use units, only : Pi, Kelvin, eV
    use fdf

    integer, intent(in)        :: N_mu
    type(ts_mu), intent(inout) :: mus(N_mu)
    real(dp), intent(in)       :: Volt

    integer :: i, j, k, c_mu
    integer :: N
    character(len=C_N_NAME_LEN), allocatable :: tmp(:), nContours(:)
    character(len=C_N_NAME_LEN) :: tmp_one
    integer :: cur
    integer, allocatable :: idx(:)
    logical :: isStart, isTail
    real(dp) :: r1, r2

    call fdf_obsolete('TS.ComplexContour.NPoles')


    ! We only allow the user to either use the old input format, or the new
    ! per-electrode input
    do i = 1 , N_mu
       if ( Eq_segs(mus(i)) == 1 ) then
          mus(i)%Eq_seg(1) = 'cont-frac-'//trim(mus(i)%name)
       else
          write(mus(i)%Eq_seg(Eq_segs(mus(i))),'(a,i0)') 'pole',i
       end if
    end do

    ! Count the number of contour segments
    N = minval(Eq_segs(mus))
    if ( N == 0 ) then
       call die('You have not assigned any contours for the chemical &
            &potentials.')
    end if
    N = sum(Eq_segs(mus))

    ! collect all equilibrium names
    allocate(tmp(N))
    tmp(:) = ' '
    N = 0
    do i = 1 , N_mu
       do j = 1 , Eq_segs(mus(i)) - 1
          N = N + 1
          tmp(N) = mus(i)%Eq_seg(j)
       end do
    end do
    ! we need to have the pole contours as the last ones
    ! to be able to easily differ them from the read in values
    do i = 1 , N_mu
       j = Eq_segs(mus(i))
       N = N + 1
       tmp(N) = mus(i)%Eq_seg(j)
    end do

    ! find all unique equilibrium names
    N = 1
    uniq_names: do i = 2 , size(tmp)
       j = 0
       do 
          j = j + 1
          if ( i <= j ) exit
          if ( tmp(i) .eq. tmp(j) ) cycle uniq_names
       end do
       N = N + 1
    end do uniq_names

    allocate(nContours(N)) ! unique names
    nContours = ' '

    ! populate unique equilibrium names
    N = 0
    uniq_names_pop: do i = 1 , size(tmp)
       do j = 1 , N
          if ( tmp(i) .eq. nContours(j) ) cycle uniq_names_pop
       end do
       N = N + 1
       if ( N > size(nContours) ) call die('Unique contour setup went wrong, contact devs')
       nContours(N) = tmp(i)
    end do uniq_names_pop
    if ( N /= size(nContours) ) then
       call die('ERROR: We have not populated enough contours')
    end if
    deallocate(tmp)

    ! Allocate space for the equilibrium contours
    nullify(Eq_io,Eq_c)
    N_Eq = N
    allocate(Eq_io(N_Eq)) ! the equilibrium contour io container
    allocate(Eq_c(N_Eq))  ! the equilibrium contour container
    Eq_io(:)%type = 'eq'  ! set the equilibrium tag

    ! Attach all the contours
    do i = 1 , N
       Eq_c(i)%c_io => Eq_io(i)
    end do

    ! read in the equilibrium contours (the last will be the poles)
    do i = 1 , N - N_mu

       ! Save the name of the contour (to be able to find the
       ! associated chemical potential)
       Eq_io(i)%name = nContours(i)

       ! Get index of associated chemical potential
       c_mu = 0
       do j = 1 , N_mu
          if ( hasC(mus(j),Eq_io(i)) ) then
             c_mu = j
             exit
          end if
       end do
       if ( c_mu == 0 ) then
          print *,Eq_io(i)%name
          call die('Error in determining associated &
            &chemical potential')
       end if
       ! Remove the name, it is needed to be ' '
       Eq_io(i)%name = ' '

       tmp_one = nContours(i)

       ! read in the contour
       if ( tmp_one(1:1) == '*' .and. &
            .not. ts_exists_contour_block('TS','',tmp_one(2:)//' ') ) then

          if ( N_mu > 2 ) then
             ! here we can have a zero bias case where 
             ! N_mu == 1 :)
             call die('When using other than 2 electrodes and &
                  &chemical potentials you are forced to setup &
                  &the contour input yourself.')
          end if

          ! this is a fake-contour, read it in
          ! if it exists, else create it...
          ! *** NOTE this is hard coded against the chem_pot code

          Eq_io(i)%name = tmp_one
          if ( tmp_one(2:2) == 'c' ) then
             ! the circle contour
             Eq_io(i)%part = 'circle'
             Eq_io(i)%N = 25
             write(Eq_io(i)%cN,'(i0)') Eq_io(i)%N
             if ( tmp_one(4:4) == 'l' ) then ! left
                Eq_io(i)%ca = '-40. eV + V/2'
                Eq_io(i)%a = - 40._dp * eV + Volt * .5_dp
                Eq_io(i)%cb = '-10 kT + V/2'
                Eq_io(i)%b = -10._dp * mus(c_mu)%kT + Volt * .5_dp
             else ! must be right
                Eq_io(i)%ca = '-40. eV - V/2'
                Eq_io(i)%a = - 40._dp * eV - Volt * .5_dp
                Eq_io(i)%cb = '-10 kT - V/2'
                Eq_io(i)%b = -10._dp * mus(c_mu)%kT - Volt * .5_dp
             end if
             Eq_io(i)%method = 'g-legendre'
             call c_io_add_opt(Eq_io(i),'right','right')

          else
             Eq_io(i)%part = 'tail'
             Eq_io(i)%N = 10
             write(Eq_io(i)%cN,'(i0)') Eq_io(i)%N
             Eq_io(i)%ca = 'prev'
             if ( tmp_one(4:4) == 'l' ) then ! left
                Eq_io(i)%a = -10._dp * mus(c_mu)%kT + Volt * .5_dp
             else ! must be right
                Eq_io(i)%a = -10._dp * mus(c_mu)%kT - Volt * .5_dp
             end if
             Eq_io(i)%cb = 'inf'
             Eq_io(i)%b = huge(1._dp)
             Eq_io(i)%method = 'g-fermi'
          end if
       else
          call ts_read_contour_block('TS','',nContours(i),Eq_io(i), &
               mus(c_mu)%kT, Volt)
       end if
       
    end do

    ! We here create the "fake" pole contours
    j = 1
    do i = N - N_mu + 1, N
       if ( Eq_segs(mus(j)) == 1 ) then
          Eq_io(i)%part = 'cont-frac'
          Eq_io(i)%method = 'continued-fraction'
          ! Currently this is just a very high number
          Eq_io(i)%a = 1.e10_dp * eV ! the continued fraction infinity point
       else
          Eq_io(i)%part = 'pole'
          Eq_io(i)%method = 'residual'
          Eq_io(i)%a = mus(j)%mu
       end if
       ! assign name to the Eq_io
       Eq_io(i)%name = nContours(i)
       Eq_io(i)%N = mus(j)%N_poles
       Eq_io(i)%b = mus(j)%kT ! Save kT in the contour
       Eq_io(i)%d = mus(j)%mu ! save the chemical potential here
       j = j + 1
           
    end do

    ! fix the read-in constants
    do i = 1 , N_mu

       cur = get_c_io_index(mus(i)%Eq_seg(1))
       if ( cur == 0 ) then
          call die('A terrible error has occured, please inform the &
               &developers')
       end if
       ! Two cases can happen,
       !  1) A continued fraction method
       !  2) The regular circle contour
       if ( leqi(Eq_io(cur)%part,'cont-frac') ) then
          
          ! The segment _has_ to be alone
          if ( Eq_segs(mus(i)) > 1 ) then
             call die('Continued fraction contours &
                  &are individual and cannot be connected.')
          end if

       else

          ! first fix the contours (not the poles)
          k = Eq_segs(mus(i)) - 1
          
          allocate(idx(k))

          ! Create index-table
          do j = 1 , k
             idx(j) = get_c_io_index(mus(i)%Eq_seg(j))
          end do
          
          ! Fix contours
          call ts_fix_contours(N_Eq, Eq_io, k, idx)

          deallocate(idx)

          ! Check equilibrium settings
          isStart = leqi(Eq_io(cur)%part,'circle') .or. &
               leqi(Eq_io(cur)%part,'square') 
          isTail = leqi(Eq_io(cur)%part,'tail')

          ! we should not check the pole
          do j = 1 , k
             cur = get_c_io_index(mus(i)%Eq_seg(j))
             call consecutive_types(Eq_io(cur),isStart,isTail, &
                  mus(i)%mu, mus(i)%kT)
          end do
          
          ! check that the last contour actually lies on the RHS
          ! of the chemical potential, at least 10 kT across!
          if ( Eq_io(cur)%b < mus(i)%mu + 10._dp * mus(i)%kT ) then
             print *,'Contour name: ',trim(Eq_io(cur)%name)
             write(*,*) 'Energies are too close: ',Eq_io(cur)%b,mus(i)%mu + 10._dp * mus(i)%kT
             call eq_die('The last contour of the chemical potential: &
                  &'//trim(Name(mus(i)))//' lies too close to the &
                  &chemical potential. It must be at least 10 kT from mu.')
          end if
          
       end if
       
    end do
    
    ! Allocate sizes of the contours
    do i = 1 , N_Eq
       ! number of different weights for each energy-point
       k = count(hasC(mus,Eq_c(i)))
       if ( k == 0 ) then
          call eq_die('No chemical potentials has been assigned this contour: '// &
               trim(Eq_c(i)%c_io%name))
       end if

       ! the id's referring to the chemical potential in the
       ! mus array.
       allocate(Eq_c(i)%ID(k))
       k = 0
       do j = 1 , N_mu
          if ( hasC(mus(j),Eq_c(i)) ) then
             k = k + 1
             Eq_c(i)%ID(k) = j

             ! We need to make sure that the number
             ! of poles is the same if they overlap
             if ( k > 1 ) then
                ! We need to ensure that the energy of the poles coincide
                ! This will happen very rarely!!!
                r1 = Pi * mus(j)%kT * 2._dp * mus(j)%N_poles
                c_mu = Eq_c(i)%ID(k-1)
                r2 = Pi * mus(c_mu)%kT * 2._dp * mus(c_mu)%N_poles

                ! Less than 5 K in difference
                if ( abs(r1 - r2) > Kelvin ) then
                   call die('If two equilibrium contours &
                        &should overlap, they should have the same &
                        &energy where the line crosses the Fermi energy!')
                end if
             end if
          end if
       end do

       if ( Eq_c(i)%c_io%N < 1 ) then
          write(*,*) 'Contour: '//trim(Eq_c(i)%c_io%Name)//' has 0 points.'
          write(*,*) 'Please ensure at least 1 point in each contour...'
          call die('Contour number of points, invalid')
       end if

       ! Allocate the different weights and initialize
       allocate(Eq_c(i)%c(Eq_c(i)%c_io%N))
       Eq_c(i)%c = 0._dp
       allocate(Eq_c(i)%w(k,Eq_c(i)%c_io%N))
       Eq_c(i)%w = 0._dp
       
    end do
    
    do i = 1 , N_mu

       call setup_Eq_contour(mus(i))
       
    end do

    ! TODO correct empty cycles i.e. when two contours share an energy-point

  contains
    
    subroutine consecutive_types(c_io,isStart,isTail,mu,kT)
      type(ts_c_io), intent(in) :: c_io
      logical, intent(inout) :: isStart, isTail
      real(dp), intent(in) :: mu, kT

      if ( (leqi(c_io%part,'circle') .or. leqi(c_io%part,'square')) &
           .and. .not. isStart ) then
         call eq_die('The circle/square contour must be the first contour type &
              &as well as connected to other circle/square contours.')
      end if
      isStart = leqi(c_io%part,'circle') .or. leqi(c_io%part,'square')

      if ( isStart ) then
         ! We need to check whether the contour lies well below the
         ! chemical potential
         if ( c_io%b > mu - 3._dp * kT ) then
            call eq_die('The circle/square contour lies too close or across the &
                 &chemical potential. This is not allowed!')
         end if
      end if

      if ( .not. leqi(c_io%part,'tail') .and. isTail ) then
         call eq_die('The tail contour must be the last contour &
              &as well as connected to other tail contours.')
      end if
      isTail = leqi(c_io%part,'tail')

    end subroutine consecutive_types

  end subroutine read_contour_eq_options


  ! This routine assures that we have setup all the 
  ! equilibrium contours for the passed electrode
  subroutine setup_Eq_contour(mu)
    use fdf, only: leqi
    type(ts_mu), intent(in) :: mu

    ! Local variables
    integer :: i, j, l, idx
    real(dp) :: a, b, R, cR, lift
    integer :: sq_prev, sq_next

    if ( Eq_segs(mu) == 0 ) &
         call die('Chemical potential: '//trim(Name(mu))//' has &
         &no equilibrium contours.')

    ! retrieve circle bounds for the electrode
    call mu_circle_bounds(mu,a,b)

    ! Calculate the circle entries
    call calc_Circle(a,b,mu%N_poles,mu%kT,R,cR,lift)

    do i = 1 , Eq_segs(mu)

       ! ensures we progress the contour from start to end
       ! and not in "random" order
       idx = get_c_io_index(mu%Eq_seg(i))

       if ( leqi(Eq_c(idx)%c_io%part,'user') ) then

          call contour_file(Eq_c(idx),mu,lift)

        else if ( leqi(Eq_c(idx)%c_io%part,'circle') ) then

          call contour_Circle(Eq_c(idx),mu,R,cR)

       else if ( leqi(Eq_c(idx)%c_io%part,'square') ) then

          ! find the previous and following square
          sq_prev = idx
          sq_next = idx
          do j = 1 , Eq_segs(mu)
             l = get_c_io_index(mu%Eq_seg(j))
             if ( leqi(Eq_c(l)%c_io%part,'square') ) then
                if ( l == idx - 1 ) sq_prev = l
                if ( l == idx + 1 ) sq_next = l
             end if
          end do

          call contour_Square(Eq_c(idx),mu, lift, i==1, &
               Eq_c(sq_prev), Eq_c(sq_next))
          
       else if ( leqi(Eq_c(idx)%c_io%part,'line') ) then

          call contour_line(Eq_c(idx),mu,lift)

       else if ( leqi(Eq_c(idx)%c_io%part,'tail') ) then

          call contour_tail(Eq_c(idx),mu,lift)

       else if ( leqi(Eq_c(idx)%c_io%part,'pole') ) then

          ! the poles all have the same weight (Pi*kT*2)
          call contour_poles(Eq_c(idx),Eq_c(idx)%c_io%d,mu%kT)

       else if ( leqi(Eq_c(idx)%c_io%part,'cont-frac') ) then

          ! the poles all have the same weight (Pi*kT*2)
          call contour_continued_fraction(Eq_c(idx),mu,lift)

       else
          
          call die('Unrecognized contour type for the &
               &equilibrium part.')

       end if

    end do

  contains

    subroutine calc_Circle(a,b,N_poles,kT,R,cR,lift)
      use units, only : Pi
      real(dp), intent(in)  :: a,b ! the circle bounds
      integer, intent(in) :: N_poles ! number of poles the 'b' point is lifted
      real(dp), intent(in) :: kT
      real(dp), intent(out) :: R, cR, lift ! the radius, center(real)
      
      ! local variables
      real(dp) :: alpha

      ! The poles are @ \pi kT, 3\pi kT, 5\pi kT, ...
      ! Middle point are @ 2\pi kT, 4\pi kT, 6\pi kT
      lift = Pi * kT * 2._dp * N_poles
      ! this means that we place the line contour right in the middle of two poles!

      cR = b - a
      ! the angle between the axis and the line from the start
      ! of the circle to the end of the circle contour
      alpha = datan( lift / cR )

      ! the radius can be calculated using two triangles in the circle
      ! there is no need to use the cosine-relations
      R = 0.5_dp * cR / cos(alpha) ** 2

      ! the real-axis center
      cR = a + R

    end subroutine calc_Circle

    ! retrieve the bounds of the circle contour
    ! this allows to split the circle integral into as many parts as needed
    subroutine mu_circle_bounds(mu,a,b)
      type(ts_mu), intent(in) :: mu
      real(dp), intent(out) :: a,b
      if ( Eq_segs(mu) < 1 ) call die('Error Eq_seg CB')
      idx = get_c_io_index(mu%Eq_seg(1))
      a = Eq_c(idx)%c_io%a
      b = Eq_c(idx)%c_io%b
      do i = 1 , Eq_segs(mu)
         idx = get_c_io_index(mu%Eq_seg(i))
         if ( Eq_c(idx)%c_io%type == 'eq' .and. &
              leqi(Eq_c(idx)%c_io%part, 'circle') ) then
            b = Eq_c(idx)%c_io%b
         end if
      end do
    end subroutine Mu_circle_bounds

  end subroutine setup_Eq_contour

  subroutine ID2idx_cw(c,ID,idx)
    type(ts_cw), intent(in) :: c
    integer,     intent(in) :: ID
    integer,    intent(out) :: idx
    do idx = 1 , size(c%ID)
       if ( c%ID(idx) == ID ) return
    end do
    idx = -1
  end subroutine ID2idx_cw

  subroutine ID2idx_cidx(c,ID,idx)
    type(ts_c_idx), intent(in) :: c
    integer, intent(in) :: ID
    integer, intent(out) :: idx
    if ( c%idx(1) /= CONTOUR_EQ ) call die('Could not locate ID')
    call ID2idx_cw(Eq_c(c%idx(2)),ID,idx)
  end subroutine ID2idx_cidx

  pure subroutine c2weight_eq(c,idx,k,W,ZW)
    type(ts_c_idx), intent(in) :: c
    integer, intent(in) :: idx
    real(dp), intent(in)  :: k
    complex(dp), intent(out) :: W,ZW
    
    if ( c%idx(1) /= CONTOUR_EQ ) then
       W = 0._dp
       ZW = 0._dp
       return
    end if

    W  = k * Eq_c(c%idx(2))%w(idx,c%idx(3))
    ZW = W * Eq_c(c%idx(2))%c(c%idx(3))
    
  end subroutine c2weight_eq

  pure function c2energy(c) result(Z)
    type(ts_c_idx), intent(in) :: c
    complex(dp) :: Z
    
    Z = Eq_c(c%idx(2))%c(c%idx(3))
    
  end function c2energy

  subroutine contour_Circle(c,mu,R,cR)
    use fdf, only: leqi
    use m_integrate
    use m_gauss_quad
    type(ts_cw), intent(inout) :: c
    type(ts_mu), intent(in) :: mu
    real(dp), intent(in) :: R, cR

    ! local variables
    character(len=c_N) :: tmpC
    logical :: set_c
    complex(dp) :: ztmp
    integer :: i, idx
    real(dp) :: a,b, tmp
    real(dp), allocatable :: ce(:), cw(:)

    if ( .not. leqi(c%c_io%part,'circle') ) &
         call die('Contour is not a circle')
   
    ! notice the switch (the circle has increasing contour
    ! from right to left)
    call calc_angle(c%c_io%a,R,cR,a)
    call calc_angle(c%c_io%b,R,cR,b)
    if ( b > a ) then
       call die('Contour equilibrium circle went wrong')
    end if

    allocate(ce(c%c_io%N))
    allocate(cw(c%c_io%N))

    select case ( method(c%c_io) ) 
    case ( CC_G_LEGENDRE ) 

       if ( c_io_has_opt(c%c_io,'right') ) then
          
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call Gauss_Legendre_Rec(2*c%c_io%N,0,2._dp*a-b,b,ce,cw)

          do i = 1 ,c%c_io%N
             ce(i) = ce(i+c%c_io%N)
             cw(i) = cw(i+c%c_io%N)
          end do

       else if ( c_io_has_opt(c%c_io,'left') ) then
          
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call Gauss_Legendre_Rec(2*c%c_io%N,0,a,2._dp*b-a,ce,cw)

       else

          call Gauss_Legendre_Rec(c%c_io%N,0,a,b,ce,cw)

       end if

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

       if ( c_io_has_opt(c%c_io,'right') ) then
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call TanhSinh_Exact(c%c_io%N*2,ce,cw,2._dp*a-b,b, p=tmp)

          do i = 1 ,c%c_io%N
             ce(i) = ce(i+c%c_io%N)
             cw(i) = cw(i+c%c_io%N)
          end do

       else if ( c_io_has_opt(c%c_io,'left') ) then
          
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call TanhSinh_Exact(c%c_io%N*2,ce,cw,a,2._dp*b-a, p=tmp)

       else

          call TanhSinh_Exact(c%c_io%N,ce,cw,a,b, p=tmp)

       end if

    case ( CC_BOOLE_MIX )

       call Booles_Simpson_38_3_rule(c%c_io%N, ce, cw, a, b)

    case ( CC_SIMP_MIX )
       
       call Simpson_38_3_rule(c%c_io%N,ce,cw,a,b)

    case ( CC_MID )
       
       call Mid_Rule(c%c_io%N,ce,cw,a,b)

    case default
       call die('Unknown method for the circle integral, please correct')
    end select

    ! I know this is "bad" practice, however, zero is a well-defined quantity in FPU
    set_c = sum(abs(c%c(:))) == 0._dp

    ! get the index in the ID array (same index in w-array)
    call ID2idx(c,mu%ID,idx)

    do i = 1 , c%c_io%N
       ztmp = R * cdexp(cmplx(0._dp,ce(i), dp))

       if ( set_c ) then
          c%c(i) = cmplx(cR,0._dp, dp) + ztmp
       else
          if ( abs(c%c(i) - (cmplx(cR,0._dp, dp) + ztmp) ) > 1.e-10_dp ) then
             print *,c%c(i),cmplx(cR,0._dp,dp) + ztmp
             call die('contours does not match')
          end if
       end if

       ! Factor i, comes from Ed\theta=dE=iR e^{i\theta}
       c%w(idx,i) = cw(i) * nf((cmplx(cR,0._dp,dp)+ztmp-mu%mu)/mu%kT) &
            * cmplx(0._dp,1._dp,dp) * ztmp

    end do

    deallocate(ce,cw)
    
  contains

    subroutine calc_angle(a,R,cR,angle)
      use units, only : Pi
      real(dp), intent(in) :: a,R,cR
      real(dp), intent(out) :: angle
      real(dp) :: E

      E = abs(cR - a)

      ! the integration angles
      if ( abs(E - R) <  1.e-8_dp ) then
         ! we are at the radius of the circle
         angle = Pi
      else
         angle = asin(E/R)
         ! correct for left/right side of radius
         if ( a < cR ) then
            angle = Pi * 0.5_dp + angle
         else if ( a > cR ) then
            angle = Pi * 0.5_dp - angle
         end if
      end if

    end subroutine calc_angle

  end subroutine contour_Circle

  subroutine contour_Square(c,mu,lift, first, prev, next)
    use fdf, only: leqi
    use m_integrate
    use m_gauss_quad
    type(ts_cw), intent(inout) :: c
    type(ts_mu), intent(in) :: mu
    real(dp), intent(in) :: lift
    ! Whether this is the first square contour (starts from 0)
    logical, intent(in) :: first
    ! the previous and following square contour
    type(ts_cw), intent(inout) :: prev, next

    ! local variables
    character(len=c_N) :: tmpC
    logical :: set_c
    complex(dp) :: ztmp(2)
    ! Segments
    integer :: N_seg, Ni(3)
    real(dp) :: li(3)
    integer :: i, j, idx
    ! The previous, current and next imaginary part
    real(dp) :: im(3)
    real(dp) :: a, b
    real(dp), allocatable :: ce(:), cw(:)

    if ( .not. leqi(c%c_io%part,'square') ) &
         call die('Contour is not a square')

    if ( first .and. .not. (prev%c_io .eq. c%c_io) ) then
       call die('Error in setting up square contour.')
    end if

    ! We need to split the square quadrature up into different parts
    ! Specifically this is a requirement as the Gaussian quadratures
    ! does not have an equi spaced distance between abscissas.
    ! Thus when climbing/descending on the complex contour we find different
    ! prefactors which is sub-optimal.

    ! Thus we figure out the number of segments and divide it into those
    ! partitions

    ! Create parametrization lengths
    if ( first ) then
       im(1) = 0._dp
    else 
       im(1) = read_eta(prev, lift)
    end if
    im(2) = read_eta(c, lift)
    if ( next%c_io .eq. c%c_io ) then
       im(3) = lift
       N_seg = 3
    else
       ! this is because the following square
       ! contour steps up!
       im(3) = im(2)
       N_seg = 2
    end if

    if ( im(2) < 0._dp ) then
       print *,'final imaginary part: ',im(2)
       call die("Your square contour crosses the real axis. &
            &This is not allowed")
    end if

    ! Allocate total weights
    allocate(ce(c%c_io%N))
    allocate(cw(c%c_io%N))
    ! I know this is "bad" practice, however, zero is a well-defined quantity in FPU
    set_c = sum(abs(c%c(:))) == 0._dp
    call ID2idx(c,mu%ID,idx)

    ! the real axis start/end
    a = c%c_io%a
    b = c%c_io%b

    ! Now 'N_seg' contains number of segments
    !  dL1 = abs(im(2) - im(1))
    li(1) = im(2) - im(1)
    !  dL2 = b - a
    li(2) = b - a
    !  dL3 = abs(im(3) - im(2)) ! which may be zero...
    li(3) = im(3) - im(2)

    ! Ensure enough points for remaining lines
    Ni = 2
    call n_distribute(c%c_io%N, N_seg, li, Ni)
    if ( any(Ni < 2) ) then
       call die('Could not divide contour points appropriately. &
            &Choose a different contour method or increase contour points.')
    end if

    ! Get first segment
    call get_abscissas(Ni(1),0._dp,li(1), ce(1), cw(1))
    ! translate to correct contours
    do i = 1 , Ni(1)
       ztmp(1) = cmplx( a, im(1) + ce(i), dp)
       ztmp(2) = cmplx( 0._dp , cw(i), dp)
       call set_abscissas(c%c(i),c%w(idx,i), ztmp)
    end do

    ! Get first segment
    call get_abscissas(Ni(2),0._dp,li(2), ce(1), cw(1))
    ! translate to correct contours
    j = Ni(1)
    do i = 1 , Ni(2)
       ! Convert to parameterisation
       ztmp(1) = cmplx( a + ce(i) , im(2), dp)
       ztmp(2) = cmplx( cw(i) , 0._dp, dp)
       call set_abscissas(c%c(j+i),c%w(idx,j+i), ztmp)
    end do
    
    if ( N_seg == 3 ) then
       ! Get first segment
       call get_abscissas(Ni(3),0._dp,li(3), ce(1), cw(1))
       ! translate to correct contours
       j = Ni(1) + Ni(2)
       do i = 1 , Ni(3)
          ztmp(1) = cmplx( b, im(2) + ce(i), dp)
          ztmp(2) = cmplx( 0._dp , cw(i), dp)
          call set_abscissas(c%c(j+i),c%w(idx,j+i), ztmp)
       end do
    end if

    deallocate(ce,cw)

  contains

    subroutine get_abscissas(N,a,b,ce,cw)
      integer, intent(in) :: N
      real(dp), intent(in) :: a, b
      real(dp), intent(inout) :: ce(N), cw(N)
      real(dp), allocatable :: tce(:), tcw(:)
      real(dp) :: tmp

      select case ( method(c%c_io) ) 
      case ( CC_G_LEGENDRE ) 
         
         if ( c_io_has_opt(c%c_io,'right') ) then
            
            allocate(tce(N*2),tcw(N*2))
            
            call Gauss_Legendre_Rec(2*N,0,2._dp*a-b,b,tce,tcw)
            
            do i = 1 ,N
               ce(i) = tce(i+N)
               cw(i) = tcw(i+N)
            end do

            deallocate(tce,tcw)
            
         else if ( c_io_has_opt(c%c_io,'left') ) then
            
            allocate(tce(N*2),tcw(N*2))
            
            call Gauss_Legendre_Rec(2*N,0,a,2._dp*b-a,tce,tcw)

            do i = 1 ,N
               ce(i) = tce(i)
               cw(i) = tcw(i)
            end do

            deallocate(tce,tcw)
            
         else
            
            call Gauss_Legendre_Rec(N,0,a,b,ce,cw)

         end if
         
      case ( CC_TANH_SINH )

         ! we should also gain an option for this
         if ( c_io_has_opt(c%c_io,'precision') ) then
            tmpC = c_io_get_opt(c%c_io,'precision')
            read(tmpC,'(g20.10)') tmp
         else
            tmp = 2.e-2_dp * abs(b-a) / real(N,dp)
            write(tmpC,'(g20.10)') tmp
            call c_io_add_opt(c%c_io,'precision',tmpC)
         end if

         if ( c_io_has_opt(c%c_io,'right') ) then

            allocate(tce(N*2),tcw(N*2))

            call TanhSinh_Exact(N*2,tce,tcw,2._dp*a-b,b, p=tmp)
            
            do i = 1 ,N
               ce(i) = tce(i+N)
               cw(i) = tcw(i+N)
            end do

         else if ( c_io_has_opt(c%c_io,'left') ) then

            allocate(tce(N*2),tcw(N*2))

            call TanhSinh_Exact(N*2,tce,tcw,a,2._dp*b-a, p=tmp)

            do i = 1 ,N
               ce(i) = tce(i)
               cw(i) = tcw(i)
            end do

            deallocate(tce,tcw)

         else

            call TanhSinh_Exact(N,ce,cw,a,b, p=tmp)

         end if

      case ( CC_BOOLE_MIX )
       
         call Booles_Simpson_38_3_rule(N, ce, cw, a, b)

      case ( CC_SIMP_MIX )

         call Simpson_38_3_rule(N,ce,cw,a,b)

      case ( CC_MID )

         call Mid_Rule(N,ce,cw,a,b)

      case default
         call die('Unknown method for the square integral, please correct')
      end select

    end subroutine get_abscissas

    subroutine set_abscissas(e,w,ztmp)
      complex(dp), intent(inout) :: e, w
      complex(dp), intent(inout) :: ztmp(2)
      real(dp) :: tmp
      
      if ( set_c ) then
         e = ztmp(1)
      else
         if ( abs(e - ztmp(1) ) > 1.e-10_dp ) then
            print *,e,ztmp(1)
            call die('contours does not match')
         end if
      end if
      
      ! the parameterisation does not alter the fermi function
      ! importantly we define the band-bottom to have a vanishing
      ! fermi function quantity.
      tmp = (real(ztmp(1),dp) - mu%mu) / mu%kT
      ztmp(1) = (ztmp(1) - mu%mu) / mu%kT
      w = nf(ztmp(1)) * ztmp(2)
      
    end subroutine set_abscissas
    
    function read_eta(c, lift) result(eta)
      use parse, only: parsed_line, digest, destroy
      type(ts_cw), intent(in) :: c
      real(dp), intent(in) :: lift
      real(dp) :: eta

      type(parsed_line), pointer :: pline
      character(len=c_N) :: tmpC

      if ( c_io_has_opt(c%c_io,'eta-add') ) then
         tmpC = c_io_get_opt(c%c_io,'eta-add')
         pline => digest(tmpC)
         eta = lift + ts_c_bphysical(pline,0,'Ry')
         call destroy(pline)
      else if ( c_io_has_opt(c%c_io,'eta') ) then
         tmpC = c_io_get_opt(c%c_io,'eta')
         pline => digest(tmpC)
         eta = ts_c_bphysical(pline,0,'Ry')
         call destroy(pline)
      else
         ! We default to add 1.5 Ry to the square
         ! to keep it well of the axis
         eta = lift + 1.5_dp
      end if

    end function read_eta
    
  end subroutine contour_Square
  
  subroutine contour_line(c,mu,Eta)
    use fdf, only: leqi
    use m_integrate
    use m_gauss_quad

    type(ts_cw), intent(inout) :: c
    type(ts_mu), intent(in) :: mu
    ! The lifting into the complex plane
    real(dp), intent(in) :: Eta

    ! local variables
    character(len=c_N) :: tmpC
    logical :: set_c
    integer :: i, idx
    real(dp) :: a,b, tmp
    real(dp), allocatable :: ce(:), cw(:)

    if ( .not. leqi(c%c_io%part,'line') ) &
         call die('Contour is not a line')

    ! get bounds
    a = c%c_io%a
    b = c%c_io%b
    
    allocate(ce(c%c_io%N))
    allocate(cw(c%c_io%N))

    select case ( method(c%c_io) )
    case ( CC_MID )
       
       call Mid_Rule(c%c_io%N,ce,cw,a,b)
       
    case ( CC_SIMP_MIX )
       
       call Simpson_38_3_rule(c%c_io%N,ce,cw,a,b)
       
    case ( CC_BOOLE_MIX )
       
       call Booles_Simpson_38_3_rule(c%c_io%N,ce,cw,a,b)       
       
    case ( CC_G_LEGENDRE ) 

       if ( c_io_has_opt(c%c_io,'right') ) then
          
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call Gauss_Legendre_Rec(c%c_io%N*2,0,2._dp*a-b,b,ce,cw)

          do i = 1 ,c%c_io%N
             ce(i) = ce(i+c%c_io%N)
             cw(i) = cw(i+c%c_io%N)
          end do

       else if ( c_io_has_opt(c%c_io,'left') ) then
          
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call Gauss_Legendre_Rec(c%c_io%N*2,0,a,2._dp*b-a,ce,cw)

       else

          call Gauss_Legendre_Rec(c%c_io%N,0,a,b,ce,cw)

       end if
       
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

       if ( c_io_has_opt(c%c_io,'right') ) then
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call TanhSinh_Exact(c%c_io%N*2,ce,cw,2._dp*a-b,b, p=tmp)

          do i = 1 ,c%c_io%N
             ce(i) = ce(i+c%c_io%N)
             cw(i) = cw(i+c%c_io%N)
          end do

       else if ( c_io_has_opt(c%c_io,'left') ) then
          
          deallocate(ce,cw)
          allocate(ce(c%c_io%N*2))
          allocate(cw(c%c_io%N*2))

          call TanhSinh_Exact(c%c_io%N*2,ce,cw,a,2._dp*b-a, p=tmp)

       else

          call TanhSinh_Exact(c%c_io%N,ce,cw,a,b, p=tmp)

       end if

    case ( CC_USER )

      ! Read the file information
      call contour_file(c,mu,Eta)

    case default
       write(*,*) 'Method for contour ',trim(c%c_io%name), &
            ' could not be deciphered: ', c%c_io%method
       call die('Could not determine the line-integral')
    end select
    
    ! get the index in the ID array (same index in w-array)
    call ID2idx(c,mu%ID,idx)
    
    if ( method(c%c_io) == CC_USER ) then

      do i = 1 , c%c_io%N
        c%w(idx,i) = c%w(idx,i) * nf((real(c%c(i), dp) - mu%mu) / mu%kT)
      end do
      
    else
      
      ! I know this is "bad" practice, however, zero is a well-defined quantity.
      set_c = sum(abs(c%c(:))) == 0._dp

      do i = 1 , c%c_io%N
        if ( set_c ) then
          c%c(i) = cmplx(ce(i),Eta, dp)
        else
          if ( abs(c%c(i) - cmplx(ce(i),Eta, dp)) > 1.e-10_dp ) then
            call die('contour_line: Error on contour match')
          end if
        end if

        c%w(idx,i) = cw(i) * nf((ce(i) - mu%mu) / mu%kT)

      end do
      
    end if
    
    deallocate(ce,cw)
    
  end subroutine contour_line

  subroutine contour_tail(c,mu,Eta)
    use fdf, only: leqi
    use m_gauss_fermi_inf
    use m_gauss_fermi_30
    use m_gauss_fermi_28
    use m_gauss_fermi_26
    use m_gauss_fermi_24
    use m_gauss_fermi_22
    use m_gauss_fermi_20
    use m_gauss_fermi_19
    use m_gauss_fermi_18
    use m_gauss_fermi_17

    type(ts_cw), intent(inout) :: c
    type(ts_mu), intent(in) :: mu
    ! This describes the lifting of the tail integral into the complex plane
    real(dp), intent(in) :: Eta

    ! local variables
    integer :: idx, offset, infinity
    real(dp) :: a,b
    real(dp), allocatable :: ce(:), cw(:)

    if ( .not. leqi(c%c_io%part,'tail') ) &
         call die('Contour is not a tail contour')

    ! get bounds
    a = c%c_io%a
    b = c%c_io%b

    allocate(ce(c%c_io%N))
    allocate(cw(c%c_io%N))
    
    select case ( method(c%c_io) )
    case ( CC_G_NF_MIN:CC_G_NF_MAX )

       offset = nint((c%c_io%a-mu%mu)/mu%kT)

       if ( abs(offset * mu%kT - (c%c_io%a-mu%mu)) > 1.e-7_dp ) then
          call die('The integer value of the kT offset for the &
               &Gauss-Fermi integral is not valid, please check input')
       end if
       if ( c%c_io%b < 0._dp ) then
          if ( c%c_io%b < -30.5_dp * mu%kT ) then
             infinity = huge(1)
          else
             infinity = nint(abs(c%c_io%b-mu%mu)/mu%kT)
          end if
       else
          if ( c%c_io%b > 30.5_dp * mu%kT ) then
             infinity = huge(1)
          else
             infinity = nint(abs(c%c_io%b-mu%mu)/mu%kT)
          end if
       end if
       if ( infinity > 30 ) infinity = huge(1)
       
       ! calculate the offset from the energy chemical potential tail
       select case ( infinity )
       case ( huge(1) )
          call GaussFermi_inf(offset,c%c_io%N,ce,cw)
       case ( 30 )
          call GaussFermi_30(offset,c%c_io%N,ce,cw)
       case ( 28 )
          call GaussFermi_28(offset,c%c_io%N,ce,cw)
       case ( 26 ) 
          call GaussFermi_26(offset,c%c_io%N,ce,cw)
       case ( 24 ) 
          call GaussFermi_24(offset,c%c_io%N,ce,cw)
       case ( 22 )
          call GaussFermi_22(offset,c%c_io%N,ce,cw)
       case ( 20 )
          call GaussFermi_20(offset,c%c_io%N,ce,cw)
       case ( 19 )
          call GaussFermi_19(offset,c%c_io%N,ce,cw)
       case ( 18 )
          call GaussFermi_18(offset,c%c_io%N,ce,cw)
       case ( 17 )
          call GaussFermi_17(offset,c%c_io%N,ce,cw)
       case default
          call die('Unknown tail integral ending')
       end select

       ! we might as well correct the method

       ! the Gauss-Fermi quadrature is wrt. E'->E/kT
       ce = ce * mu%kT + mu%mu
       cw = cw * mu%kT

       call ID2idx(c,mu%ID,idx)

       ! move over the weights and the contour values
       c%c(:) = cmplx(ce,Eta, dp)
       c%w(idx,:) = cw

    case default

      ! we revert so that we can actually use the line-integral
      ! The tail and line are equivalent in the sense that the
      ! fermi functions are applied to the weights
       c%c_io%part = 'line'

       call contour_line(c,mu,Eta)

    end select

    deallocate(ce,cw)

  end subroutine contour_tail


  ! The residuals of the fermi-function at a real-energy
  subroutine contour_poles(c, E, kT)
    use units, only : Pi

! ***********************
! * OUTPUT variables    *
! ***********************
    type(ts_cw), intent(inout) :: c

! ***********************
! * INPUT variables     *
! ***********************
    real(dp), intent(in) :: E    ! at energy
    real(dp), intent(in) :: kT   ! temperature in Ry

! ***********************
! * LOCAL variables     *
! ***********************
    integer :: i

    ! all pole-weights have the same weight (negative due to contour method)
    c%w(:,:) =  - cmplx(0._dp, Pi * kT * 2._dp, dp)
    ! Residuals
    do i = 1 , c%c_io%N
       c%c(i) = cmplx(E , Pi * kT * (2._dp*(i-1)+1._dp), dp)
    end do

  end subroutine contour_poles


  subroutine contour_continued_fraction(c,mu,Eta)
    use fdf, only: leqi
    use units, only: Pi

    type(ts_cw), intent(inout) :: c
    type(ts_mu), intent(in) :: mu
    ! The lifting into the complex plane
    real(dp), intent(in) :: Eta

    ! local variables
    logical :: set_c
    integer :: i, idx
    complex(dp) :: cc
    real(dp), allocatable :: ce(:), cw(:)

    if ( .not. leqi(c%c_io%part,'cont-frac') ) &
         call die('Contour is not a continued fraction')

    allocate(ce(c%c_io%N-1))
    allocate(cw(c%c_io%N-1))

    select case ( method(c%c_io) )
    case ( CC_CONTINUED_FRAC )
       
       call Ozaki_residue(c%c_io%N-1,ce,cw)

    case default
       write(*,*) 'Method for contour ',trim(c%c_io%name), &
            ' could not be deciphered: ', c%c_io%method
       call die('Could not determine the pole-integral')
    end select

    set_c = sum(abs(c%c(:))) == 0._dp

    ! get the index in the ID array (same index in w-array)
    call ID2idx(c,mu%ID,idx)

    do i = 1 , c%c_io%N - 1
       
       ! Calculate current contour point
       cc = cmplx(mu%mu, ce(i) * mu%kT, dp)
       
       if ( set_c ) then
          c%c(i) = cc
       else
          if ( abs(c%c(i) - cc) > 1.e-10_dp ) then
             call die('contour_cont_frac: Error on contour match')
          end if
       end if

       ! Extra minus in implementation and Im[]
       ! We also divide the weight by Pi in the loop (and it should
       ! not exist in the continued fraction scheme)
       c%w(idx,i) = cmplx( 0._dp , 2._dp * cw(i) * mu%kT * Pi, dp)

    end do

    ! The zero'th moment lies infinitely far and is from -inf -- inf
    cc = cmplx( mu%mu , c%c_io%a, dp)
    
    if ( set_c ) then
       ! The last pole is set
       c%c(c%c_io%N) = cc
    else
       if ( abs(c%c(c%c_io%N) - cc) > 1.e-10_dp ) then
          call die('contour_cont_frac: Error on contour match')
       end if
    end if

    ! The zeroth moment (extra minus in implementation and Im[])
    ! w = iR, but from -\Im we get w = R
    ! And remove the loop division by Pi
    c%w(idx,c%c_io%N) = cmplx( 0.5_dp * c%c_io%a * Pi, 0._dp, dp)

    deallocate(ce,cw)

  contains

    subroutine Ozaki_residue(N,c,w)
      ! The number of poles
      integer, intent(in) :: N
      real(dp), intent(out) :: c(N), w(N)

      ! Diagonalization matrices
      real(dp), allocatable :: A(:,:), B(:,:), we(:), work(:)
      integer :: i, j, nn

      real(dp), external :: ddot

      ! The last weight of the residue problem is the R->infinity
      ! point. Hence we remove one point from the poles.
      nn = 2*N
      allocate(A(nn,nn), B(nn,nn))
      allocate(we(nn),work(3*nn))

      A = 0._dp
      B = 0._dp
      do i = 1 , nn - 1
         B(i,i) = 2._dp * i - 1._dp
         A(i,i+1) = -0.5_dp
         A(i+1,i) = -0.5_dp
      end do
      B(nn,nn) = 2._dp * nn - 1._dp
      
      ! Matrices have been initialized, diagonalize
      ! to find residues from eigenvalues
      call dsygv(1,'V','U',nn,A,nn,B,nn,we,work,3*nn,i)
      if ( i /= 0 ) then
         write(*,'(a)')'error in Ozaki diagonalization of weights.'
         call die('Error in diagonalization of weights')
      end if


      ! Transfer weights and eigenvalues
      j = nn
      do i = 1 , N

         ! Eigenvalue
         c(i) = 1._dp / we(j)
         
         ! Residual
         w(i) = - (A(1,j)*c(i))**2 * 0.25_dp

         j = j - 1
         
      end do

      ! Checked with Ozaki and got the same
      !print '(1000(/,2(tr2,f24.12)))',( (/c(i),w(i)/) ,i=1,N)

      deallocate(A,B,we,work)

    end subroutine Ozaki_residue
    
  end subroutine contour_continued_fraction

  ! This routine will read the contour points from a file
  subroutine contour_file(c,mu,Eta)
    use m_io, only: io_assign, io_close
    use fdf, only: fdf_convfac

    type(ts_cw), intent(inout) :: c
    type(ts_mu), intent(in) :: mu
    ! The lifting into the complex plane
    real(dp), intent(in) :: Eta

    integer :: iu, iostat, ne, idx
    logical :: exist
    complex(dp) :: E , W
    real(dp) :: rE, iE, rW, iW, conv
    character(len=512) :: file, line
    character(len=16) :: unit
    
    ! The contour type contains the file name in:
    !  c%c_io%cN (weirdly enough)
    file = c%c_io%cN
    call ID2idx(c,mu%ID,idx)

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
        c%w(idx,ne) = cmplx(rW,iW, dp) * conv
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
        c%c(ne) = cmplx(rE,iE,dp) * conv
        c%w(idx,ne) = cmplx(rW,iW,dp) * conv
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
        c%c(ne) = cmplx(rE * conv,iE,dp)
        c%w(idx,ne) = cmplx(rW,iW,dp) * conv
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
  
  function Eq_E(id,step) result(c)
    integer, intent(in) :: id
    integer, intent(in), optional :: step
    type(ts_c_idx) :: c ! the configuration of the energy-segment
    integer :: lstep, i, PN
    lstep = 1
    if ( present(step) ) lstep = step
    PN = N_Eq_E()
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
  end function Eq_E

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
    do j = 1 , N_Eq ! number of contours
       if ( iE + Eq_c(j)%c_io%N < id ) then
          iE = iE + Eq_c(j)%c_io%N
          cycle
       end if
       i = id - iE
       if ( i <= Eq_c(j)%c_io%N ) then
          c%exist = .true.
          c%e     = Eq_c(j)%c(i)
          c%idx(1) = CONTOUR_EQ ! designates the equilibrium contours
          c%idx(2) = j ! designates the index of the equilibrium contour
          c%idx(3) = i ! is the index of the equilibrium contour
          return
       end if
    end do

  end function get_c

  function N_Eq_E() result(N)
    integer :: N, i
    N = 0
    do i = 1 , N_Eq
       N = N + Eq_c(i)%c_io%N
    end do
  end function N_Eq_E

  !< Calculate the linear (total) index of a given index+local contour index
  function Eq_linear_index(idx, ie) result(lidx)
    !< The index of the equilibrium contour
    integer, intent(in) :: idx
    !< The index of the energy on the `idx` contour
    integer, intent(in) :: ie
    !< The global index of the `idx,ie` contour point
    integer :: lidx

    integer :: i

    lidx = ie
    do i = 1, idx - 1
      lidx = lidx + Eq_c(i)%c_io%N
    end do

  end function Eq_linear_index
    

  subroutine print_contour_eq_block(prefix)
    use fdf, only: leqi
    character(len=*), intent(in) :: prefix
    integer :: i

    do i = 1 , N_Eq
       if ( leqi(Eq_io(i)%part,'pole') ) cycle
       if ( leqi(Eq_io(i)%part,'cont-frac') ) cycle
       call ts_print_contour_block(trim(prefix)//'.Contour.',Eq_io(i))
    end do
  end subroutine print_contour_eq_block

  subroutine print_contour_eq_options(prefix)

    use fdf, only: leqi
    use parallel, only : IONode
    use units, only : Pi

    use m_ts_io_contour

    character(len=*), intent(in) :: prefix
    character(len=200) :: chars
    real(dp) :: tmp
    integer :: i, N
    type(ts_c_opt_ll), pointer :: opt

    if ( .not. IONode ) return
    
    write(*,opt_n) ' ----------------- Contour ----------------- '

    N = 0
    write(*,opt_n) '           >> Residual contour << '
    do i = 1 , N_eq
       chars = trim(eq_io(i)%part)
       if ( .not. (leqi(chars,'pole') .or. &
            leqi(chars,'cont-frac')) ) cycle
       N = N + 1
       if ( leqi(chars,'cont-frac') ) then
          call write_e('Continued fraction chemical potential',eq_io(i)%d)
       else
          call write_e('Pole chemical potential',eq_io(i)%d)
       end if
       call write_e('   Chemical potential temperature',eq_io(i)%b, &
            unit = 'K')
       write(*,opt_int) '   Number of poles',eq_io(i)%N
       if ( .not. leqi(chars,'cont-frac') ) then
          ! Calculate energy of middle pole-point
          tmp = Pi * eq_io(i)%b * 2._dp * eq_io(i)%N
          call write_e('   Top energy point',tmp)
       end if
    end do

    if ( N < N_Eq ) &
         write(*,opt_n) '         >> Equilibrium contour << '
    do i = 1 , N_eq
       if ( leqi(eq_io(i)%part,'pole') ) cycle
       if ( leqi(eq_io(i)%part,'cont-frac') ) cycle
       chars = '  '//trim(eq_io(i)%part)
       ! the starred contours are "fakes"
       if ( eq_io(i)%name(1:1) == '*' ) then
          write(*,opt_c) 'Contour name',trim(prefix)//'.Contour.'//trim(eq_io(i)%name(2:))
       else
          write(*,opt_c) 'Contour name',trim(prefix)//'.Contour.'//trim(eq_io(i)%name)
       end if

       call write_e(trim(chars)//' contour E_min',eq_io(i)%a)
       call write_e(trim(chars)//' contour E_max',eq_io(i)%b)
       write(*,opt_int) trim(chars)//' contour points',eq_io(i)%N
       write(*,opt_c) trim(chars)//' contour method', &
            trim(longmethod2str(eq_io(i)))
       opt => eq_io(i)%opt
       do while ( associated(opt) )
          if ( len_trim(opt%val) > 0 ) then
             write(*,opt_cc) '   Option for contour method',trim(opt%opt),trim(opt%val)
          else
             write(*,opt_c)  '   Option for contour method',trim(opt%opt)
          end if
          opt => opt%next
       end do
    end do
    
  end subroutine print_contour_eq_options

  function get_c_io_index(Name) result(idx)
    character(len=*), intent(in) :: Name
    integer :: idx
    do idx = 1 , N_Eq
       if ( trim(name) == trim(Eq_io(idx)%name) ) return
    end do
    idx = 0
  end function get_c_io_index

  subroutine io_contour_eq(mus,slabel,suffix)
    type(ts_mu), intent(in) :: mus(:)
    character(len=*), intent(in) :: slabel
    character(len=*), intent(in), optional :: suffix

    integer :: i

    do i = 1 , size(mus)
       call io_contour_eq_mu(mus(i),slabel,suffix)
    end do

  end subroutine io_contour_eq

  subroutine io_contour_eq_mu(mu,slabel,suffix)
    use parallel, only : IONode
    use fdf, only : leqi
    use units, only: eV, Kelvin
    type(ts_mu), intent(in) :: mu
    character(len=*), intent(in) :: slabel
    character(len=*), intent(in), optional :: suffix
    
! *********************
! * LOCAL variables   *
! *********************
    character(len=256) :: fname
    integer :: i, unit, idx
    type(ts_c_idx) :: cidx
    
    if ( .not. IONode ) return

    if ( present(suffix) ) then
       fname = trim(slabel)//trim(suffix)
    else
       fname = trim(slabel)//'.TSCCEQ-'//trim(Name(mu))
    end if

    call io_assign( unit )
    open( unit, file=fname, status='unknown' )
    write(unit,'(a)') '# Contour path for the equilibrium contour segment.'
    write(unit,'(a)') '# This segment belongs to the chemical potential: '//trim(Name(mu))
    write(unit,'(a)') '# Chemical potential:'
    if ( mu%mu < 0._dp ) then
      write(unit, '(a,g10.4,a)')'# - ', -mu%mu / eV, ' eV'
    else
      write(unit, '(a,g10.4,a)')'# + ', mu%mu / eV, ' eV'
    end if
    write(unit,'(a)') '# Electronic temperature:'
    write(unit, '(a,g10.4,a)')'# ', mu%kT / Kelvin, ' K'
    write(unit,'(a,a24,3(tr1,a25))') '#','Re(c) [eV]','Im(c) [eV]','Re(w) [eV]','Im(w) [eV]'

    cidx%idx(1) = CONTOUR_EQ
    do i = 1 , Eq_segs(mu)

       cidx%idx(2) = get_c_io_index(mu%Eq_seg(i))
       if ( cidx%idx(2) < 1 ) call die('io_contour_eq_mu: Error in code, C-ID')
       call ID2idx(Eq_c(cidx%idx(2)),mu%ID,idx)
       if ( idx < 1 ) call die('io_contour_eq_mu: Error in code')

       !print *,trim(Name(mu)),mu%Eq_seg(i)

       ! we have now retrieved the electrode index in the contour part
       ! write it out
       call io_contour_c(unit,cidx,idx)

    end do

    call io_close( unit )

  end subroutine io_contour_eq_mu

! Write out the contour to a contour file
  subroutine io_contour_c(unit,cidx,idx)
    use parallel, only : IONode
    use units, only : Pi, eV
    use fdf, only: leqi
    integer, intent(in) :: unit
    type(ts_c_idx), intent(inout) :: cidx
    integer, intent(in) :: idx

! *********************
! * LOCAL variables   *
! *********************
    integer :: i
    logical :: is_cont_frac
    type(ts_cw), pointer :: c
    complex(dp) :: W, ZW

    if ( .not. IONode ) return
    c => Eq_c(cidx%idx(2))

    is_cont_frac = leqi(c%c_io%part,'cont-frac')
    
    do i = 1 , size(c%c)
       cidx%e = c%c(i)
       cidx%idx(3) = i
       call c2weight_eq(cidx,idx,1._dp,W,ZW)
       if ( is_cont_frac ) then
          write(unit,'(4(e25.17,tr1))') c%c(i)/eV, W/(eV*Pi)
       else
          write(unit,'(4(e25.17,tr1))') c%c(i)/eV, W/eV
       end if
    end do
    
  end subroutine io_contour_c

  subroutine eq_die(msg)
    character(len=*), intent(in) :: msg
    write(*,*) 'Killing... printing out so-far gathered information'
    call print_contour_eq_options('TS')
    call print_contour_eq_block('TS')
    call die(msg)
  end subroutine eq_die

end module m_ts_contour_eq
