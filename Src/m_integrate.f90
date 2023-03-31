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


! We neede a generalized module for integration with specific methods

! We implement the following closed form rules
!  1) rectangle (mid-point rule)
!  2) Simpons's rule
!  3) Simpson's 3/8 rule
!  4) Boole's rule

! What is nice about the methods is that we can call them recursively to
! easily construct the composite rules.
! When using composite rules it may occur that the number of points does not
! coincide with the requested number of points.
! In that case we do a tree expansion of the rules, i.e.
!  Any high-moment rules tries the lower rule for the end-points

module m_integrate

  use precision, only : dp

  implicit none

  real(dp), parameter :: Boole_Simp38_split = 3._dp * .2_dp

  public

contains

  subroutine Mid_Rule(N,x,w,a,b)
    integer, intent(in) :: N
    real(dp), intent(out) :: x(N), w(N)
    real(dp), intent(in) :: a, b

    real(dp) :: h
    integer :: i

    h = (b-a)/real(N,dp)

    do i = 1 , N
       x(i) = a + h * (real(i,dp) -.5_dp)
       w(i) = h
    end do
    
  end subroutine Mid_Rule

  ! We use a combination of the Booles, Simpson 3/8 and Simpson 3-point rules
  ! to generate arbitrary N rules
  subroutine Booles_Simpson_38_3_rule(N,x,w,a,b)
    integer,  intent(in)  :: N
    real(dp), intent(out) :: x(N), w(N)
    real(dp), intent(in)  :: a, b

    ! S_38_E = 6480._dp ** .2_dp
    real(dp), parameter :: S_38_E = 5.785155024015764_dp
    ! S_3_E  = 2880._dp ** .2_dp
    real(dp), parameter :: S_3_E  = 4.919018971698727_dp

    real(dp) :: D, mD, curA, curB
    integer :: i, j, mBo, m38

    x(:) = 0._dp
    w(:) = 0._dp

    if ( N <= 4 ) then
       call Simpson_38_3_rule(N,x,w,a,b)
       return
    else if ( N == 5 ) then
       call Booles_5(x,w,a,b)
    else
       select case ( N )
       case ( 6, 7, 10 )
          call Simpson_38_3_rule(N,x,w,a,b)
          return
       end select
    end if
    
    ! We need to calculate the best separation.
    ! This is not trivially done as the error term depends
    ! on different differentials of the function we wish to integrate
    ! Hence we need to do empirical integration.
    ! In this case we choose a fixed number.
    ! We choose the number specified by Boole_Simp38_split
    ! (which currently is 3/5)

    i = 1 
    curA = a
    
    if ( mod(N,4) == 1 ) then

       ! Everything can be in the Boole method
       mBo = (N-1) / 4
       m38 = 0

       ! We do it in here to not confuse the logic of the routine
       ! First we generate the step size of each block
       mD = (b-a)/real(mBo,dp)

       curB = a + mD

    else if ( mod(N,4) == 0 ) then

       ! We can have everything in Booles and
       ! one 3/8 rule.
       mBo = N / 4 - 1
       m38 = 1

       ! This is the actual step size for all the Boole regions
       mD = (b-a) / (real(mBo,dp) + Boole_Simp38_split)

       curB = a + mD

    else if ( mod(N,4) == 2 ) then

       ! We can have everything in Booles and one 3 rule and one 3/8 rules
       mBo = ( N - 6 ) / 4
       m38 = 1

       ! This is the percentage of the step size
       D = S_38_E / ( S_3_E + S_38_E )

       ! This is the actual step size for all the Boole regions
       mD = ( b - a ) / (real(mBo,dp) + (2._dp - D) * Boole_Simp38_split)

       ! We start with one Simpson 3
       ! First we need to calculate the step size of this one
       D = ( b - a - mBo * mD ) * ( 1._dp - D )
       curB = a + D
       call Simpson_3_3(x(i),w(i),curA,curB)
       i = i + 2

       ! enter the first Boole rule
       curA = curB
       curB = curA + mD

    else if ( mod(N,4) == 3 ) then

       ! We can have everything in Booles and 2 3/8 rules
       mBo = ( N - 3 ) / 4 - 1
       m38 = 2

       mD = ( b - a ) / (real(mBo,dp) + 2._dp * Boole_Simp38_split)

       ! We start with a Simpson 3/8
       i = 1
       ! First we need to calculate the step size of this one
       curB = ( b - a - mBo * mD ) * .5_dp
       curA = a
       curB = a + curB
       call Simpson_38_4(x(i),w(i),curA,curB)
       i = i + 3

       ! enter the first Boole rule
       curA = curB
       curB = curA + mD
       

    end if

    do j = 1 , mBo
       call Booles_5(x(i),w(i),curA,curB)
       curA = curB
       curB = curB + mD
       i = i + 4
    end do

    if ( 0 < m38 ) then
       ! We need to populate the last sequent (it must extend to b)
       call Simpson_38_4(x(i),w(i),curA,b) 
       i = i + 3
    end if

  end subroutine Booles_Simpson_38_3_rule

  ! We use a combination of the Simpson 3/8 and Simpson 3-point rule
  ! to generate arbitrary N rules
  subroutine Simpson_38_3_rule(N,x,w,a,b)
    integer, intent(in) :: N
    real(dp), intent(out) :: x(N), w(N)
    real(dp), intent(in) :: a, b

    ! S_38_E = 6480._dp ** .2_dp
    real(dp), parameter :: S_38_E = 5.785155024015764_dp
    ! S_3_E  = 2880._dp ** .2_dp
    real(dp), parameter :: S_3_E  = 4.919018971698727_dp
    real(dp) :: D, mD, curA, curB
    integer :: i, j, m38, m3

    if ( N < 1 ) &
         call die('Cannot integrate a negative &
         &number of points.')

    ! First initalize (we are doing additions)
    x(:) = 0._dp
    w(:) = 0._dp

    ! We need to ensure that we can do something (otherwise revert to the
    ! easy cases)
    if ( N <= 5 ) then
       if      ( N == 1 ) then
          x(1) = a + (b-a)*.5_dp
          w(1) = b-a
       else if ( N == 2 ) then
          call Trapez_2(x,w,a,b)
       else if ( N == 3 ) then
          call Simpson_3_3(x,w,a,b)
       else if ( N == 4 ) then
          call Simpson_38_4(x,w,a,b)
       else if ( N == 5 ) then
          call Simpson_3_3(x(1),w(1),a,(a+b)*.5_dp)
          call Simpson_3_3(x(3),w(3),(a+b)*.5_dp,b)
       end if

       return

    end if

    ! We need to calculate the best separation, we here explain the procedure
    ! for one 3/8 and one 3
    ! We can increase the precision
    ! by utilizing the known error term in the 
    ! methods.
    ! Simpson     : -(b-a)^5/2880
    ! Simpson 3/8 : -(b-a)^5/6480
    ! Best separation should be 2880^{1/5}/(2880^{1/5}+6480^{1/5}) ~ 0.45954
    ! So we take Simpson on (b-a) * 0.45954
    ! So we take Simpson 3/8 on (b-a) * 0.54046
    ! This is easily generaziled for N 3/8 separations and one 3-rule

    i = 1 
    curA = a
    
    if ( mod(N,3) == 1 ) then

       ! Everything can be in the m38 method
       m38 = (N - 1) / 3
       m3  = 0

       ! We do it in here to not confuse the logic of the routine
       ! First we generate the step size of each block
       mD = (b-a)/real(m38,dp)

       curB = a + mD

    else if ( mod(N,3) == 0 ) then

       ! We can have everything in m38 number of 3/8 rules and
       ! one 3 rule.
       m38 = N / 3 - 1
       m3  = 1

       ! This is the percentage of the step size
       D = S_38_E / ( S_3_E + S_38_E * m38 )
       ! This is the actual step size for all the m38 regions
       mD = (b-a) * D

       curB = a + mD

    else if ( mod(N,3) == 2 ) then

       ! We can have everything in m38 number of 3/8 rules and
       ! two 3 rule.
       m38 = (N - 5) / 3
       m3  = 2

       ! This is the percentage of the step size
       D = S_38_E * .5_dp / ( S_3_E + .5_dp * S_38_E * m38 )
       ! This is the actual step size for all the m38 regions
       mD = (b-a) * D

       ! We start with a Simpson 3
       ! First we need to calculate the step size of this one
       curB = (b - a - m38 * mD) * .5_dp
       curB = curA + curB
       call Simpson_3_3(x(i),w(i),curA,curB)
       i = i + 2

       ! enter the first 3/8 rule
       curA = curB
       curB = curA + mD

    end if

    do j = 1 , m38
       call Simpson_38_4(x(i),w(i),curA,curB)
       curA = curB
       curB = curB + mD
       i = i + 3
    end do

    if ( m3 > 0 ) then
       ! We need to populate the last sequent (it must extend to b)
       call Simpson_3_3(x(i),w(i),curA,b) 
       i = i + 2
    end if


  end subroutine Simpson_38_3_rule

  subroutine Booles_5(x,w,a,b)
    real(dp), intent(in out) :: x(5), w(5)
    real(dp), intent(in) :: a, b
    real(dp) :: h

    h = ( b - a ) * .25_dp
    x(1) = a
    x(2) = x(1) + h
    x(3) = x(2) + h
    x(4) = x(3) + h
    x(5) = x(4) + h

    h = ( b - a ) / 90._dp
    w(1) = w(1) + h *  7._dp
    w(2) = w(2) + h * 32._dp
    w(3) = w(3) + h * 12._dp
    w(4) = w(4) + h * 32._dp
    w(5) = w(5) + h *  7._dp

  end subroutine Booles_5

  subroutine Simpson_38_4(x,w,a,b)
    real(dp), intent(in out) :: x(4), w(4)
    real(dp), intent(in) :: a, b
    real(dp) :: h

    h = ( b - a ) / 3._dp
    x(1) = a
    x(2) = x(1) + h
    x(3) = x(2) + h
    x(4) = x(3) + h

    h = ( b - a ) * .125_dp    
    w(1) = w(1) + h
    w(2) = w(2) + h * 3._dp
    w(3) = w(3) + h * 3._dp
    w(4) = w(4) + h

  end subroutine Simpson_38_4

  subroutine Simpson_3_3(x,w,a,b)
    real(dp), intent(in out) :: x(3), w(3)
    real(dp), intent(in) :: a, b
    real(dp) :: h

    h = ( b - a ) * .5_dp
    x(1) = a
    x(2) = x(1) + h
    x(3) = x(2) + h

    h = ( b - a ) / 6._dp    
    w(1) = w(1) + h
    w(2) = w(2) + h * 4._dp
    w(3) = w(3) + h

  end subroutine Simpson_3_3

  subroutine Trapez_2(x,w,a,b)
    real(dp), intent(in out) :: x(2), w(2)
    real(dp), intent(in) :: a, b
    real(dp) :: h
    
    x(1) = a
    x(2) = b

    h = (b-a) * .5_dp    
    w(1) = w(1) + h
    w(2) = w(2) + h

  end subroutine Trapez_2

  ! Distribute 'N' points on segments composed
  ! of:
  !   len(1) -- len(2)
  !   len(2) -- len(3)
  !    ....
  ! Do this so that each segment is as evenly distributed
  ! if there are not enough points to have one on each segment
  ! a zero will be returned for all segments.
  subroutine n_distribute(N,N_seg,len,Ni)
    ! Total number of points
    integer, intent(in) :: N
    ! Number of segments
    integer, intent(in) :: N_seg
    ! Length specifications of each segment
    real(dp), intent(in) :: len(N_seg+1)
    ! Number of returned segments
    integer, intent(out) :: Ni(N_seg)

    real(dp) :: tot
    integer :: i

    ! Ensure that there are at least one point per
    ! segment
    if ( N < N_seg ) then
       Ni = 0
       return
    end if

    ! first initialize
    tot = 0._dp
    do i = 1 , N_seg
       tot = tot + abs(len(i+1) - len(i))
    end do

    ! Estimate initial distribution
    do i = 1 , N_seg
       Ni(i) = nint(abs(len(i+1) - len(i)) / tot)
    end do

    ! Even out the number of points

    do while ( sum(Ni) > N )
       i = maxloc(Ni,1)
       Ni(i) = Ni(i) - 1
    end do
    
    do while ( sum(Ni) < N )
       i = minloc(Ni,1)
       Ni(i) = Ni(i) + 1
    end do
    
  end subroutine n_distribute

end module m_integrate
    
    
    
  

  
