! Module for regular interpolation...

! Fully created by Nick Papior Andersen


module m_interpolate

  ! Generic module for interpolation
  integer, parameter :: dp = selected_real_kind(14,100)

  private

  public :: crt_pivot
  public :: interp_linear
  public :: interp_spline
  public :: prep_spline

contains

  subroutine crt_pivot(N,x,ipvt)
    integer, intent(in) :: N
    real(dp), intent(in) :: x(N)
    integer, intent(out) :: ipvt(N)
    integer :: i, j, ip, i_cur
    integer :: c_low, c_high
    real(dp) :: x_cur, xd

    if ( N < 1 ) then
      return
    else if ( N == 1 ) then
      ipvt(1) = 1
      return
    else if ( N == 2 ) then
      if ( x(1) > x(2) ) then
        ipvt(1) = 2
        ipvt(2) = 1
      else
        ipvt(1) = 1
        ipvt(2) = 2
      end if
      return
    end if

    ! Assume it is sorted
    do i = 1 , N
      ipvt(i) = i
    end do

    ! Now rotate indices until they match
    do i_cur = 2, N
      ! Current pivoting index
      x_cur = x(ipvt(i_cur))
      ! rotate and swap if possible
      i = i_cur - 1
      do while ( x(ipvt(i)) > x_cur )
        ! swap back
        ip = ipvt(i+1)
        ipvt(i+1) = ipvt(i)
        ipvt(i) = ip
        i = i - 1
        if ( i <= 0 ) exit
      end do
    end do

  end subroutine crt_pivot

  ! A simple linear interpolation algorithm
  subroutine interp_linear(N,x,y,x0,y0)
    integer, intent(in) :: N
    real(dp), intent(in) :: x(N), y(N)
    real(dp), intent(in) :: x0
    real(dp), intent(out) :: y0
    
    ! Local variables
    integer :: i, ipvt(N)
    real(dp) :: b

    if ( N <= 0 ) then
      y0 = 0._dp
      return
    else if ( N == 1 ) then
      y0 = y(1)
      return
    end if

    ! Create pivoting array for correct sorting
    ! in case the user did not provide a consecutive
    ! ordering...
    call crt_pivot(N,x,ipvt)

    ! Do extrapolation
    !  1. lower x
    !  2. higher x
    if ( x0 <= x(ipvt(1)) ) then
       b  = (y(ipvt(2)) - y(ipvt(1)))/(x(ipvt(2)) - x(ipvt(1)))
       y0 = y(ipvt(1)) + b * (x0 - x(ipvt(1)))
       return
    else if ( x(ipvt(N)) < x0 ) then
       b  = (y(ipvt(N)) - y(ipvt(N-1)))/(x(ipvt(N)) - x(ipvt(N-1)))
       y0 = y(ipvt(N)) + b * (x0 - x(ipvt(N)))
       return
    end if

    do i = 1 , N - 1
       if ( x0 <= x(ipvt(i+1)) ) then
          b  = (y(ipvt(i+1)) - y(ipvt(i)))/(x(ipvt(i+1)) - x(ipvt(i)))
          y0 = y(ipvt(i)) + b * (x0 - x(ipvt(i)))
          return
       end if
    end do
          
  end subroutine interp_linear

  subroutine interp_spline(N,x,y,x0,y0,z)
    integer, intent(in) :: N
    real(dp), intent(in) :: x(N), y(N)
    real(dp), intent(in) :: x0
    real(dp), intent(out) :: y0
    real(dp), intent(in), optional :: z(N-2) ! z as calculated by prep_spline
    
    ! Local variables
    integer :: i, ipvt(N), idx, seg
    real(dp) :: b
    real(dp), allocatable :: zz(:)

    ! If we have 2 points (or less) we do linear interpolation
    if ( N <= 2 ) then
       call interp_linear(N,x,y,x0,y0)
       return
    end if

    ! Create pivoting array for correct sorting
    ! in case the user did not provide a consecutive
    ! ordering...
    call crt_pivot(N,x,ipvt)

    ! If we are out-of-bounds we do extrapolation (from a linear
    ! perspective)
    if ( x0 <= x(ipvt(1)) ) then
       b = (y(ipvt(2)) - y(ipvt(1)))/(x(ipvt(2)) - x(ipvt(1)))
       y0 = y(ipvt(1)) + b * (x0 - x(ipvt(1)))
       return
    else if ( x(ipvt(N)) <= x0 ) then
       b = (y(ipvt(N)) - y(ipvt(N-1)))/(x(ipvt(N)) - x(ipvt(N-1)))
       y0 = y(ipvt(N)) + b * (x0 - x(ipvt(N)))
       return
    end if

    if ( present(z) ) then
       ! We calculate from already calculated parameters...
       do i = 1 , N - 1 ! z(1:N-2)
          if ( x0 <= x(ipvt(i+1)) ) then
             idx = i
             exit
          end if
       end do
       if ( idx == 1 ) then
          seg = -1
       else if ( idx == N - 1 ) then
          seg = 1
       else
          seg = 0
       end if
       i = idx
       call interp(seg,x(ipvt(i)), y(ipvt(i)), &
            x(ipvt(i+1)), y(ipvt(i+1)), z(max(i-1,1)), x0, y0)
    else

       ! We now need to calculate it manually...
       do i = 1 , N - 1 ! z(1:N-1)
          if ( x0 <= x(ipvt(i+1)) ) then
             ! We have found the position
             idx = i
             exit
          end if
       end do

       ! Allocate zz 
       i = N - idx
       allocate(zz(max(i,2)))
       zz(1:2) = 0._dp

       ! only calculate the necessary z's
       call loc_prep_spline(N,x,y,ipvt,max(i,2),zz)

       if ( idx == 1 ) then
          zz(2) = zz(1)
          zz(1) = 0._dp
       else if ( idx == N - 1 ) then
          zz(1) = zz(2)
          zz(2) = 0._dp
       end if
       i = idx
       call interp(0,x(ipvt(i)), y(ipvt(i)), &
            x(ipvt(i+1)), y(ipvt(i+1)), zz(1), x0, y0)

       deallocate(zz)

    end if

  contains

    subroutine interp(seg,x0,y0,x1,y1,z,x,y)
      integer, intent(in) :: seg ! -1 (start), 0 (middle), 1 (end)
      real(dp), intent(in) :: x0, y0, x1, y1, z(2), x
      real(dp), intent(out) :: y
      real(dp) :: h, xd, b

      ! We have found the point...
      h  = x1 - x0
      xd = x  - x0
         
      if ( seg == 1 ) then
         ! third order polynomial
         y = - z(1) / ( 6._dp * h ) * xd
         ! second order polynomial
         y = ( y + z(1) * .5_dp ) * xd
         ! first order polynomial
         b = - h / 3._dp * z(1) + ( y1 - y0 ) / h
         y = ( y + b ) * xd
      else if ( seg == 0 ) then
         ! third order polynomial
         y = ( z(2) - z(1) ) / ( 6._dp * h ) * xd
         ! second order polynomial
         y = ( y + z(1) * .5_dp ) * xd
         ! first order polynomial
         b = - h / 6._dp * z(2)
         b = b - h / 3._dp * z(1) + ( y1 - y0 ) / h
         y = ( y + b ) * xd
      else if ( seg == -1 ) then
         ! third order polynomial
         y = z(1) / ( 6._dp * h ) * xd
         ! second order polynomial
         y = y * xd
         ! first order polynomial
         b = - h / 6._dp * z(1)
         b = b + ( y1 - y0 ) / h
         y = ( y + b ) * xd
      end if
      ! zeroth order
      y = y + y0
         
    end subroutine interp

  end subroutine interp_spline

  subroutine loc_prep_spline(N,x,y,ipvt,NZ,z)
    integer, intent(in) :: N
    real(dp), intent(in) :: x(N), y(N)
    integer, intent(in) :: ipvt(N)
    integer, intent(in):: NZ
    real(dp), intent(inout) :: z(NZ)

    ! Local variables
    integer :: i, iz
    real(dp) :: u(N-2), v(N-2)
    real(dp) :: b(2), h(2)

    ! Calculate u_i/v_i
    h(1) = x(ipvt(2)) - x(ipvt(1))
    b(1) = (y(ipvt(2)) - y(ipvt(1))) / h(1)
    h(2) = x(ipvt(3)) - x(ipvt(2))
    b(2) = (y(ipvt(3)) - y(ipvt(2))) / h(2)
    u(1) = 2._dp * ( h(1) + h(2) )
    v(1) = 6._dp * ( b(2) - b(1) )
    do i = 2 , N - 2
       ! Transfer h,b
       h(1) = h(2)
       b(1) = b(2)
       h(2) = x(ipvt(i+2)) - x(ipvt(i+1))
       b(2) = ( y(ipvt(i+2)) - y(ipvt(i+1)) ) / h(2)

       ! Calculate u_i,v_i
       u(i) = 2._dp * ( h(1) + h(2) ) - h(1) ** 2 / u(i-1)
       v(i) = 6._dp * ( b(2) - b(1) ) - h(1) * v(i-1) / u(i-1)

    end do

    ! for the case of a limited number of z's, we simply
    ! calculate them like this:
    ! Calculate z_i
    iz = min(NZ,N-2)
    z(iz) = v(N-2) / u(N-2)
    do i = N - 3 , 1 , -1
       iz = iz - 1
       if ( iz == 0 ) return
       h(1) = x(ipvt(i+2)) - x(ipvt(i+1))
       z(iz) = ( v(i) - h(1) * z(iz+1) ) / u(i)
    end do

  end subroutine loc_prep_spline
    

  subroutine prep_spline(N,x,y,z,NZ)
    integer, intent(in) :: N
    real(dp), intent(in) :: x(N), y(N)
    real(dp), intent(out) :: z(:)
    integer, intent(in), optional :: NZ
    
    ! Local variables
    integer :: ipvt(N)

    if ( N <= 2 ) then
       ! We cannot do anything...
       return
    end if

    ! Create pivoting array for correct sorting
    ! in case the user did not provide a consecutive
    ! ordering...
    call crt_pivot(N,x,ipvt)
    
    if ( present(NZ) ) then
       call loc_prep_spline(N,x,y,ipvt,NZ,z)
    else
       call loc_prep_spline(N,x,y,ipvt,N-2,z)
    end if

  end subroutine prep_spline

#ifdef __M_INTERPOLATE_DEBUG
  subroutine pivot_extract(N, orig, ipvt, out)
    integer, intent(in) :: N
    real(dp), intent(in) :: orig(N)
    integer, intent(out) :: ipvt(N)
    real(dp), intent(out) :: out(N)

    integer :: i

    do i = 1, N
      out(i) = orig(ipvt(i))
    end do

  end subroutine pivot_extract
#endif

end module m_interpolate


#ifdef __M_INTERPOLATE_DEBUG
program test_interpolate_m

  use m_interpolate

  integer, parameter :: dp = selected_real_kind(14,100)

  integer :: N
  real(dp), allocatable :: x(:), y(:)
  real(dp) :: x0, y0

  N = 0
  call my_alloc(N)

  x0 = 1._dp
  call interp_spline(N, x, y, x0, y0)
  print *, 'Success N = ', N, y0 == 0._dp

  N = 1
  call my_alloc(N)
  call interp_spline(N, x, y, x0, y0)
  print *, 'Success N = ', N, y0 == 0._dp
  
  N = 2
  call my_alloc(N)
  x(1) = 0._dp
  y(1) = 3._dp
  x(2) = 1._dp
  y(2) = 4._dp
  x(3) = 2._dp
  y(3) = 5._dp
  x0 = 1.5_dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 123 = ', x0, ' y0[e/i] = ',4.5_dp, y0, 4.5_dp == y0
  x0 = 3._dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 123 = ', x0, ' y0[e/i] = ',6._dp, y0, 6._dp == y0
  x0 = -1._dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 123 = ', x0, ' y0[e/i] = ',2._dp, y0, 2._dp == y0

  x(1) = 0._dp
  y(1) = 3._dp
  x(3) = 1._dp
  y(3) = 4._dp
  x(2) = 2._dp
  y(2) = 5._dp
  x0 = 1.5_dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 132 = ', x0, ' y0[e/i] = ',4.5_dp, y0, 4.5_dp == y0
  x0 = 3._dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 132 = ', x0, ' y0[e/i] = ',6._dp, y0, 6._dp == y0
  x0 = -1._dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 132 = ', x0, ' y0[e/i] = ',2._dp, y0, 2._dp == y0

  x(2) = 0._dp
  y(2) = 3._dp
  x(1) = 1._dp
  y(1) = 4._dp
  x(3) = 2._dp
  y(3) = 5._dp
  x0 = 1.5_dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 213 = ', x0, ' y0[e/i] = ',4.5_dp, y0, 4.5_dp == y0
  x0 = 3._dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 213 = ', x0, ' y0[e/i] = ',6._dp, y0, 6._dp == y0
  x0 = -1._dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 213 = ', x0, ' y0[e/i] = ',2._dp, y0, 2._dp == y0
  
  x(2) = 0._dp
  y(2) = 3._dp
  x(3) = 1._dp
  y(3) = 4._dp
  x(1) = 2._dp
  y(1) = 5._dp
  x0 = 1.5_dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 231 = ', x0, ' y0[e/i] = ',4.5_dp, y0, 4.5_dp == y0
  x0 = 3._dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 231 = ', x0, ' y0[e/i] = ',6._dp, y0, 6._dp == y0
  x0 = -1._dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 231 = ', x0, ' y0[e/i] = ',2._dp, y0, 2._dp == y0

  x(3) = 0._dp
  y(3) = 3._dp
  x(1) = 1._dp
  y(1) = 4._dp
  x(2) = 2._dp
  y(2) = 5._dp
  x0 = 1.5_dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 312 = ', x0, ' y0[e/i] = ',4.5_dp, y0, 4.5_dp == y0
  x0 = 3._dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 312 = ', x0, ' y0[e/i] = ',6._dp, y0, 6._dp == y0
  x0 = -1._dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 312 = ', x0, ' y0[e/i] = ',2._dp, y0, 2._dp == y0

  x(3) = 0._dp
  y(3) = 3._dp
  x(2) = 1._dp
  y(2) = 4._dp
  x(1) = 2._dp
  y(1) = 5._dp
  x0 = 1.5_dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 321 = ', x0, ' y0[e/i] = ',4.5_dp, y0, 4.5_dp == y0
  x0 = 3._dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 321 = ', x0, ' y0[e/i] = ',6._dp, y0, 6._dp == y0
  x0 = -1._dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 321 = ', x0, ' y0[e/i] = ',2._dp, y0, 2._dp == y0

  N = 3
  call my_alloc(N)
  x(1) = 171.48439838412170_dp
  y(1) = -0.50797873167809759_dp
  x(2) = 171.67036059152491_dp
  y(2) = -0.50595438288526884_dp
  x(3) = 172.00273057484208_dp
  y(3) = -0.50236599173884533_dp
  x0 = 171.99999999999986_dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 = ', x0, ' y0[e/i] = ',y0

  N = 5
  call my_alloc(N)
  x(1) = 2.7305732180000000E-003_dp
  y(1) = -0.50236593563039289
  x(2) = 2.4257984650000002E-003_dp
  y(2) = -0.49724798813741200_dp
  x(3) = 1.6736149160000001E-004_dp
  y(3) = -0.45274740272530828_dp
  x(4) = 9.8428249089999996E-006_dp
  y(4) = -0.44301466374634341_dp
  x(5) = -4.8751318840000004E-003_dp
  y(5) = -0.44076526944391359_dp
  x0 = 0._dp
  call interp_spline(N, x, y, x0, y0)
  print '(a,es10.3,a,2(tr1,es10.3),tr1l)', 'x0 = ', x0, ' y0[e/i] = ',y0

contains

  subroutine my_alloc(N)
    integer, intent(in) :: N
    if ( allocated(x) ) deallocate(x,y)
    allocate(x(N),y(N))
  end subroutine my_alloc
  
end program test_interpolate_m

#endif



