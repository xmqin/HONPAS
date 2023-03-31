!!@LICENSE
!
!******************************************************************************
! MODULE array
! Provides utilities for copying data from one array shape to another
! array shape.
! Written by Nick R. Papior, Feb.2017
!******************************************************************************
!
!   PUBLIC procedures available from this module:
! array_copy    : Copies an array from one dimension to another array with another dimension
! array_add     : Adds an array from one dimension to another array with another dimension
!
!   PUBLIC parameters, types, and variables available from this module:
! none
!
! The intent of this module is to bypass the use of the intrinsic reshape
! function.
! Sadly many compilers implement the reshape with a guard on pointers such
! that a temporary array is created upon reshaping. For very large segments
! this turns out to be a substantial memory consumption on the heap which
! may easily be circumvented by doing explicit copying.
! These routines are extremely simple and does nothing but copy data from
! one shape to another.
!
! Currently only these copies and adds are implemented:
!
!   1D -> 2D [integer, real(sp), real(dp)]
!   1D -> 3D [integer, real(sp), real(dp)]
!   2D -> 1D [integer, real(sp), real(dp)]
!   3D -> 1D [integer, real(sp), real(dp)]

module m_array

  use precision, only: sp, dp
  use sys,       only: die       ! Termination routine

  private

  public :: array_copy
  public :: array_add

  interface array_copy
     module procedure ac_1d_2d_ip, ac_1d_2d_sp, ac_1d_2d_dp
     module procedure ac_1d_3d_ip, ac_1d_3d_sp, ac_1d_3d_dp
     module procedure ac_2d_1d_ip, ac_2d_1d_sp, ac_2d_1d_dp
     module procedure ac_3d_1d_ip, ac_3d_1d_sp, ac_3d_1d_dp
     module procedure ac_4d_1d_ip, ac_4d_1d_sp, ac_4d_1d_dp
  end interface array_copy

  interface array_add
     module procedure aa_1d_2d_ip, aa_1d_2d_sp, aa_1d_2d_dp
     module procedure aa_1d_3d_ip, aa_1d_3d_sp, aa_1d_3d_dp
     module procedure aa_1d_4d_ip, aa_1d_4d_sp, aa_1d_4d_dp
     module procedure aa_2d_1d_ip, aa_2d_1d_sp, aa_2d_1d_dp
     module procedure aa_3d_1d_ip, aa_3d_1d_sp, aa_3d_1d_dp
     module procedure aa_4d_1d_ip, aa_4d_1d_sp, aa_4d_1d_dp
  end interface array_add

contains

  ! Copies a 1D array to a 2D array
  subroutine ac_1d_2d_ip(from_1D, to_1D, in1D, from_2D, to_2D, out2D)
    integer, intent(in) :: from_1D, to_1D
    integer, intent(in) :: in1D(:)
    integer, intent(in) :: from_2D(2), to_2D(2)
    integer, intent(inout) :: out2D(:,:)

    ! Local variables for copying
    integer :: i1D, i2D, j2D

    i2D = from_2D(1)
    j2D = from_2D(2)
    do i1D = from_1D, to_1D
       out2D(i2D, j2D) = in1D(i1D)
       i2D = i2D + 1
       if ( i2D > to_2D(1) ) then
          i2D = from_2D(1)
          j2D = j2D + 1
       end if
    end do

    if ( i2D /= from_2D(1) ) &
         call die("integer: 1D->2D failed (i)")
    if ( j2D <= to_2D(2) ) &
         call die("integer: 1D->2D failed (j)")

  end subroutine ac_1d_2d_ip
  subroutine ac_1d_2d_sp(from_1D, to_1D, in1D, from_2D, to_2D, out2D)
    integer, intent(in) :: from_1D, to_1D
    real(sp), intent(in) :: in1D(:)
    integer, intent(in) :: from_2D(2), to_2D(2)
    real(sp), intent(inout) :: out2D(:,:)

    ! Local variables for copying
    integer :: i1D, i2D, j2D

    i2D = from_2D(1)
    j2D = from_2D(2)
    do i1D = from_1D, to_1D
       out2D(i2D, j2D) = in1D(i1D)
       i2D = i2D + 1
       if ( i2D > to_2D(1) ) then
          i2D = from_2D(1)
          j2D = j2D + 1
       end if
    end do

    if ( i2D /= from_2D(1) ) &
         call die("real: 1D->2D failed (i)")
    if ( j2D <= to_2D(2) ) &
         call die("real: 1D->2D failed (j)")

  end subroutine ac_1d_2d_sp
  subroutine ac_1d_2d_dp(from_1D, to_1D, in1D, from_2D, to_2D, out2D)
    integer, intent(in) :: from_1D, to_1D
    real(dp), intent(in) :: in1D(:)
    integer, intent(in) :: from_2D(2), to_2D(2)
    real(dp), intent(inout) :: out2D(:,:)

    ! Local variables for copying
    integer :: i1D, i2D, j2D

    i2D = from_2D(1)
    j2D = from_2D(2)
    do i1D = from_1D, to_1D
       out2D(i2D, j2D) = in1D(i1D)
       i2D = i2D + 1
       if ( i2D > to_2D(1) ) then
          i2D = from_2D(1)
          j2D = j2D + 1
       end if
    end do

    if ( i2D /= from_2D(1) ) &
         call die("double: 1D->2D failed (i)")
    if ( j2D <= to_2D(2) ) &
         call die("double: 1D->2D failed (j)")

  end subroutine ac_1d_2d_dp

  ! Copies a 1D array to a 3D array
  subroutine ac_1d_3d_ip(from_1D, to_1D, in1D, from_3D, to_3D, out3D)
    integer, intent(in) :: from_1D, to_1D
    integer, intent(in) :: in1D(:)
    integer, intent(in) :: from_3D(3), to_3D(3)
    integer, intent(inout) :: out3D(:,:,:)

    ! Local variables for copying
    integer :: i1D, i3D, j3D, k3D

    i3D = from_3D(1)
    j3D = from_3D(2)
    k3D = from_3D(3)
    do i1D = from_1D, to_1D
       out3D(i3D, j3D, k3D) = in1D(i1D)
       i3D = i3D + 1
       if ( i3D > to_3D(1) ) then
          i3D = from_3D(1)
          j3D = j3D + 1
       end if
       if ( j3D > to_3D(2) ) then
          j3D = from_3D(2)
          k3D = k3D + 1
       end if
    end do

    if ( i3D /= from_3D(1) ) &
         call die("integer: 1D->3D failed (i)")
    if ( j3D /= from_3D(2) ) &
         call die("integer: 1D->3D failed (j)")
    if ( k3D <= to_3D(3) ) &
         call die("integer: 1D->3D failed (k)")

  end subroutine ac_1d_3d_ip
  subroutine ac_1d_3d_sp(from_1D, to_1D, in1D, from_3D, to_3D, out3D)
    integer, intent(in) :: from_1D, to_1D
    real(sp), intent(in) :: in1D(:)
    integer, intent(in) :: from_3D(3), to_3D(3)
    real(sp), intent(inout) :: out3D(:,:,:)

    ! Local variables for copying
    integer :: i1D, i3D, j3D, k3D

    i3D = from_3D(1)
    j3D = from_3D(2)
    k3D = from_3D(3)
    do i1D = from_1D, to_1D
       out3D(i3D, j3D, k3D) = in1D(i1D)
       i3D = i3D + 1
       if ( i3D > to_3D(1) ) then
          i3D = from_3D(1)
          j3D = j3D + 1
       end if
       if ( j3D > to_3D(2) ) then
          j3D = from_3D(2)
          k3D = k3D + 1
       end if
    end do

    if ( i3D /= from_3D(1) ) &
         call die("real: 1D->3D failed (i)")
    if ( j3D /= from_3D(2) ) &
         call die("real: 1D->3D failed (j)")
    if ( k3D <= to_3D(3) ) &
         call die("real: 1D->3D failed (k)")

  end subroutine ac_1d_3d_sp
  subroutine ac_1d_3d_dp(from_1D, to_1D, in1D, from_3D, to_3D, out3D)
    integer, intent(in) :: from_1D, to_1D
    real(dp), intent(in) :: in1D(:)
    integer, intent(in) :: from_3D(3), to_3D(3)
    real(dp), intent(inout) :: out3D(:,:,:)

    ! Local variables for copying
    integer :: i1D, i3D, j3D, k3D

    i3D = from_3D(1)
    j3D = from_3D(2)
    k3D = from_3D(3)
    do i1D = from_1D, to_1D
       out3D(i3D, j3D, k3D) = in1D(i1D)
       i3D = i3D + 1
       if ( i3D > to_3D(1) ) then
          i3D = from_3D(1)
          j3D = j3D + 1
       end if
       if ( j3D > to_3D(2) ) then
          j3D = from_3D(2)
          k3D = k3D + 1
       end if
    end do

    if ( i3D /= from_3D(1) ) &
         call die("double: 1D->3D failed (i)")
    if ( j3D /= from_3D(2) ) &
         call die("double: 1D->3D failed (j)")
    if ( k3D <= to_3D(3) ) &
         call die("double: 1D->3D failed (k)")

  end subroutine ac_1d_3d_dp

  ! Copies a 2D array to a 1D array
  subroutine ac_2d_1d_ip(from_2D, to_2D, in2D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_2D(2), to_2D(2)
    integer, intent(in) :: in2D(:,:)
    integer, intent(in) :: from_1D, to_1D
    integer, intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i2D, j2D, i1D

    i1D = from_1D
    do j2D = from_2D(2), to_2D(2)
       do i2D = from_2D(1), to_2D(1)
          out1D(i1D) = in2D(i2D, j2D)
          i1D = i1D + 1
       end do
    end do

    if ( i1D <= to_1D ) &
         call die("integer: 2D->1D failed")

  end subroutine ac_2d_1d_ip
  subroutine ac_2d_1d_sp(from_2D, to_2D, in2D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_2D(2), to_2D(2)
    real(sp), intent(in) :: in2D(:,:)
    integer, intent(in) :: from_1D, to_1D
    real(sp), intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i2D, j2D, i1D

    i1D = from_1D
    do j2D = from_2D(2), to_2D(2)
       do i2D = from_2D(1), to_2D(1)
          out1D(i1D) = in2D(i2D, j2D)
          i1D = i1D + 1
       end do
    end do

    if ( i1D <= to_1D ) &
         call die("real: 2D->1D failed")

  end subroutine ac_2d_1d_sp
  subroutine ac_2d_1d_dp(from_2D, to_2D, in2D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_2D(2), to_2D(2)
    real(dp), intent(in) :: in2D(:,:)
    integer, intent(in) :: from_1D, to_1D
    real(dp), intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i2D, j2D, i1D

    i1D = from_1D
    do j2D = from_2D(2), to_2D(2)
       do i2D = from_2D(1), to_2D(1)
          out1D(i1D) = in2D(i2D, j2D)
          i1D = i1D + 1
       end do
    end do

    if ( i1D <= to_1D ) &
         call die("double: 2D->1D failed")

  end subroutine ac_2d_1d_dp
  
  ! Copies a 3D array to a 1D array
  subroutine ac_3d_1d_ip(from_3D, to_3D, in3D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_3D(3), to_3D(3)
    integer, intent(in) :: in3D(:,:,:)
    integer, intent(in) :: from_1D, to_1D
    integer, intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i3D, j3D, k3D, i1D

    i1D = from_1D
    do k3D = from_3D(3), to_3D(3)
       do j3D = from_3D(2), to_3D(2)
          do i3D = from_3D(1), to_3D(1)
             out1D(i1D) = in3D(i3D, j3D, k3D)
             i1D = i1D + 1
          end do
       end do
    end do

    if ( i1D <= to_1D ) &
         call die("integer: 3D->1D failed")

  end subroutine ac_3d_1d_ip
  subroutine ac_3d_1d_sp(from_3D, to_3D, in3D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_3D(3), to_3D(3)
    real(sp), intent(in) :: in3D(:,:,:)
    integer, intent(in) :: from_1D, to_1D
    real(sp), intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i3D, j3D, k3D, i1D

    i1D = from_1D
    do k3D = from_3D(3), to_3D(3)
       do j3D = from_3D(2), to_3D(2)
          do i3D = from_3D(1), to_3D(1)
             out1D(i1D) = in3D(i3D, j3D, k3D)
             i1D = i1D + 1
          end do
       end do
    end do

    if ( i1D <= to_1D ) &
         call die("real: 3D->1D failed")

  end subroutine ac_3d_1d_sp
  subroutine ac_3d_1d_dp(from_3D, to_3D, in3D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_3D(3), to_3D(3)
    real(dp), intent(in) :: in3D(:,:,:)
    integer, intent(in) :: from_1D, to_1D
    real(dp), intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i3D, j3D, k3D, i1D

    i1D = from_1D
    do k3D = from_3D(3), to_3D(3)
       do j3D = from_3D(2), to_3D(2)
          do i3D = from_3D(1), to_3D(1)
             out1D(i1D) = in3D(i3D, j3D, k3D)
             i1D = i1D + 1
          end do
       end do
    end do

    if ( i1D <= to_1D ) &
         call die("double: 3D->1D failed")

  end subroutine ac_3d_1d_dp

  ! Copies a 4D array to a 1D array
  subroutine ac_4d_1d_ip(from_4D, to_4D, in4D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_4D(4), to_4D(4)
    integer, intent(in) :: in4D(:,:,:,:)
    integer, intent(in) :: from_1D, to_1D
    integer, intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i4D, j4D, k4D, m4D, i1D

    i1D = from_1D
    do m4D = from_4D(4), to_4D(4)
       do k4D = from_4D(3), to_4D(3)
          do j4D = from_4D(2), to_4D(2)
             do i4D = from_4D(1), to_4D(1)
                out1D(i1D) = in4D(i4D, j4D, k4D, m4D)
                i1D = i1D + 1
             end do
          end do
       end do
    end do

    if ( i1D <= to_1D ) &
         call die("integer: 4D->1D failed")

  end subroutine ac_4d_1d_ip
  subroutine ac_4d_1d_sp(from_4D, to_4D, in4D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_4D(4), to_4D(4)
    real(sp), intent(in) :: in4D(:,:,:,:)
    integer, intent(in) :: from_1D, to_1D
    real(sp), intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i4D, j4D, k4D, m4D, i1D

    i1D = from_1D
    do m4D = from_4D(4), to_4D(4)
       do k4D = from_4D(3), to_4D(3)
          do j4D = from_4D(2), to_4D(2)
             do i4D = from_4D(1), to_4D(1)
                out1D(i1D) = in4D(i4D, j4D, k4D, m4D)
                i1D = i1D + 1
             end do
          end do
       end do
    end do
    
    if ( i1D <= to_1D ) &
         call die("real: 4D->1D failed")

  end subroutine ac_4d_1d_sp
  subroutine ac_4d_1d_dp(from_4D, to_4D, in4D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_4D(4), to_4D(4)
    real(dp), intent(in) :: in4D(:,:,:,:)
    integer, intent(in) :: from_1D, to_1D
    real(dp), intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i4D, j4D, k4D, m4D, i1D

    i1D = from_1D
    do m4D = from_4D(4), to_4D(4)
       do k4D = from_4D(3), to_4D(3)
          do j4D = from_4D(2), to_4D(2)
             do i4D = from_4D(1), to_4D(1)
                out1D(i1D) = in4D(i4D, j4D, k4D, m4D)
                i1D = i1D + 1
             end do
          end do
       end do
    end do
    
    if ( i1D <= to_1D ) &
         call die("double: 4D->1D failed")

  end subroutine ac_4d_1d_dp
  

  ! Adds a 1D array to a 2D array
  subroutine aa_1d_2d_ip(from_1D, to_1D, a1D, from_2D, to_2D, out2D)
    integer, intent(in) :: from_1D, to_1D
    integer, intent(in) :: a1D(:)
    integer, intent(in) :: from_2D(2), to_2D(2)
    integer, intent(inout) :: out2D(:,:)

    ! Local variables for copying
    integer :: i1D, i2D, j2D

    i2D = from_2D(1)
    j2D = from_2D(2)
    do i1D = from_1D, to_1D
       out2D(i2D, j2D) = out2D(i2D, j2D) + a1D(i1D)
       i2D = i2D + 1
       if ( i2D > to_2D(1) ) then
          i2D = from_2D(1)
          j2D = j2D + 1
       end if
    end do

    if ( i2D /= from_2D(1) ) &
         call die("integer: 1D+>2D failed (i)")
    if ( j2D <= to_2D(2) ) &
         call die("integer: 1D+>2D failed (j)")

  end subroutine aa_1d_2d_ip
  subroutine aa_1d_2d_sp(from_1D, to_1D, a1D, from_2D, to_2D, out2D)
    integer, intent(in) :: from_1D, to_1D
    real(sp), intent(in) :: a1D(:)
    integer, intent(in) :: from_2D(2), to_2D(2)
    real(sp), intent(inout) :: out2D(:,:)

    ! Local variables for copying
    integer :: i1D, i2D, j2D

    i2D = from_2D(1)
    j2D = from_2D(2)
    do i1D = from_1D, to_1D
       out2D(i2D, j2D) = out2D(i2D, j2D) + a1D(i1D)
       i2D = i2D + 1
       if ( i2D > to_2D(1) ) then
          i2D = from_2D(1)
          j2D = j2D + 1
       end if
    end do

    if ( i2D /= from_2D(1) ) &
         call die("real: 1D+>2D failed (i)")
    if ( j2D <= to_2D(2) ) &
         call die("real: 1D+>2D failed (j)")

  end subroutine aa_1d_2d_sp
  subroutine aa_1d_2d_dp(from_1D, to_1D, a1D, from_2D, to_2D, out2D)
    integer, intent(in) :: from_1D, to_1D
    real(dp), intent(in) :: a1D(:)
    integer, intent(in) :: from_2D(2), to_2D(2)
    real(dp), intent(inout) :: out2D(:,:)

    ! Local variables for copying
    integer :: i1D, i2D, j2D

    i2D = from_2D(1)
    j2D = from_2D(2)
    do i1D = from_1D, to_1D
       out2D(i2D, j2D) = out2D(i2D, j2D) + a1D(i1D)
       i2D = i2D + 1
       if ( i2D > to_2D(1) ) then
          i2D = from_2D(1)
          j2D = j2D + 1
       end if
    end do

    if ( i2D /= from_2D(1) ) &
         call die("double: 1D+>2D failed (i)")
    if ( j2D <= to_2D(2) ) &
         call die("double: 1D+>2D failed (j)")

  end subroutine aa_1d_2d_dp

  ! Adds a 1D array to a 3D array
  subroutine aa_1d_3d_ip(from_1D, to_1D, a1D, from_3D, to_3D, out3D)
    integer, intent(in) :: from_1D, to_1D
    integer, intent(in) :: a1D(:)
    integer, intent(in) :: from_3D(3), to_3D(3)
    integer, intent(inout) :: out3D(:,:,:)

    ! Local variables for copying
    integer :: i1D, i3D, j3D, k3D

    i3D = from_3D(1)
    j3D = from_3D(2)
    k3D = from_3D(3)
    do i1D = from_1D, to_1D
       out3D(i3D, j3D, k3D) = out3D(i3D, j3D, k3D) + a1D(i1D)
       i3D = i3D + 1
       if ( i3D > to_3D(1) ) then
          i3D = from_3D(1)
          j3D = j3D + 1
       end if
       if ( j3D > to_3D(2) ) then
          j3D = from_3D(2)
          k3D = k3D + 1
       end if
    end do

    if ( i3D /= from_3D(1) ) &
         call die("integer: 1D+>3D failed (i)")
    if ( j3D /= from_3D(2) ) &
         call die("integer: 1D+>3D failed (j)")
    if ( k3D <= to_3D(3) ) &
         call die("integer: 1D+>3D failed (k)")

  end subroutine aa_1d_3d_ip
  subroutine aa_1d_3d_sp(from_1D, to_1D, a1D, from_3D, to_3D, out3D)
    integer, intent(in) :: from_1D, to_1D
    real(sp), intent(in) :: a1D(:)
    integer, intent(in) :: from_3D(3), to_3D(3)
    real(sp), intent(inout) :: out3D(:,:,:)

    ! Local variables for copying
    integer :: i1D, i3D, j3D, k3D

    i3D = from_3D(1)
    j3D = from_3D(2)
    k3D = from_3D(3)
    do i1D = from_1D, to_1D
       out3D(i3D, j3D, k3D) = out3D(i3D, j3D, k3D) + a1D(i1D)
       i3D = i3D + 1
       if ( i3D > to_3D(1) ) then
          i3D = from_3D(1)
          j3D = j3D + 1
       end if
       if ( j3D > to_3D(2) ) then
          j3D = from_3D(2)
          k3D = k3D + 1
       end if
    end do

    if ( i3D /= from_3D(1) ) &
         call die("real: 1D+>3D failed (i)")
    if ( j3D /= from_3D(2) ) &
         call die("real: 1D+>3D failed (j)")
    if ( k3D <= to_3D(3) ) &
         call die("real: 1D+>3D failed (k)")

  end subroutine aa_1d_3d_sp
  subroutine aa_1d_3d_dp(from_1D, to_1D, a1D, from_3D, to_3D, out3D)
    integer, intent(in) :: from_1D, to_1D
    real(dp), intent(in) :: a1D(:)
    integer, intent(in) :: from_3D(3), to_3D(3)
    real(dp), intent(inout) :: out3D(:,:,:)

    ! Local variables for copying
    integer :: i1D, i3D, j3D, k3D

    i3D = from_3D(1)
    j3D = from_3D(2)
    k3D = from_3D(3)
    do i1D = from_1D, to_1D
       out3D(i3D, j3D, k3D) = out3D(i3D, j3D, k3D) + a1D(i1D)
       i3D = i3D + 1
       if ( i3D > to_3D(1) ) then
          i3D = from_3D(1)
          j3D = j3D + 1
       end if
       if ( j3D > to_3D(2) ) then
          j3D = from_3D(2)
          k3D = k3D + 1
       end if
    end do

    if ( i3D /= from_3D(1) ) &
         call die("double: 1D+>3D failed (i)")
    if ( j3D /= from_3D(2) ) &
         call die("double: 1D+>3D failed (j)")
    if ( k3D <= to_3D(3) ) &
         call die("double: 1D+>3D failed (k)")

  end subroutine aa_1d_3d_dp

  ! Adds a 1D array to a 4D array
  subroutine aa_1d_4d_ip(from_1D, to_1D, a1D, from_4D, to_4D, out4D)
    integer, intent(in) :: from_1D, to_1D
    integer, intent(in) :: a1D(:)
    integer, intent(in) :: from_4D(4), to_4D(4)
    integer, intent(inout) :: out4D(:,:,:,:)

    ! Local variables for copying
    integer :: i1D, i4D, j4D, k4D, m4D

    i4D = from_4D(1)
    j4D = from_4D(2)
    k4D = from_4D(3)
    m4D = from_4D(4)
    do i1D = from_1D, to_1D
       out4D(i4D, j4D, k4D, m4D) = out4D(i4D, j4D, k4D, m4D) + a1D(i1D)
       i4D = i4D + 1
       if ( i4D > to_4D(1) ) then
          i4D = from_4D(1)
          j4D = j4D + 1
       end if
       if ( j4D > to_4D(2) ) then
          j4D = from_4D(2)
          k4D = k4D + 1
       end if
       if ( k4D > to_4D(3) ) then
          k4D = from_4D(3)
          m4D = m4D + 1
       end if
    end do

    if ( i4D /= from_4D(1) ) &
         call die("integer: 1D+>4D failed (i)")
    if ( j4D /= from_4D(2) ) &
         call die("integer: 1D+>4D failed (j)")
    if ( k4D /= from_4D(3) ) &
         call die("integer: 1D+>4D failed (k)")
    if ( m4D <= to_4D(4) ) &
         call die("integer: 1D+>4D failed (m)")

  end subroutine aa_1d_4d_ip
  subroutine aa_1d_4d_sp(from_1D, to_1D, a1D, from_4D, to_4D, out4D)
    integer, intent(in) :: from_1D, to_1D
    real(sp), intent(in) :: a1D(:)
    integer, intent(in) :: from_4D(4), to_4D(4)
    real(sp), intent(inout) :: out4D(:,:,:,:)

    ! Local variables for copying
    integer :: i1D, i4D, j4D, k4D, m4D

    i4D = from_4D(1)
    j4D = from_4D(2)
    k4D = from_4D(3)
    m4D = from_4D(4)
    do i1D = from_1D, to_1D
       out4D(i4D, j4D, k4D, m4D) = out4D(i4D, j4D, k4D, m4D) + a1D(i1D)
       i4D = i4D + 1
       if ( i4D > to_4D(1) ) then
          i4D = from_4D(1)
          j4D = j4D + 1
       end if
       if ( j4D > to_4D(2) ) then
          j4D = from_4D(2)
          k4D = k4D + 1
       end if
       if ( k4D > to_4D(3) ) then
          k4D = from_4D(3)
          m4D = m4D + 1
       end if
    end do

    if ( i4D /= from_4D(1) ) &
         call die("real: 1D+>4D failed (i)")
    if ( j4D /= from_4D(2) ) &
         call die("real: 1D+>4D failed (j)")
    if ( k4D /= from_4D(3) ) &
         call die("real: 1D+>4D failed (k)")
    if ( m4D <= to_4D(4) ) &
         call die("real: 1D+>4D failed (m)")

  end subroutine aa_1d_4d_sp
  subroutine aa_1d_4d_dp(from_1D, to_1D, a1D, from_4D, to_4D, out4D)
    integer, intent(in) :: from_1D, to_1D
    real(dp), intent(in) :: a1D(:)
    integer, intent(in) :: from_4D(4), to_4D(4)
    real(dp), intent(inout) :: out4D(:,:,:,:)

    ! Local variables for copying
    integer :: i1D, i4D, j4D, k4D, m4D

    i4D = from_4D(1)
    j4D = from_4D(2)
    k4D = from_4D(3)
    m4D = from_4D(4)
    do i1D = from_1D, to_1D
       out4D(i4D, j4D, k4D, m4D) = out4D(i4D, j4D, k4D, m4D) + a1D(i1D)
       i4D = i4D + 1
       if ( i4D > to_4D(1) ) then
          i4D = from_4D(1)
          j4D = j4D + 1
       end if
       if ( j4D > to_4D(2) ) then
          j4D = from_4D(2)
          k4D = k4D + 1
       end if
       if ( k4D > to_4D(3) ) then
          k4D = from_4D(3)
          m4D = m4D + 1
       end if
    end do

    if ( i4D /= from_4D(1) ) &
         call die("double: 1D+>4D failed (i)")
    if ( j4D /= from_4D(2) ) &
         call die("double: 1D+>4D failed (j)")
    if ( k4D /= from_4D(3) ) &
         call die("double: 1D+>4D failed (k)")
    if ( m4D <= to_4D(4) ) &
         call die("double: 1D+>4D failed (m)")

  end subroutine aa_1d_4d_dp
  
  ! Adds a 2D array to a 1D array
  subroutine aa_2d_1d_ip(from_2D, to_2D, in2D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_2D(2), to_2D(2)
    integer, intent(in) :: in2D(:,:)
    integer, intent(in) :: from_1D, to_1D
    integer, intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i2D, j2D, i1D

    i1D = from_1D
    do j2D = from_2D(2), to_2D(2)
       do i2D = from_2D(1), to_2D(1)
          out1D(i1D) = out1D(i1D) + in2D(i2D, j2D)
          i1D = i1D + 1
       end do
    end do

    if ( i1D <= to_1D ) &
         call die("integer: 2D+>1D failed")

  end subroutine aa_2d_1d_ip
  subroutine aa_2d_1d_sp(from_2D, to_2D, in2D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_2D(2), to_2D(2)
    real(sp), intent(in) :: in2D(:,:)
    integer, intent(in) :: from_1D, to_1D
    real(sp), intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i2D, j2D, i1D

    i1D = from_1D
    do j2D = from_2D(2), to_2D(2)
       do i2D = from_2D(1), to_2D(1)
          out1D(i1D) = out1D(i1D) + in2D(i2D, j2D)
          i1D = i1D + 1
       end do
    end do

    if ( i1D <= to_1D ) &
         call die("real: 2D+>1D failed")

  end subroutine aa_2d_1d_sp
  subroutine aa_2d_1d_dp(from_2D, to_2D, in2D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_2D(2), to_2D(2)
    real(dp), intent(in) :: in2D(:,:)
    integer, intent(in) :: from_1D, to_1D
    real(dp), intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i2D, j2D, i1D

    i1D = from_1D
    do j2D = from_2D(2), to_2D(2)
       do i2D = from_2D(1), to_2D(1)
          out1D(i1D) = out1D(i1D) + in2D(i2D, j2D)
          i1D = i1D + 1
       end do
    end do

    if ( i1D <= to_1D ) &
         call die("double: 2D+>1D failed")

  end subroutine aa_2d_1d_dp
  
  ! Adds a 3D array to a 1D array
  subroutine aa_3d_1d_ip(from_3D, to_3D, in3D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_3D(3), to_3D(3)
    integer, intent(in) :: in3D(:,:,:)
    integer, intent(in) :: from_1D, to_1D
    integer, intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i3D, j3D, k3D, i1D

    i1D = from_1D
    do k3D = from_3D(3), to_3D(3)
       do j3D = from_3D(2), to_3D(2)
          do i3D = from_3D(1), to_3D(1)
             out1D(i1D) = out1D(i1D) + in3D(i3D, j3D, k3D)
             i1D = i1D + 1
          end do
       end do
    end do

    if ( i1D <= to_1D ) &
         call die("integer: 3D+>1D failed")

  end subroutine aa_3d_1d_ip
  subroutine aa_3d_1d_sp(from_3D, to_3D, in3D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_3D(3), to_3D(3)
    real(sp), intent(in) :: in3D(:,:,:)
    integer, intent(in) :: from_1D, to_1D
    real(sp), intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i3D, j3D, k3D, i1D

    i1D = from_1D
    do k3D = from_3D(3), to_3D(3)
       do j3D = from_3D(2), to_3D(2)
          do i3D = from_3D(1), to_3D(1)
             out1D(i1D) = out1D(i1D) + in3D(i3D, j3D, k3D)
             i1D = i1D + 1
          end do
       end do
    end do

    if ( i1D <= to_1D ) &
         call die("real: 3D+>1D failed")

  end subroutine aa_3d_1d_sp
  subroutine aa_3d_1d_dp(from_3D, to_3D, in3D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_3D(3), to_3D(3)
    real(dp), intent(in) :: in3D(:,:,:)
    integer, intent(in) :: from_1D, to_1D
    real(dp), intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i3D, j3D, k3D, i1D

    i1D = from_1D
    do k3D = from_3D(3), to_3D(3)
       do j3D = from_3D(2), to_3D(2)
          do i3D = from_3D(1), to_3D(1)
             out1D(i1D) = out1D(i1D) + in3D(i3D, j3D, k3D)
             i1D = i1D + 1
          end do
       end do
    end do

    if ( i1D <= to_1D ) &
         call die("double: 3D+>1D failed")

  end subroutine aa_3d_1d_dp

  ! Adds a 4D array to a 1D array
  subroutine aa_4d_1d_ip(from_4D, to_4D, in4D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_4D(4), to_4D(4)
    integer, intent(in) :: in4D(:,:,:,:)
    integer, intent(in) :: from_1D, to_1D
    integer, intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i4D, j4D, k4D, m4D, i1D

    i1D = from_1D
    do m4D = from_4D(4), to_4D(4)
       do k4D = from_4D(3), to_4D(3)
          do j4D = from_4D(2), to_4D(2)
             do i4D = from_4D(1), to_4D(1)
                out1D(i1D) = out1D(i1D) + in4D(i4D, j4D, k4D, m4D)
                i1D = i1D + 1
             end do
          end do
       end do
    end do
    
    if ( i1D <= to_1D ) &
         call die("integer: 4D+>1D failed")

  end subroutine aa_4d_1d_ip
  subroutine aa_4d_1d_sp(from_4D, to_4D, in4D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_4D(4), to_4D(4)
    real(sp), intent(in) :: in4D(:,:,:,:)
    integer, intent(in) :: from_1D, to_1D
    real(sp), intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i4D, j4D, k4D, m4D, i1D

    i1D = from_1D
    do m4D = from_4D(4), to_4D(4)
       do k4D = from_4D(3), to_4D(3)
          do j4D = from_4D(2), to_4D(2)
             do i4D = from_4D(1), to_4D(1)
                out1D(i1D) = out1D(i1D) + in4D(i4D, j4D, k4D, m4D)
                i1D = i1D + 1
             end do
          end do
       end do
    end do
    
    if ( i1D <= to_1D ) &
         call die("real: 4D+>1D failed")

  end subroutine aa_4d_1d_sp
  subroutine aa_4d_1d_dp(from_4D, to_4D, in4D, from_1D, to_1D, out1D)
    integer, intent(in) :: from_4D(4), to_4D(4)
    real(dp), intent(in) :: in4D(:,:,:,:)
    integer, intent(in) :: from_1D, to_1D
    real(dp), intent(inout) :: out1D(:)

    ! Local variables for copying
    integer :: i4D, j4D, k4D, m4D, i1D

    i1D = from_1D
    do m4D = from_4D(4), to_4D(4)
       do k4D = from_4D(3), to_4D(3)
          do j4D = from_4D(2), to_4D(2)
             do i4D = from_4D(1), to_4D(1)
                out1D(i1D) = out1D(i1D) + in4D(i4D, j4D, k4D, m4D)
                i1D = i1D + 1
             end do
          end do
       end do
    end do

    if ( i1D <= to_1D ) &
         call die("double: 4D+>1D failed")

  end subroutine aa_4d_1d_dp

end module m_array
    
