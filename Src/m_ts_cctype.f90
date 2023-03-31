! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
module m_ts_cctype
! Module for containing the complex contour (will also handle the "close"
! to the real axis transport contour).
! The type can easily be streamlined and is containing all relevant information
! in the type.
!
! The important information available is:
!   1. Complex contour point
!   2. Weight of point (for transport contour this is 0.0)
!   3. Which part of the contour
!   4. The type of contour (i.e. which method is used to create it)
! 

  use m_ts_io_ctype, only : ts_c_io, CC_METHOD_LEN => c_N
  use precision, only : dp

  use m_gauss_fermi_inf, only : G_NF_MIN_kT, G_NF_MAX_kT

  implicit none

  private :: dp
  public

  ! Create a type to contain the contour information
  type :: ts_cw
     type(ts_c_io), pointer   :: c_io => null()
     integer,     allocatable :: ID(:)
     complex(dp), allocatable :: c(:), w(:,:)
  end type ts_cw

  type :: ts_c_idx
     sequence
     logical     :: exist = .false.
     logical     :: fake  = .false.
     complex(dp) :: e ! the energy for the curve
     ! contains:
     ! (1) : whether this is eq,noneq
     ! (2) : index of the io contour
     ! (3) : index of the point in the io(2)'th contour-array
     integer     :: idx(3) 
  end type ts_c_idx
  
!  ! maximum length of the string that returns the type
!  integer, parameter :: CC_METHOD_LEN = 17

  ! The following Fermi Gauss-Quadratures MUST be in success
  ! I.e. NF_<x>kT == 10+x
  integer, parameter :: CC_G_NF_MIN         = 4000 ! means G_NF_MIN_kT kT
  integer, parameter :: CC_G_NF_MAX         = G_NF_MAX_kT - G_NF_MIN_kT + CC_G_NF_MIN ! means G_NF_MAX_kT kT
  integer, parameter :: CC_G_NF_0kT         = CC_G_NF_MIN - G_NF_MIN_kT ! means 0 kT
  integer, parameter :: CC_G_LEGENDRE       = 100
  integer, parameter :: CC_TANH_SINH        = 101
  integer, parameter :: CC_SIMP_MIX         = 102
  integer, parameter :: CC_BOOLE_MIX        = 103
  integer, parameter :: CC_MID              = 104
  integer, parameter :: CC_CONTINUED_FRAC   = 105
  integer, parameter :: CC_USER             = 106


  ! Converts a method to a string format
  interface method2str
     module procedure method2str_int
     module procedure method2str_ts_c_io
  end interface method2str
  private :: method2str_int, method2str_ts_c_io

  ! Converts a method to a string format
  interface longmethod2str
     module procedure longmethod2str_int
     module procedure longmethod2str_ts_c_io
  end interface longmethod2str
  private :: longmethod2str_int, longmethod2str_ts_c_io

  ! Converts a method to the equivalent string format that is required as input
  interface method2input
     module procedure method2input_int
     module procedure method2input_ts_c_io
  end interface method2input
  private :: method2input_int, method2input_ts_c_io

  ! converts a string to the integer method
  interface method
     module procedure method_str
     module procedure method_ts_c_io
  end interface method
  private :: method_str, method_ts_c_io

contains

  function method2str_ts_c_io(c) result(str)
    type(ts_c_io), intent(in) :: c
    character(len=CC_METHOD_LEN) :: str
    str = method2str(method(c))
  end function method2str_ts_c_io
  
  function method2str_int(method) result(str)
    integer, intent(in) :: method
    character(len=CC_METHOD_LEN) :: str
    select case ( method )
    case ( CC_G_NF_MIN:CC_G_NF_MAX )
       write(str,'(a,i0)') 'Gauss-Fermi_',method-CC_G_NF_0kT
    case ( CC_G_LEGENDRE )
       str = 'Gauss-Legendre'
    case ( CC_TANH_SINH )
       str = 'Tanh-Sinh'
    case ( CC_SIMP_MIX )
       str = 'Simpson 3/8-3'
    case ( CC_BOOLE_MIX )
       str = 'Boole-Simpson 3/8'
    case ( CC_MID )
       str = 'Mid-rule'
    case ( CC_CONTINUED_FRAC )
       str = 'Continued-fraction'
    case ( CC_USER )
       str = 'User'
    case default
       call die('Unknown method for the contour')
    end select
  end function method2str_int

  function method2input_int(method) result(str)
    integer, intent(in) :: method
    character(len=CC_METHOD_LEN) :: str
    select case ( method )
    case ( CC_G_NF_MIN:CC_G_NF_MAX )
       str = 'gauss-fermi'
    case ( CC_G_LEGENDRE )
       str = 'gauss-Legendre'
    case ( CC_TANH_SINH )
       str = 'Tanh-Sinh'
    case ( CC_SIMP_MIX )
       str = 'Simpson-mix'
    case ( CC_BOOLE_MIX )
       str = 'Boole-mix'
    case ( CC_MID )
       str = 'Mid-rule'
    case ( CC_CONTINUED_FRAC )
       str = 'Continued-fraction'
    case ( CC_USER )
       str = 'User defined'
    case default
       call die('Unknown method for the contour')
    end select
  end function method2input_int

  function method2input_ts_c_io(c) result(str)
    type(ts_c_io), intent(in) :: c
    character(len=CC_METHOD_LEN) :: str
    str = method2input(method(c))
  end function method2input_ts_c_io

  function longmethod2str_ts_c_io(c) result(str)
    type(ts_c_io), intent(in) :: c
    character(len=CC_METHOD_LEN*2) :: str
    str = longmethod2str(method(c))
  end function longmethod2str_ts_c_io

  function longmethod2str_int(method) result(str)
    integer, intent(in) :: method
    character(len=CC_METHOD_LEN*2) :: str
    select case ( method ) 
    case ( CC_G_NF_MIN:CC_G_NF_MAX )
       write(str,'(a,'' ('',i0,''kT)'')') 'Gauss-Fermi',method-CC_G_NF_0kT
    case ( CC_G_LEGENDRE )
       str = 'Gauss-Legendre'
    case ( CC_TANH_SINH )
       str = 'Tanh-Sinh'
    case ( CC_SIMP_MIX )
       str = 'Simpson-mix'
    case ( CC_BOOLE_MIX )
       str = 'Boole-mix'
    case ( CC_MID )
       str = 'Mid-rule'
    case ( CC_CONTINUED_FRAC )
       str = 'Continued fraction'
    case ( CC_USER )
       str = 'User defined'
    case default
       call die('Unknown method for the contour')
    end select
  end function longmethod2str_int

  function method_str(str) result(method)
    use fdf, only : leqi
    character(len=*), intent(in) :: str
    integer :: method
    character(len=20) :: tmp
    integer :: i
    if ( leqi(str,'g-legendre') .or. leqi(str, 'gauss-legendre') ) then
       method = CC_G_LEGENDRE
    else if ( leqi(str,'tanh-sinh') ) then
       method = CC_TANH_SINH
    else if ( leqi(str,'simpson-mix') .or. &
         leqi(str,'simpson') ) then
       method = CC_SIMP_MIX
    else if ( leqi(str,'boole-mix') .or. &
         leqi(str,'boole') ) then
       method = CC_BOOLE_MIX
    else if ( leqi(str,'mid-rule') .or. &
         leqi(str,'mid') ) then
       method = CC_MID
    else if ( leqi(str,'ozaki') .or. &
         leqi(str,'continued-fraction') .or. &
         leqi(str,'cont-frac') ) then
       method = CC_CONTINUED_FRAC
    else if ( leqi(str,'file') .or. &
         leqi(str,'user') ) then
       method = CC_USER
    else if ( leqi(str,'g-fermi') .or. leqi(str, 'gauss-fermi') ) then
       method = CC_G_NF_0kT
       do i = G_NF_MIN_kT , G_NF_MAX_kT
          write(tmp,'(a,i0,a)') 'g-fermi(',i,')'
          if ( leqi(str,tmp) ) then
             method = CC_G_NF_0kT + i
             return
          end if
       end do
    else
       call die('Unknown method for the contour: '//trim(str))
    end if
  end function method_str

  function method_ts_c_io(c) result(ret)
    type(ts_c_io), intent(in) :: c
    integer :: ret
    ret = method(c%method)
  end function method_ts_c_io

end module m_ts_cctype
