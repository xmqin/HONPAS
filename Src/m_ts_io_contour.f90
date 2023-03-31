!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! This code segment has been fully created by:
! Nick Papior Andersen, 2012, nickpapior@gmail.com
!

module m_ts_io_contour
!
! Routines that are used to read in and print out the contour for integration of the GFs
! 
  use precision, only : dp

  implicit none

  character(len=*), parameter :: OPT_N   = '(''ts: '',a)'
  character(len=*), parameter :: OPT_C   = '(''ts: '',a,t53,''=    '',a)'
  character(len=*), parameter :: OPT_CC  = '(''ts: '',a,t53,''=    '',a,tr2,a)'
  character(len=*), parameter :: OPT_F   = '(''ts: '',a,t53,''='',f10.4)'
  character(len=*), parameter :: OPT_INT = '(''ts: '',a,t53,''='',i5)'
  character(len=*), parameter :: OPT_F_U = '(''ts: '',a,t53,''='',f10.4,tr1,a)'
  character(len=*), parameter :: OPT_G_U = '(''ts: '',a,t53,''='',g11.4,tr1,a)'

contains

  subroutine write_e(str,val,unit)

    use units, only : eV, Kelvin

    character(len=*), intent(in) :: str
    real(dp), intent(in) :: val
    character(len=*), intent(in), optional :: unit

    ! This is our definition of infinity....
    if ( abs(val) > 10000._dp ) then
       if ( val < 0._dp ) then
          write(*,opt_c) trim(str),'-Infinity'
       else
          write(*,opt_c) trim(str),' Infinity'
       end if
    else if ( .not. present(unit) ) then
       write(*,opt_f_u) trim(str),val / eV,'eV'
    else
       if ( unit == 'eV' ) then
          write(*,opt_f_u) trim(str),val / eV,'eV'
       else if ( unit == 'Ry' ) then
          write(*,opt_f_u) trim(str),val,'Ry'
       else if ( unit == 'K' ) then
          write(*,opt_f_u) trim(str),val/Kelvin,'K'
       else
          call die('Programming error, unknown unit')
       end if
    end if

  end subroutine write_e

end module m_ts_io_contour
