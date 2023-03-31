! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!!@LICENSE
c
c     MODULE m_io
c
c Copyright Alberto Garcia, 1996, 1997, 1998
c
c This module implements an interface to the FORTRAN logical unit
c system. Based on code by Richard Maine.
c
c Alberto Garcia, December 30, 1996
c Rewritten as a single subroutine 
c with multiple entry points, March 7, 1998
c Now hybrid to comply with Siesta "die" interface.
c Converted to a module by J.M.Soler. Aug. 2009
c---------------------------------------------------------------
c
      MODULE m_io
c
c-----------------------------------------------------------------
c
c     Used module procedures
c
      USE sys, only: die   ! Termination routine

      implicit none
c
c-----------------------------------------------------------------
c
c     Public procedures provided by this module
c
      PUBLIC::
     $  io_seterr,  ! Set standard error unit
     $  io_setout,  ! Set standard output unit
     $  io_geterr,  ! Get standard error unit
     $  io_getout,  ! Get standard output unit
     $  io_assign,  ! Get some available IO unit and reserve it
     $  io_reserve, ! Reserve a specific IO unit
     $  io_close,   ! Close and free a given IO unit
     $  io_status   ! Print all used IO units

      PRIVATE ! Nothing is declared public below this point
c
c----------------------------------------------------------------
c
c     Module variables
c
c     Logical unit management. Units 0 to min_lun-1 are "reserved",
c     since most of the "typical" files (output, etc) use them.
c
c     Logical units min_lun to min_max are managed by this module.
c
      integer, parameter:: min_lun = 10
      integer, parameter:: max_lun = 99
      integer, parameter:: nunits = max_lun-min_lun+1
      integer, save:: stdout = 6
      integer, save:: stderr = 0
      logical, save:: lun_is_free(min_lun:max_lun) = .true.

      CONTAINS
c
c-----------------------------------------------------------------
c
c     Simple interfaces to modify standard units
c
      subroutine io_seterr(unit)
      integer,intent(in):: unit
      stderr = unit
      end subroutine io_seterr
c
c-----------------------------------------------------------------
c
      subroutine io_setout(unit)
      integer,intent(in):: unit
      stdout = unit
      end subroutine io_setout
c
c-----------------------------------------------------------------
c
      subroutine io_geterr(unit)
      integer,intent(out):: unit
      unit = stderr
      end subroutine io_geterr
c
c-----------------------------------------------------------------
c
      subroutine io_getout(unit)
      integer,intent(out):: unit
      unit = stdout
      end subroutine io_getout
c
c------------------------------------------------------------------     
c
c     Logical unit management
c
      subroutine io_assign(lun)
      integer,intent(out):: lun
      integer :: iostat
      logical :: used
c
c     Looks for a free unit and assigns it to lun
c
      do lun= min_lun, max_lun
         if (lun_is_free(lun)) then
            inquire(unit=lun, opened=used, iostat=iostat)
            if (iostat .ne. 0) used = .true.
            lun_is_free(lun) = .false.
            if (.not. used) return
         endif
      enddo
      call die('No luns available in io_assign')

      end subroutine io_assign
c
c------------------------------------------------------------------     
c
      subroutine io_reserve(lun)
      integer,intent(in):: lun
      integer :: iostat
      logical :: used
c
c     Useful to specify that one needs to use a particular unit number
c
c     For example, assume some legacy code expects to work with unit 15:
c
c     call io_reserve(15)   ! this call at the beginning of the program
c     ...
c     open(15,....)
c
      inquire(unit=lun, opened=used, iostat=iostat)
      if (iostat .ne. 0) used = .true.
      if (used) call die('Cannot reserve unit. Already connected')
      if (lun .ge. min_lun .and. lun .le. max_lun)
     $                      lun_is_free(lun) = .false.

      end subroutine io_reserve
c
c------------------------------------------------------------------     
c
      subroutine io_close(lun)
      integer,intent(in):: lun
c
c     Use this routine instead of a simple close!!
c
      close(lun)
      if (lun .ge. min_lun .and. lun .le. max_lun)
     $                     lun_is_free(lun) = .true.

      end subroutine io_close
c
c------------------------------------------------------------------     
c
      subroutine io_status()
c
c     Prints a list of the connected logical units and the names of
c     the associated files
c
      character(len=128) :: filename
      character(len=32) :: form
      logical :: opened, named
      integer :: i, iostat

      write(stdout,'(a)') '******** io_status ********'
      do i = 0, max_lun
         inquire(i,opened=opened,named=named,name=filename,
     $           form=form,iostat=iostat)
         if (iostat .eq. 0) then
            if (opened) then
               if (named) then
                  write(stdout,9000) i, form, filename
               else
                  write(stdout,9000) i, form, 'No name available'
               endif
            endif
         else
            write(stdout,9000) i, 'Iostat error'
         endif
      enddo
      write(stdout,'(a)') '********           ********'

 9000 format(i4,5x,a,5x,a)
      end subroutine io_status

      END MODULE m_io

