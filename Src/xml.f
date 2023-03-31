! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module xml
!
! A simple set of utilities for *writing* xml files
!

      use precision

      private
      public str, xml_dump_attribute, xml_dump_element

!-------------------------------
      interface str
        module procedure
     $     str_integer,
     $     str_double,
     $     str_string,
     $     str_logical
      end interface

      CONTAINS

      subroutine xml_dump_attribute(lun,name,value)
      integer, intent(in)  :: lun
      character(len=*), intent(in)  :: name
      character(len=*), intent(in)  :: value

      write(lun,*) trim(name),'="', trim(value), '" '
      end subroutine xml_dump_attribute

      subroutine xml_dump_element(lun,name,value)
      integer, intent(in)  :: lun
      character(len=*), intent(in)  :: name
      character(len=*), intent(in)  :: value

      write(lun,'(3a)',advance='no') '<',trim(name),'> '
      write(lun,'(a)',advance='no') trim(value)
      write(lun,'(3a)') ' </',trim(name),'> '
      end subroutine xml_dump_element
!--------------------------------------------------------
      function str_integer(int)
      integer, intent(in)  :: int
      character(len=25) str_integer

      character(len=25) dummy
      write(dummy,'(i25)')  int
      str_integer = trim(dummy)
      end function str_integer

      function str_double(num)
      real(dp), intent(in)  :: num
      character(len=25) str_double

      character(len=25) dummy
      write(dummy,'(g22.12)')  num
      str_double = trim(dummy)
      end function str_double

      function str_logical(v)
      logical, intent(in)  :: v
      character(len=5) str_logical

      if ( v ) then
         str_logical = 'true'
      else
         str_logical = 'false'
      end if
      end function str_logical

      function str_string(num)
      character(len=*), intent(in)  :: num
      character(len=50) str_string

      str_string = trim(num)
      end function str_string

      end module xml








