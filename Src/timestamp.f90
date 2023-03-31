! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_timestamp

  implicit none

  private

  public :: datestring
  public :: timestamp
  private :: iso_date

CONTAINS
  
!--------------------------------------------------
  function datestring() result (strdate)
  character(len=19) :: strdate

    integer :: values(8)

    call date_and_time(values=values)
    strdate=iso_date(values)

  end function datestring

!--------------------------------------------------
  subroutine timestamp (str)
    
    implicit none
    
    character(len=*), intent(in) :: str
    
    integer         :: sec, min, hour, day, month, year
    character(len=1), parameter :: dash = "-"
    character(len=1), parameter :: colon = ":"
    character(len=3), parameter :: prefix = ">> "
    character(len=3)  :: month_str(12) =     &
      (/'JAN','FEB','MAR','APR','MAY','JUN', &
        'JUL','AUG','SEP','OCT','NOV','DEC'/)

    integer :: values(8)

    call date_and_time(values=values)
    year = values(1)
    month = values(2)
    day = values(3)
    hour = values(5)
    min = values(6)
    sec = values(7)
    
    write(6,1000) prefix, trim(str), colon,                &
                  day, dash, month_str(month), dash, year, &
                  hour, colon, min, colon, sec
    
1000 format(2a,a1,2x,i2,a1,a3,a1,i4,2x,i2,a1,i2.2,a1,i2.2)
    
  end subroutine timestamp

!--------------------------------------------------
  function iso_date(values) result(fdate)
    character(len=19) :: fdate
    integer, intent(in) :: values(8)

    !Turn the output of the date_and_time intrinsic
    !into an ISO8601-compliant representation.
    
    write(fdate, &
      "(i4.4,'-',i2.2,'-',i2.2,'T',i2.2,'-',i2.2,'-',i2.2)") &
      values(1:3), values(5:7)
  end function iso_date
  
  
end module m_timestamp
