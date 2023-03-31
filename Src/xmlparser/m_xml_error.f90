module m_xml_error
!
! Error handling
!
use m_elstack
private

type, public  :: xml_error_t
      character(len=100)  :: message
      integer             :: line
      integer             :: column
      type(elstack_t)     :: stack
      integer             :: severity
end type xml_error_t

integer, public             ::  xml_stderr = 0    ! Unit for error info
integer, public, parameter  ::  SEVERE_ERROR_CODE=0, WARNING_CODE=1

public  :: build_error_info, default_error_handler, general_error
public  :: set_xml_stderr

CONTAINS

!-------------------------------------------------------------------------
subroutine build_error_info(error_info,message,line,column,stack,severity)
type(xml_error_t), intent(out)        :: error_info
integer, intent(in)                   :: line, column
character(len=*), intent(in)          :: message
type(elstack_t), intent(in)           :: stack
integer, intent(in)                   :: severity

error_info%message = message
error_info%line = line
error_info%column = column
error_info%stack = stack
error_info%severity = severity

end subroutine build_error_info

!--------------------------------------------------

subroutine default_error_handler(error_info)
type(xml_error_t), intent(in)            :: error_info
!
! Default error handling
!
if (error_info%severity == SEVERE_ERROR_CODE)  then
   write(unit=xml_stderr,fmt="(a)") "*** XML parsing Error:"
else if (error_info%severity == WARNING_CODE)  then
   write(unit=xml_stderr,fmt="(a)") "*** XML parsing Warning:"
endif
write(unit=xml_stderr,fmt="(a)") trim(error_info%message)
write(unit=xml_stderr,fmt="(a,i8,a,i4)") "Line: ", &
                                         error_info%line, &
                                         " Column: ", &
                                         error_info%column
write(unit=xml_stderr,fmt="(a)") "Element traceback:"
call print_elstack(error_info%stack,unit=xml_stderr)
!
!   If there is a severe error the program should stop...
!
if (error_info%severity == SEVERE_ERROR_CODE)  then
      STOP
else if (error_info%severity == WARNING_CODE)  then
   write(unit=xml_stderr,fmt="(a)") "*** Continuing after Warning..."
endif

end subroutine default_error_handler

subroutine general_error(msg,code)
character(len=*), intent(in)     :: msg
integer, intent(in)              :: code
!

if (code == SEVERE_ERROR_CODE)  then
   write(unit=xml_stderr,fmt="(2a)") "** Error: ", msg
   STOP
else if (code == WARNING_CODE)  then
   write(unit=xml_stderr,fmt="(2a)") "** Warning: ", msg
else
   stop "wrong error code"
endif

end subroutine general_error

!-------------------------------------------------------------------------
subroutine set_xml_stderr(unit)
integer, intent(in)  :: unit

xml_stderr  = unit

end subroutine set_xml_stderr

end module m_xml_error






