module m_wxml_core

use m_wxml_buffer
use m_wxml_array_str, only: assign_array_to_str
use m_wxml_array_str, only: assign_str_to_array
use m_wxml_escape, only: check_Name
use m_wxml_elstack
use m_wxml_dictionary

implicit none

logical, private, save  :: pcdata_advance_line_default = .false.
logical, private, save  :: pcdata_advance_space_default = .false.

integer, private, parameter ::  sp = selected_real_kind(6,30)
integer, private, parameter ::  dp = selected_real_kind(14,100)

private

type, public :: xmlf_t
   character, pointer      :: filename(:)
   integer                 :: lun
   type(buffer_t)          :: buffer
   type(elstack_t)         :: stack
   type(wxml_dictionary_t) :: dict
   logical                 :: start_tag_closed
   logical                 :: root_element_output
   logical                 :: indenting_requested
end type xmlf_t

public :: xml_OpenFile, xml_NewElement, xml_EndElement, xml_Close
public :: xml_AddXMLDeclaration
public :: xml_AddXMLStylesheet
public :: xml_AddXMLPI
public :: xml_AddComment, xml_AddCdataSection
public :: xml_AddPcdata, xml_AddAttribute
interface  xml_AddPcdata
   module procedure xml_AddPcdata_Ch
end interface
!
interface  xml_AddAttribute
   module procedure xml_AddAttribute_Ch
end interface
!
public :: xml_AddArray
interface xml_AddArray
   module procedure  xml_AddArray_integer,  &
                xml_AddArray_real_dp, xml_AddArray_real_sp
end interface
private :: xml_AddArray_integer,  xml_AddArray_real_dp, xml_AddArray_real_sp

private :: get_unit
private :: add_eol
private :: write_attributes

!overload error handlers to allow file info
interface wxml_warning
  module procedure wxml_warning_xf
end interface
interface wxml_error
  module procedure wxml_error_xf
end interface
interface wxml_fatal
  module procedure wxml_fatal_xf
end interface

!
! Heuristic (approximate) target for justification of output
! Large unbroken pcdatas will go beyond this limit
! 
integer, private, parameter  :: COLUMNS = 80

! TOHW - This is the longest string that may be output without
! a newline. The buffer must not be larger than this, but its size 
! can be tuned for performance.
integer, private, parameter  :: xml_recl = 4096

CONTAINS

!-------------------------------------------------------------------
subroutine xml_OpenFile(filename, xf, indent)
character(len=*), intent(in)  :: filename
type(xmlf_t), intent(inout)   :: xf
logical, intent(in), optional :: indent

integer :: iostat

allocate(xf%filename(len(filename)))
call assign_str_to_array(xf%filename,filename)

call get_unit(xf%lun,iostat)
if (iostat /= 0) call wxml_fatal(xf, "cannot open file")
!
! Use large I/O buffer in case the O.S./Compiler combination
! has hard-limits by default (i.e., NAGWare f95's 1024 byte limit)
! This is related to the maximum size of the buffer.
! TOHW - This is the longest string that may be output without
! a newline. The buffer must not be larger than this, but its size 
! can be tuned for performance.

open(unit=xf%lun, file=filename, form="formatted", status="replace", &
     action="write", position="rewind", recl=xml_recl)

call init_elstack(xf%stack)

call init_dict(xf%dict)
call reset_buffer(xf%buffer)

xf%start_tag_closed = .true.
xf%root_element_output = .false.

xf%indenting_requested = .false.
if (present(indent)) then
   xf%indenting_requested = indent
endif
end subroutine xml_OpenFile

!-------------------------------------------------------------------
subroutine xml_AddXMLDeclaration(xf,encoding)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in), optional :: encoding

if (present(encoding)) then
   call add_to_buffer("<?xml version=""1.0"" encoding=""" &
                     // trim(encoding) // """ ?>", xf%buffer)
else
   call add_to_buffer("<?xml version=""1.0"" ?>", xf%buffer)
endif
end subroutine xml_AddXMLDeclaration

!-------------------------------------------------------------------
subroutine xml_AddXMLStylesheet(xf, href, type, title, media, charset, alternate)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in) :: href
character(len=*), intent(in) :: type
character(len=*), intent(in), optional :: title
character(len=*), intent(in), optional :: media
character(len=*), intent(in), optional :: charset
logical,          intent(in), optional :: alternate

call add_eol(xf)
call add_to_buffer("<?xml-stylesheet href=""" //trim(href)// &
     """ type=""" //trim(type)// """", xf%buffer)

if (present(title)) call add_to_buffer(" title="""//trim(title)// """", xf%buffer)
if (present(media)) call add_to_buffer(" media="""//trim(media)// """", xf%buffer)
if (present(charset)) call add_to_buffer(" charset="""//trim(charset)// """", xf%buffer)
if (present(alternate)) then
   if (alternate) then
      call add_to_buffer(" alternate=""yes""", xf%buffer)
   else
      call add_to_buffer(" alternate=""no""", xf%buffer)
   endif
endif
call add_to_buffer(" ?>", xf%buffer)

end subroutine xml_AddXMLStylesheet

!-------------------------------------------------------------------
subroutine xml_AddXMLPI(xf, name, data)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in) :: name
character(len=*), intent(in), optional :: data

call add_eol(xf)
call add_to_buffer("<?" // trim(name) // " ", xf%buffer)
if(present(data)) call add_to_buffer(data, xf%buffer)
call add_to_buffer(" ?>", xf%buffer)

end subroutine xml_AddXMLPI


!-------------------------------------------------------------------
subroutine xml_AddComment(xf,comment)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in)  :: comment

call close_start_tag(xf,">")
call add_eol(xf)
call add_to_buffer("<!--", xf%buffer)
call add_to_buffer(comment, xf%buffer)
call add_to_buffer("-->", xf%buffer)
end subroutine xml_AddComment

!-------------------------------------------------------------------
subroutine xml_AddCdataSection(xf,cdata)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in)  :: cdata

call close_start_tag(xf,">")
call add_to_buffer("<![CDATA[", xf%buffer)
call add_to_buffer(cdata, xf%buffer)
call add_to_buffer("]]>", xf%buffer)
end subroutine xml_AddCdataSection

!-------------------------------------------------------------------
subroutine xml_NewElement(xf,name)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in)  :: name

if (is_empty(xf%stack)) then
   if (xf%root_element_output) call wxml_error(xf, "two root elements")
   xf%root_element_output = .true.
endif

if (.not.check_Name(name)) then
  call wxml_warning(xf, 'attribute name '//name//' is not valid')
endif


call close_start_tag(xf,">")
call push_elstack(name,xf%stack)
call add_eol(xf)
call add_to_buffer("<" // trim(name),xf%buffer)
xf%start_tag_closed = .false.
call reset_dict(xf%dict)

end subroutine xml_NewElement
!-------------------------------------------------------------------
subroutine xml_AddPcdata_Ch(xf,pcdata,space,line_feed)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in)  :: pcdata
logical, intent(in), optional  :: space
logical, intent(in), optional  :: line_feed

logical :: advance_line , advance_space

advance_line = pcdata_advance_line_default 
if (present(line_feed)) then
   advance_line = line_feed
endif

advance_space = pcdata_advance_space_default 
if (present(space)) then
   advance_space = space
endif

if (is_empty(xf%stack)) then
   call wxml_error(xf, "pcdata outside element content")
endif

call close_start_tag(xf,">")

if (advance_line) then
   call add_eol(xf)
   advance_space = .false.
else
   if (xf%indenting_requested) then
      if ((len(xf%buffer) + len_trim(pcdata) + 1) > COLUMNS ) then
         call add_eol(xf)
         advance_space = .false.
      endif
   endif
endif
if (advance_space) call add_to_buffer(" ",xf%buffer)

call add_to_buffer_escaping_markup(pcdata, xf%buffer)

end subroutine xml_AddPcdata_Ch

!-------------------------------------------------------------------
subroutine xml_AddAttribute_Ch(xf,name,value)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in)  :: name
character(len=*), intent(in)  :: value

if (is_empty(xf%stack)) then
   call wxml_error(xf, "attributes outside element content")
endif

if (xf%start_tag_closed)  then
   call wxml_error(xf, "attributes outside start tag")
endif

if (has_key(xf%dict,name)) then
   call wxml_error(xf, "duplicate att name")
endif

call add_item_to_dict(xf%dict, name, value)

end subroutine xml_AddAttribute_Ch

!-----------------------------------------------------------
subroutine xml_EndElement(xf,name)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in)  :: name

character(len=2000)  :: current

if (is_empty(xf%stack)) then
   call wxml_fatal(xf, "Out of elements to close")
endif

call get_top_elstack(xf%stack,current)
if (current /= name) then
   call wxml_fatal(xf, 'Trying to close '//Trim(name)//' but '//Trim(current)//' is open.') 
endif
if (.not. xf%start_tag_closed)  then                ! Empty element
   if (len(xf%dict) > 0) call write_attributes(xf)
   call add_to_buffer(" />",xf%buffer)
   xf%start_tag_closed = .true.
else
   call add_eol(xf)
   call add_to_buffer("</" // trim(name) // ">", xf%buffer)
endif
call pop_elstack(xf%stack,current)

end subroutine xml_EndElement

!----------------------------------------------------------------

subroutine xml_Close(xf)
type(xmlf_t), intent(inout)   :: xf

character(len=2000) :: name

do
   if (is_empty(xf%stack)) exit
   call get_top_elstack(xf%stack,name)
   call xml_EndElement(xf,trim(name))
enddo

write(unit=xf%lun,fmt="(a)") char(xf%buffer)
close(unit=xf%lun)

deallocate(xf%filename)

end subroutine xml_Close

!==================================================================
!-------------------------------------------------------------------
subroutine get_unit(lun,iostat)

! Get an available Fortran unit number

integer, intent(out)  :: lun
integer, intent(out)  :: iostat

integer :: i
logical :: unit_used

do i = 10, 99
   lun = i
   inquire(unit=lun,opened=unit_used)
   if (.not. unit_used) then
      iostat = 0
      return
   endif
enddo
iostat = -1
lun = -1
end subroutine get_unit

!----------------------------------------------------------
subroutine add_eol(xf)
type(xmlf_t), intent(inout)   :: xf

integer :: indent_level

! In case we still have a zero-length stack, we must make
! sure indent_level is not less than zero.
indent_level = max(len(xf%stack) - 1, 0)

!We must flush here (rather than just adding an eol character)
!since we don't know what the eol character is on this system.
!Flushing with a linefeed will get it automatically, though.
write(unit=xf%lun,fmt="(a)") char(xf%buffer)
call reset_buffer(xf%buffer)

if (xf%indenting_requested) &
   call add_to_buffer(repeat(' ',indent_level),xf%buffer)

end subroutine add_eol
!------------------------------------------------------------
subroutine dump_buffer(xf,lf)
type(xmlf_t), intent(inout)   :: xf
logical, intent(in), optional :: lf

if (present(lf)) then
   if (lf) then
      write(unit=xf%lun,fmt="(a)",advance="yes") char(xf%buffer)
   else
      write(unit=xf%lun,fmt="(a)",advance="no") char(xf%buffer)
   endif
else
   write(unit=xf%lun,fmt="(a)",advance="no") char(xf%buffer)
endif
call reset_buffer(xf%buffer)

end subroutine dump_buffer

!------------------------------------------------------------
subroutine close_start_tag(xf,s)
type(xmlf_t), intent(inout)   :: xf
character(len=*), intent(in)  :: s

if (.not. xf%start_tag_closed)  then
   if (len(xf%dict) > 0)  call write_attributes(xf)
   call add_to_buffer(s, xf%buffer)
   xf%start_tag_closed = .true.
endif

end subroutine close_start_tag

!-------------------------------------------------------------
subroutine write_attributes(xf)
type(xmlf_t), intent(inout)   :: xf

integer  :: i, status, size, key_len, value_len
character(len=2000)  :: key, value

do i = 1, len(xf%dict)
   call get_key(xf%dict,i,key,key_len,status)
   call get_value(xf%dict,i,value,value_len,status)
   size = key_len + value_len + 4
   if ((len(xf%buffer) + size) > COLUMNS) call add_eol(xf)
   call add_to_buffer(" ", xf%buffer)
   call add_to_buffer(key(:key_len), xf%buffer)
   call add_to_buffer("=", xf%buffer)
   call add_to_buffer("""",xf%buffer)
   call add_to_buffer_escaping_markup(value(:value_len), xf%buffer)
   call add_to_buffer("""", xf%buffer)
enddo

end subroutine write_attributes

!---------------------------------------------------------------
    subroutine xml_AddArray_integer(xf,a,format)
      type(xmlf_t), intent(inout)         :: xf
      integer, intent(in), dimension(:)   :: a
      character(len=*), intent(in), optional  :: format

      call close_start_tag(xf,">")
      if (len(xf%buffer) > 0) call dump_buffer(xf,lf=.true.)
      if (present(format)) then
         write(xf%lun,format) a
      else
         write(xf%lun,"(6(i12))") a
      endif
    end subroutine xml_AddArray_integer

!-------------------------------------------------------------------
    subroutine xml_AddArray_real_dp(xf,a,format)
      type(xmlf_t), intent(inout)         :: xf
      real(kind=dp), intent(in), dimension(:)   :: a
      character(len=*), intent(in), optional  :: format

      call close_start_tag(xf,">")
      if (len(xf%buffer) > 0) call dump_buffer(xf,lf=.true.)
      if (present(format)) then
         write(xf%lun,format) a
      else
         write(xf%lun,"(4(es20.12))") a
      endif
    end subroutine xml_AddArray_real_dp

!------------------------------------------------------------------
    subroutine xml_AddArray_real_sp(xf,a,format)
      type(xmlf_t), intent(inout)         :: xf
      real(kind=sp), intent(in), dimension(:)   :: a
      character(len=*), intent(in), optional  :: format

      call close_start_tag(xf,">")
      if (len(xf%buffer) > 0) call dump_buffer(xf,lf=.true.)
      if (present(format)) then
         write(xf%lun,format) a
      else
         write(xf%lun,"(4(es20.7))") a
      endif
    end subroutine xml_AddArray_real_sp

!---------------------------------------------------------
! Error handling/trapping routines:

    subroutine wxml_warning_xf(xf, msg)
      ! Emit warning, but carry on.
      type(xmlf_t), intent(in) :: xf
      character(len=*), intent(in) :: msg

      write(6,'(a)') 'WARNING(wxml) in writing to file ', xmlf_name(xf)
      write(6,'(a)')  msg

    end subroutine wxml_warning_xf

    subroutine wxml_error_xf(xf, msg)
      ! Emit error message, clean up file and stop.
      type(xmlf_t), intent(inout) :: xf
      character(len=*), intent(in) :: msg

      write(6,'(a)') 'ERROR(wxml) in writing to file ', xmlf_name(xf)
      write(6,'(a)')  msg

      call xml_Close(xf)
      stop

    end subroutine wxml_error_xf

    subroutine wxml_fatal_xf(xf, msg)
      !Emit error message and abort with coredump. Does not try to
      !close file, so should be used from anything xml_Close might
      !itself call (to avoid infinite recursion!)

      type(xmlf_t), intent(in) :: xf
      character(len=*), intent(in) :: msg

      write(6,'(a)') 'ERROR(wxml) in writing to file ', xmlf_name(xf)
      write(6,'(a)')  msg

      call pxfabort
      stop

    end subroutine wxml_fatal_xf

    function xmlf_name(xf) result(fn)
      Type (xmlf_t), intent(in) :: xf
      character(len=size(xf%filename)) :: fn
      call assign_array_to_str(fn,xf%filename)
    end function xmlf_name
      

end module m_wxml_core

