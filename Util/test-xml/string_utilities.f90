module string_utilities
use compare_tol_m, only: sp, dp
use m_strings

logical :: string_utilities_debug = .false.


public :: clean_string, only_numbers, compare_only_numbers, compare_alpha
private
contains

function clean_string(str) result(cstr)
  
  type(string),intent(in) :: str
  
  type(string) ::cstr
  
  if (len(str) == 0)then
     cstr=""
  else
     !Remove the new lines
     cstr = str

     if (new_lines(cstr)) then
        if (string_utilities_debug) print *,"         clean_string: Removing new lines!"
        call remove_new_lines(cstr)
     endif

     cstr = adjustl(cstr)
     cstr = trim(cstr)
  endif

end function clean_string

!---------------------------------------------------

function only_numbers(str)
  type(string),intent(in)::str
  logical :: only_numbers

  character(len=53)::letters
  type(string) :: c_str
  integer :: length,position

  !Remove leading/trailing white spaces
  !print *,"original:","|",char(str),"|"
  
  letters="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ*"
  !Copy, we don't want to modify the xml file.
  c_str = str

  only_numbers = .false.
  
  if (len(c_str) >0)then
     c_str = clean_string(c_str)
     if (string_utilities_debug) print *,"           only_numbers in string:|",char(c_str),"|?"
	
     length = len_trim(c_str)
     position = scan(char(c_str),letters)
	
     if (position == 0)then
        only_numbers = .true.
     else
        only_numbers = .false.
     endif
	
     if (string_utilities_debug) print *,"           only_numbers:",only_numbers
  endif
 
end function only_numbers
!---------------------------------------------------

function new_lines(str)
  use m_strings
  type(string),intent(in)::str
  logical :: new_lines
  
  integer :: position
  character :: nl
  
  nl = achar(10)

  if (len(str) > 0)then
     position = index(str,nl)
  else
     position = 0
  endif

  if (position == 0) then
     new_lines = .false.
  else
     new_lines = .true.
  endif

end function new_lines

!---------------------------------------------------
recursive subroutine remove_new_lines(str)
  use m_strings
  type(string),intent(inout)::str

  integer :: position
  character(len=1) :: nl
  
  nl = achar(10)
 
  do
     position = index(str,nl)
     if (position == 0)then
        exit
     else
        str = remove(str,position,position+1)
        call remove_new_lines(str)
     endif
  enddo

end subroutine remove_new_lines

!---------------------------------------------------
function compare_only_numbers(c_ref,c_mod,label,tol)
  type(string), intent(in) :: c_ref,c_mod
  character(len=*)         :: label
  real(dp), intent(in)         :: tol
  logical compare_only_numbers
  
  character :: ws
  integer   :: position
  real(dp)     :: r_ref,r_mod 
  type(string) :: ref,mod


  !Remove the new lines
  ref = c_ref
  mod = c_mod

  if (len(ref) > 0 .and. len(mod) >0)then
	  ref = clean_string(ref)
	  mod = clean_string(mod)
  endif

  !Check if there are whitespaces:
  ws = achar(32)
  position = index(ref,ws)
  
  if (position == 0) then !There are no white spaces
     r_ref = str_to_scalar(ref)
     r_mod = str_to_scalar(mod)
     compare_only_numbers = compare_numbers(r_ref,r_mod,tol)
  else
     !There are whitespaces: compare array
     compare_only_numbers = compare_array(ref,mod,label,tol)     
  endif
end function compare_only_numbers

!---------------------------------------------------

function compare_alpha(c_ref,c_mod)
  type(string), intent(in) :: c_ref,c_mod
  logical compare_alpha
  
  
  !TODO: maybe something more elaborate?

  if (char(c_ref) /= char(c_mod)) then

     !print *,"  Alphanumeric strings don't match"
     !print *,"     Reference:",char(c_ref)
     !print *,"     Modified:",char(c_mod)
     !print *, "ref:",char(c_ref)
     !print *, "mod:",char(c_mod)
     compare_alpha = .true.
  else
     compare_alpha = .true.
  endif

end function compare_alpha

!---------------------------------------------------

function str_to_scalar(str) result(scal)
type(string),intent(in) :: str
real(dp) :: scal
integer  :: int

character(len=100) :: c_str

c_str = char(str)
!print*, "string:","|",char(str),"|"
if (is_integer(str)) then
   if (string_utilities_debug) print *,"          str_to_scalar: string is integer:",c_str
   read(c_str,"(i3)") int
   scal = real(int)
else
   read(c_str,"(f16.6)") scal
endif
if (string_utilities_debug) print *,"          str_to_scalar: Scalar:|",scal,"|"
end function str_to_scalar

!---------------------------------------------------

function is_integer(str)
  type(string), intent(in) :: str
  logical is_integer

  integer :: position
  character(len=2) :: po
  po = ".*"
  position = scan(str,po)
  if (position == 0)then
     is_integer = .true.
  else
     is_integer = .false.
  endif

end function is_integer

!---------------------------------------------------

function compare_numbers(ref,mod,tol)
  real(dp), intent(in) :: ref,mod,tol
  logical compare_numbers


  if (abs(abs(ref) - abs(mod)) .gt. tol)then
     compare_numbers = .false.
  else
     compare_numbers = .true.
  endif

end function compare_numbers
!---------------------------------------------------

function compare_array(ref,mod,label,tol)
  type(string), intent(in) :: ref, mod
  character(len=*),intent(in) ::label
  real(dp) , intent(in) :: tol
  logical :: compare_array

  integer :: i, length,j
  real(dp),pointer,dimension(:)    :: a_ref,a_mod

  call find_array(ref,mod, a_ref,a_mod)

  length = size(a_ref)

  !print *, "length=",length

  compare_array = .true.

  if (length > 0)then
     do i=1,length
        
        if (abs(abs(a_ref(i))-abs(a_mod(i))) > tol) then 
           print *, "          Compare array: error in array: ",trim(label)
           print *, "            ",trim(ref), " | ", trim(mod)
           !print *, "            parsed arrays:"
           !print *, "                ref=",(a_ref(j),j=1,size(a_ref))
           !print *, "                mod=",(a_mod(j),j=1,size(a_mod))               
           print "(a,i4,2f10.6)","             Index, values: ",i, abs(a_ref(i)),abs(a_mod(i))
           print *,"              Diff,tol: ", abs(a_ref(i))-abs(a_mod(i)),tol
           compare_array = .false.
           exit
        endif
     enddo
  endif

end function compare_array

!---------------------------------------------------

subroutine find_array(ref,mod,a_ref,a_mod)
  type(string), intent(in) :: ref, mod
  real(dp),pointer,dimension(:) :: a_ref,a_mod

  integer      :: length = 1
  type(string) :: l_ref,l_mod,sub_ref,sub_mod
  integer :: pos_ref, pos_mod
  character :: ws
  character(len=2) ::ws2
  real(dp)         :: r_ref,r_mod

  l_ref = ref
  l_mod = mod

  allocate(a_ref(1),a_mod(1))
  
  ws = achar(32) !\n
  pos_ref = index(l_ref,ws)
  pos_mod = index(l_mod,ws)

  length = 0

  ws2="  "
  
  if (string_utilities_debug) then
        print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
  	print*,"           find_array: parsing:"
  	print*,"             ",char(l_ref)
  	print*,"             ",char(l_mod)
  endif
  
  do
     pos_ref = index(l_ref,ws2)
     !print *,"pos_ref 2sp:",pos_ref
     if (pos_ref == 0)then
        exit
     else
        l_ref = remove(l_ref,pos_ref,pos_ref)
        !print *,"New ref without 2ws:",char(l_ref)
     endif
  enddo
 
  do
     pos_mod = index(l_mod,ws2)
     !print *,"pos_mod 2sp:",pos_mod
     if (pos_mod == 0)then
        exit
     else
        l_mod = remove(l_mod,pos_mod,pos_mod)
        !print *,"New mod without 2ws:",char(l_mod)
     endif
  enddo
  
  pos_ref = index(l_ref,ws)
  pos_mod = index(l_mod,ws)

  do

     if (pos_ref == 0) exit

     sub_ref = extract(l_ref,1,pos_ref-1)
     sub_mod = extract(l_mod,1,pos_mod-1)
     
     !Get the scalars from the string.
     r_ref = str_to_scalar(sub_ref)
     r_mod = str_to_Scalar(sub_mod)
     
     if (string_utilities_debug) print*, "            Strings ref|mod:",char(sub_ref),"|", &
          char(sub_mod)
     if (string_utilities_debug) print*, "            New scalars ref|mod:",r_ref,"|",r_mod

     !Remove the previous string
     l_ref = remove(l_ref,1,pos_ref)
     l_mod = remove(l_mod,1,pos_mod)
     
     !Store the values
     length = length + 1

     !Resize the array
     call resize(a_ref,length)
     call resize(a_mod,length)
      
     !store the new values
     a_ref(length) = r_ref
     a_mod(length) = r_mod

     !Look for the next ws.
     pos_ref = index(l_ref,ws)
     pos_mod = index(l_mod,ws)

     if (pos_ref > 0)then
        sub_ref = extract(l_ref,1,pos_ref+1)
        sub_mod = extract(l_mod,1,pos_mod+1)
     endif

  enddo
  

  if (len(sub_ref) > 0)then
     !print *, "one missing"
     r_ref = str_to_scalar(sub_ref)
     r_mod = str_to_Scalar(sub_mod)
     length = length + 1
     call resize(a_ref,length)
     call resize(a_mod,length)
      
     !store the new values
     a_ref(length) = r_ref
     a_mod(length) = r_mod
  endif
  
  !print *,"Final ref:",(a_ref(i),"|",i=1,size(a_ref))
  !print *,"Final mod:",(a_mod(i),"|",i=1,size(a_mod))
  !print *,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
end subroutine find_array

!---------------------------------------------------
subroutine resize(a,new_length) 

  real(dp), dimension(:), pointer, intent(inout) :: a
  integer, intent(in)            :: new_length

  !Internal vars
  real(dp), dimension(:), pointer :: old_a
  integer                         :: length,i
  
  !print *,"resize"
  if (.not. associated (a)) then
     allocate(a(1:new_length))
     a = 0.0
  else
     length = size(a)
     if (length > 0) then
        allocate(old_a(1:length))
        old_a = a
        deallocate(a)
        allocate(a(1:new_length))
        a = 0.0
        do i=1,size(old_a)
           a(i) = old_a(i)
        enddo
        !print *,"before old_a"
        deallocate(old_a)
        !print *,"after old_a"
     endif
  endif

end subroutine resize
!---------------------------------------------------
end module string_utilities
