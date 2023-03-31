module m_char

  implicit none

  public :: lcase, ucase, ccase

  character(len=26), parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(len=26), parameter :: lower = 'abcdefghijklmnopqrstuvwxyz'

contains

  ! Lower-case a string
  pure function lcase(str) 
    character(len=*), intent(in) :: str
    character(len=len(str)) :: lcase
    integer :: ic, i

    ! Capitalize each letter if it is lowecase
    lcase = str
    i = scan(lcase,upper)
    do while ( i > 0 ) 
       ! Get the conversion index
       ic = index(upper,lcase(i:i))
       lcase(i:i) = lower(ic:ic)
       ic = scan(lcase(i+1:),upper)
       if ( ic > 0 ) then
          i = i + ic
       else
          i = 0
       end if
    end do

  end function lcase

  ! Upper-case a string
  pure function ucase(str) 
    character(len=*), intent(in) :: str
    character(len=len(str)) :: ucase
    integer :: ic, i

    ! Capitalize each letter if it is lowecase
    ucase = str
    i = scan(ucase,lower)
    do while ( i > 0 ) 
       ! Get the conversion index
       ic = index(lower,ucase(i:i))
       ucase(i:i) = upper(ic:ic)
       ic = scan(ucase(i+1:),lower)
       if ( ic > 0 ) then
          i = i + ic
       else
          i = 0
       end if
    end do

  end function ucase

  ! convert a string based on input
  pure function ccase(str,case) 
    character(len=*), intent(in) :: str
    character(len=1), intent(in), optional :: case
    character(len=len(str)) :: ccase

    if ( present(case) ) then
       if ( scan(case,'Uu') > 0 ) then
          ccase = ucase(str)
       else if ( scan(case,'Ll') > 0 ) then
          ccase = lcase(str)
       else
          ccase = str
       end if
    else
       ccase = str
    end if
  end function ccase

end module m_char
  
