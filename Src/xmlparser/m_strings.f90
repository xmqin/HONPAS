!*******************************************************************************
! module STRINGS
! Mart Rentmeester, Mart.Rentmeester@sci.kun.nl
! http://nn-online.sci.kun.nl/fortran
! Version 1.0
!*******************************************************************************

      module m_strings

      private

      type string
          private
          integer                 :: len = 0
          integer                 :: size = 0

          character, pointer      :: chars(:) => null()

      end type string

      character, parameter :: blank = ' '

!     GENERIC PROCEDURE INTERFACE DEFINITIONS

!---- LEN interface
      interface len
          module procedure len_s
      end interface

!---- Conversion (to CHAR) procedure interfaces
      interface char
          module procedure s_to_c, &! string to character
                           s_to_slc       ! string to specified length character
      end interface

!---- ASSIGNMENT interfaces
      interface assignment(=)
          module procedure assign_s_to_s, &! string = string
                           assign_s_to_c, &! character = string
                           assign_c_to_s        ! string = character
      end interface

!---- // operator interfaces
      interface operator(//)
          module procedure s_concat_s, &! string // string
                           s_concat_c, &! string // character
                           c_concat_s     ! character // string
      end interface

!---- INSERT_IN_STRING interface
      interface insert_in_string
          module procedure insert_in_string_c, insert_in_string_s
      end interface

!---- PREPEND_TO_STRING interface
      interface prepend_to_string
          module procedure prepend_to_string_c, prepend_to_string_s
      end interface

!---- APPEND_TO_STRING interface
      interface append_to_string
          module procedure append_to_string_c, append_to_string_s
      end interface

!---- REPLACE_IN_STRING interface
      interface replace_in_string
          module procedure replace_in_string_sc_s,   replace_in_string_ss_s,   &
                           replace_in_string_sc_sf,  replace_in_string_ss_sf,  &
                           replace_in_string_scc,    replace_in_string_ssc,    &
                           replace_in_string_scs,    replace_in_string_sss,    &
                           replace_in_string_scc_f,  replace_in_string_ssc_f,  &
                           replace_in_string_scs_f,  replace_in_string_sss_f
      end interface


!---- REPEAT interface
      interface repeat
          module procedure repeat_s
      end interface

!---- ==  .eq. comparison operator interfaces
      interface operator(==)
          module procedure s_eq_s, &! string == string
                           s_eq_c, &! string == character
                           c_eq_s         ! character == string
      end interface

!---- /=  .ne. comparison operator interfaces
      interface operator(/=)
          module procedure s_ne_s, &! string /= string
                           s_ne_c, &! string /= character
                           c_ne_s         ! character /= string
      end interface

!---- <  .lt. comparison operator interfaces
      interface operator(<)
          module procedure s_lt_s, &! string < string
                           s_lt_c, &! string < character
                           c_lt_s         ! character < string
      end interface

!---- <=  .le. comparison operator interfaces
      interface operator(<=)
          module procedure s_le_s, &! string <= string
                           s_le_c, &! string <= character
                           c_le_s         ! character <= string
      end interface

!---- >=  .ge. comparison operator interfaces
      interface operator(>=)
          module procedure s_ge_s, &! string >= string
                           s_ge_c, &! string >= character
                           c_ge_s         ! character >= string
      end interface

!---- >  .gt. comparison operator interfaces
      interface operator(>)
          module procedure s_gt_s, &! string > string
                           s_gt_c, &! string > character
                           c_gt_s         ! character > string
      end interface

!---- .aeq. comparison operator interfaces
      interface operator(.aeq.)
          module procedure a_eq_a, &! array == array
                           a_eq_c, &! array == character
                           c_eq_a         ! character == array
      end interface

!---- .ane. comparison operator interfaces
      interface operator(.ane.)
          module procedure a_ne_a, &! array /= array
                           a_ne_c, &! array /= character
                           c_ne_a         ! character /= array
      end interface

!---- .alt. comparison operator interfaces
      interface operator(.alt.)
          module procedure a_lt_a, &! array < array
                           a_lt_c, &! array < character
                           c_lt_a         ! character < array
      end interface

!---- .ale. comparison operator interfaces
      interface operator(.ale.)
          module procedure a_le_a, &! array <= array
                           a_le_c, &! array <= character
                           c_le_a         ! character <= array
      end interface

!---- .age. comparison operator interfaces
      interface operator(.age.)
          module procedure a_ge_a, &! array >= array
                           a_ge_c, &! array >= character
                           c_ge_a         ! character >= array
      end interface

!---- .agt. comparison operator interfaces
      interface operator(.agt.)
          module procedure a_gt_a, &! array > array
                           a_gt_c, &! array > character
                           c_gt_a         ! character > array
      end interface

!---- LLT comparison function interfaces
      interface llt
          module procedure s_llt_s, &! llt(string,string)
                           s_llt_c, &! llt(string,character)
                           c_llt_s        ! llt(character,string)
      end interface

!---- LLE comparison function interfaces
      interface lle
          module procedure s_lle_s, &! lle(string,string)
                           s_lle_c, &! lle(string,character)
                           c_lle_s        ! lle(character,string)
      end interface

!---- LGE comparison function interfaces
      interface lge
          module procedure s_lge_s, &! lge(string,string)
                           s_lge_c, &! lge(string,character)
                           c_lge_s        ! lge(character,string)
      end interface

!---- LGT comparison function interfaces
      interface lgt
          module procedure s_lgt_s, &! lgt(string,string)
                           s_lgt_c, &! lgt(string,character)
                           c_lgt_s        ! lgt(character,string)
      end interface

!---- ALLT comparison function interfaces
      interface allt
          module procedure a_allt_a, &! allt(array,array)
                           a_allt_c, &! allt(array,character)
                           c_allt_a        ! allt(character,array)
      end interface

!---- ALLE comparison function interfaces
      interface alle
          module procedure a_alle_a, &! alle(array,array)
                           a_alle_c, &! alle(array,character)
                           c_alle_a        ! alle(character,array)
      end interface

!---- ALGE comparison function interfaces
      interface alge
          module procedure a_alge_a, &! alge(array,array)
                           a_alge_c, &! alge(array,character)
                           c_alge_a        ! alge(character,array)
      end interface

!---- ALGT comparison function interfaces
      interface algt
          module procedure a_algt_a, &! algt(array,array)
                           a_algt_c, &! algt(array,character)
                           c_algt_a        ! algt(character,array)
      end interface

!---- INDEX procedure
      interface index
          module procedure index_ss, index_sc, index_cs
      end interface

!---- AINDEX procedure
      interface aindex
          module procedure aindex_aa, aindex_ac, aindex_ca
      end interface

!---- SCAN procedure
      interface scan
          module procedure scan_ss, scan_sc, scan_cs
      end interface

!---- ASCAN procedure
      interface ascan
          module procedure ascan_aa, ascan_ac, ascan_ca
      end interface

!---- VERIFY procedure
      interface verify
          module procedure verify_ss, verify_sc, verify_cs
      end interface

!---- AVERIFY procedure
      interface averify
          module procedure averify_aa, averify_ac, averify_ca
      end interface

!---- TRIM interface
      interface len_trim
          module procedure len_trim_s
      end interface

!---- LEN_TRIM interface
      interface trim
          module procedure trim_s
      end interface

!---- IACHAR interface
      interface iachar
          module procedure iachar_s
      end interface

!---- ICHAR interface
      interface ichar
          module procedure ichar_s
      end interface

!---- ADJUSTL interface
      interface adjustl
          module procedure adjustl_s
      end interface

!---- ADJUSTR interface
      interface adjustr
          module procedure adjustr_s
      end interface

!---- LEN_STRIP interface
      interface len_strip
          module procedure len_strip_c, len_strip_s
      end interface

!---- STRIP interface
      interface strip
          module procedure strip_c, strip_s
      end interface

!---- UPPERCASE interface
      interface uppercase
          module procedure uppercase_s, uppercase_c
      end interface

!---- TO_UPPERCASE interface
      interface to_uppercase
          module procedure to_uppercase_s, to_uppercase_c
      end interface

!---- LOWERCASE interface
      interface lowercase
          module procedure lowercase_s, lowercase_c
      end interface

!---- TO_LOWERCASE interface
      interface to_lowercase
          module procedure to_lowercase_s, to_lowercase_c
      end interface

!---- EXTRACT interface
      interface extract
          module procedure extract_s, extract_c
      end interface

!---- SUBSTRING interface
      interface substring
          module procedure extract_s, extract_c
      end interface

!---- REMOVE interface
      interface remove
          module procedure remove_s, remove_c
      end interface

!---- INSERT interface
      interface insert
          module procedure insert_ss, insert_cs, insert_sc, insert_cc
      end interface

!---- REPLACE interface
      interface replace
          module procedure replace_cc_s,   replace_cs_s,   &
                           replace_sc_s,   replace_ss_s,   &
                           replace_cc_sf,  replace_cs_sf,  &
                           replace_sc_sf,  replace_ss_sf,  &
                           replace_ccc,    replace_csc,    &
                           replace_ccs,    replace_css,    &
                           replace_scc,    replace_ssc,    &
                           replace_scs,    replace_sss,    &
                           replace_ccc_f,  replace_csc_f,  &
                           replace_ccs_f,  replace_css_f,  &
                           replace_scc_f,  replace_ssc_f,  &
                           replace_scs_f,  replace_sss_f
      end interface

!---- SORT interface
      interface sort
          module procedure sort_c, sort_s
      end interface

!---- LSORT interface
      interface lsort
          module procedure lsort_c, lsort_s
      end interface

!---- RANK interface
      interface rank
          module procedure rank_c, rank_s
      end interface

!---- LRANK interface
      interface lrank
          module procedure lrank_c, lrank_s
      end interface



!---- Publically accessible entities
      public :: string
      public :: assignment(=),unstring
      public :: insert,replace,remove,extract,substring
      public :: repeat,index,scan,verify
      public :: operator(//)
      public :: operator(==),operator(/=)
      public :: operator(<),operator(<=)
      public :: operator(>),operator(>=)
      public :: llt,lle,lge,lgt
      public :: char,len,len_trim,trim,iachar,ichar,adjustl,adjustr
      public :: lowercase,to_lowercase,uppercase,to_uppercase
      public :: strip,len_strip
      public :: sort,rank,lsort,lrank

      public :: resize_string,string_size,swap_strings
      public :: trim_string,strip_string
      public :: adjustl_string,adjustr_string
      public :: insert_in_string,remove_from_string
      public :: prepend_to_string,append_to_string
      public :: replace_in_string




      contains

!*******************************************************************************
!     LEN
!*******************************************************************************

      elemental function len_s(s)

      implicit none
      type(string), intent(in)  :: s
      integer                   :: len_s


      len_s = s%len

      end function len_s

!*******************************************************************************
!     STRING_SIZE
!*******************************************************************************

      elemental function string_size(s)

      implicit none
      type(string), intent(in)  :: s
      integer                   :: string_size


      string_size = s%size

      end function string_size

!*******************************************************************************
!     CHAR
!*******************************************************************************
!     Returns the characters of string as an automatically sized character

      pure function s_to_c(s)

      implicit none
      type(string),intent(in)   :: s
      character(len(s))         :: s_to_c


      s_to_c = transfer(s%chars(1:len(s)),s_to_c)

      end function s_to_c

!*******************************************************************************
!     Returns the character of fixed length, length, containing the characters
!     of string either padded with blanks or truncated on the right to fit

      pure function s_to_slc(s,length)

      implicit none
      type(string),intent(in)  :: s
      integer, intent(in)      :: length
      character(length)        :: s_to_slc
      integer                  :: i,lc


      lc = min(len(s),length)
      s_to_slc(1:lc) = transfer(s%chars(1:lc),s_to_slc)

!     Result longer than string: padding needed
      if (lc < length) s_to_slc(lc+1:length) = blank

      end function s_to_slc

!*******************************************************************************
! Assign a string value to a string variable overriding default assignement.
! Reallocates string variable to size of string value and copies characters.

      elemental subroutine assign_s_to_s(var,expr)

      implicit none
      type(string), intent(out)  :: var
      type(string), intent(in)   :: expr



      if (associated(var%chars,expr%chars)) then
!         Identity assignment: nothing to be done
          continue
      else
          if (associated(var%chars)) deallocate(var%chars)

          var%size = expr%size
          var%len = expr%len
!AG
          if (associated(expr%chars)) then
             allocate(var%chars(1:var%size))
             var%chars(1:var%len) = expr%chars(1:var%len)
          endif
      endif


      end subroutine assign_s_to_s

!*******************************************************************************
! Assign a string value to a character variable.
! If the string is longer than the character truncate the string on the right.
! If the string is shorter the character is blank padded on the right.

      elemental subroutine assign_s_to_c(var,expr)

      implicit none
      character(*), intent(out)  :: var
      type(string), intent(in)   :: expr
      integer                    :: i,lc,ls


      lc = len(var);
      ls = min(len(expr),lc)

      var(1:ls) = transfer(expr%chars(1:ls),var(1:ls))

      do i=ls+1,lc
          var(i:i) = blank
      enddo

      end subroutine assign_s_to_c

!*******************************************************************************
!     Assign a character value to a string variable.
!     Disassociates the string variable from its current value, allocates new
!     space to hold the characters and copies them from the character value
!     into this space.

      elemental subroutine assign_c_to_s(var,expr)

      implicit none
      type(string), intent(out)  :: var
      character(*), intent(in)   :: expr
      integer                    :: i,lc



      if (associated(var%chars)) deallocate(var%chars)


      lc = len(expr)
      var%len = lc
      var%size = lc
      allocate(var%chars(1:lc))
!!AG: NAG compiler uses temporaries here:
      var%chars(:) = (/ (expr(i:i), i=1,lc) /)

      endsubroutine assign_c_to_s

!*******************************************************************************
!     RESIZE_STRING procedure
!*******************************************************************************

!*** return code
!*** n < 0  --> deallocate?

!     pure subroutine resize_string(s,newsize,status)
      pure subroutine resize_string(s,newsize)

      implicit none
      type(string), intent(inout)     :: s
      integer, intent(in)             :: newsize
!     integer, intent(out), optional  :: status

      character, pointer              :: c(:)

      integer                         :: i


      if (newsize <= 0) return


      if (associated(s%chars)) then

          i = min(newsize,s%len)
          allocate(c(i))
          c(:) = s%chars(1:i)
          deallocate(s%chars)

          s%chars => c

          s%len = i
          s%size = newsize
      else
          s%size = newsize
          s%len = 0
          allocate(s%chars(s%size))
      endif

      end subroutine resize_string

!*******************************************************************************
!     SWAP_STRINGS
!*******************************************************************************
      subroutine swap_strings(s1,s2)


      implicit none
      type(string), intent(inout)  :: s1,s2
      integer                      :: l,s
      character, pointer           :: c(:)


      l = s1%len
      s = s1%size
      c => s1%chars
      s1%len = s2%len
      s1%size = s2%size
      s1%chars => s2%chars
      s2%len = l
      s2%size = s
      s2%chars => c

      end subroutine swap_strings

!*******************************************************************************
!     TRIM_STRINGSIZE
!*******************************************************************************

      subroutine trim_stringsize(s)

      implicit none
      type(string), intent(inout)  :: s


      call resize_string(s,len(s))

      end subroutine trim_stringsize

!*******************************************************************************
!     TRIM_STRING
!*******************************************************************************

      subroutine trim_string(s)

      implicit none
      type(string), intent(inout)  :: s


      s%len = len_trim(s)

      end subroutine trim_string

!*******************************************************************************
!    STRIP
!*******************************************************************************

     pure subroutine strip_string(s)

     implicit none
     type(string), intent(inout)  :: s
     integer                      :: i,i1,i2


     do i1=1,len(s)
         if (s%chars(i1) /= blank) exit
     enddo
     do i2=len(s),1,-1
         if (s%chars(i2) /= blank) exit
     enddo
     do i=i1,i2
         s%chars(i-i1+1) = s%chars(i)
     enddo
     s%len = i2 - i1 + 1

     end subroutine strip_string

!*******************************************************************************
!     ADJUSTL_STRING
!*******************************************************************************
! Returns as a character variable the string adjusted to the left,
! removing leading blanks and inserting trailing blanks.

      pure subroutine adjustl_string(s)

      implicit none
      type(string), intent(inout)  :: s
      integer                      :: i,j


      do i=1,len(s)
          if (s%chars(i) /= blank) exit
      enddo
      do j=i,len(s)
          s%chars(j-i:j-i) = s%chars(j)
      enddo
      s%chars(j+1:) = blank

      end subroutine adjustl_string

!*******************************************************************************
!     ADJUSTR_STRING
!*******************************************************************************
! Returns as a character variable the string adjusted to the right,
! removing trailing blanks and inserting leading blanks.

      pure subroutine adjustr_string(s)

      implicit none
      type(string), intent(inout)  :: s
      integer                      :: i,j,l,lt


      l = len(s)
      lt = len_trim(s)

      i = l - lt

      do j=1,lt
          s%chars(j+i:j+i) = s%chars(j)
      enddo
      s%chars(1:i) = blank


      end subroutine adjustr_string

!*******************************************************************************
!     PREPEND_TO_STRING
!*******************************************************************************

      pure subroutine prepend_to_string_s(s1,s2)

      implicit none
      type(string), intent(inout)  :: s1
      type(string), intent(in)     :: s2
      integer                      :: i,ls1,ls2

      character, pointer           :: ss(:)


      ls1 = len(s1)
      ls2 = len(s2)
      if (ls2 == 0) return
      if (ls1+ls2 > string_size(s1)) then
          allocate(ss(ls1+ls2))
          do i=1,ls2
              ss(i) = s2%chars(i)
          enddo
          do i=1,ls1
              ss(ls2+i) = s1%chars(i)
          enddo
          deallocate(s1%chars)

          s1%chars => ss

          s1%len = ls1 + ls2
          s1%size = s1%len
      else
          do i=ls1,1,-1
              s1%chars(ls2+i) = s1%chars(i)
          enddo
          do i=1,ls2
              s1%chars(i) = s2%chars(i)
          enddo
          s1%len = ls1 + ls2
      endif

      end subroutine prepend_to_string_s

!*******************************************************************************

      pure subroutine prepend_to_string_c(s,c)

      implicit none
      type(string), intent(inout)  :: s
      character(*), intent(in)     :: c
      integer                      :: i,ls,lc

      character, pointer           :: ss(:)



      ls = len(s)
      lc = len(c)
      if (lc == 0) return
      if (ls+lc > string_size(s)) then
          allocate(ss(ls+lc))
          do i=1,lc
              ss(i) = c(i:i)
          enddo
          do i=1,ls
              ss(lc+i) = s%chars(i)
          enddo
          deallocate(s%chars)

          s%chars => ss

          s%len = ls + lc
          s%size = s%len
      else
          do i=ls,1,-1
              s%chars(lc+i) = s%chars(i)
          enddo
          do i=1,lc
              s%chars(i) = c(i:i)
          enddo
          s%len = ls + lc
      endif

      end subroutine prepend_to_string_c

!*******************************************************************************
!     APPEND_TO_STRING
!*******************************************************************************

      pure subroutine append_to_string_s(s1,s2)

      implicit none
      type(string), intent(inout)  :: s1
      type(string), intent(in)     :: s2
      integer                      :: i,ls1,ls2

      character, pointer           :: ss(:)


      ls1 = len(s1)
      ls2 = len(s2)
      if (ls2 == 0) return
      if (ls1+ls2 > string_size(s1)) then
          allocate(ss(ls1+ls2))
          do i=1,ls1
              ss(i) = s1%chars(i)
          enddo
          do i=ls1+1,ls1+ls2
              ss(i) = s2%chars(i-ls1)
          enddo
          deallocate(s1%chars)

          s1%chars => ss

          s1%len = ls1 + ls2
          s1%size = s1%len
      else
          do i=ls1+1,ls1+ls2
              s1%chars(i) = s2%chars(i-ls1)
          enddo
          s1%len = ls1 + ls2
      endif

      end subroutine append_to_string_s

!*******************************************************************************

      pure subroutine append_to_string_c(s,c)

      implicit none
      type(string), intent(inout)  :: s
      character(*), intent(in)     :: c
      integer                      :: i,ls,lc

      character, pointer           :: ss(:)



      ls = len(s)
      lc = len(c)
      if (lc == 0) return
      if (ls+lc > string_size(s)) then
          allocate(ss(ls+lc))
          do i=1,ls
              ss(i) = s%chars(i)
          enddo
          do i=ls+1,ls+lc
              ss(i) = c(i-ls:i-ls)
          enddo
          deallocate(s%chars)

          s%chars => ss

          s%len = ls + lc
          s%size = s%len
      else
          do i=ls+1,ls+lc
              s%chars(i) = c(i-ls:i-ls)
          enddo
          s%len = ls + lc
      endif

      end subroutine append_to_string_c

!*******************************************************************************
!     INSERT_IN_STRING
!*******************************************************************************

      pure subroutine insert_in_string_s(s1,start,s2)

      implicit none
      type(string), intent(inout)  :: s1
      type(string), intent(in)     :: s2
      integer, intent(in)          :: start
      integer                      :: i,ip,is,ls1,ls2

      character, pointer           :: ss(:)


      ls1 = len(s1)
      ls2 = len(s2)
      if (ls2 == 0) return
      if (ls1+ls2 > string_size(s1)) then
          allocate(ss(ls1+ls2))
          is = max(start,1)
          ip = min(ls1+1,is)
          do i=1,ip-1
              ss(i) = s1%chars(i)
          enddo
          do i=ip,ip+ls2-1
              ss(i) = s2%chars(i-ip+1)
          enddo
          do i=ip+ls2,ls1+ls2
              ss(i) = s1%chars(i-ls2)
          enddo
          deallocate(s1%chars)

          s1%chars => ss

          s1%len = ls1 + ls2
          s1%size = s1%len
      else
          is = max(start,1)
          ip = min(ls1+1,is)
          do i=ls1+ls2,ip+ls2,-1
              s1%chars(i) = s1%chars(i-ls2)
          enddo
          do i=ip,ip+ls2-1
              s1%chars(i) = s2%chars(i-ip+1)
          enddo
          s1%len = ls1 + ls2
      endif

      end subroutine insert_in_string_s

!*******************************************************************************

      pure subroutine insert_in_string_c(s,start,c)

      implicit none
      type(string), intent(inout)  :: s
      character(*), intent(in)     :: c
      integer, intent(in)          :: start
      integer                      :: i,ip,is,ls,lc

      character, pointer           :: ss(:)



      ls = len(s)
      lc = len(c)
      if (lc == 0) return
      if (ls+lc > string_size(s)) then
          allocate(ss(ls+lc))
          is = max(start,1)
          ip = min(ls+1,is)
          do i=1,ip-1
              ss(i) = s%chars(i)
          enddo
          do i=ip,ip+lc-1
              ss(i) = c(i-ip+1:i-ip+1)
          enddo
          do i=ip+lc,ls+lc
              ss(i) = s%chars(i-lc)
          enddo
          deallocate(s%chars)

          s%chars => ss

          s%len = ls + lc
          s%size = s%len
      else
          is = max(start,1)
          ip = min(ls+1,is)
          do i=ls+lc,ip+lc,-1
              s%chars(i) = s%chars(i-lc)
          enddo
          do i=ip,ip+lc-1
              s%chars(i) = c(i-ip+1:i-ip+1)
          enddo
          s%len = ls + lc
      endif

      end subroutine insert_in_string_c

!*******************************************************************************
!     REPLACE_IN_STRING
!*******************************************************************************
!     pure subroutine replace_in_string_ss_s(s,start,ss)
!
!     implicit none
!     type(string), intent(inout)  :: s
!     type(string), intent(in)     :: ss
!     integer, intent(in)          :: start
!
!
!     call replace_in_string_sc_s(s,start,char(ss))
!
!     end subroutine replace_in_string_ss_s
!*******************************************************************************

!*******************************************************************************

      pure subroutine replace_in_string_ss_s(s,start,ss)

      implicit none
      type(string), intent(inout)  :: s
      type(string), intent(in)     :: ss
      integer, intent(in)          :: start
      integer                      :: i,ip,is,lr,lss,ls
      character, pointer           :: rs(:)
      logical                      :: new


      lr = lr_ss_s(s,start,ss)
      lss = len(ss)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)

      new = lr > string_size(s)

      if (new) then
          allocate(rs(lr))
      else
          rs => s%chars
      endif

      do i=lr,ip+lss,-1
          rs(i) = s%chars(i)
      enddo
      do i=lss,1,-1
          rs(ip-1+i) = ss%chars(i)
      enddo
      if (new) then
          do i=1,ip-1
              rs(i) = s%chars(i)
          enddo
      endif

      if (new) then
          deallocate(s%chars)
          s%chars => rs
          s%size = lr
      else
          nullify(rs)
      endif
      s%len = lr

      end subroutine replace_in_string_ss_s

!*******************************************************************************
!     pure subroutine replace_in_string_ss_sf(s,start,finish,ss)
!
!     implicit none
!     type(string), intent(inout)  :: s
!     type(string), intent(in)     :: ss
!     integer, intent(in)          :: start,finish
!
!
!     call replace_in_string_sc_sf(s,start,finish,char(ss))
!
!     end subroutine replace_in_string_ss_sf
!*******************************************************************************

!*******************************************************************************

      pure subroutine replace_in_string_ss_sf(s,start,finish,ss)

      implicit none
      type(string), intent(inout)  :: s
      type(string), intent(in)     :: ss
      integer, intent(in)          :: start,finish
      integer                      :: i,if,ip,is,lr,ls,lss
      character, pointer           :: rs(:)
      logical                      :: new


      lr = lr_ss_sf(s,start,finish,ss)
      lss = len(ss)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)
      if = max(ip-1,min(finish,ls))

      new = lr > string_size(s)

      if (new) then
          allocate(rs(lr))
      else
          rs => s%chars
      endif

      do i=1,lr-ip-lss+1
          rs(i+ip+lss-1) = s%chars(if+i)
      enddo
      do i=lss,1,-1
          rs(i+ip-1) = ss%chars(i)
      enddo
      if (new) then
          do i=1,ip-1
              rs(i) = s%chars(i)
          enddo
      endif

      if (new) then
          deallocate(s%chars)
          s%chars => rs
          s%size = lr
      else
          nullify(rs)
      endif
      s%len = lr

      end subroutine replace_in_string_ss_sf

!*******************************************************************************

!*******************************************************************************

      pure subroutine replace_in_string_sc_s(s,start,c)

      implicit none
      type(string), intent(inout)  :: s
      character(*), intent(in)     :: c
      integer, intent(in)          :: start
      integer                      :: i,ip,is,lc,lr,ls
      character, pointer           :: rs(:)
      logical                      :: new


      lr = lr_sc_s(s,start,c)
      lc = len(c)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)

      new = lr > string_size(s)

      if (new) then
          allocate(rs(lr))
      else
          rs => s%chars
      endif

      do i=lr,ip+lc,-1
          rs(i) = s%chars(i)
      enddo
      do i=lc,1,-1
          rs(ip-1+i) = c(i:i)
      enddo
      if (new) then
          do i=1,ip-1
              rs(i) = s%chars(i)
          enddo
      endif

      if (new) then
          deallocate(s%chars)
          s%chars => rs
          s%size = lr
      else
          nullify(rs)
      endif
      s%len = lr

      end subroutine replace_in_string_sc_s

!*******************************************************************************

!*******************************************************************************

      pure subroutine replace_in_string_sc_sf(s,start,finish,c)

      implicit none
      type(string), intent(inout)  :: s
      character(*), intent(in)     :: c
      integer, intent(in)          :: start,finish
      integer                      :: i,if,ip,is,lc,lr,ls
      character, pointer           :: rs(:)
      logical                      :: new


      lr = lr_sc_sf(s,start,finish,c)
      lc = len(c)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)
      if = max(ip-1,min(finish,ls))

      new = lr > string_size(s)

      if (new) then
          allocate(rs(lr))
      else
          rs => s%chars
      endif

      do i=1,lr-ip-lc+1
          rs(i+ip+lc-1) = s%chars(if+i)
      enddo
      do i=lc,1,-1
          rs(i+ip-1) = c(i:i)
      enddo
      if (new) then
          do i=1,ip-1
              rs(i) = s%chars(i)
          enddo
      endif

      if (new) then
          deallocate(s%chars)
          s%chars => rs
          s%size = lr
      else
          nullify(rs)
      endif
      s%len = lr

      end subroutine replace_in_string_sc_sf

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************

      pure subroutine replace_in_string_scc(s,target,ss)

      implicit none
      type(string), intent(inout)  :: s
      character(*), intent(in)     :: target,ss


      call x_replace_in_string_scc(s,target,ss,'first')


      end subroutine replace_in_string_scc

!*******************************************************************************

      pure subroutine replace_in_string_scc_f(s,target,ss,action)

      implicit none
      type(string), intent(inout)  :: s
      character(*), intent(in)     :: target,ss,action


      call x_replace_in_string_scc(s,target,ss,action)

      end subroutine replace_in_string_scc_f

!*******************************************************************************

      pure subroutine x_replace_in_string_scc(s,target,ss,action)

      implicit none
      type(string), intent(inout)  :: s
      character(*), intent(in)     :: target,ss,action
      logical                      :: every,back
      integer                      :: lr,ls,lt,lss
      integer                      :: i,i1,i2,k1,k2,m1,m2

      character, pointer           :: rs(:)



      lr = lr_scc(s,target,ss,action)
      ls = len(s)
      lt = len(target)
      lss = len(ss)

      if (lt == 0) then
          if (ls == 0) then
              do i=1,lss
                  s%chars(i) = ss(i:i)
              enddo
              s%len = lss
          endif
          return
      endif

      select case(uppercase(action))
      case('FIRST')
          back = .false.
          every = .false.
      case('LAST')
          back = .true.
          every = .false.
      case('ALL')
          back = .false.
          every = .true.
      case default
          back = .false.
          every = .false.
      end select

      allocate(rs(lr))

      if (back) then
!         Backwards search

!         k2 points to the absolute position one before the target in string
          k2 = ls
          m2 = lr
          do
!             find the next occurrence of target
              i1 = aindex(s%chars(:k2),target,back)
              if (i1 == 0) then
!                 fill up to the end
                  rs(:m2) = s%chars(:k2)
                  exit
              endif
!             i1 points to the absolute position of the first
!             letter of target in string
!             i2 points to the absolute position of the last
!             letter of target in string
              i2 = i1 + lt - 1

!             copy the unaffected text string chunk after it
!             k1 points to the absolute position one after target in string
              k1 = i2 + 1
              m1 = m2 + k1 - k2
              rs(m1:m2) = s%chars(k1:k2)
              m2 = m1 - 1
              m1 = m2 - lss + 1
!             copy the replacement substring for target
              do i=1,lss
                  rs(m1+i-1) = ss(i:i)
              enddo

!             k2 points to the absolute position one before the target in string
              k2 = i1 - 1
              m2 = m1 - 1
              if (.not.every) then
                  rs(:m2) = s%chars(:k2)
                  exit
              endif
          enddo
      else
!         Forward search

!         k1 points to the absolute position one after target in string
          k1 = 1
          m1 = 1
          do
!             find the next occurrence of target
              i1 = aindex(s%chars(k1:),target)
              if (i1 == 0) then
!                 fill up to the end
                  rs(m1:lr) = s%chars(k1:ls)
                  exit
              endif
!             i1 points to the absolute position of the first
!             letter of target in string
              i1 = k1 + (i1 - 1)
!             i2 points to the absolute position of the last
!             letter of target in string
              i2 = i1 + lt - 1

!             copy the unaffected text string chunk before it
!             k2 points to the absolute position one before the target in string
              k2 = i1 - 1
              m2 = m1 + k2 - k1
              rs(m1:m2) = s%chars(k1:k2)
              m1 = m2 + 1
              m2 = m1 + lss - 1
!             copy the replacement substring for target
              do i=1,lss
                  rs(m1+i-1) = ss(i:i)
              enddo

!             k1 points to the absolute position one after target in string
              k1 = i2 + 1
              m1 = m2 + 1
              if (.not.every) then
                  rs(m1:lr) = s%chars(k1:ls)
                  exit
              endif
          enddo
      endif


      if (associated(s%chars)) deallocate(s%chars)
      s%chars => rs

      s%len = lr
      s%size = size(s%chars)

      end subroutine x_replace_in_string_scc

!*******************************************************************************

      pure subroutine replace_in_string_ssc(s,target,ss)

      implicit none
      type(string), intent(inout)  :: s
      type(string), intent(in)     :: target
      character(*), intent(in)     :: ss


      call x_replace_in_string_scc(s,char(target),ss,'first')

      end subroutine replace_in_string_ssc

!*******************************************************************************

      pure subroutine replace_in_string_ssc_f(s,target,ss,action)

      implicit none
      type(string), intent(inout)  :: s
      type(string), intent(in)     :: target
      character(*), intent(in)     :: ss,action


      call x_replace_in_string_scc(s,char(target),ss,action)

      end subroutine replace_in_string_ssc_f

!*******************************************************************************

      pure subroutine replace_in_string_scs(s,target,ss)

      implicit none
      type(string), intent(inout)  :: s
      type(string), intent(in)     :: ss
      character(*), intent(in)     :: target


      call x_replace_in_string_scc(s,target,char(ss),'first')

      end subroutine replace_in_string_scs

!*******************************************************************************

      pure subroutine replace_in_string_scs_f(s,target,ss,action)

      implicit none
      type(string), intent(inout)  :: s
      type(string), intent(in)     :: ss
      character(*), intent(in)     :: target,action


      call x_replace_in_string_scc(s,target,char(ss),action)

      end subroutine replace_in_string_scs_f

!*******************************************************************************

      pure subroutine replace_in_string_sss(s,target,ss)

      implicit none
      type(string), intent(inout)  :: s
      type(string), intent(in)     :: ss,target


      call x_replace_in_string_scc(s,char(target),char(ss),'first')

      end subroutine replace_in_string_sss

!*******************************************************************************

      pure subroutine replace_in_string_sss_f(s,target,ss,action)

      implicit none
      type(string), intent(inout)  :: s
      type(string), intent(in)     :: ss,target
      character(*), intent(in)     :: action


      call x_replace_in_string_scc(s,char(target),char(ss),action)

      end subroutine replace_in_string_sss_f

!*******************************************************************************
!     REMOVE_FROM_STRING
!*******************************************************************************

      pure subroutine remove_from_string(s,start,finish)

      implicit none
      type(string), intent(inout)                      :: s
      integer, intent(in)                              :: start,finish
      integer                                          :: i,if,is,le,ls


      is = max(1,start)
      ls = len(s)
      if = min(ls,finish)
      if (if < is) return

      le = if - is + 1            ! = len_extract
      do i=if+1,ls
          s%chars(i-le) = s%chars(i)
      enddo
      s%len = s%len - le

      end subroutine remove_from_string

!*******************************************************************************
!     UNSTRING procedure
!*******************************************************************************
!     Deallocate the chars in the string to avoid leaking of memory
!     Use this in functions and subroutines on locally declared variables
!     of type string after their use. (I.e. garbage collecting).

      elemental subroutine unstring(s)

      implicit none
      type(string), intent(inout)  :: s



      if (associated(s%chars)) deallocate(s%chars)
      nullify(s%chars)

      s%size = 0
      s%len = 0

      end subroutine unstring

!*******************************************************************************
!     //
!*******************************************************************************
!     string // string

      pure function s_concat_s(s1,s2)

      implicit none
      type(string), intent(in)    :: s1,s2
      character(len(s1)+len(s2))  :: s_concat_s
      integer                     :: l1,l2


      l1 = len(s1)
      l2 = len(s2)
      s_concat_s(1:l1) = s1
      s_concat_s(1+l1:l1+l2) = s2

      end function s_concat_s

!*******************************************************************************
!     string // character

      pure function s_concat_c(s,c)

      implicit none
      type(string), intent(in)  :: s
      character(*), intent(in)  :: c
      character(len(s)+len(c))  :: s_concat_c
      integer                   :: ls,lc


      ls = len(s)
      lc = len(c)
      s_concat_c(1:ls) = s
      s_concat_c(1+ls:ls+lc) = c

      end function s_concat_c

!*******************************************************************************
!     character // string

      pure function c_concat_s(c,s)

      implicit none
      character(*), intent(in)  :: c
      type(string), intent(in)  :: s
      character(len(s)+len(c))  :: c_concat_s
      integer                   :: lc,ls


      lc = len(c)
      ls = len(s)
      c_concat_s(1:lc) = c
      c_concat_s(1+lc:lc+ls) = s

      end function c_concat_s

!*******************************************************************************
!     REPEAT
!*******************************************************************************

      function repeat_s(s,ncopies)

      implicit none
      type(string), intent(in)   :: s
      integer, intent(in)        :: ncopies
      character(ncopies*len(s))  :: repeat_s


      if (ncopies < 0) stop 'Negative ncopies requested in REPEAT'
      repeat_s = repeat(char(s),ncopies)

      end function repeat_s

!*******************************************************************************
!     LEN_TRIM
!*******************************************************************************

      elemental function len_trim_s(s)

      implicit none
      type(string), intent(in)  :: s
      integer                   :: len_trim_s

      if (len(s) == 0) then
        len_trim_s = 0
        return
      endif
      do len_trim_s = len(s),1,-1
          if (s%chars(len_trim_s) /= blank) return
      end do

      end function len_trim_s

!*******************************************************************************
!     TRIM
!*******************************************************************************

      pure function trim_s(s)

      implicit none
      type(string), intent(in)  :: s
      character(len_trim(s))    :: trim_s
      integer                   :: i


      do i=1,len(trim_s)
          trim_s(i:i) = s%chars(i)
      enddo

      end function trim_s

!*******************************************************************************
!     IACHAR
!*******************************************************************************
! Returns the position of the character string in the ISO 646 collating
! sequence. String must be of length one, otherwise result is as for
! intrinsic iachar.

      elemental function iachar_s(s)

      implicit none
      type(string), intent(in) :: s
      integer                  :: iachar_s


      iachar_s = iachar(char(s))

      end function iachar_s

!*******************************************************************************
!     ICHAR
!*******************************************************************************
! Returns the position of character from string in the processor collating
! sequence. String must be of length one, otherwise it will behave as the
! intrinsic ichar with the equivalent character string.

      elemental function ichar_s(s)

      implicit none
      type(string), intent(in) ::  s
      integer                  ::  ichar_s


      ichar_s = ichar(char(s))

      end function ichar_s

!*******************************************************************************
!     ADJUSTL
!*******************************************************************************
! Returns as a character variable the string adjusted to the left,
! removing leading blanks and inserting trailing blanks.

      pure function adjustl_s(s)

      implicit none
      type(string), intent(in)  :: s
      character(len(s))         :: adjustl_s


      adjustl_s = adjustl(char(s))

      end function adjustl_s

!*******************************************************************************
!     ADJUSTR
!*******************************************************************************
! Returns as a character variable the string adjusted to the right,
! removing trailing blanks and inserting leading blanks.

      pure function adjustr_s(s)

      implicit none
      type(string), intent(in)  :: s
      character(len(s))         :: adjustr_s


      adjustr_s = adjustr(char(s))

      end function adjustr_s

!*******************************************************************************
!    LEN_STRIP
!*******************************************************************************

     elemental function len_strip_s(s)

     implicit none
     type(string), intent(in) :: s
     integer                  :: len_strip_s
     integer                  :: i1,i2


     do i1=1,len(s)
         if (s%chars(i1) /= blank) exit
     enddo
     do i2=len(s),1,-1
         if (s%chars(i2) /= blank) exit
     enddo
     len_strip_s = max(0,i2-i1+1)

     end function len_strip_s

!*******************************************************************************
!    STRIP
!*******************************************************************************

     pure function strip_s(s)

     implicit none
     type(string), intent(in)  :: s
     character(len_strip(s))   :: strip_s
     integer                   :: i,j


     do i=1,len(s)
         if (s%chars(i) /= blank) exit
     enddo
     do j=1,len(strip_s)
         strip_s(j:j) = s%chars(i+j-1)
     enddo

     end function strip_s

!*******************************************************************************

     elemental function len_strip_c(c)

     implicit none
     character(*), intent(in) :: c
     integer                  :: len_strip_c
     integer                  :: i1,i2


     do i1=1,len(c)
         if (c(i1:i1) /= blank) exit
     enddo
     i2 = len_trim(c)
     len_strip_c = max(0,i2-i1+1)

     end function len_strip_c

!*******************************************************************************

     pure function strip_c(c)

     implicit none
     character(*), intent(in)  :: c
     character(len_strip(c))   :: strip_c
     integer                   :: i


     do i=1,len(c)
         if (c(i:i) /= blank) exit
     enddo
     strip_c(1:) = c(i:)

     end function strip_c

!*******************************************************************************
!     EXTRACT
!*******************************************************************************

  elemental function len_extract_s(s,start,finish)

      implicit none
      type(string), intent(in)  :: s
      integer, intent(in)       :: start,finish
      integer                   :: len_extract_s
      integer                   :: is,if


      is = max(1,start)
      if = min(len(s),finish)
      if (if < is) then
          len_extract_s = 0
      else
          len_extract_s = max(0,if-is) + 1
      endif

      end function len_extract_s

      !*******************************************************************************


      elemental function len_extract_c(c,start,finish)

      implicit none
      character(*), intent(in)  :: c
      integer, intent(in)       :: start,finish
      integer                   :: len_extract_c
      integer                   :: is,if


      is = max(1,start)
      if = min(len(c),finish)
      if (if < is) then
          len_extract_c = 0
      else
          len_extract_c = max(0,if-is) + 1
      endif

      end function len_extract_c


!*******************************************************************************

      pure function extract_s(s,start,finish)

      implicit none
      type(string), intent(in)                  :: s
      integer, intent(in)                       :: start,finish
      character(len_extract_s(s,start,finish))  :: extract_s
      integer                                   :: i,is,if


      is = max(1,start)
      if = min(len(s),finish)
      if (if < is) then
          extract_s = ''
      else
          do i=1,max(0,if-is+1)
              extract_s(i:i) = s%chars(is+i-1)
          enddo
      endif

      end function extract_s

!*******************************************************************************

      pure function extract_c(c,start,finish)

      implicit none
      character(*), intent(in)                  :: c
      integer, intent(in)                       :: start,finish
      character(len_extract_c(c,start,finish))  :: extract_c
      integer                                   :: is,if


      is = max(1,start)
      if = min(len(c),finish)
      if (if < is) then
          extract_c = ''
      else
          extract_c(1:if-is+1) = c(is:if)
      endif

      end function extract_c


 
!*******************************************************************************
!     INSERT
!*******************************************************************************

      pure function insert_ss(s1,start,s2)

      implicit none
      type(string), intent(in)    :: s1,s2
      integer, intent(in)         :: start
      character(len(s1)+len(s2))  :: insert_ss
      integer                     :: i,ip,is,ls1,ls2


      ls1 = len(s1)
      ls2 = len(s2)
      is = max(start,1)
      ip = min(ls1+1,is)
      do i=1,ip-1
          insert_ss(i:i) = s1%chars(i)
      enddo
      do i=ip,ip+ls2-1
          insert_ss(i:i) = s2%chars(i-ip+1)
      enddo
      do i=ip+ls2,ls1+ls2
          insert_ss(i:i) = s1%chars(i-ls2)
      enddo

      end function insert_ss

!*******************************************************************************

      pure function insert_sc(s1,start,c2)

      implicit none
      type(string), intent(in)    :: s1
      character(*), intent(in)    :: c2
      integer, intent(in)         :: start
      character(len(s1)+len(c2))  :: insert_sc
      integer                     :: i,ip,is,ls1,ls2


      ls1 = len(s1)
      ls2 = len(c2)
      is = max(start,1)
      ip = min(ls1+1,is)
      do i=1,ip-1
          insert_sc(i:i) = s1%chars(i)
      enddo
      insert_sc(ip:ip+ls2-1) = c2
      do i=ip+ls2,ls1+ls2
          insert_sc(i:i) = s1%chars(i-ls2)
      enddo

      end function insert_sc

!*******************************************************************************

      pure function insert_cs(c1,start,s2)

      implicit none
      character(*), intent(in)    :: c1
      type(string), intent(in)    :: s2
      integer, intent(in)         :: start
      character(len(c1)+len(s2))  :: insert_cs
      integer                     :: i,ip,is,ls1,ls2


      ls1 = len(c1)
      ls2 = len(s2)
      is = max(start,1)
      ip = min(ls1+1,is)
      insert_cs(1:ip-1) = c1(1:ip-1)
      do i=ip,ip+ls2-1
          insert_cs(i:i) = s2%chars(i-ip+1)
      enddo
      insert_cs(ip+ls2:ls1+ls2) = c1(ip:ls1)

      end function insert_cs

!*******************************************************************************

      pure function insert_cc(c1,start,c2)

      implicit none
      character(*), intent(in)    :: c1,c2
      integer, intent(in)         :: start
      character(len(c1)+len(c2))  :: insert_cc
      integer                     :: ip,is,ls1,ls2


      ls1 = len(c1)
      ls2 = len(c2)
      is = max(start,1)
      ip = min(ls1+1,is)
      insert_cc(1:ip-1) = c1(1:ip-1)
      insert_cc(ip:ip+ls2-1) = c2
      insert_cc(ip+ls2:ls1+ls2) = c1(ip:ls1)

      end function insert_cc

!*******************************************************************************
!     REMOVE
!*******************************************************************************

      pure function remove_c(c,start,finish)

      implicit none
      character(*), intent(in)                         :: c
      integer, intent(in)                              :: start,finish
      character(len(c)-len_extract_c(c,start,finish))  :: remove_c
      integer                                          :: if,is,ls


      is = max(1,start)
      ls = len(c)
      if = min(ls,finish)
      if (if < is) then
          remove_c = c
      else
          remove_c = c(1:is-1) // c(if+1:)
      endif

      end function remove_c

!*******************************************************************************

      pure function remove_s(s,start,finish)

      implicit none
      type(string), intent(in)                         :: s
      integer, intent(in)                              :: start,finish
      character(len(s)-len_extract_s(s,start,finish))  :: remove_s
      integer                                          :: i,if,is,le,ls


      is = max(1,start)
      ls = len(s)
      if = min(ls,finish)
      if (if < is) then
           remove_s = s
      else
          le = if - is + 1
          do i=1,is-1
              remove_s(i:i) = s%chars(i)
          enddo
          do i=if+1,ls
              remove_s(i-le:i-le) = s%chars(i)
          enddo
      endif

      end function remove_s

!*******************************************************************************
!     REPLACE
!*******************************************************************************

      pure function lr_cc_s(s,start,ss) result(l)

      implicit none
      character(*), intent(in)  :: s,ss
      integer, intent(in)       :: start
      integer                   :: l
      integer                   :: ip,is,ls,lss


      l = max(len(s),min(len(s)+1,max(start,1)+len(ss)-1))

      end function lr_cc_s

!*******************************************************************************
!  Calculate the result string by the following actions:
!  Insert the characters from substring SS into string STR beginning
!  at position START replacing the following LEN(SUBSTRING) characters of
!  the string and enlarging string if necessary. If START is greater than
!  LEN(STRING) substring is simply appended to string by concatenation.
!  If START is less than 1, substring replaces characters in string
!  starting at 1

      function replace_cc_s(s,start,ss) result(r)

      implicit none
      character(*), intent(in)        :: s,ss
      integer, intent(in)             :: start
      character(lr_cc_s(s,start,ss))  :: r
      integer                         :: ip,is,l,lss,ls


      lss = len(ss)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)
      l = len(r)

      r(1:ip-1) = s(1:ip-1)
      r(ip:ip+lss-1) = ss
      r(ip+lss:l) = s(ip+lss:ls)

      end function replace_cc_s

!*******************************************************************************

      pure function lr_cc_sf(s,start,finish,ss) result(l)

      implicit none
      character(*), intent(in)  :: s,ss
      integer, intent(in)       :: start,finish
      integer                   :: l
      integer                   :: if,ip,is,ls,lss


      lss = len(ss)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)
      if = max(ip-1,min(finish,ls))
      l = lss + ls - if+ip-1

      end function lr_cc_sf

!*******************************************************************************
!  Calculates the result string by the following actions:
!  Insert the substring SS into string STR beginning at position
!  START replacing the following FINISH-START+1 characters of the string
!  and enlarging or shrinking the string if necessary.
!  If start is greater than LEN(STRING) substring is simply appended to
!  string by concatenation. If START is less than 1, START = 1 is used.
!  If FINISH is greater than LEN(STRING), FINISH = LEN(STRING) is used.
!  If FINISH is less than START, substring is inserted before START.

      function replace_cc_sf(s,start,finish,ss) result(r)

      implicit none
      character(*), intent(in)                :: s,ss
      integer, intent(in)                     :: start,finish
      character(lr_cc_sf(s,start,finish,ss))  :: r
      integer                                 :: i,if,ip,is,l,ls,lss


      lss = len(ss)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)
      if = max(ip-1,min(finish,ls))
      l = len(r)

      r(1:ip-1) = s(1:ip-1)
      do i=1,lss
          r(i+ip-1:i+ip-1) = ss(i:i)
      enddo
      do i=1,l-ip-lss+1
          r(i+ip+lss-1:i+ip+lss-1) = s(if+i:if+i)
      enddo

      end function replace_cc_sf

!*******************************************************************************

      pure function lr_cs_s(s,start,ss) result(l)

      implicit none
      character(*), intent(in)  :: s
      type(string), intent(in)  :: ss
      integer, intent(in)       :: start
      integer                   :: l
      integer                   :: ip,is,ls,lss


      l = max(len(s),min(len(s)+1,max(start,1)+len(ss)-1))

      end function lr_cs_s

!*******************************************************************************
!  Calculate the result string by the following actions:
!  Insert the characters from substring SS into string STR beginning
!  at position START replacing the following LEN(SUBSTRING) characters of
!  the string and enlarging string if necessary. If START is greater than
!  LEN(STRING) substring is simply appended to string by concatenation.
!  If START is less than 1, substring replaces characters in string
!  starting at 1

      function replace_cs_s(s,start,ss) result(r)

      implicit none
      character(*), intent(in)        :: s
      type(string), intent(in)        :: ss
      integer, intent(in)             :: start
      character(lr_cs_s(s,start,ss))  :: r
      integer                         :: i,ip,is,l,lss,ls


      lss = len(ss)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)
      l = len(r)

      r(1:ip-1) = s(1:ip-1)
      r(ip:ip+lss-1) = transfer(ss%chars(1:lss),r(1:lss))
      r(ip+lss:l) = s(ip+lss:ls)

      end function replace_cs_s

!*******************************************************************************

      pure function lr_cs_sf(s,start,finish,ss) result(l)

      implicit none
      character(*), intent(in)  :: s
      type(string), intent(in)  :: ss
      integer, intent(in)       :: start,finish
      integer                   :: l
      integer                   :: if,ip,is,ls,lss


      lss = len(ss)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)
      if = max(ip-1,min(finish,ls))
      l = lss + ls - if+ip-1

      end function lr_cs_sf

!*******************************************************************************
!  Calculates the result string by the following actions:
!  Insert the substring SS into string STR beginning at position
!  START replacing the following FINISH-START+1 characters of the string
!  and enlarging or shrinking the string if necessary.
!  If start is greater than LEN(STRING) substring is simply appended to
!  string by concatenation. If START is less than 1, START = 1 is used.
!  If FINISH is greater than LEN(STRING), FINISH = LEN(STRING) is used.
!  If FINISH is less than START, substring is inserted before START.

      function replace_cs_sf(s,start,finish,ss) result(r)

      implicit none
      character(*), intent(in)                :: s
      type(string), intent(in)                :: ss
      integer, intent(in)                     :: start,finish
      character(lr_cs_sf(s,start,finish,ss))  :: r
      integer                                 :: i,if,ip,is,l,ls,lss


      lss = len(ss)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)
      if = max(ip-1,min(finish,ls))
      l = len(r)

      r(1:ip-1) = s(1:ip-1)

      r(i+ip:i+ip+lss-1) = transfer(ss%chars(1:lss),r(1:lss))

      do i=1,lss
          r(i+ip-1:i+ip-1) = ss%chars(i)
      enddo

      do i=1,l-ip-lss+1
          r(i+ip+lss-1:i+ip+lss-1) = s(if+i:if+i)
      enddo

      end function replace_cs_sf

!*******************************************************************************

      pure function lr_sc_s(s,start,ss) result(l)

      implicit none
      type(string), intent(in)  :: s
      character(*), intent(in)  :: ss
      integer, intent(in)       :: start
      integer                   :: l
      integer                   :: ip,is,ls,lss


      l = max(len(s),min(len(s)+1,max(start,1)+len(ss)-1))

      end function lr_sc_s

!*******************************************************************************
!  Calculate the result string by the following actions:
!  Insert the characters from substring SS into string STR beginning
!  at position START replacing the following LEN(SUBSTRING) characters of
!  the string and enlarging string if necessary. If START is greater than
!  LEN(STRING) substring is simply appended to string by concatenation.
!  If START is less than 1, substring replaces characters in string
!  starting at 1

      function replace_sc_s(s,start,ss) result(r)

      implicit none
      type(string), intent(in)        :: s
      character(*), intent(in)        :: ss
      integer, intent(in)             :: start
      character(lr_sc_s(s,start,ss))  :: r
      integer                         :: i,ip,is,l,lss,ls


      lss = len(ss)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)
      l = len(r)

      do i=1,ip-1
          r(i:i) = s%chars(i)
      enddo

      do i=1,lss
          r(i+ip-1:i+ip-1) = ss(i:i)
      enddo

      do i=ip+lss,l
          r(i:i) = s%chars(i)
      enddo

      end function replace_sc_s

!*******************************************************************************

      pure function lr_sc_sf(s,start,finish,ss) result(l)

      implicit none
      type(string), intent(in)  :: s
      character(*), intent(in)  :: ss
      integer, intent(in)       :: start,finish
      integer                   :: l
      integer                   :: if,ip,is,ls,lss


      lss = len(ss)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)
      if = max(ip-1,min(finish,ls))
      l = lss + ls - if+ip-1

      end function lr_sc_sf

!*******************************************************************************
!  Calculates the result string by the following actions:
!  Insert the substring SS into string STR beginning at position
!  START replacing the following FINISH-START+1 characters of the string
!  and enlarging or shrinking the string if necessary.
!  If start is greater than LEN(STRING) substring is simply appended to
!  string by concatenation. If START is less than 1, START = 1 is used.
!  If FINISH is greater than LEN(STRING), FINISH = LEN(STRING) is used.
!  If FINISH is less than START, substring is inserted before START.

      function replace_sc_sf(s,start,finish,ss) result(r)

      implicit none
      type(string), intent(in)                :: s
      character(*), intent(in)                :: ss
      integer, intent(in)                     :: start,finish
      character(lr_sc_sf(s,start,finish,ss))  :: r
      integer                                 :: i,if,ip,is,l,ls,lss


      lss = len(ss)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)
      if = max(ip-1,min(finish,ls))
      l = len(r)

      do i=1,ip-1
          r(i:i) = s%chars(i)
      enddo

      r(ip:ip+lss-1) = ss

      do i=1,l-ip-lss+1
          r(i+ip+lss-1:i+ip+lss-1) = s%chars(if+i)
      enddo

      end function replace_sc_sf

!*******************************************************************************

      pure function lr_ss_s(s,start,ss) result(l)

      implicit none
      type(string), intent(in)  :: s,ss
      integer, intent(in)       :: start
      integer                   :: l
      integer                   :: ip,is,ls,lss


      l = max(len(s),min(len(s)+1,max(start,1)+len(ss)-1))

      end function lr_ss_s

!*******************************************************************************
!  Calculate the result string by the following actions:
!  Insert the characters from substring SS into string STR beginning
!  at position START replacing the following LEN(SUBSTRING) characters of
!  the string and enlarging string if necessary. If START is greater than
!  LEN(STRING) substring is simply appended to string by concatenation.
!  If START is less than 1, substring replaces characters in string
!  starting at 1

      function replace_ss_s(s,start,ss) result(r)

      implicit none
      type(string), intent(in)        :: s,ss
      integer, intent(in)             :: start
      character(lr_ss_s(s,start,ss))  :: r
      integer                         :: i,ip,is,l,lss,ls


      lss = len(ss)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)
      l = len(r)

      do i=1,ip-1
          r(i:i) = s%chars(i)
      enddo

      do i=1,lss
          r(ip-1+i:ip-1+i) = ss%chars(i)
      enddo

      do i=ip+lss,l
          r(i:i) = s%chars(i)
      enddo

      end function replace_ss_s

!*******************************************************************************

      pure function lr_ss_sf(s,start,finish,ss) result(l)

      implicit none
      type(string), intent(in)  :: s,ss
      integer, intent(in)       :: start,finish
      integer                   :: l
      integer                   :: if,ip,is,ls,lss


      lss = len(ss)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)
      if = max(ip-1,min(finish,ls))
      l = lss + ls - if+ip-1

      end function lr_ss_sf

!*******************************************************************************
!  Calculates the result string by the following actions:
!  Insert the substring SS into string STR beginning at position
!  START replacing the following FINISH-START+1 characters of the string
!  and enlarging or shrinking the string if necessary.
!  If start is greater than LEN(STRING) substring is simply appended to
!  string by concatenation. If START is less than 1, START = 1 is used.
!  If FINISH is greater than LEN(STRING), FINISH = LEN(STRING) is used.
!  If FINISH is less than START, substring is inserted before START.

      function replace_ss_sf(s,start,finish,ss) result(r)

      implicit none
      type(string), intent(in)                :: s,ss
      integer, intent(in)                     :: start,finish
      character(lr_ss_sf(s,start,finish,ss))  :: r
      integer                                 :: i,if,ip,is,l,ls,lss


      lss = len(ss)
      ls = len(s)
      is = max(start,1)
      ip = min(ls+1,is)
      if = max(ip-1,min(finish,ls))
      l = len(r)

      do i=1,ip-1
          r(i:i) = s%chars(i)
      enddo

      do i=1,lss
          r(i+ip-1:i+ip-1) = ss%chars(i)
      enddo

      do i=1,l-ip-lss+1
          r(i+ip+lss-1:i+ip+lss-1) = s%chars(if+i)
      enddo

      end function replace_ss_sf

!*******************************************************************************

      pure function lr_ccc(s,target,ss,action) result(l)

      implicit none
      character(*), intent(in)       :: s,target,ss,action
      integer                        :: l
      logical                        :: every,back
      integer                        :: ls,lt,lss,ipos,nr


      ls = len(s)
      lt = len(target)
      lss = len(ss)

      if (lt == 0) then
          if (ls == 0) then
              l = lss
          else
              l = ls
          endif
          return
      endif

      if (lt == lss) then
          l = ls
          return
      endif

      select case(uppercase(action))
      case('FIRST')
          back = .false.
          every = .false.
      case('LAST')
          back = .true.
          every = .false.
      case('ALL')
          back = .false.
          every = .true.
      case default
          back = .false.
          every = .false.
      end select

      nr = 0
      if (back) then
          ipos = ls
          do while (ipos > 0)
              ipos = index(s(:ipos),target,back)
              if (ipos == 0) exit
              nr = nr + 1
              if (.not. every) exit
              ipos = ipos - 1
          enddo
      else
          ipos = 1
          do while (ipos <= ls-lt+1)
              l = index(s(ipos:),target)
              if (l == 0) exit
              nr = nr + 1
              if (.not. every) exit
              ipos = ipos + l + 1
              ipos = ipos + 1
          enddo
      endif
      l = ls + nr*(lss-lt)

      end function lr_ccc

!*******************************************************************************

      function replace_ccc(s,target,ss) result(r)

      implicit none
      character(*), intent(in)                :: s,target,ss
      character(lr_ccc(s,target,ss,'first'))  :: r


      call x_replace_ccc(s,target,ss,'first',r)

      end function replace_ccc

!*******************************************************************************

      function replace_ccc_f(s,target,ss,action) result(r)

      implicit none
      character(*), intent(in)               :: s,target,ss,action
      character(lr_ccc(s,target,ss,action))  :: r


      call x_replace_ccc(s,target,ss,action,r)

      end function replace_ccc_f

!*******************************************************************************
!  Calculate the result string by the following actions:
!  Search for occurences of TARGET in string S, and replaces these with
!  substring SS.  If BACK present with value true search is backward otherwise
!  search is done forward.  If EVERY present with value true all accurences
!  of TARGET in S are replaced, otherwise only the first found is
!  replaced.  If TARGET is not found the result is the same as S.

      subroutine x_replace_ccc(s,target,ss,action,r)

      implicit none
      character(*), intent(in)               :: s,target,ss,action
      character(*), intent(inout)            :: r
      logical                                :: every,back
      integer                                :: lr,ls,lt,lss
      integer                                :: i1,i2,k1,k2,m1,m2


      lr = len(r)
      ls = len(s)
      lt = len(target)
      lss = len(ss)

      if (lt == 0) then
          if (ls == 0) then
              r = ss
          else
              r = s
          endif
          return
      endif

      select case(uppercase(action))
      case('FIRST')
          back = .false.
          every = .false.
      case('LAST')
          back = .true.
          every = .false.
      case('ALL')
          back = .false.
          every = .true.
      case default
          back = .false.
          every = .false.
      end select

      if (back) then
          k2 = ls
          m2 = lr
          do
              i1 = index(s(:k2),target,back)
              if (i1 == 0) then
                  r(:m2) = s(:k2)
                  return
              endif
              i2 = i1 + lt - 1
              k1 = i2 + 1
              m1 = m2 + k1 - k2
              r(m1:m2) = s(k1:k2)
              m2 = m1 - 1
              m1 = m2 - lss + 1
              r(m1:m2) = ss
              k2 = i1 - 1
              m2 = m1 - 1
              if (.not. every) then
                  r(:m2) = s(:k2)
                  return
              endif
          enddo
      else
          k1 = 1
          m1 = 1
          do
              i1 = index(s(k1:),target)
              if (i1 == 0) then
                  r(m1:) = s(k1:)
                  return
              endif
              i1 = k1 + (i1 - 1)
              i2 = i1 + lt - 1
              k2 = i1 - 1
              m2 = m1 + k2 - k1
              r(m1:m2) = s(k1:k2)
              m1 = m2 + 1
              m2 = m1 + lss - 1
              r(m1:m2) = ss
              k1 = i2 + 1
              m1 = m2 + 1
              if (.not. every) then
                  r(m1:) = s(k1:)
                  return
              endif
          enddo
      endif

      end subroutine x_replace_ccc

!*******************************************************************************

      function replace_csc(s,target,ss) result(r)

      implicit none
      character(*), intent(in)                      :: s,ss
      type(string), intent(in)                      :: target
      character(lr_ccc(s,char(target),ss,'first'))  :: r


      call x_replace_ccc(s,char(target),ss,'first',r)

      end function replace_csc

!*******************************************************************************

      function replace_csc_f(s,target,ss,action) result(r)

      implicit none
      character(*), intent(in)                     :: s,ss,action
      type(string), intent(in)                     :: target
      character(lr_ccc(s,char(target),ss,action))  :: r


      call x_replace_ccc(s,char(target),ss,action,r)

      end function replace_csc_f

!*******************************************************************************
!*******************************************************************************

      function replace_ccs(s,target,ss) result(r)

      implicit none
      character(*), intent(in)                      :: s,target
      type(string), intent(in)                      :: ss
      character(lr_ccc(s,target,char(ss),'first'))  :: r


      call x_replace_ccc(s,target,char(ss),'first',r)

      end function replace_ccs

!*******************************************************************************

      function replace_ccs_f(s,target,ss,action) result(r)

      implicit none
      character(*), intent(in)                     :: s,target,action
      type(string), intent(in)                     :: ss
      character(lr_ccc(s,target,char(ss),action))  :: r


      call x_replace_ccc(s,target,char(ss),action,r)

      end function replace_ccs_f

!*******************************************************************************
!*******************************************************************************

      function replace_css(s,target,ss) result(r)

      implicit none
      character(*), intent(in)                            :: s
      type(string), intent(in)                            :: ss,target
      character(lr_ccc(s,char(target),char(ss),'first'))  :: r


      call x_replace_ccc(s,char(target),char(ss),'first',r)

      end function replace_css

!*******************************************************************************

      function replace_css_f(s,target,ss,action) result(r)

      implicit none
      character(*), intent(in)                           :: s,action
      type(string), intent(in)                           :: ss,target
      character(lr_ccc(s,char(target),char(ss),action))  :: r


      call x_replace_ccc(s,char(target),char(ss),action,r)

      end function replace_css_f

!*******************************************************************************
!*******************************************************************************
      pure function lr_scc(s,target,ss,action) result(l)

      implicit none
      type(string), intent(in)       :: s
      character(*), intent(in)       :: target,ss,action
      integer                        :: l
      logical                        :: every,back
      integer                        :: ls,lt,lss,ipos,nr


      ls = len(s)
      lt = len(target)
      lss = len(ss)

      if (lt == 0) then
          if (ls == 0) then
              l = lss
          else
              l = ls
          endif
          return
      endif
      if (lt == lss) then
          l = ls
          return
      endif

      select case(uppercase(action))
      case('FIRST')
          back = .false.
          every = .false.
      case('LAST')
          back = .true.
          every = .false.
      case('ALL')
          back = .false.
          every = .true.
      case default
          back = .false.
          every = .false.
      end select

      nr = 0
      if (back) then
          ipos = ls
          do while (ipos > 0)
              ipos = aindex(s%chars(:ipos),target,back)
              if (ipos == 0) exit
              nr = nr + 1
              if (.not. every) exit
              ipos = ipos - 1
          enddo

      else
          ipos = 1
          do while (ipos <= ls-lt+1)
              l = aindex(s%chars(ipos:),target)
              if (l == 0) exit
              nr = nr + 1
              if (.not. every) exit
              ipos = ipos + l + 1
          enddo
      endif
      l = ls + nr*(lss-lt)

      end function lr_scc

!*******************************************************************************

      function replace_scc(s,target,ss) result(r)

      implicit none
      type(string), intent(in)                :: s
      character(*), intent(in)                :: target,ss
      character(lr_scc(s,target,ss,'first'))  :: r


      call x_replace_scc(s,target,ss,'first',r)


      end function replace_scc

!*******************************************************************************

      function replace_scc_f(s,target,ss,action) result(r)

      implicit none
      type(string), intent(in)               :: s
      character(*), intent(in)               :: target,ss,action
      character(lr_scc(s,target,ss,action))  :: r


      call x_replace_scc(s,target,ss,action,r)

      end function replace_scc_f

!*******************************************************************************
!  Calculate the result string by the following actions:
!  Search for occurences of TARGET in string S, and replaces these with
!  substring SS.  If BACK present with value true search is backward otherwise
!  search is done forward.  If EVERY present with value true all accurences
!  of TARGET in S are replaced, otherwise only the first found is
!  replaced.  If TARGET is not found the result is the same as S.

      subroutine x_replace_scc(s,target,ss,action,r)

      implicit none
      type(string), intent(in)               :: s
      character(*), intent(in)               :: target,ss,action
      character(*), intent(inout)            :: r
      logical                                :: every,back
      integer                                :: lr,ls,lt,lss
      integer                                :: i1,i2,k1,k2,m1,m2


      lr = len(r)
      ls = len(s)
      lt = len(target)
      lss = len(ss)

      if (lt == 0) then
          if (ls == 0) then
              r = ss
          else
              r = s
          endif
          return
      endif

      select case(uppercase(action))
      case('FIRST')
          back = .false.
          every = .false.
      case('LAST')
          back = .true.
          every = .false.
      case('ALL')
          back = .false.
          every = .true.
      case default
          back = .false.
          every = .false.
      end select

      if (back) then
          k2 = ls
          m2 = lr
          do
              i1 = aindex(s%chars(:k2),target,back)
              if (i1 == 0) then
                  r(:m2) = transfer(s%chars(:k2),r(:m2))
                  return
              endif
              i2 = i1 + lt - 1
              k1 = i2 + 1
              m1 = m2 + k1 - k2
              r(m1:m2) = transfer(s%chars(k1:k2),r(m1:m2))
              m2 = m1 - 1
              m1 = m2 - lss + 1
              r(m1:m2) = ss
              k2 = i1 - 1
              m2 = m1 - 1
              if (.not.every) then
                  r(:m2) = transfer(s%chars(:k2),r(:m2))
                  return
              endif
          enddo
      else
          k1 = 1
          m1 = 1
          do
              i1 = aindex(s%chars(k1:),target)
              if (i1 == 0) then
                  r(m1:) = transfer(s%chars(k1:),r(m1:))
                  return
              endif
              i1 = k1 + (i1 - 1)
              i2 = i1 + lt - 1
              k2 = i1 - 1
              m2 = m1 + k2 - k1
              r(m1:m2) = transfer(s%chars(k1:k2),r(m1:m2))
              m1 = m2 + 1
              m2 = m1 + lss - 1
              r(m1:m2) = ss
              k1 = i2 + 1
              m1 = m2 + 1
              if (.not.every) then
                  r(m1:) = transfer(s%chars(k1:),r(m1:))
                  return
              endif
          enddo
      endif

      end subroutine x_replace_scc

!*******************************************************************************

      function replace_ssc(s,target,ss) result(r)

      implicit none
      type(string), intent(in)                      :: s,target
      character(*), intent(in)                      :: ss
      character(lr_scc(s,char(target),ss,'first'))  :: r


      call x_replace_scc(s,char(target),ss,'first',r)


      end function replace_ssc

!*******************************************************************************

      function replace_ssc_f(s,target,ss,action) result(r)

      implicit none
      type(string), intent(in)                     :: s,target
      character(*), intent(in)                     :: ss,action
      character(lr_scc(s,char(target),ss,action))  :: r


      call x_replace_scc(s,char(target),ss,action,r)

      end function replace_ssc_f

!*******************************************************************************

      function replace_scs(s,target,ss) result(r)

      implicit none
      type(string), intent(in)                      :: s,ss
      character(*), intent(in)                      :: target
      character(lr_scc(s,target,char(ss),'first'))  :: r


      call x_replace_scc(s,target,char(ss),'first',r)

      end function replace_scs

!*******************************************************************************

      function replace_scs_f(s,target,ss,action) result(r)

      implicit none
      type(string), intent(in)                     :: s,ss
      character(*), intent(in)                     :: target,action
      character(lr_scc(s,target,char(ss),action))  :: r


      call x_replace_scc(s,target,char(ss),action,r)

      end function replace_scs_f

!*******************************************************************************

      function replace_sss(s,target,ss) result(r)

      implicit none
      type(string), intent(in)                            :: s,ss,target
      character(lr_scc(s,char(target),char(ss),'first'))  :: r


      call x_replace_scc(s,char(target),char(ss),'first',r)

      end function replace_sss

!*******************************************************************************

      function replace_sss_f(s,target,ss,action) result(r)

      implicit none
      type(string), intent(in)                           :: s,ss,target
      character(*), intent(in)                           :: action
      character(lr_scc(s,char(target),char(ss),action))  :: r


      call x_replace_scc(s,char(target),char(ss),action,r)

      end function replace_sss_f

!*******************************************************************************
!     SORT, LSORT
!*******************************************************************************
!*******************************************************************************
! Sorts A into ascending order, from A(1) to A(N).
! Reference: Richard C. Singleton, Algorithm 347, SORT.
! Comm. ACM 3, 321 (March 1969).
! Algorithm is Copyright 1969 Association of Computing Machinery,
!*******************************************************************************

      subroutine sort_c(a)

      implicit none
      character(*), intent(inout)  :: a(:)
      character(len(a))            :: t,s
      integer                      :: p,i,j,k,l,m
      integer                      :: is(0:63)


      m = 0
      i = 1
      j = size(a)

    5 continue
      if (i >= j) goto 70

   10 continue
      p = (i + j)/2
      t = a(p)
      if (a(i) > t) then
          a(p) = a(i)
          a(i) = t
          t = a(p)
      endif
      if (a(j) < t) then
          a(p) = a(j)
          a(j) = t
          t = a(p)
          if (a(i) > t) then
              a(p) = a(i)
              a(i) = t
              t = a(p)
          endif
      endif

      k = i
      l = j
      do
          do
              l = l - 1
              if (a(l) <= t) exit
          enddo
          s = a(l)
          do
              k = k + 1
              if (a(k) >= t) exit
          enddo
          if (k > l) exit
          a(l) = a(k)
          a(k) = s
      enddo

      if (l-i > j-k) then
          is(m) = i
          m = m + 1
          is(m) = l
          m = m + 1
          i = k
      else
          is(m) = k
          m = m + 1
          is(m) = j
          m = m + 1
          j = l
      endif
      goto 80

   70 continue
      if (m == 0) return
      m = m - 1
      j = is(m)
      m = m - 1
      i = is(m)

   80 continue
      if (j-i >= 11) goto 10
      if (i == 1) goto 5
      i = i - 1

      do
          i = i + 1
          if (i == j) goto 70
          t = a(i+1)
          if (a(i) <= t) cycle
          k = i
          do
              a(k+1) = a(k)
              k = k - 1
              if (t >= a(k)) exit
          enddo
          a(k+1) = t
      enddo

      end subroutine sort_c

!*******************************************************************************
! Sorts A into ascending order, from A(1) to A(N).
! Reference: Richard C. Singleton, Algorithm 347, SORT.
! Comm. ACM 3, 321 (March 1969).
! Algorithm is Copyright 1969 Association of Computing Machinery,
!*******************************************************************************

      subroutine sort_s(a)

      implicit none
      type(string), intent(inout)  :: a(:)
      type(string)                 :: s,t
      integer                      :: p,i,j,k,l,m
      integer                      :: is(0:63)


      m = 0
      i = 1
      j = size(a)

    5 continue
      if (i >= j) goto 70

   10 continue
      p = (i + j)/2
      call pstring(t,a(p))
      if (a(i) > t) then
          call pstring(a(p),a(i))
          call pstring(a(i),t)
          call pstring(t,a(p))
      endif
      if (a(j) < t) then
          call pstring(a(p),a(j))
          call pstring(a(j),t)
          call pstring(t,a(p))
          if (a(i) > t) then
              call pstring(a(p),a(i))
              call pstring(a(i),t)
              call pstring(t,a(p))
          endif
      endif

      k = i
      l = j
      do
          do
              l = l - 1
              if (a(l) <= t) exit
          enddo
          call pstring(s,a(l))
          do
              k = k + 1
              if (a(k) >= t) exit
          enddo
          if (k > l) exit
          call pstring(a(l),a(k))
          call pstring(a(k),s)
      enddo

      if (l-i > j-k) then
          is(m) = i
          m = m + 1
          is(m) = l
          m = m + 1
          i = k
      else
          is(m) = k
          m = m + 1
          is(m) = j
          m = m + 1
          j = l
      endif
      goto 80

   70 continue
      if (m == 0) return
      m = m - 1
      j = is(m)
      m = m - 1
      i = is(m)

   80 continue
      if (j-i >= 11) goto 10
      if (i == 1) goto 5
      i = i - 1

      do
          i = i + 1
          if (i == j) goto 70
          call pstring(t,a(i+1))
          if (a(i) <= t) cycle
          k = i
          do
              call pstring(a(k+1),a(k))
              k = k - 1
              if (t >= a(k)) exit
          enddo
          call pstring(a(k+1),t)
      enddo

      contains

!-------------------------------------------------------------------------------
      subroutine pstring(p,t)

      implicit none
      type(string), intent(inout)  :: p
      type(string), intent(in)     :: t


      p%len = t%len
      p%size = t%size
      p%chars => t%chars


      end subroutine pstring
!-------------------------------------------------------------------------------

      end subroutine sort_s

!*******************************************************************************
! Sorts A into ascending order, from A(1) to A(N).
! Reference: Richard C. Singleton, Algorithm 347, SORT.
! Comm. ACM 3, 321 (March 1969).
! Algorithm is Copyright 1969 Association of Computing Machinery,
! reproduced with permission.
!*******************************************************************************

      subroutine lsort_c(a)

      implicit none
      character(*), intent(inout)  :: a(:)
      character(len(a))            :: t,s
      integer                      :: p,i,j,k,l,m
      integer                      :: is(0:63)


      m = 0
      i = 1
      j = size(a)

    5 continue
      if (i >= j) goto 70

   10 continue
      p = (i + j)/2
      t = a(p)
      if (lgt(a(i),t)) then
          a(p) = a(i)
          a(i) = t
          t = a(p)
      endif
      if (llt(a(j),t)) then
          a(p) = a(j)
          a(j) = t
          t = a(p)
          if (lgt(a(i),t)) then
              a(p) = a(i)
              a(i) = t
              t = a(p)
          endif
      endif

      k = i
      l = j
      do
          do
              l = l - 1
              if (lle(a(l),t)) exit
          enddo
          s = a(l)
          do
              k = k + 1
              if (lge(a(k),t)) exit
          enddo
          if (k > l) exit
          a(l) = a(k)
          a(k) = s
      enddo

      if (l-i > j-k) then
          is(m) = i
          m = m + 1
          is(m) = l
          m = m + 1
          i = k
      else
          is(m) = k
          m = m + 1
          is(m) = j
          m = m + 1
          j = l
      endif
      goto 80

   70 continue
      if (m == 0) return
      m = m - 1
      j = is(m)
      m = m - 1
      i = is(m)

   80 continue
      if (j-i >= 11) goto 10
      if (i == 1) goto 5
      i = i - 1

      do
          i = i + 1
          if (i == j) goto 70
          t = a(i+1)
          if (lle(a(i),t)) cycle
          k = i
          do
              a(k+1) = a(k)
              k = k - 1
              if (lge(t,a(k))) exit
          enddo
          a(k+1) = t
      enddo

      end subroutine lsort_c

!*******************************************************************************
! Sorts A into ascending order, from A(1) to A(N).
! Reference: Richard C. Singleton, Algorithm 347, SORT.
! Comm. ACM 3, 321 (March 1969).
! Algorithm is Copyright 1969 Association of Computing Machinery,
!*******************************************************************************

      subroutine lsort_s(a)

      implicit none
      type(string), intent(inout)  :: a(:)
      type(string)                 :: s,t
      integer                      :: p,i,j,k,l,m
      integer                      :: is(0:63)


      m = 0
      i = 1
      j = size(a)

    5 continue
      if (i >= j) goto 70

   10 continue
      p = (i + j)/2
      call pstring(t,a(p))
      if (lgt(a(i),t)) then
          call pstring(a(p),a(i))
          call pstring(a(i),t)
          call pstring(t,a(p))
      endif
      if (llt(a(j),t)) then
          call pstring(a(p),a(j))
          call pstring(a(j),t)
          call pstring(t,a(p))
          if (lgt(a(i),t)) then
              call pstring(a(p),a(i))
              call pstring(a(i),t)
              call pstring(t,a(p))
          endif
      endif

      k = i
      l = j
      do
          do
              l = l - 1
              if (lle(a(l),t)) exit
          enddo
          call pstring(s,a(l))
          do
              k = k + 1
              if (lge(a(k),t)) exit
          enddo
          if (k > l) exit
          call pstring(a(l),a(k))
          call pstring(a(k),s)
      enddo

      if (l-i > j-k) then
          is(m) = i
          m = m + 1
          is(m) = l
          m = m + 1
          i = k
      else
          is(m) = k
          m = m + 1
          is(m) = j
          m = m + 1
          j = l
      endif
      goto 80

   70 continue
      if (m == 0) return
      m = m - 1
      j = is(m)
      m = m - 1
      i = is(m)

   80 continue
      if (j-i >= 11) goto 10
      if (i == 1) goto 5
      i = i - 1

      do
          i = i + 1
          if (i == j) goto 70
          call pstring(t,a(i+1))
          if (lle(a(i),t)) cycle
          k = i
          do
              call pstring(a(k+1),a(k))
              k = k - 1
              if (lge(t,a(k))) exit
          enddo
          call pstring(a(k+1),t)
      enddo

      contains

!-------------------------------------------------------------------------------
      subroutine pstring(p,t)

      implicit none
      type(string), intent(inout)  :: p
      type(string), intent(in)     :: t


      p%len = t%len
      p%size = t%size
      p%chars => t%chars


      end subroutine pstring
!-------------------------------------------------------------------------------

      end subroutine lsort_s

!*******************************************************************************
!     RANK, LRANK
!*******************************************************************************
!*******************************************************************************
! Sorts A into ascending order, from A(1) to A(N).
! Reference: Richard C. Singleton, Algorithm 347, SORT.
! Comm. ACM 3, 321 (March 1969).
! Algorithm is Copyright 1969 Association of Computing Machinery,
! reproduced with permission.
!*******************************************************************************

      subroutine rank_c(a,r)

      implicit none
      character(*), intent(in)  :: a(:)
      integer, intent(out)      :: r(size(a))
      character(len(a))         :: t
      integer                   :: i,j,k,l,m,n,p,rs,rt
      integer                   :: is(0:63)


      n = size(a)
      r(:) = (/ (i, i=1,n) /)
      m = 0
      i = 1
      j = n

    5 continue
      if (i >= j) goto 70

   10 continue
      p = (j+i)/2
      rt = r(p)
      t = a(rt)
      if (a(r(i)) > t) then
          r(p) = r(i)
          r(i) = rt
          rt = r(p)
          t = a(rt)
      endif
      if (a(r(j)) < t) then
          r(p) = r(j)
          r(j) = rt
          rt = r(p)
          t = a(rt)
          if (a(r(i)) > t) then
              r(p) = r(i)
              r(i) = rt
              rt = r(p)
              t = a(rt)
          endif
      endif

      k = i
      l = j
      do
          do
              l = l - 1
              if (a(r(l)) <= t) exit
          enddo
          rs = r(l)
          do
              k = k + 1
              if (a(r(k)) >= t) exit
          enddo
          if (k > l) exit
          r(l) = r(k)
          r(k) = rs
      enddo

      if (l-i > j-k) then
          is(m) = i
          m = m + 1
          is(m) = l
          m = m + 1
          i = k
      else
          is(m) = k
          m = m + 1
          is(m) = j
          m = m + 1
          j = l
      endif
      goto 80

   70 continue
      if (m == 0) return
      m = m - 1
      j = is(m)
      m = m - 1
      i = is(m)

   80 continue
      if (j-i >= 11) goto 10
      if (i == 1) goto 5
      i = i - 1

      do
          i = i + 1
          if (i == j) goto 70
          rt = r(i+1)
          t = a(rt)
          if (a(r(i)) <= t) cycle
          k = i
          do
              r(k+1) = r(k)
              k = k - 1
              if (t >= a(r(k))) exit
          enddo
          r(k+1) = rt
      enddo

      end subroutine rank_c

!*******************************************************************************
! Sorts A into ascending order, from A(1) to A(N).
! Reference: Richard C. Singleton, Algorithm 347, SORT.
! Comm. ACM 3, 321 (March 1969).
! Algorithm is Copyright 1969 Association of Computing Machinery,
!*******************************************************************************

      subroutine rank_s(a,r)

      implicit none
      type(string), intent(in)  :: a(:)
      integer, intent(out)      :: r(size(a))
      type(string)              :: t
      integer                   :: i,j,k,l,m,n,p,rs,rt
      integer                   :: is(0:63)


      n = size(a)
      r(:) = (/ (i, i=1,n) /)
      m = 0
      i = 1
      j = n

    5 continue
      if (i >= j) goto 70

   10 continue
      p = (j+i)/2
      rt = r(p)
      call pstring(t,a(rt))
      if (a(r(i)) > t) then
          r(p) = r(i)
          r(i) = rt
          rt = r(p)
          call pstring(t,a(rt))
      endif
      if (a(r(j)) < t) then
          r(p) = r(j)
          r(j) = rt
          rt = r(p)
          call pstring(t,a(rt))
          if (a(r(i)) > t) then
              r(p) = r(i)
              r(i) = rt
              rt = r(p)
              call pstring(t,a(rt))
          endif
      endif

      k = i
      l = j
      do
          do
              l = l - 1
              if (a(r(l)) <= t) exit
          enddo
          rs = r(l)
          do
              k = k + 1
              if (a(r(k)) >= t) exit
          enddo
          if (k > l) exit
          r(l) = r(k)
          r(k) = rs
      enddo

      if (l-i > j-k) then
          is(m) = i
          m = m + 1
          is(m) = l
          m = m + 1
          i = k
      else
          is(m) = k
          m = m + 1
          is(m) = j
          m = m + 1
          j = l
      endif
      goto 80

   70 continue
      if (m == 0) return
      m = m - 1
      j = is(m)
      m = m - 1
      i = is(m)

   80 continue
      if (j-i >= 11) goto 10
      if (i == 1) goto 5
      i = i - 1

      do
          i = i + 1
          if (i == j) goto 70
          rt = r(i+1)
          call pstring(t,a(rt))
          if (a(r(i)) <= t) cycle
          k = i
          do
              r(k+1) = r(k)
              k = k - 1
              if (t >= a(r(k))) exit
          enddo
          r(k+1) = rt
      enddo

      contains

!-------------------------------------------------------------------------------
      subroutine pstring(p,t)

      implicit none
      type(string), intent(inout)  :: p
      type(string), intent(in)     :: t


      p%len = t%len
      p%size = t%size
      p%chars => t%chars


      end subroutine pstring
!-------------------------------------------------------------------------------

      end subroutine rank_s

!*******************************************************************************
! Sorts A into ascending order, from A(1) to A(N).
! Reference: Richard C. Singleton, Algorithm 347, SORT.
! Comm. ACM 3, 321 (March 1969).
! Algorithm is Copyright 1969 Association of Computing Machinery,
!*******************************************************************************

      subroutine lrank_c(a,r)

      implicit none
      character(*), intent(in)  :: a(:)
      integer, intent(out)      :: r(size(a))
      character(len(a))         :: t
      integer                   :: i,j,k,l,m,n,p,rs,rt
      integer                   :: is(0:63)


      n = size(a)
      r(:) = (/ (i, i=1,n) /)
      m = 0
      i = 1
      j = n

    5 continue
      if (i >= j) goto 70

   10 continue
      p = (j+i)/2
      rt = r(p)
      t = a(rt)
      if (lgt(a(r(i)),t)) then
          r(p) = r(i)
          r(i) = rt
          rt = r(p)
          t = a(rt)
      endif
      if (llt(a(r(j)),t)) then
          r(p) = r(j)
          r(j) = rt
          rt = r(p)
          t = a(rt)
          if (llt(a(r(i)),t)) then
              r(p) = r(i)
              r(i) = rt
              rt = r(p)
              t = a(rt)
          endif
      endif

      k = i
      l = j
      do
          do
              l = l - 1
              if (lle(a(r(l)),t)) exit
          enddo
          rs = r(l)
          do
              k = k + 1
              if (lge(a(r(k)),t)) exit
          enddo
          if (k > l) exit
          r(l) = r(k)
          r(k) = rs
      enddo

      if (l-i > j-k) then
          is(m) = i
          m = m + 1
          is(m) = l
          m = m + 1
          i = k
      else
          is(m) = k
          m = m + 1
          is(m) = j
          m = m + 1
          j = l
      endif
      goto 80

   70 continue
      if (m == 0) return
      m = m - 1
      j = is(m)
      m = m - 1
      i = is(m)

   80 continue
      if (j-i >= 11) goto 10
      if (i == 1) goto 5
      i = i - 1

      do
          i = i + 1
          if (i == j) goto 70
          rt = r(i+1)
          t = a(rt)
          if (lle(a(r(i)),t)) cycle
          k = i
          do
              r(k+1) = r(k)
              k = k - 1
              if (lge(t,a(r(k)))) exit
          enddo
          r(k+1) = rt
      enddo

      end subroutine lrank_c

!*******************************************************************************
! Sorts A into ascending order, from A(1) to A(N).
! Reference: Richard C. Singleton, Algorithm 347, SORT.
! Comm. ACM 3, 321 (March 1969).
! Algorithm is Copyright 1969 Association of Computing Machinery,
!*******************************************************************************

      subroutine lrank_s(a,r)

      implicit none
      type(string), intent(in)  :: a(:)
      integer, intent(out)      :: r(size(a))
      type(string)              :: t
      integer                   :: i,j,k,l,m,n,p,rs,rt
      integer                   :: is(0:63)


      n = size(a)
      r(:) = (/ (i, i=1,n) /)
      m = 0
      i = 1
      j = n

    5 continue
      if (i >= j) goto 70

   10 continue
      p = (j+i)/2
      rt = r(p)
      call pstring(t,a(rt))
      if (lgt(a(r(i)),t)) then
          r(p) = r(i)
          r(i) = rt
          rt = r(p)
          call pstring(t,a(rt))
      endif
      if (llt(a(r(j)),t)) then
          r(p) = r(j)
          r(j) = rt
          rt = r(p)
          call pstring(t,a(rt))
          if (lgt(a(r(i)),t)) then
              r(p) = r(i)
              r(i) = rt
              rt = r(p)
              call pstring(t,a(rt))
          endif
      endif

      k = i
      l = j
      do
          do
              l = l - 1
              if (lle(a(r(l)),t)) exit
          enddo
          rs = r(l)
          do
              k = k + 1
              if (lge(a(r(k)),t)) exit
          enddo
          if (k > l) exit
          r(l) = r(k)
          r(k) = rs
      enddo

      if (l-i > j-k) then
          is(m) = i
          m = m + 1
          is(m) = l
          m = m + 1
          i = k
      else
          is(m) = k
          m = m + 1
          is(m) = j
          m = m + 1
          j = l
      endif
      goto 80

   70 continue
      if (m == 0) return
      m = m - 1
      j = is(m)
      m = m - 1
      i = is(m)

   80 continue
      if (j-i >= 11) goto 10
      if (i == 1) goto 5
      i = i - 1

      do
          i = i + 1
          if (i == j) goto 70
          rt = r(i+1)
          call pstring(t,a(rt))
          if (lle(a(r(i)),t)) cycle
          k = i
          do
              r(k+1) = r(k)
              k = k - 1
              if (lge(t,a(r(k)))) exit
          enddo
          r(k+1) = rt
      enddo

      contains

!-------------------------------------------------------------------------------
      subroutine pstring(p,t)

      implicit none
      type(string), intent(inout)  :: p
      type(string), intent(in)     :: t


      p%len = t%len
      p%size = t%size
      p%chars => t%chars


      end subroutine pstring
!-------------------------------------------------------------------------------

      end subroutine lrank_s

!*******************************************************************************
!     COMPARE, LCOMPARE, ACOMPARE, ALCOMPARE
!*******************************************************************************
!*******************************************************************************

      elemental function compare_ss(s1,s2) result(css)

      implicit none
      type(string), intent(in)  :: s1,s2
      character(2)              :: css
      integer                   :: i,l1,l2


      l1 = len(s1)
      l2 = len(s2)
      do i=1,min(l1,l2)
          if (s1%chars(i) < s2%chars(i)) then
              css = 'LT'
              return
          elseif (s1%chars(i) > s2%chars(i)) then
              css = 'GT'
              return
          endif
      enddo
      if (l1 < l2) then
          do i=l1+1,l2
              if (blank < s2%chars(i)) then
                  css = 'LT'
                  return
              elseif (blank > s2%chars(i)) then
                  css = 'GT'
                  return
              endif
          enddo
      elseif (l1 > l2) then
          do i=l2+1,l1
              if (s1%chars(i) < blank) then
                  css = 'LT'
                  return
              elseif (s1%chars(i) > blank) then
                  css = 'GT'
                  return
              endif
          enddo
      endif
      css = 'EQ'

      end function compare_ss

!*******************************************************************************

      elemental function compare_cs(c,s) result(css)

      implicit none
      character(*), intent(in)  :: c
      type(string), intent(in)  :: s
      character(2)              :: css
      integer                   :: i,lc,ls


      lc = len(c)
      ls = len(s)
      do i=1,min(lc,ls)
          if (c(i:i) < s%chars(i)) then
              css = 'LT'
              return
          elseif (c(i:i) > s%chars(i)) then
              css = 'GT'
              return
          endif
      enddo
      if (lc < ls) then
          do i=lc+1,ls
              if (blank < s%chars(i)) then
                  css = 'LT'
                  return
              elseif (blank > s%chars(i)) then
                  css = 'GT'
                  return
              endif
          enddo
      elseif (lc > ls) then
          do i=ls+1,lc
              if (c(i:i) < blank) then
                  css = 'LT'
                  return
              elseif (c(i:i) > blank) then
                  css = 'GT'
                  return
              endif
          enddo
      endif
      css = 'EQ'

      end function compare_cs

!*******************************************************************************
!     ==
!*******************************************************************************
! string == string

      elemental function s_eq_s(s1,s2)

      implicit none
      type(string), intent(in)  :: s1,s2
      logical                   :: s_eq_s
      integer                   :: l1,l2


      l1 = len(s1)
      l2 = len(s2)
      if (l1 > l2) then
          s_eq_s = all(s1%chars(1:l2) == s2%chars) .and.  &
                   all(s1%chars(l2+1:l1) == blank)
      elseif (l1 < l2) then
          s_eq_s = all(s1%chars == s2%chars(1:l1)) .and.  &
                   all(blank == s2%chars(l1+1:l2))
      else
          s_eq_s = all(s1%chars == s2%chars)
      endif

      end function s_eq_s

!*******************************************************************************
! string == character

      elemental function s_eq_c(s,c)

      implicit none
      type(string), intent(in)  :: s
      character(*), intent(in)  :: c
      logical                   :: s_eq_c
      integer                   :: i,ls,lc


      ls = len(s)
      lc = len(c)
      do i=1,min(ls,lc)
          if (s%chars(i) /= c(i:i)) then
              s_eq_c = .false.
              return
          endif
      enddo
      if ((ls > lc) .and. any(s%chars(lc+1:ls) /= blank)) then
          s_eq_c = .false.
      elseif ((ls < lc) .and. (blank /= c(ls+1:lc))) then
          s_eq_c = .false.
      else
          s_eq_c = .true.
      endif

      end function s_eq_c

!*******************************************************************************
! character == string

      elemental function c_eq_s(c,s)

      implicit none
      character(*), intent(in)  :: c
      type(string), intent(in)  :: s
      logical                   :: c_eq_s
      integer                   :: i,lc,ls


      lc = len(c)
      ls = len(s)
      do i=1,min(lc,ls)
          if (c(i:i) /= s%chars(i)) then
              c_eq_s = .false.
              return
          endif
      enddo
      if ((lc > ls) .and. (c(ls+1:lc) /= blank)) then
          c_eq_s = .false.
      elseif ((lc < ls) .and. any(blank /= s%chars(lc+1:ls) ) )then
          c_eq_s = .false.
      else
          c_eq_s = .true.
      endif

      end function c_eq_s

!*******************************************************************************
!     /=
!*******************************************************************************
! string /= string

      elemental function s_ne_s(s1,s2)

      implicit none
      type(string), intent(in)  :: s1,s2
      logical                   :: s_ne_s
      integer                   :: l1,l2


      l1 = len(s1)
      l2 = len(s2)
      if (l1 > l2) then
          s_ne_s = any(s1%chars(1:l2) /= s2%chars) .or.  &
                   any(s1%chars(l2+1:l1) /= blank)
      elseif (l1 < l2) then
          s_ne_s = any(s1%chars /= s2%chars(1:l1)) .or. &
                   any(blank /= s2%chars(l1+1:l2))
      else
          s_ne_s = any(s1%chars /= s2%chars)
      endif

      end function s_ne_s

!*******************************************************************************
! string /= character

      elemental function s_ne_c(s,c)

      implicit none
      type(string), intent(in)  :: s
      character(*), intent(in)  :: c
      logical                   :: s_ne_c
      integer                   :: i,ls,lc


      ls = len(s)
      lc = len(c)
      do i=1,min(ls,lc)
          if (s%chars(i) /= c(i:i) )then
              s_ne_c = .true.
              return
          endif
      enddo
      if ((ls > lc) .and. any(s%chars(ls+1:lc) /= blank)) then
          s_ne_c = .true.
      elseif ((ls < lc) .and. blank /= c(ls+1:lc)) then
          s_ne_c = .true.
      else
          s_ne_c = .false.
      endif

      end function s_ne_c

!*******************************************************************************
! character /= string

      elemental function c_ne_s(c,s)

      implicit none
      character(*), intent(in)  :: c
      type(string), intent(in)  :: s
      logical                   :: c_ne_s
      integer                   :: i,lc,ls


      lc = len(c)
      ls = len(s)
      do i=1,min(lc,ls)
          if (c(i:i) /= s%chars(i)) then
              c_ne_s = .true.
              return
          endif
      enddo
      if ((lc > ls) .and. c(ls+1:lc) /= blank) then
          c_ne_s = .true.
      elseif ((lc < ls) .and. any(blank /= s%chars(lc+1:ls))) then
          c_ne_s = .true.
      else
          c_ne_s = .false.
      endif

      end function c_ne_s

!*******************************************************************************
!     < operators
!*******************************************************************************
! string < string

      elemental function s_lt_s(s1,s2)

      implicit none
      type(string), intent(in)  :: s1,s2
      logical                   :: s_lt_s


      s_lt_s = compare_ss(s1,s2) == 'LT'

      end function s_lt_s

!*******************************************************************************
! string < character

      elemental function s_lt_c(s,c)

      implicit none
      type(string), intent(in)  :: s
      character(*), intent(in)  :: c
      logical                   :: s_lt_c


      s_lt_c = compare_cs(c,s) == 'GT'

      end function s_lt_c

!*******************************************************************************
! character < string

      elemental function c_lt_s(c,s)

      implicit none
      character(*), intent(in)  :: c
      type(string), intent(in)  :: s
      logical                   :: c_lt_s


      c_lt_s = compare_cs(c,s) == 'LT'

      end function c_lt_s

!*******************************************************************************
!     <=  operators
!*******************************************************************************
! string <= string

      elemental function s_le_s(s1,s2)

      implicit none
      type(string), intent(in)  :: s1,s2
      logical                   :: s_le_s


      s_le_s = compare_ss(s1,s2) /= 'GT'

      end function s_le_s

!*******************************************************************************
! string <= character

      elemental function s_le_c(s,c)

      implicit none
      type(string), intent(in)  :: s
      character(*), intent(in)  :: c
      logical                   :: s_le_c


      s_le_c = compare_cs(c,s) /= 'LT'

      end function s_le_c

!*******************************************************************************
! character <= string

      elemental function c_le_s(c,s)

      implicit none
      character(*), intent(in)  :: c
      type(string), intent(in)  :: s
      logical                   :: c_le_s


      c_le_s = compare_cs(c,s) /= 'GT'

      end function c_le_s

!*******************************************************************************
!     >=  operators
!*******************************************************************************
! string >= string

      elemental function s_ge_s(s1,s2)

      implicit none
      type(string), intent(in) :: s1,s2
      logical                  :: s_ge_s


      s_ge_s = compare_ss(s1,s2) /= 'LT'

      end function s_ge_s

!*******************************************************************************
! string >= character

      elemental function s_ge_c(s,c)

      implicit none
      type(string), intent(in)  :: s
      character(*), intent(in)  :: c
      logical                   :: s_ge_c


      s_ge_c = compare_cs(c,s) /= 'GT'

      end function s_ge_c

!*******************************************************************************
! character >= string

      elemental function c_ge_s(c,s)

      implicit none
      character(*), intent(in)  :: c
      type(string), intent(in)  :: s
      logical                   :: c_ge_s


      c_ge_s = compare_cs(c,s) /= 'LT'

      end function c_ge_s

!*******************************************************************************
!     >  operators
!*******************************************************************************
! string > string

      elemental function s_gt_s(s1,s2)

      implicit none
      type(string), intent(in) :: s1,s2
      logical                  :: s_gt_s


      s_gt_s = compare_ss(s1,s2) == 'GT'

      end function s_gt_s

!*******************************************************************************
! string > character

      elemental function s_gt_c(s,c)

      implicit none
      type(string), intent(in)  :: s
      character(*), intent(in)  :: c
      logical                   :: s_gt_c


      s_gt_c = compare_cs(c,s) == 'LT'

      end function s_gt_c

!*******************************************************************************
! character > string

      elemental function c_gt_s(c,s)

      implicit none
      character(*), intent(in)  :: c
      type(string), intent(in)  :: s
      logical                   :: c_gt_s


      c_gt_s = compare_cs(c,s) == 'GT'

      end function c_gt_s

!*******************************************************************************

      elemental function lcompare_ss(s1,s2) result(css)

      implicit none
      type(string), intent(in)  :: s1,s2
      character(2)              :: css
      integer                   :: i,l1,l2


      l1 = len(s1)
      l2 = len(s2)
      do i=1,min(l1,l2)
          if (llt(s1%chars(i),s2%chars(i))) then
              css = 'LT'
              return
          elseif (lgt(s1%chars(i),s2%chars(i))) then
              css = 'GT'
              return
          endif
      enddo
      if (l1 < l2) then
          do i=l1+1,l2
              if (llt(blank,s2%chars(i))) then
                  css = 'LT'
                  return
              elseif (lgt(blank,s2%chars(i))) then
                  css = 'GT'
                  return
              endif
          enddo
      elseif (l1 > l2) then
          do i=l2+1,l1
              if (llt(s1%chars(i),blank)) then
                  css = 'LT'
                  return
              elseif (lgt(s1%chars(i),blank)) then
                  css = 'GT'
                  return
              endif
          enddo
      endif
      css = 'EQ'

      end function lcompare_ss

!*******************************************************************************

      elemental function lcompare_cs(c,s) result(css)

      implicit none
      character(*), intent(in)  :: c
      type(string), intent(in)  :: s
      character(2)              :: css
      integer                   :: i,lc,ls


      lc = len(c)
      ls = len(s)
      do i=1,min(lc,ls)
          if (llt(c(i:i),s%chars(i))) then
              css = 'LT'
              return
          elseif (lgt(c(i:i),s%chars(i))) then
              css = 'GT'
              return
          endif
      enddo
      if (lc < ls) then
          do i=lc+1,ls
              if (llt(blank,s%chars(i))) then
                  css = 'LT'
                  return
              elseif (lgt(blank,s%chars(i))) then
                  css = 'GT'
                  return
              endif
          enddo
      elseif (lc > ls) then
          do i=ls+1,lc
              if (llt(c(i:i),blank)) then
                  css = 'LT'
                  return
              elseif (lgt(c(i:i),blank)) then
                  css = 'GT'
                  return
              endif
          enddo
      endif
      css = 'EQ'

      end function lcompare_cs

!*******************************************************************************
!     LLT function
!*******************************************************************************
!     llt(string,string)

      elemental function s_llt_s(s1,s2)

      implicit none
      type(string), intent(in)  :: s1,s2
      logical                   :: s_llt_s

      s_llt_s = (lcompare_ss(s1,s2) == 'LT')

      end function s_llt_s

!*******************************************************************************
!     llt(string,character)

      elemental function s_llt_c(s1,c2)

      implicit none
      type(string), intent(in)  :: s1
      character(*), intent(in)  :: c2
      logical                   :: s_llt_c

      s_llt_c = (lcompare_cs(c2,s1) == 'GT')

      end function s_llt_c

!*******************************************************************************
!     llt(character,string)

      elemental function c_llt_s(c1,s2)

      implicit none
      type(string), intent(in)  :: s2
      character(*), intent(in)  :: c1
      logical                   :: c_llt_s

      c_llt_s = (lcompare_cs(c1,s2) == 'LT')

      end function c_llt_s

!*******************************************************************************
!     LGT function
!*******************************************************************************
!     lgt(string,string)

      elemental function s_lgt_s(s1,s2)

      implicit none
      type(string), intent(in)  :: s1,s2
      logical                   :: s_lgt_s

      s_lgt_s = (lcompare_ss(s1,s2) == 'GT')

      end function s_lgt_s

!*******************************************************************************
!     lgt(string,character)

      elemental function s_lgt_c(s1,c2)

      implicit none
      type(string), intent(in)  :: s1
      character(*), intent(in)  :: c2
      logical                   :: s_lgt_c

      s_lgt_c = (lcompare_cs(c2,s1) == 'LT')

      end function s_lgt_c

!*******************************************************************************
!     lgt(character,string)

      elemental function c_lgt_s(c1,s2)

      implicit none
      type(string), intent(in)  :: s2
      character(*), intent(in)  :: c1
      logical                   :: c_lgt_s

      c_lgt_s = (lcompare_cs(c1,s2) == 'GT')

      end function c_lgt_s

!*******************************************************************************
!     LGE function
!*******************************************************************************
!     lge(string,string)

      elemental function s_lge_s(s1,s2)

      implicit none
      type(string), intent(in)  :: s1,s2
      logical                   :: s_lge_s

      s_lge_s = (lcompare_ss(s1,s2) /= 'LT')

      end function s_lge_s

!*******************************************************************************
!     lge(string,character)

      elemental function s_lge_c(s1,c2)

      implicit none
      type(string), intent(in)  :: s1
      character(*), intent(in)  :: c2
      logical                   :: s_lge_c

      s_lge_c = (lcompare_cs(c2,s1) /= 'GT')

      end function s_lge_c

!*******************************************************************************
!     lge(character,string)

      elemental function c_lge_s(c1,s2)

      implicit none
      type(string), intent(in)  :: s2
      character(*), intent(in)  :: c1
      logical                   :: c_lge_s

      c_lge_s = (lcompare_cs(c1,s2) /= 'LT')

      end function c_lge_s

!*******************************************************************************
!     LLE function
!*******************************************************************************
!     lle(string,string)

      elemental function s_lle_s(s1,s2)

      implicit none
      type(string), intent(in)  :: s1,s2
      logical                   :: s_lle_s

      s_lle_s = (lcompare_ss(s1,s2) /= 'GT')

      end function s_lle_s

!*******************************************************************************
!     lle(string,character)

      elemental function s_lle_c(s1,c2)

      implicit none
      type(string), intent(in)  :: s1
      character(*), intent(in)  :: c2
      logical                   :: s_lle_c

      s_lle_c = (lcompare_cs(c2,s1) /= 'LT')

      end function s_lle_c

!*******************************************************************************
!     lle(character,string)

      elemental function c_lle_s(c1,s2)

      implicit none
      type(string), intent(in)  :: s2
      character(*), intent(in)  :: c1
      logical                   :: c_lle_s

      c_lle_s = (lcompare_cs(c1,s2) /= 'GT')

      end function c_lle_s

!*******************************************************************************

      pure function acompare_aa(a1,a2) result(caa)

      implicit none
      character, intent(in)  :: a1(:),a2(:)
      character(2)           :: caa
      integer                :: i,l1,l2


      l1 = size(a1)
      l2 = size(a2)
      do i=1,min(l1,l2)
          if (a1(i) < a2(i)) then
              caa = 'LT'
              return
          elseif (a1(i) > a2(i)) then
              caa = 'GT'
              return
          endif
      enddo
      if (l1 < l2) then
          do i=l1+1,l2
              if (blank < a2(i)) then
                  caa = 'LT'
                  return
              elseif (blank > a2(i)) then
                  caa = 'GT'
                  return
              endif
          enddo
      elseif (l1 > l2) then
          do i=l2+1,l1
              if (a1(i) < blank) then
                  caa = 'LT'
                  return
              elseif (a1(i) > blank) then
                  caa = 'GT'
                  return
              endif
          enddo
      endif
      caa = 'EQ'

      end function acompare_aa

!*******************************************************************************

      pure function acompare_ca(c,a) result(cca)

      implicit none
      character(*), intent(in)  :: c
      character, intent(in)     :: a(:)
      character(2)              :: cca
      integer                   :: i,lc,la


      lc = len(c)
      la = size(a)
      do i=1,min(lc,la)
          if (c(i:i) < a(i)) then
              cca = 'LT'
              return
          elseif (c(i:i) > a(i)) then
              cca = 'GT'
              return
          endif
      enddo
      if (lc < la) then
          do i=lc+1,la
              if (blank < a(i)) then
                  cca = 'LT'
                  return
              elseif (blank > a(i)) then
                  cca = 'GT'
                  return
              endif
          enddo
      elseif (lc > la) then
          do i=la+1,lc
              if (c(i:i) < blank) then
                  cca = 'LT'
                  return
              elseif (c(i:i) > blank) then
                  cca = 'GT'
                  return
              endif
          enddo
      endif
      cca = 'EQ'

      end function acompare_ca

!*******************************************************************************
!     ==
!*******************************************************************************
! array == array

      pure function a_eq_a(a1,a2)

      implicit none
      character, intent(in)  :: a1(:),a2(:)
      logical                :: a_eq_a
      integer                :: l1,l2


      l1 = size(a1)
      l2 = size(a2)
      if (l1 > l2) then
          a_eq_a = all(a1(1:l2) == a2) .and.  &
                   all(a1(l2+1:l1) == blank)
      elseif (l1 < l2) then
          a_eq_a = all(a1 == a2(1:l1)) .and.  &
                   all(blank == a2(l1+1:l2))
      else
          a_eq_a = all(a1 == a2)
      endif

      end function a_eq_a

!*******************************************************************************
! array == character

      pure function a_eq_c(a,c)

      implicit none
      character, intent(in)     :: a(:)
      character(*), intent(in)  :: c
      logical                   :: a_eq_c
      integer                   :: i,la,lc


      la = len(a)
      lc = len(c)
      do i=1,min(la,lc)
          if (a(i) /= c(i:i)) then
              a_eq_c = .false.
              return
          endif
      enddo
      if ((la > lc) .and. any(a(lc+1:la) /= blank)) then
          a_eq_c = .false.
      elseif ((la < lc) .and. (blank /= c(la+1:lc))) then
          a_eq_c = .false.
      else
          a_eq_c = .true.
      endif

      end function a_eq_c

!*******************************************************************************
! character == array

      pure function c_eq_a(c,a)

      implicit none
      character(*), intent(in)  :: c
      character, intent(in)     :: a(:)
      logical                   :: c_eq_a


      c_eq_a = a_eq_c(a,c)

      end function c_eq_a

!*******************************************************************************
!     /=
!*******************************************************************************
! array /= array

      pure function a_ne_a(a1,a2)

      implicit none
      character, intent(in)  :: a1(:),a2(:)
      logical                :: a_ne_a
      integer                :: l1,l2


      l1 = size(a1)
      l2 = size(a2)
      if (l1 > l2) then
          a_ne_a = any(a1(1:l2) /= a2) .or.  &
                   any(a1(l2+1:l1) /= blank)
      elseif (l1 < l2) then
          a_ne_a = any(a1 /= a2(1:l1)) .or. &
                   any(blank /= a2(l1+1:l2))
      else
          a_ne_a = any(a1 /= a2)
      endif

      end function a_ne_a

!*******************************************************************************
! array /= character

      pure function a_ne_c(a,c)

      implicit none
      character, intent(in)     :: a(:)
      character(*), intent(in)  :: c
      logical                   :: a_ne_c
      integer                   :: i,la,lc


      la = size(a)
      lc = len(c)
      do i=1,min(la,lc)
          if (a(i) /= c(i:i) )then
              a_ne_c = .true.
              return
          endif
      enddo
      if ((la > lc) .and. any(a(la+1:lc) /= blank)) then
          a_ne_c = .true.
      elseif ((la < lc) .and. blank /= c(la+1:lc)) then
          a_ne_c = .true.
      else
          a_ne_c = .false.
      endif

      end function a_ne_c

!*******************************************************************************
! character /= array

      pure function c_ne_a(c,a)

      implicit none
      character(*), intent(in)  :: c
      character, intent(in)     :: a(:)
      logical                   :: c_ne_a


      c_ne_a = acompare_ca(c,a) /= 'EQ'

      end function c_ne_a

!*******************************************************************************
!     < operators
!*******************************************************************************
! array < array

      pure function a_lt_a(a1,a2)

      implicit none
      character, intent(in)  :: a1(:),a2(:)
      logical                :: a_lt_a


      a_lt_a = acompare_aa(a1,a2) == 'LT'

      end function a_lt_a

!*******************************************************************************
! array < character

      pure function a_lt_c(a,c)

      implicit none
      character, intent(in)     :: a(:)
      character(*), intent(in)  :: c
      logical                   :: a_lt_c


      a_lt_c = acompare_ca(c,a) == 'GT'

      end function a_lt_c

!*******************************************************************************
! character < array

      pure function c_lt_a(c,a)

      implicit none
      character(*), intent(in)  :: c
      character, intent(in)     :: a(:)
      logical                   :: c_lt_a


      c_lt_a = acompare_ca(c,a) == 'LT'

      end function c_lt_a

!*******************************************************************************
!     <=  operators
!*******************************************************************************
! array <= array

      pure function a_le_a(a1,a2)

      implicit none
      character, intent(in)  :: a1(:),a2(:)
      logical                :: a_le_a


      a_le_a = acompare_aa(a1,a2) /= 'GT'

      end function a_le_a

!*******************************************************************************
! array <= character

      pure function a_le_c(a,c)

      implicit none
      character, intent(in)     :: a(:)
      character(*), intent(in)  :: c
      logical                   :: a_le_c


      a_le_c = acompare_ca(c,a) /= 'LT'

      end function a_le_c

!*******************************************************************************
! character <= array

      pure function c_le_a(c,a)

      implicit none
      character(*), intent(in)  :: c
      character, intent(in)     :: a(:)
      logical                   :: c_le_a


      c_le_a = acompare_ca(c,a) /= 'GT'

      end function c_le_a

!*******************************************************************************
!     >=  operators
!*******************************************************************************
! array >= array

      pure function a_ge_a(a1,a2)

      implicit none
      character, intent(in)  :: a1(:),a2(:)
      logical                :: a_ge_a


      a_ge_a = acompare_aa(a1,a2) /= 'LT'

      end function a_ge_a

!*******************************************************************************
! array >= character

      pure function a_ge_c(a,c)

      implicit none
      character, intent(in)     :: a(:)
      character(*), intent(in)  :: c
      logical                   :: a_ge_c


      a_ge_c = acompare_ca(c,a) /= 'GT'

      end function a_ge_c

!*******************************************************************************
! character >= array

      pure function c_ge_a(c,a)

      implicit none
      character(*), intent(in)  :: c
      character, intent(in)     :: a(:)
      logical                   :: c_ge_a


      c_ge_a = acompare_ca(c,a) /= 'LT'

      end function c_ge_a

!*******************************************************************************
!     >  operators
!*******************************************************************************
! array > array

      pure function a_gt_a(a1,a2)

      implicit none
      character, intent(in)  :: a1(:),a2(:)
      logical                :: a_gt_a


      a_gt_a = acompare_aa(a1,a2) == 'GT'

      end function a_gt_a

!*******************************************************************************
! array > character

      pure function a_gt_c(a,c)

      implicit none
      character, intent(in)     :: a(:)
      character(*), intent(in)  :: c
      logical                   :: a_gt_c


      a_gt_c = acompare_ca(c,a) == 'LT'

      end function a_gt_c

!*******************************************************************************
! character > array

      pure function c_gt_a(c,a)

      implicit none
      character(*), intent(in)  :: c
      character, intent(in)     :: a(:)
      logical                   :: c_gt_a


      c_gt_a = acompare_ca(c,a) == 'GT'

      end function c_gt_a

!*******************************************************************************

      pure function alcompare_aa(a1,a2) result(caa)

      implicit none
      character, intent(in)  :: a1(:),a2(:)
      character(2)           :: caa
      integer                :: i,l1,l2


      l1 = size(a1)
      l2 = size(a2)
      do i=1,min(l1,l2)
          if (llt(a1(i),a2(i))) then
              caa = 'LT'
              return
          elseif (lgt(a1(i),a2(i))) then
              caa = 'GT'
              return
          endif
      enddo
      if (l1 < l2) then
          do i=l1+1,l2
              if (llt(blank,a2(i))) then
                  caa = 'LT'
                  return
              elseif (lgt(blank,a2(i))) then
                  caa = 'GT'
                  return
              endif
          enddo
      elseif (l1 > l2) then
          do i=l2+1,l1
              if (llt(a1(i),blank)) then
                  caa = 'LT'
                  return
              elseif (lgt(a1(i),blank)) then
                  caa = 'GT'
                  return
              endif
          enddo
      endif
      caa = 'EQ'

      end function alcompare_aa

!*******************************************************************************

      pure function alcompare_ca(c,a) result(cca)

      implicit none
      character(*), intent(in)  :: c
      character, intent(in)     :: a(:)
      character(2)              :: cca
      integer                   :: i,lc,la


      lc = len(c)
      la = size(a)
      do i=1,min(lc,la)
          if (llt(c(i:i),a(i))) then
              cca = 'LT'
              return
          elseif (lgt(c(i:i),a(i))) then
              cca = 'GT'
              return
          endif
      enddo
      if (lc < la) then
          do i=lc+1,la
              if (llt(blank,a(i))) then
                  cca = 'LT'
                  return
              elseif (lgt(blank,a(i))) then
                  cca = 'GT'
                  return
              endif
          enddo
      elseif (lc > la) then
          do i=la+1,lc
              if (llt(c(i:i),blank)) then
                  cca = 'LT'
                  return
              elseif (lgt(c(i:i),blank)) then
                  cca = 'GT'
                  return
              endif
          enddo
      endif
      cca = 'EQ'

      end function alcompare_ca

!*******************************************************************************
!     LLT operators
!*******************************************************************************
! array < array

      pure function a_allt_a(a1,a2)

      implicit none
      character, intent(in)  :: a1(:),a2(:)
      logical                :: a_allt_a


      a_allt_a = alcompare_aa(a1,a2) == 'LT'

      end function a_allt_a

!*******************************************************************************
! array < character

      pure function a_allt_c(a1,c2)

      implicit none
      character, intent(in)     :: a1(:)
      character(*), intent(in)  :: c2
      logical                   :: a_allt_c


      a_allt_c = alcompare_ca(c2,a1) == 'GT'

      end function a_allt_c

!*******************************************************************************
! character < array

      pure function c_allt_a(c1,a2)

      implicit none
      character(*), intent(in)  :: c1
      character, intent(in)     :: a2(:)
      logical                   :: c_allt_a


      c_allt_a = alcompare_ca(c1,a2) == 'LT'

      end function c_allt_a

!*******************************************************************************
!     LLE  operators
!*******************************************************************************
! array <= array

      pure function a_alle_a(a1,a2)

      implicit none
      character, intent(in)  :: a1(:),a2(:)
      logical                :: a_alle_a


      a_alle_a = alcompare_aa(a1,a2) /= 'GT'

      end function a_alle_a

!*******************************************************************************
! array <= character

      pure function a_alle_c(a1,c2)

      implicit none
      character, intent(in)     :: a1(:)
      character(*), intent(in)  :: c2
      logical                   :: a_alle_c


      a_alle_c = alcompare_ca(c2,a1) /= 'LT'

      end function a_alle_c

!*******************************************************************************
! character <= array

      pure function c_alle_a(c1,a2)

      implicit none
      character(*), intent(in)  :: c1
      character, intent(in)     :: a2(:)
      logical                   :: c_alle_a


      c_alle_a = alcompare_ca(c1,a2) /= 'GT'

      end function c_alle_a

!*******************************************************************************
!     LGE  operators
!*******************************************************************************
! array >= array

      pure function a_alge_a(a1,a2)

      implicit none
      character, intent(in)  :: a1(:),a2(:)
      logical                :: a_alge_a


      a_alge_a = alcompare_aa(a1,a2) /= 'LT'

      end function a_alge_a

!*******************************************************************************
! array >= character

      pure function a_alge_c(a1,c2)

      implicit none
      character, intent(in)     :: a1(:)
      character(*), intent(in)  :: c2
      logical                   :: a_alge_c


      a_alge_c = alcompare_ca(c2,a1) /= 'GT'

      end function a_alge_c

!*******************************************************************************
! character >= array

      pure function c_alge_a(c1,a2)

      implicit none
      character(*), intent(in)  :: c1
      character, intent(in)     :: a2(:)
      logical                   :: c_alge_a


      c_alge_a = alcompare_ca(c1,a2) /= 'LT'

      end function c_alge_a

!*******************************************************************************
!     LGT  operators
!*******************************************************************************
! array > array

      pure function a_algt_a(a1,a2)

      implicit none
      character, intent(in)  :: a1(:),a2(:)
      logical                :: a_algt_a


      a_algt_a = alcompare_aa(a1,a2) == 'GT'

      end function a_algt_a

!*******************************************************************************
! array > character

      pure function a_algt_c(a1,c2)

      implicit none
      character, intent(in)     :: a1(:)
      character(*), intent(in)  :: c2
      logical                   :: a_algt_c


      a_algt_c = alcompare_ca(c2,a1) == 'LT'

      end function a_algt_c

!*******************************************************************************
! character > array

      pure function c_algt_a(c1,a2)

      implicit none
      character(*), intent(in)  :: c1
      character, intent(in)     :: a2(:)
      logical                   :: c_algt_a


      c_algt_a = alcompare_ca(c1,a2) == 'GT'

      end function c_algt_a

!*******************************************************************************
!     INDEX
!*******************************************************************************

      elemental function index_ss(s,sub,back)

      implicit none
      type(string), intent(in)       :: s,sub
      logical, intent(in), optional  :: back
      integer                        :: index_ss
      logical                        :: dir_switch
      integer                        :: ls,lsub


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = len(s)
      lsub = len(sub)
      index_ss = aindex(s%chars(:ls),sub%chars(:lsub),dir_switch)

      end function index_ss

!*******************************************************************************

      elemental function index_sc(s,sub,back)

      implicit none
      type(string), intent(in)       :: s
      character(*), intent(in)       :: sub
      logical, intent(in), optional  :: back
      integer                        :: index_sc
      logical                        :: dir_switch
      integer                        :: ls


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = len(s)
      index_sc = aindex(s%chars(:ls),sub,dir_switch)

      end function index_sc

!*******************************************************************************

      elemental function index_cs(s,sub,back)

      implicit none
      character(*), intent(in)       :: s
      type(string), intent(in)       :: sub
      logical, intent(in), optional  :: back
      integer                        :: index_cs
      logical                        :: dir_switch


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      index_cs = index(s,char(sub),dir_switch)

      end function index_cs

!*******************************************************************************
!     AINDEX
!*******************************************************************************

      pure function aindex_aa(s,sub,back) result(index_aa)

      implicit none
      character, intent(in)          :: s(:),sub(:)
      logical, intent(in), optional  :: back
      integer                        :: index_aa
      logical                        :: dir_switch
      integer                        :: i,ls,lss


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = size(s)
      lss = size(sub)

      if (lss == 0) then
          if (dir_switch) then
              index_aa = ls + 1
          else
              index_aa = 1
          endif
          return
      endif

      if (dir_switch) then
!         backwards search
          do i=ls-lss+1,1,-1
              if (all(s(i:i+lss-1) == sub(1:lss))) then
                  index_aa = i
                  return
              endif
          enddo
          index_aa = 0
      else
!         forward search
          do i=1,ls-lss+1
              if (all(s(i:i+lss-1) == sub(1:lss))) then
                  index_aa = i
                  return
              endif
          enddo
          index_aa = 0
      endif

      end function aindex_aa

!*******************************************************************************

      pure function aindex_ac(s,sub,back) result(index_ac)

      implicit none
      character, intent(in)          :: s(:)
      character(*), intent(in)       :: sub
      logical, intent(in), optional  :: back
      integer                        :: index_ac
      logical                        :: dir_switch,matched
      integer                        :: i,j,ls,lss


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = size(s)
      lss = len(sub)

      if (lss == 0) then
          if (dir_switch) then
              index_ac = ls + 1
          else
              index_ac = 1
          endif
          return
      endif

      if (dir_switch) then
          index_ac = 0
          do i=ls-lss+1,1,-1
              matched = all(s(i:i+lss-1) == (/ (sub(j:j), j=1,lss) /))
              if (matched) then
                  index_ac = i
                  return
              endif
          enddo
      else
          index_ac = 0
          do i=1,ls-lss+1
              matched = all(s(i:i+lss-1) == (/ (sub(j:j), j=1,lss) /))
              if (matched) then
                  index_ac = i
                  return
              endif
          enddo
      endif

      end function aindex_ac

!*******************************************************************************

      pure function aindex_ca(s,sub,back) result(index_ca)

      implicit none
      character(*), intent(in)       :: s
      character, intent(in)          :: sub(:)
      logical, intent(in), optional  :: back
      integer                        :: index_ca
      logical                        :: dir_switch,matched
      integer                        :: i,j,ls,lss


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = len(s)
      lss = size(sub)

      if (lss == 0) then
          if (dir_switch) then
              index_ca = ls + 1
          else
              index_ca = 1
          endif
          return
      endif

      if (dir_switch) then
          do i=ls-lss+1,1,-1
              matched = .true.
              do j=1,lss
                  if (s(i+j-1:i+j-1) /= sub(j)) then
                      matched = .false.
                      exit
                  endif
              enddo
              if (matched) then
                  index_ca = i
                  return
              endif
          enddo
          index_ca = 0
      else
          do i=1,ls-lss+1
              matched = .true.
              do j=1,lss
                  if (s(i+j-1:i+j-1) /= sub(j)) then
                      matched = .false.
                      exit
                  endif
              enddo
              if (matched) then
                  index_ca = i
                  return
              endif
          enddo
          index_ca = 0
      endif

      end function aindex_ca

!*******************************************************************************
!     SCAN
!*******************************************************************************

      elemental function scan_ss(s,set,back)

      implicit none
      type(string), intent(in)       :: s,set
      logical, intent(in), optional  :: back
      integer                        :: scan_ss
      logical                        :: dir_switch
      integer                        :: ls,lset


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = len(s)
      lset = len(set)
      scan_ss = ascan_aa(s%chars(1:ls),set%chars(1:lset),dir_switch)

      end function scan_ss

!*******************************************************************************

      elemental function scan_sc(s,set,back)

      implicit none
      type(string), intent(in)       :: s
      character(*), intent(in)       :: set
      logical, intent(in), optional  :: back
      integer                        :: scan_sc
      logical                        :: dir_switch
      integer                        :: ls


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = len(s)
      scan_sc = ascan_ac(s%chars(1:ls),set,dir_switch)

      end function scan_sc

!*******************************************************************************

      elemental function scan_cs(s,set,back)

      implicit none
      character(*), intent(in)       :: s
      type(string), intent(in)       :: set
      logical, intent(in), optional  :: back
      integer                        :: scan_cs
      logical                        :: dir_switch
      integer                        :: lset


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      lset = len(set)
      scan_cs = ascan_ca(s,set%chars(1:lset),dir_switch)

      end function scan_cs
!*******************************************************************************
!     ASCAN
!*******************************************************************************

      pure function ascan_aa(s,set,back)

      implicit none
      character, intent(in)          :: s(:),set(:)
      logical, intent(in), optional  :: back
      integer                        :: ascan_aa
      logical                        :: dir_switch
      integer                        :: i,ls,lset


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = size(s)
      lset = size(set)
      if (dir_switch) then
!         backwards search
          do i=ls,1,-1
              if (any(set(1:lset) == s(i))) then
                  ascan_aa = i
                  return
              endif
          enddo
          ascan_aa = 0
      else
!         forward search
          do i=1,ls
              if (any(set(1:lset) == s(i))) then
                  ascan_aa = i
                  return
              endif
          enddo
          ascan_aa = 0
      endif

      end function ascan_aa

!*******************************************************************************

      pure function ascan_ac(s,set,back)

      implicit none
      character, intent(in)          :: s(:)
      character(*), intent(in)       :: set
      logical, intent(in), optional  :: back
      integer                        :: ascan_ac
      logical                        :: dir_switch,matched
      integer                        :: i,j,ls,lset


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = size(s)
      lset = len(set)
      if (dir_switch) then
!         backwards search
          do i=ls,1,-1
              matched = .false.
              do j=1,lset
                  if (s(i) == set(j:j)) then
                      matched = .true.
                      exit
                  endif
              enddo
              if (matched) then
                  ascan_ac = i
                  return
              endif
          enddo
          ascan_ac = 0
      else
!         forward search
          do i=1,ls
              matched = .false.
              do j=1,lset
                  if (s(i) == set(j:j)) then
                      matched = .true.
                      exit
                  endif
              enddo
              if (matched) then
                  ascan_ac = i
                  return
              endif
          enddo
          ascan_ac = 0
      endif

      end function ascan_ac

!*******************************************************************************

      pure function ascan_ca(s,set,back)

      implicit none
      character(*), intent(in)       :: s
      character, intent(in)          :: set(:)
      logical, intent(in), optional  :: back
      integer                        :: ascan_ca
      logical                        :: dir_switch,matched
      integer                        :: i,j,ls,lset


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = len(s)
      lset = size(set)
      if (dir_switch) then
!         backwards search
          do i=ls,1,-1
              matched = .false.
              do j=1,lset
                  if (s(i:i) == set(j)) then
                      matched = .true.
                      exit
                  endif
              enddo
              if (matched) then
                  ascan_ca = i
                  return
              endif
          enddo
          ascan_ca = 0
      else
!         forward search
          do i=1,ls
              matched = .false.
              do j=1,lset
                  if (s(i:i) == set(j)) then
                      matched = .true.
                      exit
                  endif
              enddo
              if (matched) then
                  ascan_ca = i
                  return
              endif
          enddo
          ascan_ca = 0
      endif

      end function ascan_ca

!*******************************************************************************
!     VERIFY
!*******************************************************************************

      elemental function verify_ss(s,set,back)

      implicit none
      type(string), intent(in)       :: s,set
      logical, intent(in), optional  :: back
      integer                        :: verify_ss
      logical                        :: dir_switch
      integer                        :: ls,lset


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = len(s)
      lset = len(set)
      verify_ss = averify_aa(s%chars(1:ls),set%chars(1:lset),dir_switch)

      end function verify_ss

!*******************************************************************************

      elemental function verify_sc(s,set,back)

      implicit none
      type(string), intent(in)       :: s
      character(*), intent(in)       :: set
      logical, intent(in), optional  :: back
      integer                        :: verify_sc
      logical                        :: dir_switch
      integer                        :: ls


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = len(s)
      verify_sc = averify_ac(s%chars(1:ls),set,dir_switch)

      end function verify_sc

!*******************************************************************************

      elemental function verify_cs(s,set,back)

      implicit none
      character(*), intent(in)       :: s
      type(string), intent(in)       :: set
      logical, intent(in), optional  :: back
      integer                        :: verify_cs
      logical                        :: dir_switch
      integer                        :: lset


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      lset = len(set)
      verify_cs = averify_ca(s,set%chars(1:lset),dir_switch)

      end function verify_cs

!*******************************************************************************
!     AVERIFY
!*******************************************************************************

      pure function averify_aa(s,set,back)

      implicit none
      character, intent(in)          :: s(:),set(:)
      logical, intent(in), optional  :: back
      integer                        :: averify_aa
      logical                        :: dir_switch
      integer                        :: i,ls,lset


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = size(s)
      lset = size(set)
      if (dir_switch) then
!         backwards search
          do i=ls,1,-1
              if (.not.(any(set(1:lset) == s(i)))) then
                  averify_aa = i
                  return
              endif
          enddo
          averify_aa = 0
      else
!         forward search
          do i=1,ls
              if (.not.(any(set(1:lset) == s(i)))) then
                  averify_aa = i
                  return
              endif
          enddo
          averify_aa = 0
      endif

      end function averify_aa

!*******************************************************************************

      pure function averify_ac(s,set,back)

      implicit none
      character, intent(in)          :: s(:)
      character(*), intent(in)       :: set
      logical, intent(in), optional  :: back
      integer                        :: averify_ac
      logical                        :: dir_switch
      integer                        :: i,j,ls,lset


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = size(s)
      lset = len(set)
      if (dir_switch) then
!         backwards search
b:        do i=ls,1,-1
              do j=1,lset
                  if (s(i) == set(j:j)) cycle b
              enddo
              averify_ac = i
              return
          enddo b
          averify_ac = 0
      else
!         forward search
f:        do i=1,ls
              do j=1,lset
                  if (s(i) == set(j:j)) cycle f
              enddo
              averify_ac = i
              return
          enddo f
          averify_ac = 0
      endif

      end function averify_ac

!*******************************************************************************

      pure function averify_ca(s,set,back)

      implicit none
      character(*), intent(in)       :: s
      character, intent(in)          :: set(:)
      logical, intent(in), optional  :: back
      integer                        :: averify_ca
      logical                        :: dir_switch
      integer                        :: i,j,ls,lset


      if (present(back)) then
          dir_switch = back
      else
          dir_switch = .false.
      endif

      ls = len(s)
      lset = size(set)
      if (dir_switch) then
!         backwards search
b:        do i=ls,1,-1
              do j=1,lset
                  if (s(i:i) == set(j)) cycle b
              enddo
              averify_ca = i
              return
          enddo b
          averify_ca = 0
      else
!         forward search
f:        do i=1,ls
              do j=1,lset
                  if (s(i:i) == set(j)) cycle f
              enddo
              averify_ca = i
              return
          enddo f
          averify_ca = 0
      endif

      end function averify_ca

!*******************************************************************************
!     UPPERCASE
!*******************************************************************************

      pure function uppercase_s(s,begin,end)

      implicit none
      type(string), intent(in)       :: s
      integer, intent(in), optional  :: begin,end
      character(len(s))              :: uppercase_s
      integer                        :: i,i1,i2,j


      i1 = 1
      if (present(begin)) i1 = max(i1,begin)
      i2 = len(s)
      if (present(end)) i2 = min(i2,end)

      do i=1,i1-1
          uppercase_s(i:i) = s%chars(i)
      enddo
      do i=i1,i2
          j = iachar(s%chars(i))
          select case(j)
          case(97:122)
              uppercase_s(i:i) = achar(j-32)
          case default
              uppercase_s(1:i) = s%chars(i)
          end select
      enddo
      do i=i2+1,len(s)
          uppercase_s(i:i) = s%chars(i)
      enddo

      end function uppercase_s

!*******************************************************************************

      pure function uppercase_c(c,begin,end)

      implicit none
      character(*), intent(in)       :: c
      integer, intent(in), optional  :: begin,end
      character(len(c))              :: uppercase_c
      integer                        :: i,i1,i2,j


      i1 = 1
      if (present(begin)) i1 = max(i1,begin)
      i2 = len(c)
      if (present(end)) i2 = min(i2,end)

      uppercase_c(:i1-1) = c(:i1-1)
      do i=i1,i2
          j = iachar(c(i:i))
          select case(j)
          case(97:122)
              uppercase_c(i:i) = achar(j-32)
          case default
              uppercase_c(i:i) = c(i:i)
          end select
      enddo
      uppercase_c(i2+1:) = c(i2+1:)

      end function uppercase_c

!*******************************************************************************

      elemental subroutine to_uppercase_s(s,begin,end)

      implicit none
      type(string), intent(inout)    :: s
      integer, intent(in), optional  :: begin,end
      integer                        :: i,i1,i2,j


      i1 = 1
      if (present(begin)) i1 = max(i1,begin)
      i2 = len(s)
      if (present(end)) i2 = min(i2,end)

      do i=i1,i2
          j = iachar(s%chars(i))
          select case(j)
          case(97:122)
              s%chars(i) = achar(j-32)
          case default
              continue
          end select
      enddo

      end subroutine to_uppercase_s

!*******************************************************************************

      elemental subroutine to_uppercase_c(c,begin,end)

      implicit none
      character(*), intent(inout)    :: c
      integer, intent(in), optional  :: begin,end
      integer                        :: i,i1,i2,j


      i1 = 1
      if (present(begin)) i1 = max(i1,begin)
      i2 = len(c)
      if (present(end)) i2 = min(i2,end)

      do i=i1,i2
          j = iachar(c(i:i))
          select case(j)
          case(97:122)
              c(i:i) = achar(j-32)
          case default
              continue
          end select
      enddo

      end subroutine to_uppercase_c

!*******************************************************************************
!     LOWERCASE
!*******************************************************************************

      pure function lowercase_s(s,begin,end)

      implicit none
      type(string), intent(in)       :: s
      integer, intent(in), optional  :: begin,end
      character(len(s))              :: lowercase_s
      integer                        :: i,i1,i2,j


      i1 = 1
      if (present(begin)) i1 = max(i1,begin)
      i2 = len(s)
      if (present(end)) i2 = min(i2,end)

      do i=1,i1-1
          lowercase_s(i:i) = s%chars(i)
      enddo
      do i=i1,i2
          j = iachar(s%chars(i))
          select case(j)
          case(65:90)
              lowercase_s(i:i) = achar(j+32)
          case default
              lowercase_s(i:i) = s%chars(i)
          end select
      enddo
      do i=i2+1,len(s)
          lowercase_s(i:i) = s%chars(i)
      enddo

      end function lowercase_s

!*******************************************************************************

      pure function lowercase_c(c,begin,end)

      implicit none
      character(*), intent(in)       :: c
      integer, intent(in), optional  :: begin,end
      character(len(c))              :: lowercase_c
      integer                        :: i,i1,i2,j


      i1 = 1
      if (present(begin)) i1 = max(i1,begin)
      i2 = len(c)
      if (present(end)) i2 = min(i2,end)

      lowercase_c(:i1-1) = c(:i1-1)
      do i=i1,i2
          j = iachar(c(i:i))
          select case(j)
          case(65:90)
              lowercase_c(i:i) = achar(j+32)
          case default
              lowercase_c(i:i) = c(i:i)
          end select
      enddo
      lowercase_c(i2+1:) = c(i2+1:)

      end function lowercase_c

!*******************************************************************************

      elemental subroutine to_lowercase_s(s,begin,end)

      implicit none
      type(string), intent(inout)    :: s
      integer, intent(in), optional  :: begin,end
      integer                        :: i,i1,i2,j


      i1 = 1
      if (present(begin)) i1 = max(i1,begin)
      i2 = len(s)
      if (present(end)) i2 = min(i2,end)

      do i=i1,i2
          j = iachar(s%chars(i))
          select case(j)
          case(65:90)
              s%chars(i) = achar(j+32)
          case default
              continue
          end select
      enddo

      end subroutine to_lowercase_s

!*******************************************************************************

      elemental subroutine to_lowercase_c(c,begin,end)

      implicit none
      character(*), intent(inout)    :: c
      integer, intent(in), optional  :: begin,end
      integer                        :: i,i1,i2,j


      i1 = 1
      if (present(begin)) i1 = max(i1,begin)
      i2 = len(c)
      if (present(end)) i2 = min(i2,end)

      do i=i1,i2
          j = iachar(c(i:i))
          select case(j)
          case(65:90)
              c(i:i) = achar(j+32)
          case default
              continue
          end select
      enddo

      end subroutine to_lowercase_c

!*******************************************************************************

!*******************************************************************************

      end module m_strings
