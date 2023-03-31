! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
program kind_explorer

!
! Simple program to detect the kinds of real numbers
! Standard output: one or two integers, corresponding
!                  to the "single" and "double precision"
!
integer, parameter  ::  sp = selected_real_kind(6,34)
integer, parameter  ::  dp = selected_real_kind(15,300)

integer, parameter   ::  default_p = kind(1.0)

if (sp > 0) then
   if (dp > 0) then
      if (sp /= dp) then
         write(unit=*,fmt=*) sp, dp
         write(unit=0,fmt=*) "Single: ", sp, " Double: ", dp
      else
         write(unit=*,fmt=*) dp
         write(unit=0,fmt=*) "In your computer sp=dp"
      endif
   else
      write(unit=*,fmt=*) sp
      write(unit=0,fmt=*) "Your computer does not have double prec reals!"
   endif
else
   if (dp > 0) then
      write(unit=*,fmt=*) dp
      write(unit=0,fmt=*) "Your computer does not have single prec reals!"
   else
      write(unit=0,fmt=*) "Your computer's real number system is weird"
   endif
endif
write(unit=0,fmt=*) "Your computer's default precision real kind is ", default_p
      

end program kind_explorer

