! 
! This file is part of the SIESTA package.
!
! Copyright (c) Fundacion General Universidad Autonoma de Madrid:
! E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
! and J.M.Soler, 1996- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
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

