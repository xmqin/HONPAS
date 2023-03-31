! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
program int_explorer

!
! Standard output: one or two integers, corresponding
!                  to the "standard" and "8-byte" integers
!
integer, parameter  ::  i8b = selected_int_kind(12)
integer, parameter  ::  idef = kind(1)

if (idef > 0) then
   if (i8b > 0) then
      if (idef /= i8b) then
         write(unit=*,fmt=*) idef, i8b
         write(unit=0,fmt=*) "int def: ", idef, " 8-byte int kind: ", i8b
      else
         write(unit=*,fmt=*) idef
         write(unit=0,fmt=*) "In your computer the standard int is 8-byte"
      endif
   else
      write(unit=*,fmt=*) idef
      write(unit=0,fmt=*) "Your computer does not have 8-byte ints!"
      STOP
   endif
else
      write(unit=0,fmt=*) "Your computer's int number system is weird..."
      STOP
endif

write(unit=0,fmt=*) "Your computer's default precision int kind is: ", idef
      

end program int_explorer

