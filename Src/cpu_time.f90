! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
! Emulation of the F95 standard cputime, for those compilers
! that do not support it (such as pgf90)
!
! This version for BSD-style system call
!
subroutine cpu_time(t)
real, intent(out)   :: t

real tarray(2)
external etime
real etime

t = etime(tarray)

end subroutine cpu_time

