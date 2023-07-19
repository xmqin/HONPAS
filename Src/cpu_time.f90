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

