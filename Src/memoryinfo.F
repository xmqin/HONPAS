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
      module memoryinfo

C
C WordsInteger = number of words of memory for integer variables
C WordsSP      = number of words of memory for single precision variables
C WordsDP      = number of words of memory for double precision variables
C WordsLogical = number of words of memory for logical variables
C WordsSC      = number of words of memory for single complex variables
C WordsDC      = number of words of memory for double complex variables
C PeakMemory   = maximum amount of dynamic memory used
C PeakRoutine  = routine name where memory peak was reached
C ByteSize     = array of values specifying the size of a word for
C                each data type
C
      implicit none

      integer, save ::
     .  WordsInteger, WordsSP, WordsDP, WordsLogical, WordsSC, WordsDC, 
     .  WordsString,
     .  PeakMemory, CurrentMemory

      character(len=30), save ::
     .  PeakRoutine

      data
     .  WordsInteger / 0 /,
     .  WordsSP / 0 /,
     .  WordsDP / 0 /,
     .  WordsLogical / 0 /,
     .  WordsSC / 0 /,
     .  WordsDC / 0 /,
     .  WordsString / 0 /,
     .  PeakMemory / 0 /,
     .  CurrentMemory / 0 /
      
      integer, parameter, dimension(7) :: ByteSize = 
#ifdef OLD_CRAY
     .     (/8,8,8,8,8,16,8/)
#else
     .     (/4,4,8,4,8,16,4/)
#endif

      end module memoryinfo

      module diagmemory
C
C  Stores the factor used to scale the default memory in rdiag/cdiag
C  By increasing this value it is possible to avoid failure to
C  converge eigenvalues.
C
C  real*8  MemoryFactor      : factor by which memory is scaled
C
      use precision, only : dp

      implicit none

      real(dp), save :: MemoryFactor

      end module diagmemory
