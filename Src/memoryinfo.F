! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module memoryinfo

      use precision, only: i8b

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

      integer(i8b), save ::
     .  WordsInteger, WordsSP, WordsDP, WordsLogical, WordsSC, WordsDC, 
     .  WordsString, WordsLongInt,
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
     .  WordsLongInt / 0 /,
     .  PeakMemory / 0 /,
     .  CurrentMemory / 0 /
      
      integer, parameter, dimension(8) :: ByteSize = 
#ifdef OLD_CRAY
     .     (/8,8,8,8,8,16,8,8/)
#else
     .     (/4,4,8,4,8,16,4,8/)
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
