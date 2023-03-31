! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
      module m_memory

      integer, public :: mem_stat    ! For use in calls to allocate 
                                     ! and deallocate

      public :: memory
      private

      CONTAINS

      subroutine memory( Task, Type, NElements, CallingRoutine,
     $                   stat,id)
C 
C This subroutine keeps track of information relating to the use 
C of dynamic memory
C
C Input :
C
C character*1 Task  : job type = 'A' -> allocate
C                   :            'D' -> deallocate
C character*1 Type  : type of variable = 'I' = integer
C                   :                    'S' = single precision real
C                   :                    'D' = double precision real
C                   :                    'X' = grid precision real
C                   :                    'L' = logical
C                   :                    'C' = single precision complex
C                   :                    'Z' = double precision complex
C                   :                    'S' = character data (we assume takes one word)
C                   :                    'E' = double precision integer
C integer NElements : number of array elements being 
C                   : allocated/deallocated
C character         :
C   CallingRoutine  : string containing the name of the calling routine
C
C Created by J.D. Gale, October 1999
C
C Stat and ID keyword arguments added by Alberto Garcia, 2005

      use memoryinfo
      use sys, only: die
      use parallel, only: Node

      implicit none

      integer, intent(in)                 :: NElements
      character(len=1), intent(in)        :: Task, Type
      character(len=*), intent(in)        :: CallingRoutine
      integer, intent(in), optional       :: stat
      character(len=*), intent(in), optional    :: id

C Local variables
      integer    Sign
C
C CAG
C       if (Node.eq.0) then
C         write(66,'(a,a,a,a,i10)') "---memory: ", CallingRoutine // " ",
C      $                          task // " ",
C      $                          type // " ", nelements
C         call flush(66)
C       endif
      if (present(stat)) then
         if (stat .ne. 0) then
            if (present(id)) then
               call die(Task // "-llocation failed in " //
     $                  CallingRoutine // id)
            else
               call die(Task // "-llocation failed in " //
     $                  CallingRoutine)
            endif
         endif
      endif

C Assign sign of operation based on task
      if (Task .eq. 'A' .or. Task .eq. 'a') then
        Sign = +1
      else
        Sign = -1
      endif

C Add number of elements to appropriate pointer
      if (Type .eq. 'I' .or. Type .eq. 'i') then
        WordsInteger = WordsInteger + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(1)
      elseif (Type .eq. 'S' .or. Type .eq. 's') then
        WordsSP = WordsSP + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(2)
      elseif (Type .eq. 'D' .or. Type .eq. 'd') then
        WordsDP = WordsDP + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(3)
!
!     Support for flexible precision in grid-related arrays
!
      elseif (Type .eq. 'X' .or. Type .eq. 'x') then
#ifdef GRID_DP
        WordsDP = WordsDP + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(3)
#else
        WordsSP = WordsSP + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(2)
#endif

      elseif (Type .eq. 'L' .or. Type .eq. 'l') then
        WordsLogical = WordsLogical + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(4)
      elseif (Type .eq. 'C' .or. Type .eq. 'c') then
        WordsSC = WordsSC + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(5)
      elseif (Type .eq. 'Z' .or. Type .eq. 'z') then
        WordsDC = WordsDC + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(6)
      elseif (Type .eq. 'S' .or. Type .eq. 's') then
        WordsString = WordsString + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(7)
      elseif (Type .eq. 'E' .or. Type .eq. 'e') then
        WordsLongInt  = WordsLongInt + Sign*NElements
        CurrentMemory = CurrentMemory + Sign*NElements*ByteSize(8)
      endif

C Check whether memory is greater than peak memory so far
      if (CurrentMemory .gt. PeakMemory) then
        PeakMemory = CurrentMemory
        PeakRoutine = CallingRoutine
      endif

      end subroutine memory

      end module m_memory
