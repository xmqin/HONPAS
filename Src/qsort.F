! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

C ==================================================================
C Reorder a vector of numbers in decreasing order using the
C quicksort algorithm. A permutation vector is also created.
C
C Written by Rogeli Grima (BSC) Dec.2007
C
C ==================================================================
C SUBROUTINE myQsort( N, A, perm )
C
C INPUT:
C integer   N       : Number of elements
C
C OUTPUT:
C integer   A(N)    : Vector to be reordered
C integer   perm(N) : Permutation vector
C
C BEHAVIOR:
C This is a recursive algorithm. At every step:
C   - Choose a pivot element (the first of the array)
C   - Move elements bigger than pivot to the begining of the array
C     and smaller ones at the end. Variable 'iq' marks the frontier
C     between both parts
C   - Reorder every part calling recursively myQsort
C
C ==================================================================
      recursive subroutine myQsort( N, A, perm )
      implicit none
C     Input variables
      integer       :: N, A(N), perm(N)
C     Local variables
      integer       :: iq

!------------------------------------------------------------------------- BEGIN
      if (N.gt.1) then
        call Partition( N, A, perm, iq )
        call myQsort( iq-1, A, perm )
        call myQsort( N-iq+1, A(iq), perm(iq) )
      endif
!--------------------------------------------------------------------------- END
      end subroutine myQsort

C ==================================================================
C Set the first element of an array as a pivot and moves values
C greater than or equal to the pivot element to the begining of the
C array. It also returns the position of the intersection between
C both groups. The array perm suffer the same reordering.
C ==================================================================
C SUBROUTINE Partition( N, A, perm, marker )
C
C INPUT:
C integer   N       : Number of elements
C
C OUTPUT:
C integer   A(N)    : Vector to be reordered
C integer   perm(N) : Permutation vector
C integer   marker  : frontier between "big" values and "small" values
C
C BEHAVIOR:
C Choose the first element of the array as pivot.
C At every iteration:
C   - Find a value bigger than the pivot starting from the end
C     of the array.
C   - Find a value smaller than the pivot starting from the begining
c     of the array
C   - Swap both elements
C Quit once we have checked all values.
C
C NOTE: Values equal to pivot can be at any side.
C
C ==================================================================
      subroutine Partition( N, A, perm, marker )
C     Input variables
      integer       :: N, A(N),  perm(N), marker
C     Local variables
      integer       :: i, j, temp, x

!------------------------------------------------------------------------- BEGIN
      x = A(1)
      i = 0
      j = N + 1
      do
        j = j-1
        do
          if (A(j) >= x) exit
          j = j-1
        end do
        i = i+1
        do
          if (A(i) <= x) exit
          i = i+1
        end do
        if (i < j) then
C         exchange A(i) and A(j)
          temp    = A(i)
          A(i)    = A(j)
          A(j)    = temp
          temp    = perm(i)
          perm(i) = perm(j)
          perm(j) = temp
        elseif (i == j) then
          marker = i+1
          return
        else
          marker = i
          return
        endif
      end do
!--------------------------------------------------------------------------- END
      end subroutine Partition
