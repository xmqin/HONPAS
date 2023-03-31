! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
! Random number generator based on the subtractive method proposed
! by D. Knuth as explained for the routine ran3 in 'Numerical Recipes, 
! The Art of Scientific Computing' by W.H. Press, S.A. Teukolsky, 
! W.T. Veterling and B.P. Flannery, Cambridge U.P. 1987 and 1992.

function randomg(ipar) result(ran3)

integer, parameter :: dp = selected_real_kind(10,100)

integer, intent(inout)  :: ipar
real(dp)                :: ran3

! --------- internal
integer,  parameter :: tt = 0, bignum = 1000000000
integer,  parameter :: seedl=161803398, nseq=55
real(dp), parameter :: smallnum = 1./bignum
integer,  save      :: first=0, count1, count2, lista(nseq)
integer             :: i1, i2, k, j, kk

! --------- smallnum should be 1._dp/bignum, kept as is for backward check

! --------- initialization (seed)

if ((ipar < 0) .or. (first == 0)) then
   j = seedl - iabs(ipar)
   j = mod(j,bignum)
   lista(nseq) = j
   first = 1

   kk = 1
   do i1 = 1, nseq - 1
      i2 = mod(21*i1,nseq)
      lista(i2) = kk
      kk = j - kk
      j = lista(i2)
      if (kk < tt) kk = kk + bignum 
   end do 

   do k = 1, 4
      do i1 = 1, nseq
         lista(i1) = lista(i1) - lista(1 + mod(i1+30,nseq))
         if (lista(i1) < tt) lista(i1) = lista(i1) + bignum 
      end do  ! i1
   end do  ! k

   ipar = 1
   count1 = 0
   count2 = 31
endif

! --------- generation

count1 = count1 + 1; count2 = count2 + 1
if (count1 == (nseq + 1)) count1 = 1
if (count2 == (nseq + 1)) count2 = 1

j = lista(count1) - lista(count2)
if (j < tt) j = j + bignum
lista(count1) = j

ran3 = j*smallnum

end function randomg
