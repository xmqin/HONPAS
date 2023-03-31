! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      subroutine matvect(matrix,vector,product)
*
*     This subroutine calculates the vector result of the product of a 
*    matrix (3,3) and a vector.
*
      use precision
      implicit none
      integer i,j
      real(dp) matrix(3,3)
      real(dp) vector(3),product(3)

      do i = 1,3
	 product(i) = 0.0_dp
	 do j = 1,3
	    product(i) = product(i) + matrix(i,j) * vector(j)
         enddo
      enddo

      return
      end

