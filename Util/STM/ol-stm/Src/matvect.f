      subroutine matvect(matrix,vector,product)
*
*     This subroutine calculates the vector result of the product of a 
*    matrix (3,3) and a vector.
*
      implicit none
      integer i,j
      real*8 matrix(3,3)
      real*8 vector(3),product(3)

      do i = 1,3
	 product(i) = 0.d0
	 do j = 1,3
	    product(i) = product(i) + matrix(i,j) * vector(j)
         enddo
      enddo

      return
      end

