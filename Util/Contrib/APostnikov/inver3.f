      subroutine inver3(a,b)
C
C     Inverts a 3x3 matrix

      implicit none

      integer           i
      double precision  a(3,3), b(3,3), c

      b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      b(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      b(1,3) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
      b(2,1) = a(2,3)*a(3,1) - a(3,3)*a(2,1)
      b(2,2) = a(3,3)*a(1,1) - a(1,3)*a(3,1)
      b(2,3) = a(1,3)*a(2,1) - a(2,3)*a(1,1)
      b(3,1) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
      b(3,2) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
      b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
      do i = 1, 3
         c=1.d0/(a(1,i)*b(i,1) + a(2,i)*b(i,2) + a(3,i)*b(i,3) )
         b(i,1)=b(i,1)*c
         b(i,2)=b(i,2)*c
         b(i,3)=b(i,3)*c
      enddo
      end
