cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     calculate square root of hermitian matrix 
c

      subroutine zsqrtM(n,M,SQM)
      implicit none
c IN
      integer n
      complex*16 M(n,n)

c OUT
      complex*16 SQM(n,n)

C LAPACK diagononalization:
      complex*16, dimension(:), allocatable:: zwork, smat
      complex*16, dimension(:,:), allocatable:: svect
      real*8, dimension(:), allocatable:: rwork, seig
      integer info

c helpers
      integer i,j,k
      complex*16 csum 
c      complex*16 smat(n*(n+1)/2)
c      real*8 seig(n)
c      complex*16 svect(n,n)

      
      

c BEGIN

      allocate(zwork(2*n-1))
      allocate(smat(n*(n+1)/2))
      allocate(rwork(3*n-2))
      allocate(svect(n,n))
      allocate(seig(n))




c upper triangular form
      do j=1,n
         do i=1,n
            smat(i + ((j-1)*j/2)) = M(i,j)
         end do                 !i
      end do                    !j
      
c     Diagonalize 
      call zhpev('V','U',n,smat,seig,svect,n,zwork,rwork,info)
      if(info.ne.0) then
         write(6,*) 'INFO = ',info,' when diagonalizing in sqrtM'
      end if
      
c      write(6,*) 'Eigenvalues of M = '
c      write(6,*) (seig(i),i=1,n)
      

c
c     There are faster ways of doing this, but let's play safe for
c     now.  Form the sqM matrix
c

      do i = 1,n
         if(seig(i) .lt. 0d0) seig(i) = 0d0
         seig(i) = dsqrt(seig(i))
      end do

      do i = 1,n
         do j = 1,i
            csum = dcmplx(0d0,0d0)
            do k = 1,n
               csum = csum + seig(k)*svect(i,k)*dconjg(svect(j,k))
            end do
            sqM(i,j) = csum
            sqM(j,i) = dconjg(sqM(i,j))
         end do
      end do


      deallocate(zwork)
      deallocate(rwork)
      deallocate(smat)
      deallocate(svect)
      deallocate(seig)

      RETURN
      END



