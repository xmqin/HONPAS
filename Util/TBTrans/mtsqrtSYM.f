      PROGRAM matrixroot
      implicit none

      integer N
      parameter(N=3)
      integer i,j,k
      complex*16 zwork(2*N-1),smat(N*(N+1)/2)
      real*8 rwork(3*N-2)
      complex*16 M(N,N), SQM(N,N)

      complex*16 V(N,N),csum
      
      do i=1,N
         do k=1,N
            V(i,k) = dcmplx(0d0,1d0*k)
            if(i.eq.k) V(i,i) = dcmplx(1d0*i, 0d0)
            
         end do 
      end do


      do i=1,N
         do j=1,N
            csum=dcmplx(0d0,0d0)
            do k=1,N
               csum = csum + conjg(V(i,k))*V(j,k)
            end do              !k
            M(i,j)=csum
         end do
      end do 
      

      write(6,*) 'M:'
      do i=1,N
         write(6,*) (M(i,j),j=1,N)
      end do 

      call sqrtM(N,M,smat,zwork,rwork,SQM)

      write(6,*) 'sqrt(M):'
      do i=1,N
         write(6,*) (SQM(i,j),j=1,N)
      end do 


      CALL zgemm('N','N',N,N,N,dcmplx(1d0,0d0),
     &     SQM,N,SQM,N,dcmplx(0d0,0d0),M,N)
      


      write(6,*) 'sqrt(M)*sqrt(M):'
      do i=1,N
         write(6,*) (M(i,j),j=1,N)
      end do 


      

      END


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     calculate square root of hermitian matrix 
c

      subroutine sqrtM(n,M,smat,zwork,rwork,SQM)
      implicit none
c IN
      integer n
      complex*16 M(n,n),eig(n)

C LAPACK diagononalization:
      complex*16 zwork(2*n-1)
      real*8 rwork(3*n-2)
      integer info
c helpers
      integer i,j,k
      complex*16 csum 
      complex*16 smat(n*(n+1)/2)
      real*8 seig(n)
      complex*16 svect(n,n)

      
c return
      complex*16 SQM(n,n)
      

c BEGIN


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
      
      write(6,*) 'Eigenvalues of M = '
      write(6,*) (seig(i),i=1,n)
      

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
      

      RETURN
      END




