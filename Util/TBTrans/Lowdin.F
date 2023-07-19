        
C ##################################################################
C ##                Make Löwdin transformation                    ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##################################################################

      subroutine Lowdin(flag,n,M,T)
c     iop=0  T is actually the S matrix itself:
c
      
      implicit none

c INPUT
      integer n
      logical flag

c INPUT/OUTPUT     

      complex*16 M(n,n)         ! on return: T.M.T where T=S**(-1/2)

c flag = T: T = S and on return T=S**(-1/2)
c flag = F: T = S**(-1/2) and is not changed

      complex*16 T(n,n)     

c helpers
      complex*16, dimension(:,:), allocatable :: X
      complex*16 a,b

      integer i,j
      
      allocate(X(n,n))

      a = dcmplx(1d0,0d0)
      b = dcmplx(0d0,0d0)


      if(flag) then
         do j=1,n     
            do i=1,n
               X(i,j)=T(i,j)
            end do
         end do
         call zinvsqrtM(n,X,T)
      end if

      call zgemm('N','N',n,n,n,a,T,n,M,n,b,X,n)      
      call zgemm('N','N',n,n,n,a,X,n,T,n,b,M,n)
      
      deallocate(X)

      return
      end



C ##################################################################
C ##         Calculate square root of hermitian matrix            ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##################################################################


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



C ##################################################################
C ##     Calculate inverse square root of hermitian matrix        ##
C ##                            By                                ##
C ##              Mads Brandbyge, mbr@mic.dtu.dk                  ##
C ##################################################################


      subroutine zinvsqrtM(n,M,SQM)
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
      real*8 eps
      parameter(eps=1.e-10)
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
c      write(6,*)      

c
c     There are faster ways of doing this, but let's play safe for
c     now.  Form the sqM matrix
c

      do i = 1,n
         if(seig(i) .lt. eps) then 
            seig(i)=0d0
         else 
            seig(i) = 1.d0/dsqrt(seig(i))
         end if
      end do
      
      do i = 1,n
         do j = 1,i
            csum = dcmplx(0d0,0d0)
            do k = 1,n
               csum = csum + seig(k)*svect(i,k)*dconjg(svect(j,k))
            end do
            SQM(i,j) = csum
            SQM(j,i) = dconjg(SQM(i,j))
         end do
      end do


      deallocate(zwork)
      deallocate(rwork)
      deallocate(smat)
      deallocate(svect)
      deallocate(seig)

      RETURN
      END

