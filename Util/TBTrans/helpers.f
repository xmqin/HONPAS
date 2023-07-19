C ====================================================================
C Write out a real matrix in a format so that it can be plotted
C directly with Mathematica using e.g. ListDensityPlot
C ====================================================================


      subroutine writermat(N,name,X)
      implicit none
      
      integer N
      real*8 X(N,N)
      character*4 name ! Use 4 chars as a name

      integer i,j,unit1
      external io_assign,io_close


c
c     BEGIN
c
      call io_assign(unit1)

      open (unit1,
     &     file=name,
     &     status='unknown')
      

      write(unit1,*) '   '
      write(unit1,*) name,' ={'
      do j=N,2,-1
         write(unit1,*) '{'
         write(unit1,'(F12.8,a1)') (X(i,j),',',i=1,N-1)
         write(unit1,'(F12.8,a2)') (X(i,j),'},',
     .        i=N,N)
      end do  
      j=1
      write(unit1,*) '{'
      write(unit1,'(F12.8,a1)') 
     .     (X(i,j),',',i=1,N-1)
      write(unit1,'(F12.8,a1)') 
     .     (X(i,j),'}',i=N,N)
      write(unit1,*) '};'
         


      call io_close(unit1)

      return
      end

C ====================================================================
C Write out a complex matrix in a format so that it can be plotted
C directly with Mathematica using e.g. ListDensityPlot
C ====================================================================


      subroutine writezmat(N,name,Z)
      implicit none
      
      integer N
      complex*16 Z(N,N)
      character*4 name ! Use 4 chars as a name

      integer i,j,unit1
      external io_assign,io_close


c
c     BEGIN
c
      call io_assign(unit1)

      open (unit1,
     &     file=name,
     &     status='unknown')
      

      write(unit1,*) '   '
      write(unit1,*) name,' ={'
      do j=N,2,-1
         write(unit1,*) '{'
         write(unit1,'(F12.8,a3,F12.8,a1)') 
     .        (dreal(Z(i,j)),'+I*',dimag(Z(i,j)),',',i=1,N-1)
         write(unit1,'(F12.8,a3,F12.8,a2)') (dreal(Z(i,j)),'+I*',
     .        dimag(Z(i,j)),'},',
     .        i=N,N)
      end do  
      j=1
      write(unit1,*) '{'
      write(unit1,'(F12.8,a3,F12.8,a1)') 
     .     (dreal(Z(i,j)),'+I*',dimag(Z(i,j)),',',i=1,N-1)
      write(unit1,'(F12.8,a3,F12.8,a1)') 
     .     (dreal(Z(i,j)),'+I*',dimag(Z(i,j)),'}',
     .     i=N,N)
      write(unit1,*) '};'
         


      call io_close(unit1)

      return
      end



c======================================================================
      subroutine writerealmatrixmath(n,M,Name)

      implicit none

      integer i,j,n,unit1
      complex*16 M(n,n)
      character*3 Name
      external io_assign,io_close

c
c     BEGIN
c
      call io_assign(unit1)
      OPEN(unit1,FILE='REmatrix.ma', STATUS='UNKNOWN')
      write(unit1,*) Name,'={'
      do i=1,n-1
         write(unit1,*) '{'
         do j=1,n-1
            write(unit1,'(f11.7,1a)') DREAL(M(i,j)),','
         end do                 !j
         write(unit1,'(f11.7,1a)') DREAL(M(i,n)),'},'
      end do                    !i
      write(unit1,*) '{'
      do j=1,n-1
         write(unit1,'(f11.7,1a)') DREAL(M(i,j)),','
      end do                    !j
      write(unit1,'(f11.7,1a)') DREAL(M(i,n)),'}'         
      write(unit1,*) '}'
      call io_close(unit1)
c======================================================================
      return 
      end
c======================================================================
c
c
c
c======================================================================
      subroutine writeimagmatrixmath(n,M,Name)

      implicit none

      integer i,j,n,unit1
      complex*16 M(n,n)
      character*3 Name

      call io_assign(unit1)
      OPEN(unit1,FILE='IMmatrix.ma', STATUS='UNKNOWN')
      write(unit1,*) Name,'={'
      do i=1,n-1
         write(unit1,*) '{'
         do j=1,n-1
            write(unit1,'(f11.7,1a)') DIMAG(M(i,j)),','
         end do                 !j
         write(unit1,'(f11.7,1a)') DIMAG(M(i,n)),'},'
      end do                    !i
      write(unit1,*) '{'
      do j=1,n-1
         write(unit1,'(f11.7,1a)') DIMAG(M(i,j)),','
      end do                    !j
      write(unit1,'(f11.7,1a)') DIMAG(M(i,n)),'}'         
      write(unit1,*) '}'
      call io_close(unit1)
c======================================================================
      return 
      end
c======================================================================



C ====================================================================
      SUBROUTINE writemolden(r,kind,N,maxN)
C ====================================================================
      implicit none

      INTEGER N,maxN 
      real*8 BOHR
      parameter(BOHR=0.529d0)
      real*8 r(maxN,3)
      character*2 elem
      INTEGER i, kind(maxN),unit1

      call io_assign(unit1)             
      OPEN(unit1,FILE="struc.xyz", STATUS='UNKNOWN')


      write(unit1,'(i5)') N
      write(unit1,*) '  '
      do i=1,N
         if(kind(i).EQ.79) elem='Au'
         if(kind(i).EQ.14) elem='Si'
         write(unit1,100)
     &        elem,r(i,1)*BOHR,r(i,2)*BOHR,r(i,3)*BOHR
      end do
 100  format(a2,3(f14.6))
      call io_close(unit1)
      RETURN
C ====================================================================
      END                       ! of write molden
C ====================================================================



