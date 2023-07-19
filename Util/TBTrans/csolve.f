      subroutine csolve(N,A,B,ipvt,joutfile)

c INPUT
      integer N            ! problem size
      complex*16 A(N,N)
      integer joutfile          ! file unit for write out

c INPUT/OUTPUT
      complex*16 B(N,N)

c OUTPUT
      integer ipvt(N)           !pivoting vector

c helpers, tempos ...
      integer INVINFO           ! Used by LAPACK mtr. inv.

      call  csolveg(N,N,A,B,ipvt,joutfile)

      return 
      end

      subroutine csolveg(N,NRHS,A,B,ipvt,joutfile)

      implicit none

c Finds the solution to a general complex N x N matrix equation by calling
c lapack or ESSL.
c
c     A.X = B;   returns X = Inv(A).B
c

c INPUT
      integer N,NRHS            ! problem size
      complex*16 A(N,N)
      integer joutfile          ! file unit for write out
      integer i,j

c INPUT/OUTPUT
      complex*16 B(N,NRHS)

c OUTPUT
      integer ipvt(N)           !pivoting vector

c helpers, tempos ...
      integer INVINFO           ! Used by LAPACK mtr. inv.
      
c BEGIN



c ESSL: 
c      CALL zgef(A,N,N,ipvt)     ! factorize A for inv.
c      CALL zgesm('N',A,N,N,ipvt,B,N,NRHS) ! inversion
c -----------
c LAPACK:

      INVINFO=0
      CALL zgesv(N,NRHS,A,N,ipvt,B,N,INVINFO) 
      IF(INVINFO.NE.0) THEN
         write(6,*) 'ERROR: MATRIX INVERSION FAILED'
         write(6,*) 'ERROR: LAPACK INFO = ',INVINFO
         write(joutfile,*) 'ERROR: MATRIX INVERSION FAILED'
         write(joutfile,*) 'ERROR: LAPACK INFO = ',INVINFO
c          do i =1,N 
c           do j =1,N
c              write(92,*)i,j,a(i,j)
c            enddo
c          enddo
         stop ' LAPACK zgesv error in csolveg'
      END IF
c -----------       

      return 
      end
